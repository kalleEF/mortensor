#' Master pipeline functions for cause-specific mortality forecasting
#'
#' @name pipeline
NULL

# ---- Internal helpers ---------------------------------------------------

#' Append group columns to a data.table row (generic)
#' @keywords internal
append_group_cols <- function(row_dt, groups, k) {
  for (g in names(groups))
    data.table::set(row_dt, j = g, value = groups[[g]][k])
  row_dt
}

#' Build a row filter expression matching group k
#' @keywords internal
group_filter <- function(dt, groups, k) {
  mask <- rep(TRUE, nrow(dt))
  for (g in names(groups))
    mask <- mask & dt[[g]] == groups[[g]][k]
  mask
}

# ---- Run single-cause forecast pipeline ---------------------------------

#' Run the full mortality forecast pipeline for a single cause
#'
#' Executes 17 stages: data loading, array building, smoothing, Tucker3
#' calibration, model calibration, holdout comparison, model selection /
#' ensemble, production forecast, bootstrap PIs, PI re-centering, coherence
#' enforcement, fitted values, plots, disaggregation, diagnostics, and output
#' saving.
#'
#' @param data Either a [mort_data] object (generic path) or a [mort_config]
#'   object (ONS path, backward compatible).
#' @param config Optional [mort_config] object providing pipeline settings
#'   (forecast horizon, number of bootstrap simulations, etc.). Required when
#'   `data` is a `mort_data` object. Ignored when `data` is itself a
#'   `mort_config`.
#' @return Invisibly, a list with all results (data, forecasts, PIs, tables).
#' @export
run_forecast <- function(data, config = NULL) {
  # --- Resolve data vs config ---
  if (inherits(data, "mort_config")) {
    config <- data
    # ONS path: load from files
    dt <- load_and_prepare_ons(config)
    cause_short <- config$cause_short
    ages   <- config$age_starts
    years  <- sort(unique(dt$year))
    sexes  <- c("Male", "Female")
    quints <- seq_len(5L)
    groups <- as.data.frame(data.table::CJ(sex = sexes, quintile = quints))

    mx_arr  <- build_array(dt, paste0("mx_", cause_short), ages, years, groups)
    ac_arr  <- if (config$is_cause_specific) {
      build_array(dt, "mx_all", ages, years, groups)
    }
    pop_arr <- build_array(dt, "population", ages, years, groups, floor = NULL)
  } else if (inherits(data, "mort_data")) {
    if (is.null(config)) stop("'config' is required when 'data' is a mort_data object.")
    mx_arr  <- data$mx
    pop_arr <- data$pop
    ages    <- data$ages
    years   <- data$years
    groups  <- data$groups
    ac_arr  <- NULL
    cause_short <- config$cause_short
    dt <- NULL
  } else {
    stop("'data' must be a 'mort_data' or 'mort_config' object.")
  }

  set.seed(config$seed)
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("  MORTALITY FORECAST:", config$cause_label, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  cause_short <- config$cause_short

  # ---- 2. BUILD ARRAYS ----
  cat("\n=== Building arrays ===\n")
  n_a <- length(ages); n_t <- length(years); n_g <- nrow(groups)
  pop_arr[is.na(pop_arr) | pop_arr <= 0] <- 1

  cat("  Dimensions:", n_a, "ages x", n_t, "years x", n_g, "groups\n")

  # ---- 3. PRE-SMOOTHING ----
  cat("\n=== Pre-smoothing (", config$smoothing_method, ") ===\n")
  smooth_res <- smooth_mortality(mx_arr, method = config$smoothing_method,
                                  ages = ages, years = years,
                                  covid_years = config$covid_years %||% integer(0))
  mx_smooth <- smooth_res$arr
  lambda_used <- smooth_res$lambda
  gam_fits <- smooth_res$gam_fits  # stored for GAM-Bayes PIs and GAM disagg

  if (config$is_cause_specific && !is.null(ac_arr)) {
    ac_smooth <- smooth_mortality(ac_arr, method = config$smoothing_method,
                                   ages = ages, years = years,
                                   covid_years = config$covid_years %||% integer(0))$arr
  }

  smooth_check <- check_smoothing(mx_arr, mx_smooth, years,
                                  config$covid_years)
  plot_smoothing_comparison(mx_arr, mx_smooth, years, ages,
                            config$output_dir, config$cause_short,
                            config$covid_years)

  # ---- 4. CALIBRATE TUCKER3 ----
  cat("\n=== Calibrating Tucker3 parameters ===\n")
  h_tucker3 <- length(config$holdout_test)
  min_train_tucker3 <- max(8L, length(config$holdout_train) - h_tucker3)
  calib_coda <- calibrate_tucker3(mx_smooth, years, ages,
                                  method = "coda", width = config$age_width,
                                  h_val = h_tucker3,
                                  min_train = min_train_tucker3)
  calib_direct <- calibrate_tucker3(mx_smooth, years, ages,
                                    method = "direct", width = config$age_width,
                                    h_val = h_tucker3,
                                    min_train = min_train_tucker3)

  # ---- 5. CALIBRATE ALL OTHER MODELS ----
  calib_all <- NULL
  if (isTRUE(config$calibrate_models)) {
    cat("\n=== Calibrating all models (expanding-window CV) ===\n")
    calib_all <- calibrate_all_models(mx_smooth, pop_arr, years, ages, groups,
                                      h_val = h_tucker3,
                                      n_origins = config$cv_n_origins,
                                      width = config$age_width)
  }

  # ---- 6. HOLDOUT TEST ----
  cat("\n=== Holdout test ===\n")
  train_idx <- which(years %in% config$holdout_train)
  test_idx  <- which(years %in% config$holdout_test)
  if (length(test_idx) == 0) {
    cat("  WARNING: no test years found! Using last 3 years.\n")
    test_idx  <- (n_t - 2):n_t
    train_idx <- 1:(n_t - 3)
  }

  comp_res <- compare_models(mx_smooth, pop_arr, years, ages, groups,
                             train_idx, test_idx,
                             calib_coda, calib_direct,
                             width = config$age_width,
                             mx_eval_ref = mx_arr,
                             calib_all = calib_all)
  comp_dt <- comp_res$results

  cat("\n  Holdout comparison:\n")
  print(comp_dt[order(comp_dt$mase_combined)])
  data.table::fwrite(comp_dt, file.path(config$output_dir,
                     paste0("validation_results_", cause_short, ".csv")))
  plot_model_comparison(comp_dt, config$output_dir, config$cause_label,
                        config$cause_short)

  # ---- 7. SELECT BEST (or ENSEMBLE) ----
  ens_weights <- comp_res$ens_best_weights
  ens_optimal_names <- comp_res$ens_best_names
  ens_optimal_k <- comp_res$ens_best_k

  use_ensemble <- isTRUE(config$ensemble_mode) &&
    "Ensemble" %in% comp_dt$model &&
    comp_dt[model == "Ensemble" & !is.na(mase_combined), .N] > 0

  if (use_ensemble) {
    ens_comb <- comp_dt[model == "Ensemble", mase_combined]
    best_single_eligible <- comp_dt[model != "Ensemble" &
                                      !is.na(mase_combined) &
                                      (complete == TRUE | is.na(complete))]
    if (nrow(best_single_eligible) == 0)
      best_single_eligible <- comp_dt[model != "Ensemble" &
                                        !is.na(mase_combined)]
    best_single <- best_single_eligible[which.min(mase_combined)]
    if (is.finite(ens_comb) && ens_comb <= best_single$mase_combined * 1.05) {
      best_model <- "Ensemble"
      cat("  Best model: Ensemble (combined =", round(ens_comb, 4), ")\n")
    } else {
      use_ensemble <- FALSE
    }
  }

  if (!use_ensemble) {
    eligible <- comp_dt[model != "Ensemble" & !is.na(mase_combined) &
                          (complete == TRUE | is.na(complete))]
    if (nrow(eligible) == 0)
      eligible <- comp_dt[model != "Ensemble" & !is.na(mase_combined)]
    best_row <- eligible[which.min(mase_combined)]
    best_model <- best_row$model
    cat("  Best model:", best_model, "(combined MASE =",
        round(best_row$mase_combined, 4), ")\n")
  }

  # ---- 8. FULL FORECAST ----
  cat("\n=== Full forecast (h =", config$forecast_h, ") ===\n")
  forecast_years <- (max(years) + 1):(max(years) + config$forecast_h)
  mx_raw <- mx_arr

  if (use_ensemble) {
    ens_names <- ens_optimal_names
    if (!is.null(ens_weights) && length(ens_weights) > 0) {
      ew <- as.numeric(ens_weights)
      names(ew) <- names(ens_weights)
    } else {
      valid_single <- comp_dt[model %in% ens_names & !is.na(mase_combined)]
      inv_v <- 1 / valid_single$mase_combined
      ew <- inv_v / sum(inv_v)
      names(ew) <- valid_single$model
    }
    cat("  Ensemble:", paste(ens_names, collapse = ", "),
        "\n  Weights:", paste(round(ew, 3), collapse = ", "), "\n")

    ens_fc_list <- list(); ens_fc_results <- list(); ens_fitted_list <- list()
    for (enm in ens_names) {
      cat("  Running", enm, "...\n")
      fc_res_e <- dispatch_forecast(enm, mx_smooth, pop_arr, years, ages,
                                    groups, config$forecast_h,
                                    calib_coda, calib_direct,
                                    width = config$age_width,
                                    calib_all = calib_all,
                                    return_fitted = TRUE)
      ens_fc_list[[enm]]    <- fc_res_e$mx_forecast
      ens_fc_results[[enm]] <- fc_res_e
      if (!is.null(fc_res_e$mx_fitted))
        ens_fitted_list[[enm]] <- fc_res_e$mx_fitted
    }

    # Weighted average point forecast
    n_a_fc <- dim(ens_fc_list[[1]])[1]
    h_fc <- config$forecast_h; n_g_fc <- n_g
    ens_point <- array(0, dim = c(n_a_fc, h_fc, n_g_fc))
    w_total <- 0
    for (ii in seq_along(ens_names)) {
      fa <- ens_fc_list[[ens_names[ii]]]
      if (!is.null(fa)) {
        fa[is.na(fa)] <- 0
        ens_point <- ens_point + ew[ii] * fa
        w_total   <- w_total + ew[ii]
      }
    }
    if (w_total > 0 && w_total < 1) ens_point <- ens_point / w_total

    # Ensemble fitted values
    fitted_ens <- NULL
    if (length(ens_fitted_list) > 0) {
      fitted_ens <- array(0, dim = c(n_a, n_t, n_g))
      wf_total <- 0
      for (ii in seq_along(ens_names)) {
        fa_fit <- ens_fitted_list[[ens_names[ii]]]
        if (!is.null(fa_fit) && all(dim(fa_fit) == c(n_a, n_t, n_g))) {
          fa_fit[is.na(fa_fit)] <- 0
          fitted_ens <- fitted_ens + ew[ii] * fa_fit
          wf_total   <- wf_total + ew[ii]
        }
      }
      if (wf_total > 0 && wf_total < 1) fitted_ens <- fitted_ens / wf_total
      if (wf_total == 0) fitted_ens <- NULL
    }

    t3_name <- intersect(ens_names, c("Tucker3-Direct", "Tucker3-CoDa"))[1]
    t3_fc <- if (!is.null(t3_name)) ens_fc_results[[t3_name]]$full_fc else NULL
    full_fc <- list(
      mx_forecast = ens_point, fitted_ens = fitted_ens,
      fit_pct     = if (!is.null(t3_fc)) t3_fc$fit_pct else NA,
      X_log_full  = if (!is.null(t3_fc)) t3_fc$X_log_full else NULL,
      X_clr_full  = if (!is.null(t3_fc)) t3_fc$X_clr_full else NULL,
      B           = if (!is.null(t3_fc)) t3_fc$B else NULL,
      arima_fits  = if (!is.null(t3_fc)) t3_fc$arima_fits else NULL,
      n_t         = if (!is.null(t3_fc)) t3_fc$n_t else n_t)
    prod_method <- "ensemble"
  } else {
    cat("  Running", best_model, "...\n")
    fc_res <- dispatch_forecast(best_model, mx_smooth, pop_arr, years, ages,
                                groups, config$forecast_h,
                                calib_coda, calib_direct,
                                width = config$age_width,
                                calib_all = calib_all,
                                return_fitted = TRUE)
    full_fc <- fc_res$full_fc
    full_fc$fitted_ens <- fc_res$mx_fitted
    prod_method <- fc_res$method
  }

  # All-causes forecast for constraint
  ac_fc <- NULL
  if (config$is_cause_specific) {
    ac_fc <- tryCatch(
      tucker3_direct_forecast(ac_smooth, years, config$forecast_h,
                              3, 3, 2, "ets", fast = FALSE),
      error = function(e) NULL)
  }

  # ---- 9. PREDICTION INTERVALS ----
  pi_method <- config$pi_method %||% "bootstrap"
  cat("\n=== Prediction intervals (", pi_method, ", n_sim =", config$n_sim, ") ===\n")

  if (pi_method == "gam_bayes" && !is.null(gam_fits) && use_ensemble) {
    bs <- bootstrap_gam_bayes(gam_fits, ages, years, pop_arr, groups,
                               config$forecast_h, config$n_sim,
                               ens_names, ew,
                               calib_coda, calib_direct, calib_all,
                               width = config$age_width)
  } else if (use_ensemble) {
    sim_arrays <- list()
    for (enm in ens_names) {
      cat("  Bootstrapping", enm, "...\n")
      sim_arrays[[enm]] <- dispatch_model_sims(
        enm, mx_smooth, pop_arr, years, ages, groups,
        config$forecast_h, config$n_sim, ens_fc_results[[enm]],
        width = config$age_width, calib_all = calib_all,
        mx_raw = mx_raw)
    }
    bs <- ensemble_bootstrap(ens_names, ew, sim_arrays,
                             n_a, config$forecast_h, n_g)
  } else {
    cat("  Bootstrapping", best_model, "...\n")
    sims_single <- dispatch_model_sims(
      best_model, mx_smooth, pop_arr, years, ages, groups,
      config$forecast_h, config$n_sim, fc_res,
      width = config$age_width, calib_all = calib_all,
      mx_raw = mx_raw)
    bs <- list(
      mx_lower  = apply(sims_single, c(1, 2, 3), stats::quantile,
                         0.025, na.rm = TRUE),
      mx_upper  = apply(sims_single, c(1, 2, 3), stats::quantile,
                         0.975, na.rm = TRUE),
      mx_median = apply(sims_single, c(1, 2, 3), stats::quantile,
                         0.50, na.rm = TRUE))
  }

  # PI re-center + cap
  log_fc  <- log(pmax(full_fc$mx_forecast, 1e-12))
  log_lo  <- log(pmax(bs$mx_lower, 1e-12))
  log_hi  <- log(pmax(bs$mx_upper, 1e-12))
  log_med <- log(pmax(bs$mx_median, 1e-12))
  shift   <- log_fc - log_med
  log_lo  <- log_lo + shift; log_hi <- log_hi + shift
  bs$mx_median <- full_fc$mx_forecast
  cap <- 1.5
  log_lo <- pmax(log_lo, log_fc - cap)
  log_hi <- pmin(log_hi, log_fc + cap)
  bs$mx_lower <- exp(log_lo); bs$mx_upper <- exp(log_hi)
  cat("  PI re-centered (max |shift| =", round(max(abs(shift)), 4),
      ") and capped at +/-", cap, "on log-scale.\n")

  # ---- 10. BUILD FORECAST TABLE ----
  cat("\n=== Building forecast tables ===\n")
  fc_list <- list()
  for (k in seq_len(n_g)) {
    for (j in seq_len(config$forecast_h)) {
      fc_list[[length(fc_list) + 1]] <- append_group_cols(
        data.table::data.table(
          cause = config$cause_label, year = forecast_years[j],
          age_group = config$age_groups, age_start = ages,
          mx_point  = full_fc$mx_forecast[, j, k],
          mx_median = bs$mx_median[, j, k],
          mx_lo95   = bs$mx_lower[, j, k],
          mx_hi95   = bs$mx_upper[, j, k]),
        groups, k)
    }
  }
  fc_dt <- data.table::rbindlist(fc_list)

  # Cause-specific <= all-causes cap
  if (config$is_cause_specific && !is.null(ac_fc)) {
    ac_list <- list()
    for (k in seq_len(n_g)) for (j in seq_len(config$forecast_h)) {
      ac_list[[length(ac_list) + 1]] <- append_group_cols(
        data.table::data.table(
          year = forecast_years[j],
          age_group = config$age_groups, age_start = ages,
          mx_all_fc = ac_fc$mx_forecast[, j, k]),
        groups, k)
    }
    ac_dt <- data.table::rbindlist(ac_list)
    merge_keys <- c("year", names(groups), "age_group", "age_start")
    fc_dt <- merge(fc_dt, ac_dt, by = merge_keys, all.x = TRUE)
    fc_dt[!is.na(mx_all_fc) & mx_point > 0.95 * mx_all_fc,
          mx_point := 0.95 * mx_all_fc]
    fc_dt[, mx_all_fc := NULL]
  }

  # Observed table (raw data)
  obs_list <- list()
  for (k in seq_len(n_g)) for (j in seq_along(years)) {
    obs_list[[length(obs_list) + 1]] <- append_group_cols(
      data.table::data.table(
        cause = config$cause_label, year = years[j],
        age_group = config$age_groups, age_start = ages,
        mx_obs = mx_arr[, j, k]),
      groups, k)
  }
  obs_dt <- data.table::rbindlist(obs_list)

  # Smoothed table
  smooth_list <- list()
  for (k in seq_len(n_g)) for (j in seq_along(years)) {
    smooth_list[[length(smooth_list) + 1]] <- append_group_cols(
      data.table::data.table(
        cause = config$cause_label, year = years[j],
        age_group = config$age_groups, age_start = ages,
        mx_smooth = mx_smooth[, j, k]),
      groups, k)
  }
  smooth_dt <- data.table::rbindlist(smooth_list)

  # ---- 11. COHERENCE ----
  fc_dt <- enforce_coherence(fc_dt, obs_dt,
                             sex_constraint = config$sex_constraint)

  # ---- 12. FITTED VALUES ----
  fitted_list <- list()
  has_blended <- !is.null(full_fc$fitted_ens)
  has_decomp  <- !is.null(full_fc$X_log_full) || !is.null(full_fc$X_clr_full)
  for (k in seq_len(n_g)) for (j in seq_along(years)) {
    if (has_blended) {
      mx_fit <- full_fc$fitted_ens[, j, k]
    } else if (has_decomp && prod_method == "coda" &&
               !is.null(full_fc$X_clr_full)) {
      dx_fit <- inv_clr_to_dx(full_fc$X_clr_full[, j, k])
      mx_fit <- dx_to_mx(dx_fit, ages, config$age_width)
    } else if (has_decomp && !is.null(full_fc$X_log_full)) {
      mx_fit <- exp(full_fc$X_log_full[, j, k])
    } else {
      mx_fit <- mx_raw[, j, k]
    }
    fitted_list[[length(fitted_list) + 1]] <- append_group_cols(
      data.table::data.table(
        cause = config$cause_label, year = years[j],
        age_group = config$age_groups, age_start = ages,
        mx_fitted = mx_fit),
      groups, k)
  }
  fitted_dt <- data.table::rbindlist(fitted_list)

  # ---- 13. PLOTS ----
  cat("\n=== Generating plots ===\n")
  grp_names <- names(groups)
  plot_forecast_trends(obs_dt, fitted_dt, fc_dt, config$output_dir,
                       config$cause_label, config$cause_short, forecast_years,
                       group_cols = grp_names)
  plot_obs_vs_smooth(obs_dt, smooth_dt, config$output_dir,
                     config$cause_label, config$cause_short,
                     group_cols = grp_names)
  plot_forecast_trends_by_age(obs_dt, fitted_dt, fc_dt, config$output_dir,
                              config$cause_label, config$cause_short,
                              forecast_years, group_cols = grp_names)
  plot_forecast_age(fc_dt, config$output_dir, config$cause_label,
                    config$cause_short, config$age_groups,
                    group_cols = grp_names)
  if (!is.null(full_fc$B) && !is.null(full_fc$arima_fits)) {
    plot_B_components(full_fc, forecast_years, years, config$output_dir,
                      config$cause_label, config$cause_short)
  }

  # ---- 14. DISAGGREGATE TO SINGLE YEARS ----
  disagg_method <- config$disagg_method %||% "natural_gompertz"
  cat("\n=== Disaggregating to single year of age (", disagg_method, ") ===\n")
  target_ages <- config$target_ages

  # Select disaggregation function
  disagg_fn <- switch(disagg_method,
    hyman_kannisto   = disaggregate_hyman_kannisto,
    natural_gompertz = disaggregate_natural_gompertz,
    disaggregate_natural_gompertz  # default fallback
  )

  if (disagg_method == "gam") {
    # GAM-based 2D disaggregation with parametric tail blending
    mx_fc_arr <- full_fc$mx_forecast
    disagg_res <- disaggregate_forecast_gam(mx_smooth, ages, years,
                                             mx_fc_arr, forecast_years, config)
    # Build output table from GAM disaggregation arrays
    disagg_list <- list()
    for (k in seq_len(n_g)) {
      for (j in seq_len(config$forecast_h)) {
        yr <- forecast_years[j]
        mx_1yr <- disagg_res$mx_fc_1yr[, j, k]
        # PI via ratio interpolation on grouped PI
        sub_mask <- group_filter(fc_dt, groups, k) &
          fc_dt$cause == config$cause_label & fc_dt$year == yr
        sub <- fc_dt[sub_mask][order(age_start)]
        pi_res <- disaggregate_pi_ratios(mx_1yr, sub$mx_adjusted,
                                          sub$mx_lo95_adj, sub$mx_hi95_adj,
                                          sub$mx_median_adj,
                                          target_ages, ages, config$age_width)
        row <- data.table::data.table(
          cause = config$cause_label, year = yr, age = target_ages,
          mx_1yr = mx_1yr, mx_1yr_median = pi_res$mx_1yr_median,
          mx_1yr_lo95 = pi_res$mx_1yr_lo95, mx_1yr_hi95 = pi_res$mx_1yr_hi95)
        disagg_list[[length(disagg_list) + 1]] <- append_group_cols(row, groups, k)
      }
    }
  } else {
    # Compare methods on one example (use first group)
    ex_mask <- group_filter(fc_dt, groups, 1) &
      fc_dt$year == forecast_years[1]
    ex_mx <- fc_dt[ex_mask][order(age_start)]$mx_adjusted
    check_disaggregation(ex_mx, disagg_fn(ex_mx, ages, target_ages,
                                           config$age_width),
                         ages, target_ages)
    plot_disagg_comparison(ex_mx, ages, target_ages, config$output_dir,
                           config$cause_short, config$age_width)

    # Full disaggregation with ratio-based PIs
    disagg_list <- list()
    for (k in seq_len(n_g)) {
      for (j in seq_len(config$forecast_h)) {
        yr <- forecast_years[j]
        sub_mask <- group_filter(fc_dt, groups, k) &
          fc_dt$cause == config$cause_label & fc_dt$year == yr
        sub <- fc_dt[sub_mask][order(age_start)]
        mx_1yr <- disagg_fn(sub$mx_adjusted, ages, target_ages,
                             config$age_width)
        pi_res <- disaggregate_pi_ratios(mx_1yr, sub$mx_adjusted,
                                          sub$mx_lo95_adj, sub$mx_hi95_adj,
                                          sub$mx_median_adj,
                                          target_ages, ages, config$age_width)
        row <- data.table::data.table(
          cause = config$cause_label, year = yr, age = target_ages,
          mx_1yr = mx_1yr, mx_1yr_median = pi_res$mx_1yr_median,
          mx_1yr_lo95 = pi_res$mx_1yr_lo95, mx_1yr_hi95 = pi_res$mx_1yr_hi95)
        disagg_list[[length(disagg_list) + 1]] <- append_group_cols(row, groups, k)
      }
    }
  }
  mx_1yr_dt <- data.table::rbindlist(disagg_list)
  cat("  Single-year rows:", nrow(mx_1yr_dt), "\n")

  plot_disaggregation(mx_1yr_dt, config$output_dir, config$cause_label,
                      config$cause_short, group_cols = grp_names)

  # ---- 15. DIAGNOSTICS ----
  cat("\n=== Final diagnostics ===\n")
  for (yr in forecast_years[c(1, config$forecast_h)]) {
    for (k in seq_len(min(n_g, 2L))) {
      sub_mask <- group_filter(mx_1yr_dt, groups, k) &
        mx_1yr_dt$year == yr
      sub <- mx_1yr_dt[sub_mask]
      if (nrow(sub) < 3) next
      d2 <- diff(log(sub$mx_1yr), differences = 2)
      max_j <- max(abs(d2))
      glbl <- paste(vapply(names(groups),
                           function(g) paste0(g, "=", groups[[g]][k]),
                           character(1)), collapse = ", ")
      cat("  Age smoothness (", yr, ", ", glbl, "): max|d2| =",
          round(max_j, 3), ifelse(max_j < 0.5, "OK", "WARN"), "\n")
    }
  }

  # ---- 16. SAVE OUTPUT ----
  cat("\n=== Saving output ===\n")
  data.table::fwrite(fc_dt, file.path(config$output_dir,
                     paste0("mx_forecast_grouped_", cause_short, ".csv")))
  data.table::fwrite(mx_1yr_dt, file.path(config$output_dir,
                     paste0("mx_forecast_single_year_", cause_short, ".csv")))

  # ---- 17. SUMMARY ----
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("SUMMARY\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat("  Cause:        ", config$cause_label, "\n")
  cat("  Best model:   ", best_model, "\n")
  cat("  Forecast:     ", min(forecast_years), "-",
      max(forecast_years), "\n")
  cat("  Bootstrap:    ", config$n_sim, "sims\n")
  cat("  Grouped rows: ", nrow(fc_dt), "\n")
  cat("  Single-yr:    ", nrow(mx_1yr_dt), "\n")

  if (use_ensemble) {
    cat("  Ensemble:      ", paste(ens_names, collapse = ", "), "\n")
    cat("  Weights:       ", paste(round(ew, 3), collapse = ", "), "\n")
  }
  if (!is.null(calib_all)) {
    cat("\n  Calibrated model params:\n")
    for (cn in names(calib_all))
      cat("    ", cn, ":", paste(names(calib_all[[cn]]), calib_all[[cn]],
                                 sep = "=", collapse = ", "), "\n")
  }
  cat(paste(rep("=", 60), collapse = ""), "\n")

  ensemble_info <- if (use_ensemble) {
    list(models = ens_names, weights = ew,
         mase = comp_dt[model == "Ensemble", mase])
  } else NULL

  # ---- 18. BUILD mort_forecast RETURN OBJECT ----
  # Build a 3D fitted array from fitted_ens or decomposition
  fitted_arr <- full_fc$fitted_ens
  if (is.null(fitted_arr) && !is.null(full_fc$X_log_full)) {
    fitted_arr <- exp(full_fc$X_log_full[, seq_len(n_t), , drop = FALSE])
  }

  best_mase <- NA_real_
  bm_row <- comp_dt[comp_dt$model == best_model & !is.na(comp_dt$mase_combined), ]
  if (nrow(bm_row) > 0) best_mase <- bm_row$mase_combined[1]

  train_md <- mort_data(mx = mx_smooth, ages = ages, years = years,
                         groups = groups, pop = pop_arr,
                         cause = config$cause_label, smoothed = TRUE)

  result <- mort_forecast(
    mx_point       = full_fc$mx_forecast,
    train          = train_md,
    forecast_years = forecast_years,
    model          = best_model,
    mx_fitted      = fitted_arr,
    mx_raw         = mx_raw,
    mx_lower       = bs$mx_lower,
    mx_upper       = bs$mx_upper,
    mx_median      = bs$mx_median,
    fit_info       = list(
      mase             = best_mase,
      fit_pct          = full_fc$fit_pct,
      model_comparison = comp_dt
    ),
    extras = list(
      # Pipeline-specific data.tables
      fc_dt = fc_dt, mx_1yr_dt = mx_1yr_dt,
      obs_dt = obs_dt, fitted_dt = fitted_dt, smooth_dt = smooth_dt,
      # Calibration & decomposition
      full_fc = full_fc, method = prod_method,
      calib_coda = calib_coda, calib_direct = calib_direct,
      calib_all = calib_all, config = config,
      ensemble_info = ensemble_info,
      # Raw data for bootstrap / diagnostics
      mx_smooth = mx_smooth, pop_arr = pop_arr, gam_fits = gam_fits,
      smooth_check = smooth_check, lambda = lambda_used,
      # Tucker3 components (for plot type = "components")
      A = full_fc$A, B = full_fc$B, C = full_fc$C,
      G_mat = full_fc$G_mat, arima_fits = full_fc$arima_fits,
      n_t = full_fc$n_t
    )
  )
  invisible(result)
}

# ---- Reconciled multi-disease pipeline -----------------------------------

#' Run the reconciled multi-disease forecast pipeline
#'
#' Runs independent forecasts for each disease and all-causes, then
#' reconciles via variance-weighted (WLS) adjustment so that disease
#' forecasts sum to the all-cause forecast.
#'
#' @param disease_configs Named list of `mort_config` objects (one per cause).
#' @param all_cause_config `mort_config` for all-causes.
#' @param shared Named list of shared parameters (output_dir, forecast_h,
#'   n_sim, age_width, seed, etc.).
#' @return Invisibly, a list with all reconciled results.
#' @export
run_reconciled_forecast <- function(disease_configs, all_cause_config, shared) {
  set.seed(shared$seed)
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("  RECONCILED MULTI-DISEASE FORECAST\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  disease_order <- names(disease_configs)
  K <- length(disease_order)

  # Stage 1: All-cause forecast
  cat("=== STAGE 1: All-cause level forecast ===\n")
  ac_result <- run_forecast(all_cause_config)
  if (is.null(ac_result)) stop("All-cause forecast failed!")

  years  <- ac_result$train_years; ages <- ac_result$ages; groups <- ac_result$groups
  n_a <- length(ages); n_t <- length(years); n_g <- nrow(groups)
  h <- shared$forecast_h
  forecast_years <- (max(years) + 1):(max(years) + h)
  width <- shared$age_width
  mx_all_fc <- ac_result$mx_point

  # Stage 2: Independent per-disease forecasts
  cat("\n=== STAGE 2: Independent per-disease forecasts ===\n")
  disease_results <- list()
  for (nm in disease_order) {
    cfg <- disease_configs[[nm]]
    cat("\n*** Disease:", cfg$cause_label, "***\n")
    disease_results[[nm]] <- tryCatch(run_forecast(cfg),
                                      error = function(e) {
      cat("  !!! FAILED:", conditionMessage(e), "\n"); NULL
    })
  }
  failed <- vapply(disease_results, is.null, logical(1))
  if (any(failed)) {
    disease_order <- disease_order[!failed]; K <- length(disease_order)
  }

  # Stage 3: Collect forecasts & bootstrap
  cat("\n=== STAGE 3: Collecting independent forecasts ===\n")
  mx_indep <- list(); bs_indep <- list()
  for (nm in disease_order) {
    mx_fc_nm <- disease_results[[nm]]$mx_point
    bs_nm    <- list(mx_lower  = disease_results[[nm]]$mx_lower,
                     mx_upper  = disease_results[[nm]]$mx_upper,
                     mx_median = disease_results[[nm]]$mx_median)
    if (is.null(mx_fc_nm) || length(dim(mx_fc_nm)) != 3 ||
        !all(dim(mx_fc_nm) == c(n_a, h, n_g)))
      mx_fc_nm <- array(NA_real_, dim = c(n_a, h, n_g))
    mx_fc_nm[!is.finite(mx_fc_nm)] <- 1e-12
    mx_indep[[nm]] <- mx_fc_nm
    if (is.null(bs_nm) || is.null(bs_nm$mx_lower)) {
      bs_indep[[nm]] <- list(mx_lower = mx_fc_nm * 0.8,
                             mx_upper = mx_fc_nm * 1.2,
                             mx_median = mx_fc_nm)
    } else {
      for (bc in c("mx_lower", "mx_upper", "mx_median"))
        if (!is.null(bs_nm[[bc]]))
          bs_nm[[bc]][!is.finite(bs_nm[[bc]])] <- mx_fc_nm[!is.finite(bs_nm[[bc]])]
      bs_indep[[nm]] <- bs_nm
    }
  }

  # Stage 4: WLS reconciliation
  cat("\n=== STAGE 4: WLS reconciliation ===\n")
  recon <- wls_reconcile(mx_all_fc, mx_indep, bs_indep, disease_order)
  mx_recon <- recon$mx_recon; bs_recon <- recon$bs_recon; delta <- recon$delta

  # Stage 5: Per-disease output
  cat("\n=== STAGE 5: Per-disease reconciled output ===\n")
  fc_dt_combined <- list()
  all_reconciled_results <- list()
  target_ages <- shared$target_ages

  for (nm in disease_order) {
    cfg <- disease_configs[[nm]]; res <- disease_results[[nm]]
    cat("\n  ---", cfg$cause_label, "(Reconciled) ---\n")

    fc_list <- list()
    for (k in seq_len(n_g)) for (j in seq_len(h)) {
      fc_list[[length(fc_list) + 1]] <- append_group_cols(
        data.table::data.table(
          cause = cfg$cause_label, year = forecast_years[j],
          age_group = shared$age_groups, age_start = ages,
          mx_point  = mx_recon[[nm]][, j, k],
          mx_median = bs_recon[[nm]]$mx_median[, j, k],
          mx_lo95   = bs_recon[[nm]]$mx_lower[, j, k],
          mx_hi95   = bs_recon[[nm]]$mx_upper[, j, k]),
        groups, k)
    }
    fc_dt <- data.table::rbindlist(fc_list)
    res_obs_dt    <- res$extras$obs_dt
    res_fitted_dt <- res$extras$fitted_dt
    res_mx_smooth <- res$extras$mx_smooth

    fc_dt <- enforce_coherence(fc_dt, res_obs_dt,
                               sex_constraint = shared$sex_constraint)

    rsuf <- paste0(nm, "_reconciled")
    plot_forecast_trends(res_obs_dt, res_fitted_dt, fc_dt, shared$output_dir,
                         paste0(cfg$cause_label, " (Reconciled)"), rsuf,
                         forecast_years, group_cols = names(groups))

    # Build smooth_dt
    smooth_list_r <- list()
    for (k in seq_len(n_g)) for (j in seq_along(years)) {
      smooth_list_r[[length(smooth_list_r) + 1]] <- append_group_cols(
        data.table::data.table(
          cause = cfg$cause_label, year = years[j],
          age_group = shared$age_groups, age_start = ages,
          mx_smooth = res_mx_smooth[, j, k]),
        groups, k)
    }
    smooth_dt_r <- data.table::rbindlist(smooth_list_r)
    plot_obs_vs_smooth(res_obs_dt, smooth_dt_r, shared$output_dir,
                       paste0(cfg$cause_label, " (Reconciled)"), rsuf,
                       group_cols = names(groups))
    plot_forecast_trends_by_age(res_obs_dt, res_fitted_dt, fc_dt,
                                shared$output_dir,
                                paste0(cfg$cause_label, " (Reconciled)"), rsuf,
                                forecast_years, group_cols = names(groups))
    plot_forecast_age(fc_dt, shared$output_dir,
                      paste0(cfg$cause_label, " (Reconciled)"), rsuf,
                      shared$age_groups, group_cols = names(groups))

    # Disaggregate
    disagg_list <- list()
    for (k in seq_len(n_g)) {
      for (j in seq_len(h)) {
        yr <- forecast_years[j]
        sub_mask <- group_filter(fc_dt, groups, k) &
          fc_dt$cause == cfg$cause_label & fc_dt$year == yr
        sub <- fc_dt[sub_mask][order(age_start)]
        mx_1yr <- disaggregate_natural_gompertz(sub$mx_adjusted, ages,
                                                target_ages, width)
        grp_pt  <- pmax(sub$mx_adjusted, 1e-15)
        r_lo    <- sub$mx_lo95_adj   / grp_pt
        r_hi    <- sub$mx_hi95_adj   / grp_pt
        r_med   <- sub$mx_median_adj / grp_pt
        grp_mids <- c(ages[-length(ages)] + width / 2, ages[length(ages)] + 5)
        fn_lo  <- stats::splinefun(grp_mids, r_lo,  method = "natural")
        fn_hi  <- stats::splinefun(grp_mids, r_hi,  method = "natural")
        fn_med <- stats::splinefun(grp_mids, r_med, method = "natural")
        mx_1yr_lo  <- mx_1yr * pmin(fn_lo(target_ages + 0.5),  1)
        mx_1yr_hi  <- mx_1yr * pmax(fn_hi(target_ages + 0.5),  1)
        mx_1yr_med <- mx_1yr * fn_med(target_ages + 0.5)
        mx_1yr_med <- pmin(pmax(mx_1yr_med, mx_1yr_lo), mx_1yr_hi)
        row <- data.table::data.table(
          cause = cfg$cause_label, year = yr, age = target_ages,
          mx_1yr = mx_1yr, mx_1yr_median = mx_1yr_med,
          mx_1yr_lo95 = mx_1yr_lo, mx_1yr_hi95 = mx_1yr_hi)
        disagg_list[[length(disagg_list) + 1]] <- append_group_cols(row, groups, k)
      }
    }
    mx_1yr_dt <- data.table::rbindlist(disagg_list)

    data.table::fwrite(fc_dt, file.path(shared$output_dir,
                       paste0("mx_forecast_grouped_", rsuf, ".csv")))
    data.table::fwrite(mx_1yr_dt, file.path(shared$output_dir,
                       paste0("mx_forecast_single_year_", rsuf, ".csv")))

    grp_names <- names(groups)
    sel_cols <- c("year", grp_names, "age_group", "age_start", "cause")
    fc_dt_combined[[length(fc_dt_combined) + 1]] <- fc_dt[, c(
      sel_cols, "mx_adjusted"), with = FALSE]
    data.table::setnames(fc_dt_combined[[length(fc_dt_combined)]],
                         "mx_adjusted", "mx_point")
    all_reconciled_results[[nm]] <- list(
      config = cfg, fc_grouped = fc_dt, fc_single_yr = mx_1yr_dt,
      obs_dt = res_obs_dt, fitted_dt = res_fitted_dt)
  }
  fc_dt_combined <- data.table::rbindlist(fc_dt_combined)

  # Stage 6: Consistency check
  cat("\n=== STAGE 6: Consistency check ===\n")
  ac_fc_dt_list <- list()
  for (k in seq_len(n_g)) for (j in seq_len(h)) {
    ac_fc_dt_list[[length(ac_fc_dt_list) + 1]] <- append_group_cols(
      data.table::data.table(
        year = forecast_years[j],
        age_group = shared$age_groups, age_start = ages,
        mx_point = mx_all_fc[, j, k]),
      groups, k)
  }
  ac_fc_dt <- data.table::rbindlist(ac_fc_dt_list)

  grp_keys <- c("year", names(groups), "age_group")
  disease_sum <- fc_dt_combined[, .(mx_sum = sum(mx_point, na.rm = TRUE)),
                                 by = grp_keys]
  ac_rename <- ac_fc_dt[, c(grp_keys, "mx_point"), with = FALSE]
  data.table::setnames(ac_rename, "mx_point", "mx_all")
  ac_check <- merge(disease_sum, ac_rename, by = grp_keys)
  ac_check[, ratio := mx_sum / mx_all]
  cat("  Sum/All-cause: min=", round(min(ac_check$ratio), 6),
      " max=", round(max(ac_check$ratio), 6), "\n")

  plot_consistency_check(fc_dt_combined, ac_fc_dt, shared$output_dir,
                         group_cols = names(groups))

  data.table::fwrite(fc_dt_combined,
                     file.path(shared$output_dir,
                               "mx_forecast_combined_reconciled.csv"))

  invisible(list(
    ac_result = ac_result,
    disease_results = disease_results,
    reconciled_results = all_reconciled_results,
    delta = delta,
    fc_combined = fc_dt_combined,
    disease_order = disease_order
  ))
}

# ---- Coherent (Tucker4-ILR) multi-disease pipeline ----------------------

#' Run the coherent multi-disease forecast pipeline (Tucker4-ILR)
#'
#' Jointly models cause-of-death shares via a 4D Tucker decomposition
#' on ILR-transformed shares, ensuring compositions sum to 1 by
#' construction.
#'
#' @param disease_configs Named list of `mort_config` objects (one per cause).
#' @param all_cause_config `mort_config` for all-causes.
#' @param shared Named list of shared parameters.
#' @return Invisibly, a list with all coherent results.
#' @export
run_coherent_forecast <- function(disease_configs, all_cause_config, shared) {
  set.seed(shared$seed)
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("  COHERENT MULTI-DISEASE FORECAST (Tucker4-ILR)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")

  disease_order <- names(disease_configs)
  K <- length(disease_order)

  # Stage 1: All-cause forecast
  cat("=== STAGE 1: All-cause level forecast ===\n")
  ac_result <- run_forecast(all_cause_config)
  if (is.null(ac_result)) stop("All-cause forecast failed!")

  years  <- ac_result$train_years; ages <- ac_result$ages; groups <- ac_result$groups
  n_a <- length(ages); n_t <- length(years); n_g <- nrow(groups)
  forecast_years <- (max(years) + 1):(max(years) + shared$forecast_h)
  width <- shared$age_width

  mx_all_smooth <- ac_result$extras$mx_smooth
  mx_all_raw    <- ac_result$mx_raw

  ac_fc_obj <- ac_result$extras$full_fc
  if (!is.null(ac_fc_obj$X_log_full)) {
    mx_all_fc <- array(NA_real_, dim = c(n_a, shared$forecast_h, n_g))
    for (j in seq_len(shared$forecast_h))
      for (k in seq_len(n_g))
        mx_all_fc[, j, k] <- exp(ac_fc_obj$X_log_full[, n_t + j, k])
  } else {
    mx_all_fc <- ac_fc_obj$mx_forecast
  }

  # Stage 2: Load all diseases, build 4D tensor
  cat("=== STAGE 2: Loading diseases & building 4D tensor ===\n")
  mx_disease_raw <- list(); mx_disease_smooth <- list()
  for (nm in disease_order) {
    cfg <- disease_configs[[nm]]
    if (inherits(cfg, "mort_data")) {
      mx_d <- cfg$mx
      sm_cfg <- cfg
    } else {
      dt_d <- load_and_prepare_ons(cfg)
      mx_d <- build_array(dt_d, paste0("mx_", nm), ages, years, groups)
      sm_cfg <- cfg
    }
    sm_d <- smooth_mortality(mx_d, ages, years, sm_cfg)
    mx_disease_raw[[nm]]    <- mx_d
    mx_disease_smooth[[nm]] <- sm_d$arr
  }

  shares_obs <- compute_cause_shares(mx_disease_smooth, mx_all_smooth,
                                     disease_order)
  X_ilr <- ilr_cause_transform(shares_obs)
  cat("  ILR tensor:", paste(dim(X_ilr), collapse = " x "), "\n")

  # Stage 3: Calibrate Tucker4
  cat("\n=== STAGE 3: Calibrating Tucker4 ===\n")
  cv_years_idx <- which(years %in% shared$cv_years)
  ilr_cv <- X_ilr[, cv_years_idx, , , drop = FALSE]
  mx_all_cv <- mx_all_smooth[, cv_years_idx, , drop = FALSE]
  mx_disease_cv <- lapply(mx_disease_smooth,
                          function(arr) arr[, cv_years_idx, , drop = FALSE])
  cv_yrs <- years[cv_years_idx]

  calib <- calibrate_tucker4(ilr_cv, mx_all_cv, mx_disease_cv, disease_order,
                             cv_yrs, h_val = shared$cv_h,
                             min_train = shared$cv_min_train)

  # Stage 4: Holdout comparison
  cat("\n=== STAGE 4: Holdout comparison ===\n")
  train_idx <- which(years %in% shared$holdout_train)
  test_idx  <- which(years %in% shared$holdout_test)
  holdout_comp <- compare_coherent_vs_independent(
    X_ilr, mx_all_smooth, mx_disease_smooth, mx_disease_smooth,
    ac_result$extras$pop_arr, disease_order, disease_configs,
    years, ages, groups, train_idx, test_idx, calib, width = width)
  print(holdout_comp)
  data.table::fwrite(holdout_comp,
                     file.path(shared$output_dir,
                               "validation_coherent_vs_independent.csv"))

  # Stage 5: Production Tucker4 forecast
  cat("\n=== STAGE 5: Production Tucker4 forecast ===\n")
  full_fc <- tucker4_coherent_forecast(
    X_ilr, years, shared$forecast_h,
    calib$R_A, calib$R_B, calib$R_C, calib$R_D, calib$ts_method,
    fast = FALSE, return_full = TRUE)
  cat("  Tucker4 fit:", round(full_fc$fit_pct, 1), "%\n")

  shares_insample <- array(NA_real_, dim = c(n_a, n_t, n_g, K + 1))
  for (i in seq_len(n_a)) for (j in seq_len(n_t)) for (k in seq_len(n_g))
    shares_insample[i, j, k, ] <- ilr_cause_inverse(
      full_fc$X_ilr_full[i, j, k, ])

  # Stage 6: Bootstrap PIs
  cat("\n=== STAGE 6: Bootstrap PIs ===\n")
  bs <- bootstrap_coherent_forecast(full_fc, mx_all_fc, shared$n_sim,
                                    shared$forecast_h)

  # Stage 7: Per-disease output
  cat("\n=== STAGE 7: Per-disease output ===\n")
  ac_fc_dt_list <- list()
  for (k in seq_len(n_g)) for (j in seq_len(shared$forecast_h)) {
    ac_fc_dt_list[[length(ac_fc_dt_list) + 1]] <- append_group_cols(
      data.table::data.table(
        year = forecast_years[j],
        age_group = shared$age_groups, age_start = ages,
        mx_point = mx_all_fc[, j, k]),
      groups, k)
  }
  ac_fc_dt <- data.table::rbindlist(ac_fc_dt_list)

  fc_dt_combined <- list(); all_disease_results <- list()

  for (d in seq_len(K)) {
    nm <- disease_order[d]; cfg <- disease_configs[[nm]]
    cat("\n  ---", cfg$cause_label, "(Coherent) ---\n")

    fc_list <- list()
    for (k in seq_len(n_g)) for (j in seq_len(shared$forecast_h)) {
      mx_pt  <- full_fc$shares_forecast[, j, k, d] * mx_all_fc[, j, k]
      fc_list[[length(fc_list) + 1]] <- append_group_cols(
        data.table::data.table(
          cause = cfg$cause_label, year = forecast_years[j],
          age_group = shared$age_groups, age_start = ages,
          mx_point = mx_pt, mx_median = bs$mx_median[, j, k, d],
          mx_lo95 = bs$mx_lower[, j, k, d],
          mx_hi95 = bs$mx_upper[, j, k, d]),
        groups, k)
    }
    fc_dt <- data.table::rbindlist(fc_list)

    # Observed & smoothed tables
    obs_list <- list(); smooth_list <- list(); fitted_list <- list()
    for (k in seq_len(n_g)) for (j in seq_along(years)) {
      obs_list[[length(obs_list) + 1]] <- append_group_cols(
        data.table::data.table(
          cause = cfg$cause_label, year = years[j],
          age_group = shared$age_groups, age_start = ages,
          mx_obs = mx_disease_raw[[nm]][, j, k]),
        groups, k)
      smooth_list[[length(smooth_list) + 1]] <- append_group_cols(
        data.table::data.table(
          cause = cfg$cause_label, year = years[j],
          age_group = shared$age_groups, age_start = ages,
          mx_smooth = mx_disease_smooth[[nm]][, j, k]),
        groups, k)
      fitted_list[[length(fitted_list) + 1]] <- append_group_cols(
        data.table::data.table(
          cause = cfg$cause_label, year = years[j],
          age_group = shared$age_groups, age_start = ages,
          mx_fitted = shares_insample[, j, k, d] * mx_all_raw[, j, k]),
        groups, k)
    }
    obs_dt    <- data.table::rbindlist(obs_list)
    smooth_dt <- data.table::rbindlist(smooth_list)
    fitted_dt <- data.table::rbindlist(fitted_list)

    fc_dt <- enforce_coherence(fc_dt, obs_dt,
                               sex_constraint = shared$sex_constraint)

    csuf <- paste0(nm, "_coherent")
    plot_forecast_trends(obs_dt, fitted_dt, fc_dt, shared$output_dir,
                         paste0(cfg$cause_label, " (Coherent)"), csuf,
                         forecast_years, group_cols = names(groups))
    plot_obs_vs_smooth(obs_dt, smooth_dt, shared$output_dir,
                       paste0(cfg$cause_label, " (Coherent)"), csuf,
                       group_cols = names(groups))
    plot_forecast_trends_by_age(obs_dt, fitted_dt, fc_dt, shared$output_dir,
                                paste0(cfg$cause_label, " (Coherent)"), csuf,
                                forecast_years, group_cols = names(groups))
    plot_forecast_age(fc_dt, shared$output_dir,
                      paste0(cfg$cause_label, " (Coherent)"), csuf,
                      shared$age_groups, group_cols = names(groups))

    # Disaggregate
    disagg_list <- list()
    target_ages <- shared$target_ages
    for (k in seq_len(n_g)) {
      for (j in seq_len(shared$forecast_h)) {
        yr <- forecast_years[j]
        sub <- fc_dt[cause == cfg$cause_label & year == yr &
                       group_filter(.SD, groups, k)][order(age_start)]
        mx_1yr <- disaggregate_natural_gompertz(sub$mx_adjusted, ages,
                                                target_ages, width)
        grp_pt  <- pmax(sub$mx_adjusted, 1e-15)
        r_lo    <- sub$mx_lo95_adj   / grp_pt
        r_hi    <- sub$mx_hi95_adj   / grp_pt
        r_med   <- sub$mx_median_adj / grp_pt
        grp_mids <- c(ages[-length(ages)] + width / 2, ages[length(ages)] + 5)
        fn_lo  <- stats::splinefun(grp_mids, r_lo,  method = "natural")
        fn_hi  <- stats::splinefun(grp_mids, r_hi,  method = "natural")
        fn_med <- stats::splinefun(grp_mids, r_med, method = "natural")
        mx_1yr_lo  <- mx_1yr * pmin(fn_lo(target_ages + 0.5),  1)
        mx_1yr_hi  <- mx_1yr * pmax(fn_hi(target_ages + 0.5),  1)
        mx_1yr_med <- mx_1yr * fn_med(target_ages + 0.5)
        mx_1yr_med <- pmin(pmax(mx_1yr_med, mx_1yr_lo), mx_1yr_hi)
        disagg_list[[length(disagg_list) + 1]] <- append_group_cols(
          data.table::data.table(
            cause = cfg$cause_label, year = yr, age = target_ages,
            mx_1yr = mx_1yr, mx_1yr_median = mx_1yr_med,
            mx_1yr_lo95 = mx_1yr_lo, mx_1yr_hi95 = mx_1yr_hi),
          groups, k)
      }
    }
    mx_1yr_dt <- data.table::rbindlist(disagg_list)

    plot_disaggregation(mx_1yr_dt, shared$output_dir,
                        paste0(cfg$cause_label, " (Coherent)"), csuf,
                        group_cols = names(groups))

    data.table::fwrite(fc_dt, file.path(shared$output_dir,
                       paste0("mx_forecast_grouped_", csuf, ".csv")))
    data.table::fwrite(mx_1yr_dt, file.path(shared$output_dir,
                       paste0("mx_forecast_single_year_", csuf, ".csv")))

    grp_names <- names(groups)
    fc_dt_combined[[length(fc_dt_combined) + 1]] <- fc_dt[,
      c("year", grp_names, "age_group", "cause", "mx_adjusted"),
      with = FALSE]
    data.table::setnames(
      fc_dt_combined[[length(fc_dt_combined)]], "mx_adjusted", "mx_point")
    all_disease_results[[nm]] <- list(
      config = cfg, fc_grouped = fc_dt, fc_single_yr = mx_1yr_dt,
      obs_dt = obs_dt, fitted_dt = fitted_dt)
  }
  fc_dt_combined <- data.table::rbindlist(fc_dt_combined)

  # Stage 8: Remainder & consistency
  cat("\n=== STAGE 8: Consistency check ===\n")
  rem_list <- list()
  for (k in seq_len(n_g)) for (j in seq_len(shared$forecast_h)) {
    rem_list[[length(rem_list) + 1]] <- append_group_cols(
      data.table::data.table(
        year = forecast_years[j],
        age_group = shared$age_groups, cause = "Remainder",
        mx_point = full_fc$shares_forecast[, j, k, K + 1] *
          mx_all_fc[, j, k]),
      groups, k)
  }
  fc_dt_combined <- data.table::rbindlist(
    list(fc_dt_combined, data.table::rbindlist(rem_list)), fill = TRUE)

  grp_keys <- c("year", names(groups), "age_group")
  disease_sum <- fc_dt_combined[, .(mx_sum = sum(mx_point, na.rm = TRUE)),
                                 by = grp_keys]
  ac_rename <- ac_fc_dt[, c(grp_keys, "mx_point"), with = FALSE]
  data.table::setnames(ac_rename, "mx_point", "mx_all")
  ac_check <- merge(disease_sum, ac_rename, by = grp_keys)
  ac_check[, ratio := mx_sum / mx_all]
  cat("  Sum/All-cause: min=", round(min(ac_check$ratio), 6),
      " max=", round(max(ac_check$ratio), 6), "\n")

  plot_consistency_check(fc_dt_combined, ac_fc_dt, shared$output_dir,
                         group_cols = names(groups))
  plot_share_evolution(shares_obs, full_fc$shares_forecast,
                       years, forecast_years, disease_order,
                       shared$output_dir, group_index = 1, groups = groups)

  invisible(list(
    ac_result = ac_result,
    disease_results = all_disease_results,
    holdout_comp = holdout_comp,
    calib_tucker4 = calib,
    full_fc_tucker4 = full_fc,
    shares_obs = shares_obs,
    X_ilr = X_ilr,
    fc_combined = fc_dt_combined,
    disease_order = disease_order
  ))
}
