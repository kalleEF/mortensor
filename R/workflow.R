#' Step-by-step mortality forecasting workflow
#'
#' These functions provide a modular, interactive workflow for mortality
#' forecasting as an alternative to the monolithic [run_forecast()]:
#'
#' 1. [read_mortality()] or [read_mortality_from_arrays()] -- load data
#' 2. [smooth_mortality()] -- smooth the data (already exists)
#' 3. [calibrate()] -- calibrate and compare models
#' 4. [forecast.mort_calibration()] -- produce forecasts
#' 5. [add_pi()] -- add prediction intervals
#' 6. [disaggregate()] -- disaggregate to single-year ages
#'
#' @name workflow
NULL

#' @export
forecast.mort_data <- function(object, h = NULL, ...) {
  stop("Cannot forecast a mort_data object directly.\n",
       "  Run calibrate() first, then forecast the result:\n",
       "    cal <- calibrate(sm)\n",
       "    fc  <- forecast(cal, h = 10)",
       call. = FALSE)
}

# ---- S3 class: mort_calibration ----

#' Create a mortality calibration object
#'
#' @param data A `mort_data` object (smoothed training data).
#' @param data_raw A `mort_data` object (unsmoothed raw data), or `NULL`.
#' @param models Character vector of model names evaluated.
#' @param comparison data.table of model comparison results.
#' @param best_model Character: name of the best model.
#' @param ensemble_names Character vector of ensemble component names.
#' @param ensemble_weights Named numeric vector of ensemble weights.
#' @param calib_coda Tucker3-CoDa calibration result.
#' @param calib_direct Tucker3-Direct calibration result.
#' @param calib_all List of calibration results for other models.
#' @param config A `mort_config` object (optional pipeline settings).
#' @return A `mort_calibration` object.
#' @export
mort_calibration <- function(data, data_raw = NULL, models, comparison,
                              best_model, ensemble_names = character(0),
                              ensemble_weights = NULL,
                              calib_coda = NULL, calib_direct = NULL,
                              calib_all = NULL, config = NULL) {
  stopifnot(inherits(data, "mort_data"))
  structure(
    list(
      data              = data,
      data_raw          = data_raw,
      models            = models,
      comparison        = comparison,
      best_model        = best_model,
      ensemble_names    = ensemble_names,
      ensemble_weights  = ensemble_weights,
      calib_coda        = calib_coda,
      calib_direct      = calib_direct,
      calib_all         = calib_all,
      config            = config
    ),
    class = "mort_calibration"
  )
}

#' @export
print.mort_calibration <- function(x, ...) {
  cat("<mort_calibration>\n")
  cat("  Cause:     ", x$data$cause %||% "(not set)", "\n")
  cat("  Data:      ", x$data$n_a, "ages x", x$data$n_t, "years x",
      x$data$n_g, "groups\n")
  cat("  Models:    ", paste(x$models, collapse = ", "), "\n")
  cat("  Best:      ", x$best_model, "\n")
  if (length(x$ensemble_names) > 0) {
    cat("  Ensemble:  ", paste(x$ensemble_names, collapse = ", "), "\n")
    cat("  Weights:   ", paste(round(x$ensemble_weights, 3), collapse = ", "), "\n")
  }
  if (!is.null(x$comparison))
    cat("  Use summary() for full model comparison.\n")
  invisible(x)
}

#' @export
summary.mort_calibration <- function(object, ...) {
  cat("<mort_calibration> summary\n\n")
  cat("Best model:", object$best_model, "\n\n")
  if (!is.null(object$comparison)) {
    cat("Model comparison (sorted by combined MASE):\n")
    print(object$comparison[order(object$comparison$mase_combined), ])
  }
  if (!is.null(object$calib_coda)) {
    cat("\nTucker3-CoDa ranks: R_A=", object$calib_coda$R_A,
        " R_B=", object$calib_coda$R_B, " R_C=", object$calib_coda$R_C,
        " ts=", object$calib_coda$ts_method, "\n")
  }
  invisible(object)
}


# ---- calibrate ----

#' Calibrate and compare mortality forecasting models
#'
#' Runs cross-validated calibration for Tucker3 (CoDa and Direct) and
#' optionally other models (Lee-Carter, FDM, StMoMo, ARIMA, ETS, RWD),
#' then compares all on a holdout test set and builds an ensemble.
#'
#' @param data A `mort_data` object. If not yet smoothed, smoothing is
#'   applied automatically.
#' @param models Character vector of models to evaluate. Available:
#'   `"tucker3_coda"`, `"tucker3_direct"`, `"lee_carter"`, `"fdm"`,
#'   `"stmomo"`, `"arima"`, `"ets"`, `"rwd"`.
#' @param holdout_years Integer vector of years for holdout test. If `NULL`,
#'   the last `holdout_frac` of years is used.
#' @param holdout_frac Fraction of years for holdout (default 0.2). Ignored
#'   if `holdout_years` is provided.
#' @param ensemble Logical; build an inverse-MASE ensemble? (default TRUE)
#' @param ensemble_top_k Maximum number of models in the ensemble.
#' @param rank_selection Character; rank selection strategy for Tucker3 models.
#'   `"cv"` (default) uses rolling-origin cross-validation, `"bic"` or `"aic"`
#'   selects ranks via information criteria (faster). Passed to
#'   [calibrate_tucker3()].
#' @param cv_metric Character; scoring metric for cross-validation.
#'   `"mase"` (default) uses cell-level MASE, `"mase_combined"` uses
#'   0.7 * cell + 0.3 * trend MASE, `"mase_trend"` uses trend MASE only.
#'   Applied to all models during calibration.
#' @param config Optional `mort_config` for advanced settings (smoothing
#'   method, COVID years, etc.).
#' @param ... Additional arguments passed to calibration functions.
#' @return A `mort_calibration` object.
#' @export
calibrate <- function(data,
                      models = c("tucker3_coda", "tucker3_direct",
                                 "lee_carter", "fdm",
                                 "ets", "arima", "rwd"),
                      holdout_years = NULL, holdout_frac = 0.2,
                      ensemble = TRUE, ensemble_top_k = 3L,
                      rank_selection = c("cv", "bic", "aic"),
                      cv_metric = c("mase", "mase_combined", "mase_trend"),
                      config = NULL, ...) {
  stopifnot(inherits(data, "mort_data"))
  rank_selection <- match.arg(rank_selection)
  cv_metric <- match.arg(cv_metric)
  if (is.null(config)) config <- mort_config()

  data_raw <- data
  ages   <- data$ages
  years  <- data$years
  groups <- data$groups
  n_t    <- data$n_t

  # Auto-smooth if needed
  if (!isTRUE(data$smoothed)) {
    cat("=== Smoothing data ===\n")
    data <- smooth_mortality(data,
                             method = config$smoothing_method,
                             covid_years = config$covid_years)
  }

  mx_smooth <- data$mx
  pop_arr   <- data$pop
  if (is.null(pop_arr)) {
    pop_arr <- array(1, dim = dim(mx_smooth))
  }

  # Derive holdout split
  if (is.null(holdout_years)) {
    n_test <- max(2L, round(n_t * holdout_frac))
    test_idx  <- (n_t - n_test + 1):n_t
    train_idx <- seq_len(n_t - n_test)
  } else {
    test_idx  <- which(years %in% holdout_years)
    train_idx <- which(!years %in% holdout_years)
  }
  if (length(test_idx) == 0) stop("No test years in holdout set.")
  h_val <- length(test_idx)

  # Calibrate Tucker3
  cat("=== Calibrating Tucker3 ===\n")
  calib_coda <- calib_direct <- NULL
  if ("tucker3_coda" %in% models) {
    calib_coda <- calibrate_tucker3(mx_smooth, years, ages, method = "coda",
                                     width = config$age_width, h_val = h_val,
                                     rank_selection = rank_selection,
                                     cv_metric = cv_metric)
    cat("  CoDa: R_A=", calib_coda$R_A, " R_B=", calib_coda$R_B,
        " R_C=", calib_coda$R_C, " ts=", calib_coda$ts_method,
        " CV-score=", round(calib_coda$cv_score, 4), "\n")
  }
  if ("tucker3_direct" %in% models) {
    calib_direct <- calibrate_tucker3(mx_smooth, years, ages, method = "direct",
                                       width = config$age_width, h_val = h_val,
                                       rank_selection = rank_selection,
                                       cv_metric = cv_metric)
    cat("  Direct: R_A=", calib_direct$R_A, " R_B=", calib_direct$R_B,
        " R_C=", calib_direct$R_C, " ts=", calib_direct$ts_method,
        " CV-score=", round(calib_direct$cv_score, 4), "\n")
  }

  # Calibrate other models
  calib_all <- NULL
  non_tucker <- intersect(models, c("lee_carter", "fdm", "stmomo",
                                     "arima", "ets", "rwd"))
  if (length(non_tucker) > 0 && isTRUE(config$calibrate_models)) {
    cat("=== Calibrating other models ===\n")
    calib_all <- calibrate_all_models(mx_smooth, pop_arr, years, ages, groups,
                                       h_val = h_val,
                                       n_origins = config$cv_n_origins,
                                       width = config$age_width,
                                       cv_metric = cv_metric)
  }

  # Model comparison on holdout
  cat("=== Comparing models on holdout ===\n")
  if (is.null(calib_coda))
    calib_coda <- list(R_A = 3, R_B = 3, R_C = 2, ts_method = "ets")
  if (is.null(calib_direct))
    calib_direct <- list(R_A = 3, R_B = 3, R_C = 2, ts_method = "ets")

  comp_res <- compare_models(mx_smooth, pop_arr, years, ages, groups,
                              train_idx, test_idx,
                              calib_coda, calib_direct,
                              width = config$age_width,
                              mx_eval_ref = data_raw$mx,
                              calib_all = calib_all)
  comp_dt <- comp_res$results
  cat("\n")
  print(comp_dt[order(comp_dt$mase_combined), ])

  # Select best
  ens_names <- character(0)
  ens_weights <- NULL
  use_ensemble <- ensemble && !is.null(comp_res$ens_best_names) &&
    length(comp_res$ens_best_names) >= 2

  if (use_ensemble) {
    ens_comb <- comp_dt[comp_dt$model == "Ensemble", ]$mase_combined
    best_single <- comp_dt[comp_dt$model != "Ensemble" &
                             !is.na(comp_dt$mase_combined), ]
    best_single <- best_single[which.min(best_single$mase_combined), ]
    if (is.finite(ens_comb) && ens_comb <= best_single$mase_combined * 1.05) {
      best_model <- "Ensemble"
      ens_names  <- comp_res$ens_best_names
      ens_weights <- comp_res$ens_best_weights
    } else {
      use_ensemble <- FALSE
    }
  }
  if (!use_ensemble) {
    eligible <- comp_dt[comp_dt$model != "Ensemble" &
                          !is.na(comp_dt$mase_combined), ]
    best_model <- eligible$model[which.min(eligible$mase_combined)]
  }
  cat("\nBest model:", best_model, "\n")

  mort_calibration(
    data = data, data_raw = data_raw, models = models,
    comparison = comp_dt, best_model = best_model,
    ensemble_names = ens_names, ensemble_weights = ens_weights,
    calib_coda = calib_coda, calib_direct = calib_direct,
    calib_all = calib_all, config = config)
}


# ---- forecast.mort_calibration ----

#' Produce mortality forecasts from a calibrated model
#'
#' @param object A `mort_calibration` object from [calibrate()].
#' @param h Forecast horizon in years (default: from config or 10).
#' @param ... Ignored.
#' @return A `mort_forecast` object.
#' @export
forecast.mort_calibration <- function(object, h = NULL, ...) {
  config <- object$config %||% mort_config()
  if (is.null(h)) h <- config$forecast_h
  data   <- object$data
  ages   <- data$ages
  years  <- data$years
  groups <- data$groups
  pop_arr <- data$pop
  if (is.null(pop_arr)) pop_arr <- array(1, dim = dim(data$mx))

  forecast_years <- (max(years) + 1):(max(years) + h)
  cat("=== Forecasting (h =", h, ") ===\n")

  if (object$best_model == "Ensemble" && length(object$ensemble_names) >= 2) {
    ew <- object$ensemble_weights
    ens_names <- object$ensemble_names
    n_a <- data$n_a; n_g <- data$n_g

    ens_point <- array(0, dim = c(n_a, h, n_g))
    ens_fitted <- array(0, dim = c(n_a, data$n_t, n_g))
    w_total <- 0; wf_total <- 0
    full_fc_ref <- NULL

    for (ii in seq_along(ens_names)) {
      enm <- ens_names[ii]
      fc_res <- dispatch_forecast(enm, data$mx, pop_arr, years, ages,
                                   groups, h, object$calib_coda,
                                   object$calib_direct,
                                   width = config$age_width,
                                   calib_all = object$calib_all,
                                   return_fitted = TRUE)
      if (!is.null(fc_res$mx_forecast)) {
        fa <- fc_res$mx_forecast; fa[is.na(fa)] <- 0
        ens_point <- ens_point + ew[ii] * fa
        w_total <- w_total + ew[ii]
      }
      if (!is.null(fc_res$mx_fitted) &&
          all(dim(fc_res$mx_fitted) == c(n_a, data$n_t, n_g))) {
        ff <- fc_res$mx_fitted; ff[is.na(ff)] <- 0
        ens_fitted <- ens_fitted + ew[ii] * ff
        wf_total <- wf_total + ew[ii]
      }
      if (is.null(full_fc_ref)) full_fc_ref <- fc_res
    }
    if (w_total > 0 && w_total < 1) ens_point <- ens_point / w_total
    if (wf_total > 0 && wf_total < 1) ens_fitted <- ens_fitted / wf_total
    if (wf_total == 0) ens_fitted <- NULL

    mx_forecast <- ens_point
    mx_fitted   <- ens_fitted
    extras <- if (!is.null(full_fc_ref$full_fc)) full_fc_ref$full_fc else list()
    extras$method <- full_fc_ref$method  # preserve dispatch method (coda/direct/lca/...)
    model_label <- "Ensemble"
  } else {
    fc_res <- dispatch_forecast(object$best_model, data$mx, pop_arr, years,
                                 ages, groups, h, object$calib_coda,
                                 object$calib_direct,
                                 width = config$age_width,
                                 calib_all = object$calib_all,
                                 return_fitted = TRUE)
    mx_forecast <- fc_res$mx_forecast
    mx_fitted   <- fc_res$mx_fitted
    extras <- if (!is.null(fc_res$full_fc)) fc_res$full_fc else list()
    extras$method <- fc_res$method  # preserve dispatch method (coda/direct/lca/...)
    model_label <- object$best_model
  }

  # Build mort_data for train metadata
  train_md <- data

  best_mase <- NA_real_
  comp_best <- object$comparison[object$comparison$model == object$best_model, ]
  if (nrow(comp_best) > 0) best_mase <- comp_best$mase_combined[1]

  # Store calibration context in extras so add_pi() can use it
  extras$calib_coda   <- extras$calib_coda   %||% object$calib_coda
  extras$calib_direct <- extras$calib_direct  %||% object$calib_direct
  extras$calib_all    <- extras$calib_all     %||% object$calib_all
  extras$pop_arr      <- pop_arr
  extras$mx_smooth    <- data$mx
  extras$ens_names    <- if (model_label == "Ensemble") ens_names else NULL
  extras$ens_weights  <- if (model_label == "Ensemble") ew else NULL
  extras$gam_fits     <- data$gam_fits  # from smooth_mortality(), may be NULL

  mort_forecast(
    mx_point       = mx_forecast,
    train          = train_md,
    forecast_years = forecast_years,
    model          = model_label,
    mx_fitted      = mx_fitted,
    mx_raw         = if (!is.null(object$data_raw)) object$data_raw$mx else NULL,
    fit_info       = list(
      mase             = best_mase,
      fit_pct          = extras$fit_pct,
      model_comparison = object$comparison
    ),
    extras = extras
  )
}


# ---- add_pi ----

#' Add prediction intervals to a mort_forecast
#'
#' Computes bootstrap prediction intervals and attaches them to the
#' forecast object.
#'
#' @param fc A `mort_forecast` object.
#' @param n_sim Number of bootstrap simulations (default 500).
#' @param method One of `"bootstrap"` (standard residual/ARIMA bootstrap)
#'   or `"gam_bayes"` (GAM posterior propagation). The latter requires
#'   that GAM fits are stored in the forecast extras.
#' @param cap Maximum PI width on log scale (default 1.5).
#' @param calibration Optional `mort_calibration` object (needed for
#'   `"gam_bayes"` method to access ensemble weights and GAM fits).
#' @param ... Ignored.
#' @return An updated `mort_forecast` with `mx_lower`, `mx_upper`,
#'   `mx_median` filled in.
#' @export
add_pi <- function(fc, n_sim = 500L, method = c("bootstrap", "gam_bayes"),
                   cap = 1.5, calibration = NULL, ...) {
  valid_methods <- c("bootstrap", "gam_bayes")
  method <- match.arg(method, choices = valid_methods)
  stopifnot(inherits(fc, "mort_forecast"))

  h    <- fc$h
  n_a  <- fc$n_a
  n_g  <- fc$n_g
  ages <- fc$ages

  cat("=== Computing PIs (", method, ", n_sim =", n_sim, ") ===\n")

  extras <- fc$extras

  if (method == "bootstrap") {
    # Recover smoothed mx and pop from extras (stored by forecast.mort_calibration)
    mx_smooth <- extras$mx_smooth
    if (is.null(mx_smooth)) mx_smooth <- fc$mx_fitted
    if (is.null(mx_smooth)) mx_smooth <- fc$mx_raw
    if (is.null(mx_smooth))
      mx_smooth <- array(1e-4, dim = c(n_a, length(fc$train_years), n_g))

    pop_arr <- extras$pop_arr
    if (is.null(pop_arr))
      pop_arr <- array(1, dim = c(n_a, length(fc$train_years), n_g))

    groups     <- fc$groups
    disp_meth  <- extras$method %||% "direct"
    calib_all  <- extras$calib_all
    if (!is.null(calibration)) calib_all <- calibration$calib_all

    fc_result <- list(
      mx_forecast = fc$mx_point,
      mx_fitted   = fc$mx_fitted,
      full_fc     = extras,
      method      = disp_meth
    )

    sims <- dispatch_model_sims(
      fc$model, mx_smooth, pop_arr, fc$train_years, ages, groups,
      h, n_sim, fc_result,
      width = 5L,
      calib_all = calib_all,
      mx_raw = fc$mx_raw)

    bs <- list(
      mx_lower  = apply(sims, c(1, 2, 3), stats::quantile, 0.025, na.rm = TRUE),
      mx_upper  = apply(sims, c(1, 2, 3), stats::quantile, 0.975, na.rm = TRUE),
      mx_median = apply(sims, c(1, 2, 3), stats::quantile, 0.50,  na.rm = TRUE))

  } else if (method == "gam_bayes") {
    # Recover GAM fits from extras (stored by forecast.mort_calibration)
    gam_fits <- extras$gam_fits
    if (is.null(gam_fits))
      stop("gam_bayes requires GAM fits from smooth_mortality(method = 'gam'). ",
           "No gam_fits found in forecast extras.")

    pop_arr <- extras$pop_arr
    if (is.null(pop_arr))
      pop_arr <- array(1, dim = c(n_a, length(fc$train_years), n_g))

    groups      <- fc$groups
    calib_coda  <- extras$calib_coda
    calib_direct <- extras$calib_direct
    calib_all   <- extras$calib_all
    if (!is.null(calibration)) {
      calib_coda   <- calibration$calib_coda
      calib_direct <- calibration$calib_direct
      calib_all    <- calibration$calib_all
    }

    # Determine ensemble names and weights
    ens_names   <- extras$ens_names %||% fc$model
    ens_weights <- extras$ens_weights %||% 1
    if (length(ens_names) == 1 && !is.character(ens_names))
      ens_names <- fc$model
    if (length(ens_weights) == 1 && is.null(names(ens_weights)))
      names(ens_weights) <- ens_names

    sims_raw <- bootstrap_gam_bayes(
      gam_fits, ages, fc$train_years, pop_arr, groups,
      h, n_sim, ens_names, ens_weights,
      calib_coda, calib_direct, calib_all,
      width = 5L)

    # Extract quantiles from sims_raw (returns list with mx_lower/mx_upper/mx_median
    # or a 4D array depending on implementation)
    if (is.list(sims_raw) && !is.null(sims_raw$mx_lower)) {
      bs <- sims_raw
      sims <- sims_raw$sims
    } else {
      sims <- sims_raw
      bs <- list(
        mx_lower  = apply(sims, c(1, 2, 3), stats::quantile, 0.025, na.rm = TRUE),
        mx_upper  = apply(sims, c(1, 2, 3), stats::quantile, 0.975, na.rm = TRUE),
        mx_median = apply(sims, c(1, 2, 3), stats::quantile, 0.50,  na.rm = TRUE))
    }
  }

  # Re-center and cap PIs
  log_fc  <- log(pmax(fc$mx_point, 1e-12))
  log_lo  <- log(pmax(bs$mx_lower, 1e-12))
  log_hi  <- log(pmax(bs$mx_upper, 1e-12))
  log_med <- log(pmax(bs$mx_median, 1e-12))
  shift   <- log_fc - log_med
  log_lo  <- log_lo + shift
  log_hi  <- log_hi + shift

  log_lo <- pmax(log_lo, log_fc - cap)
  log_hi <- pmin(log_hi, log_fc + cap)

  fc$mx_lower  <- exp(log_lo)
  fc$mx_upper  <- exp(log_hi)
  fc$mx_median <- fc$mx_point
  fc$sims      <- if (exists("sims", inherits = FALSE)) sims else NULL

  cat("  PI re-centered (max |shift| =", round(max(abs(shift)), 4),
      "), capped at +/-", cap, "on log scale.\n")
  fc
}


# ---- disaggregate ----

#' Disaggregate a mort_forecast to single-year ages
#'
#' Converts grouped-age forecasts to single-year age forecasts using
#' spline interpolation with parametric tail closure.
#'
#' @param fc A `mort_forecast` object with grouped ages.
#' @param target_ages Integer vector of single-year target ages
#'   (default `30:100`).
#' @param method One of `"natural_gompertz"`, `"hyman_kannisto"`, or
#'   `"gam"`.
#' @param width Width of closed age groups (default 5).
#' @param ... Additional arguments for the disaggregation function.
#' @return A new `mort_forecast` object with single-year ages.
#' @export
disaggregate <- function(fc, target_ages = 30L:100L,
                         method = c("natural_gompertz", "hyman_kannisto", "gam"),
                         width = 5L, ...) {
  method <- match.arg(method)
  stopifnot(inherits(fc, "mort_forecast"))

  ages <- fc$ages
  n_target <- length(target_ages)

  cat("=== Disaggregating to single-year ages (", method, ") ===\n")

  # Disaggregate point forecasts
  mx_1yr <- array(NA_real_, dim = c(n_target, fc$h, fc$n_g))
  disagg_fn <- switch(method,
    natural_gompertz = disaggregate_natural_gompertz,
    hyman_kannisto   = disaggregate_hyman_kannisto,
    gam = NULL
  )

  if (method == "gam") {
    stop("GAM disaggregation requires observed data. Use ",
         "disaggregate_forecast_gam() directly or run_forecast() pipeline.")
  }

  for (k in seq_len(fc$n_g)) {
    for (j in seq_len(fc$h)) {
      mx_1yr[, j, k] <- disagg_fn(fc$mx_point[, j, k], ages,
                                    target_ages, width)
    }
  }

  # Disaggregate PIs via ratio interpolation
  mx_1yr_lo <- mx_1yr_hi <- mx_1yr_med <- NULL
  has_pi <- !is.null(fc$mx_lower) && !is.null(fc$mx_upper)
  if (has_pi) {
    mx_1yr_lo  <- array(NA_real_, dim = c(n_target, fc$h, fc$n_g))
    mx_1yr_hi  <- array(NA_real_, dim = c(n_target, fc$h, fc$n_g))
    mx_1yr_med <- array(NA_real_, dim = c(n_target, fc$h, fc$n_g))

    for (k in seq_len(fc$n_g)) {
      for (j in seq_len(fc$h)) {
        pi_res <- disaggregate_pi_ratios(
          mx_1yr[, j, k], fc$mx_point[, j, k],
          fc$mx_lower[, j, k], fc$mx_upper[, j, k],
          fc$mx_median[, j, k] %||% fc$mx_point[, j, k],
          target_ages, ages, width)
        mx_1yr_lo[, j, k]  <- pi_res$mx_1yr_lo95
        mx_1yr_hi[, j, k]  <- pi_res$mx_1yr_hi95
        mx_1yr_med[, j, k] <- pi_res$mx_1yr_median
      }
    }
  }

  # Disaggregate fitted values to single-year ages (for accuracy/plotting)
  mx_fitted_1yr <- NULL
  if (!is.null(fc$mx_fitted)) {
    n_t_train <- dim(fc$mx_fitted)[2]
    mx_fitted_1yr <- array(NA_real_, dim = c(n_target, n_t_train, fc$n_g))
    for (k in seq_len(fc$n_g)) {
      for (j in seq_len(n_t_train)) {
        mx_fitted_1yr[, j, k] <- tryCatch(
          disagg_fn(fc$mx_fitted[, j, k], ages, target_ages, width),
          error = function(e) rep(NA_real_, n_target))
      }
    }
    if (all(is.na(mx_fitted_1yr))) mx_fitted_1yr <- NULL
  }

  # Disaggregate raw observations to single-year ages
  mx_raw_1yr <- NULL
  if (!is.null(fc$mx_raw)) {
    n_t_train <- dim(fc$mx_raw)[2]
    mx_raw_1yr <- array(NA_real_, dim = c(n_target, n_t_train, fc$n_g))
    for (k in seq_len(fc$n_g)) {
      for (j in seq_len(n_t_train)) {
        mx_raw_1yr[, j, k] <- tryCatch(
          disagg_fn(fc$mx_raw[, j, k], ages, target_ages, width),
          error = function(e) rep(NA_real_, n_target))
      }
    }
    if (all(is.na(mx_raw_1yr))) mx_raw_1yr <- NULL
  }

  # Build new age labels for single years
  age_labels_1yr <- as.character(target_ages)

  # Create a fake train mort_data for the constructor
  train_1yr <- structure(
    list(
      ages       = target_ages,
      years      = fc$train_years,
      groups     = fc$groups,
      cause      = fc$cause,
      age_labels = age_labels_1yr,
      n_a        = n_target,
      n_t        = length(fc$train_years),
      n_g        = fc$n_g
    ),
    class = "mort_data"
  )

  mort_forecast(
    mx_point       = mx_1yr,
    train          = train_1yr,
    forecast_years = fc$forecast_years,
    model          = paste0(fc$model, " (disagg:", method, ")"),
    mx_median      = mx_1yr_med,
    mx_lower       = mx_1yr_lo,
    mx_upper       = mx_1yr_hi,
    mx_fitted      = mx_fitted_1yr,
    mx_raw         = mx_raw_1yr,
    fit_info       = fc$fit_info,
    extras         = fc$extras
  )
}

# ---- Multi-cause reconciliation (S3 workflow) ----

#' Reconcile multiple cause-specific forecasts to sum to an all-cause forecast
#'
#' Takes independently produced `mort_forecast` objects (one per cause) and
#' an all-cause forecast, then applies WLS reconciliation so disease forecasts
#' sum to the all-cause total. Each reconciled result remains a valid
#' `mort_forecast`.
#'
#' @param forecasts Named list of `mort_forecast` objects (one per disease).
#' @param all_cause A `mort_forecast` object for all-cause mortality.
#' @param recon_cap Soft cap: max |delta_k| as fraction of point forecast
#'   (default 0.5).
#' @return A named list of reconciled `mort_forecast` objects.
#' @export
reconcile <- function(forecasts, all_cause, recon_cap = 0.5) {
  stopifnot(inherits(all_cause, "mort_forecast"))
  stopifnot(is.list(forecasts), length(forecasts) > 0)
  for (nm in names(forecasts))
    stopifnot(inherits(forecasts[[nm]], "mort_forecast"))

  disease_order <- names(forecasts)
  mx_all_fc <- all_cause$mx_point

  # Extract arrays from mort_forecast objects

  mx_indep <- list()
  bs_indep <- list()
  for (nm in disease_order) {
    fc <- forecasts[[nm]]
    mx_indep[[nm]] <- fc$mx_point
    if (!is.null(fc$mx_lower) && !is.null(fc$mx_upper)) {
      bs_indep[[nm]] <- list(
        mx_lower  = fc$mx_lower,
        mx_upper  = fc$mx_upper,
        mx_median = fc$mx_median %||% fc$mx_point)
    } else {
      # No PI: use +/- 20% as fallback variance estimate
      bs_indep[[nm]] <- list(
        mx_lower  = fc$mx_point * 0.8,
        mx_upper  = fc$mx_point * 1.2,
        mx_median = fc$mx_point)
    }
  }

  recon <- wls_reconcile(mx_all_fc, mx_indep, bs_indep, disease_order,
                          recon_cap = recon_cap)

  # Wrap reconciled arrays back into mort_forecast objects
  result <- list()
  for (nm in disease_order) {
    fc <- forecasts[[nm]]
    result[[nm]] <- mort_forecast(
      mx_point       = recon$mx_recon[[nm]],
      train          = structure(
        list(mx = fc$extras$mx_smooth %||%
               array(NA_real_, c(fc$n_a, length(fc$train_years), fc$n_g)),
             ages = fc$ages, years = fc$train_years,
             groups = fc$groups, pop = NULL,
             cause = fc$cause, age_labels = fc$age_labels,
             smoothed = TRUE, n_a = fc$n_a,
             n_t = length(fc$train_years), n_g = fc$n_g),
        class = "mort_data"),
      forecast_years = fc$forecast_years,
      model          = paste0(fc$model, " (reconciled)"),
      mx_fitted      = fc$mx_fitted,
      mx_raw         = fc$mx_raw,
      mx_lower       = recon$bs_recon[[nm]]$mx_lower,
      mx_upper       = recon$bs_recon[[nm]]$mx_upper,
      mx_median      = recon$bs_recon[[nm]]$mx_median,
      fit_info       = fc$fit_info,
      extras         = fc$extras
    )
  }
  result
}
