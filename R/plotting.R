#' Plotting functions for mortality forecasts
#'
#' @name plotting
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_ribbon
#'   geom_col geom_vline geom_hline geom_abline geom_area geom_boxplot
#'   facet_wrap facet_grid scale_color_brewer scale_fill_brewer
#'   scale_color_manual scale_fill_manual scale_fill_brewer
#'   coord_flip labs theme_minimal theme element_text ggsave
NULL

#' Standard colour palette for model comparison
#' @export
model_colors <- c(
 "Tucker3-CoDa" = "#E41A1C", "Tucker3-Direct" = "#FF7F00",
 "Lee-Carter" = "#377EB8", "FDM" = "#4DAF4A",
 "StMoMo-APC" = "#984EA3", "RWD" = "#A65628",
 "ARIMA-perSeries" = "#F781BF", "ETS-Damped" = "#999999",
 "Observed" = "black"
)

# Internal: determine color and facet columns from group_cols vector
.plot_group_setup <- function(group_cols) {
  if (length(group_cols) == 0) return(list(color_col = NULL, facet_cols = NULL))
  list(
    color_col  = group_cols[length(group_cols)],
    facet_cols = if (length(group_cols) > 1) group_cols[-length(group_cols)] else NULL
  )
}

#' Bar chart of model MASE scores
#' @param comp_dt  data.table from [compare_models()].
#' @param output_dir Output directory.
#' @param cause_label Full cause name.
#' @param cause_short Short cause key.
#' @export
plot_model_comparison <- function(comp_dt, output_dir, cause_label, cause_short) {
  p <- ggplot2::ggplot(comp_dt[!is.na(mase)],
                       ggplot2::aes(x = reorder(model, mase), y = mase, fill = model)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    ggplot2::scale_fill_manual(values = model_colors) +
    ggplot2::labs(title = paste0("Model Comparison: ", cause_label, " (MASE)"),
                  subtitle = "Dashed = naive (MASE=1)", x = NULL, y = "MASE") +
    ggplot2::coord_flip() + ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(output_dir, paste0("plot_model_mase_", cause_short, ".png")),
                  p, width = 8, height = 5, dpi = 150)
}

#' Forecast trend plot (age-averaged, by group)
#' @param obs_dt    Observed data.table with `mx_obs`.
#' @param fitted_dt Fitted data.table with `mx_fitted`.
#' @param fc_dt     Forecast data.table with `mx_adjusted`, `mx_lo95_adj`, `mx_hi95_adj`.
#' @param output_dir Output directory path.
#' @param cause_label Full cause name.
#' @param cause_short Short cause key.
#' @param forecast_years Integer vector of forecast years.
#' @param group_cols Character vector of group column names in the data.tables.
#' @param smooth_dt Optional data.table with `mx_lo95`, `mx_hi95` for Bayesian CI ribbon.
#' @export
plot_forecast_trends <- function(obs_dt, fitted_dt, fc_dt,
                                 output_dir, cause_label, cause_short,
                                 forecast_years, group_cols = character(0),
                                 smooth_dt = NULL) {
  gs <- .plot_group_setup(group_cols)
  cc <- gs$color_col
  by_cols <- c("year", group_cols)

  obs_tr <- obs_dt[, .(mx = mean(mx_obs)), by = by_cols]
  fit_tr <- fitted_dt[, .(mx = mean(mx_fitted)), by = by_cols]
  fc_tr  <- fc_dt[, .(mx    = mean(mx_adjusted),
                       mx_lo = mean(mx_lo95_adj),
                       mx_hi = mean(mx_hi95_adj)),
                  by = by_cols]

  has_ci <- !is.null(smooth_dt) && "mx_lo95" %in% names(smooth_dt) &&
            !all(is.na(smooth_dt$mx_lo95))
  if (has_ci) {
    ci_tr <- smooth_dt[, .(mx_lo_s = mean(mx_lo95, na.rm = TRUE),
                           mx_hi_s = mean(mx_hi95, na.rm = TRUE)),
                       by = by_cols]
    fit_tr <- merge(fit_tr, ci_tr, by = by_cols, all.x = TRUE)
  }

  bridge <- fit_tr[year == max(fit_tr$year)]
  bridge[, `:=`(mx_lo = mx, mx_hi = mx)]
  fc_comb <- data.table::rbindlist(list(bridge, fc_tr), fill = TRUE)

  if (!is.null(cc)) {
    obs_tr[, .group := factor(obs_tr[[cc]])]
    fit_tr[, .group := factor(fit_tr[[cc]])]
    fc_comb[, .group := factor(fc_comb[[cc]])]
  }

  p <- ggplot2::ggplot()
  if (has_ci && !is.null(cc)) {
    p <- p + ggplot2::geom_ribbon(data = fit_tr,
      ggplot2::aes(x = year, ymin = log(mx_lo_s), ymax = log(mx_hi_s),
                   fill = .group), alpha = 0.08)
  } else if (has_ci) {
    p <- p + ggplot2::geom_ribbon(data = fit_tr,
      ggplot2::aes(x = year, ymin = log(mx_lo_s), ymax = log(mx_hi_s)),
      alpha = 0.08)
  }
  if (!is.null(cc)) {
    p <- p +
      ggplot2::geom_point(data = obs_tr,
        ggplot2::aes(x = year, y = log(mx), color = .group),
        size = 0.8, alpha = 0.4) +
      ggplot2::geom_line(data = fit_tr,
        ggplot2::aes(x = year, y = log(mx), color = .group),
        linewidth = 0.7) +
      ggplot2::geom_ribbon(data = fc_comb,
        ggplot2::aes(x = year, ymin = log(mx_lo), ymax = log(mx_hi),
                     fill = .group), alpha = 0.15) +
      ggplot2::geom_line(data = fc_comb,
        ggplot2::aes(x = year, y = log(mx), color = .group),
        linewidth = 0.8, linetype = "dashed") +
      ggplot2::scale_color_brewer(palette = "RdYlGn", direction = -1, name = cc) +
      ggplot2::scale_fill_brewer(palette = "RdYlGn", direction = -1, name = cc)
  } else {
    p <- p +
      ggplot2::geom_point(data = obs_tr,
        ggplot2::aes(x = year, y = log(mx)),
        size = 0.8, alpha = 0.4) +
      ggplot2::geom_line(data = fit_tr,
        ggplot2::aes(x = year, y = log(mx)),
        linewidth = 0.7) +
      ggplot2::geom_ribbon(data = fc_comb,
        ggplot2::aes(x = year, ymin = log(mx_lo), ymax = log(mx_hi)),
        alpha = 0.15) +
      ggplot2::geom_line(data = fc_comb,
        ggplot2::aes(x = year, y = log(mx)),
        linewidth = 0.8, linetype = "dashed")
  }
  p <- p +
    ggplot2::geom_vline(xintercept = min(forecast_years) - 0.5,
                        linetype = "dotted", color = "grey50")
  if (!is.null(gs$facet_cols))
    p <- p + ggplot2::facet_wrap(stats::reformulate(gs$facet_cols),
                                 scales = "free_y")
  p <- p +
    ggplot2::labs(title = paste0(cause_label, " log(mx) Forecast Trends"),
                  subtitle = if (has_ci) "Dots=observed, solid=fitted (95% Bayes CI), dashed=forecast (95% PI)"
                             else "Dots=observed, solid=fitted, dashed=forecast (95% PI)",
                  x = "Year", y = "Mean log(mx)") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(file.path(output_dir, paste0("plot_fc_trend_", cause_short, ".png")),
                  p, width = 10, height = 5, dpi = 150)
}

#' Observed vs smoothed comparison plot
#' @param obs_dt    data.table with `mx_obs`.
#' @param smooth_dt data.table with `mx_smooth` (and optionally `mx_lo95`, `mx_hi95`).
#' @param output_dir,cause_label,cause_short Standard plot params.
#' @param group_cols Character vector of group column names.
#' @export
plot_obs_vs_smooth <- function(obs_dt, smooth_dt, output_dir,
                               cause_label, cause_short,
                               group_cols = character(0)) {
  gs <- .plot_group_setup(group_cols)
  cc <- gs$color_col
  by_cols <- c("year", group_cols)

  obs_tr <- obs_dt[, .(mx = mean(mx_obs)), by = by_cols]
  sm_tr  <- smooth_dt[, .(mx = mean(mx_smooth)), by = by_cols]

  has_ci <- "mx_lo95" %in% names(smooth_dt) &&
            !all(is.na(smooth_dt$mx_lo95))
  if (has_ci) {
    ci_tr <- smooth_dt[, .(mx_lo = mean(mx_lo95, na.rm = TRUE),
                           mx_hi = mean(mx_hi95, na.rm = TRUE)),
                       by = by_cols]
    sm_tr <- merge(sm_tr, ci_tr, by = by_cols)
  }

  if (!is.null(cc)) {
    obs_tr[, .group := factor(obs_tr[[cc]])]
    sm_tr[, .group := factor(sm_tr[[cc]])]
  }

  p <- ggplot2::ggplot()
  if (has_ci && !is.null(cc)) {
    p <- p + ggplot2::geom_ribbon(data = sm_tr,
      ggplot2::aes(x = year, ymin = log(mx_lo), ymax = log(mx_hi),
                   fill = .group), alpha = 0.12)
  } else if (has_ci) {
    p <- p + ggplot2::geom_ribbon(data = sm_tr,
      ggplot2::aes(x = year, ymin = log(mx_lo), ymax = log(mx_hi)),
      alpha = 0.12)
  }
  if (!is.null(cc)) {
    p <- p +
      ggplot2::geom_point(data = obs_tr,
        ggplot2::aes(x = year, y = log(mx), color = .group),
        size = 0.8, alpha = 0.4) +
      ggplot2::geom_line(data = sm_tr,
        ggplot2::aes(x = year, y = log(mx), color = .group),
        linewidth = 0.7) +
      ggplot2::scale_color_brewer(palette = "RdYlGn", direction = -1, name = cc) +
      ggplot2::scale_fill_brewer(palette = "RdYlGn", direction = -1, name = cc)
  } else {
    p <- p +
      ggplot2::geom_point(data = obs_tr,
        ggplot2::aes(x = year, y = log(mx)),
        size = 0.8, alpha = 0.4) +
      ggplot2::geom_line(data = sm_tr,
        ggplot2::aes(x = year, y = log(mx)),
        linewidth = 0.7)
  }
  if (!is.null(gs$facet_cols))
    p <- p + ggplot2::facet_wrap(stats::reformulate(gs$facet_cols),
                                 scales = "free_y")
  p <- p +
    ggplot2::labs(title = paste0(cause_label, " Observed vs Smoothed log(mx)"),
                  subtitle = if (has_ci) "Dots=observed, lines=GAM smooth, shading=95% Bayesian CI"
                             else "Dots=raw observed, lines=smoothed",
                  x = "Year", y = "Mean log(mx)") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(
    file.path(output_dir, paste0("plot_obs_smooth_", cause_short, ".png")),
    p, width = 10, height = 5, dpi = 150)
}

#' Forecast trends faceted by age group
#' @inheritParams plot_forecast_trends
#' @export
plot_forecast_trends_by_age <- function(obs_dt, fitted_dt, fc_dt,
                                        output_dir, cause_label, cause_short,
                                        forecast_years, group_cols = character(0),
                                        smooth_dt = NULL) {
  gs <- .plot_group_setup(group_cols)
  cc <- gs$color_col
  by_cols <- c("year", group_cols, "age_group")

  rep_ages <- c("40-44", "50-54", "60-64", "70-74", "80-84", "90+")
  rep_ages <- intersect(rep_ages, unique(obs_dt$age_group))
  if (length(rep_ages) == 0) return(invisible(NULL))

  obs_sub <- obs_dt[age_group %in% rep_ages]
  fit_sub <- fitted_dt[age_group %in% rep_ages]
  fc_sub  <- fc_dt[age_group %in% rep_ages]

  obs_ag <- obs_sub[, .(mx = mean(mx_obs)), by = by_cols]
  fit_ag <- fit_sub[, .(mx = mean(mx_fitted)), by = by_cols]
  fc_ag  <- fc_sub[, .(mx    = mean(mx_adjusted),
                        mx_lo = mean(mx_lo95_adj),
                        mx_hi = mean(mx_hi95_adj)),
                   by = by_cols]

  has_ci <- !is.null(smooth_dt) && "mx_lo95" %in% names(smooth_dt) &&
            !all(is.na(smooth_dt$mx_lo95))
  if (has_ci) {
    sm_sub <- smooth_dt[age_group %in% rep_ages]
    ci_ag <- sm_sub[, .(mx_lo_s = mean(mx_lo95, na.rm = TRUE),
                        mx_hi_s = mean(mx_hi95, na.rm = TRUE)),
                    by = by_cols]
    fit_ag <- merge(fit_ag, ci_ag, by = by_cols, all.x = TRUE)
  }

  bridge <- fit_ag[year == max(fit_ag$year)]
  bridge[, `:=`(mx_lo = mx, mx_hi = mx)]
  fc_comb <- data.table::rbindlist(list(bridge, fc_ag), fill = TRUE)

  obs_ag[, age_group := factor(age_group, levels = rep_ages)]
  fit_ag[, age_group := factor(age_group, levels = rep_ages)]
  fc_comb[, age_group := factor(age_group, levels = rep_ages)]

  if (!is.null(cc)) {
    obs_ag[, .group := factor(obs_ag[[cc]])]
    fit_ag[, .group := factor(fit_ag[[cc]])]
    fc_comb[, .group := factor(fc_comb[[cc]])]
  }

  p <- ggplot2::ggplot()
  if (has_ci && !is.null(cc)) {
    p <- p + ggplot2::geom_ribbon(data = fit_ag,
      ggplot2::aes(x = year, ymin = log(mx_lo_s), ymax = log(mx_hi_s),
                   fill = .group), alpha = 0.08)
  } else if (has_ci) {
    p <- p + ggplot2::geom_ribbon(data = fit_ag,
      ggplot2::aes(x = year, ymin = log(mx_lo_s), ymax = log(mx_hi_s)),
      alpha = 0.08)
  }
  if (!is.null(cc)) {
    p <- p +
      ggplot2::geom_point(data = obs_ag,
        ggplot2::aes(x = year, y = log(mx), color = .group),
        size = 0.5, alpha = 0.35) +
      ggplot2::geom_line(data = fit_ag,
        ggplot2::aes(x = year, y = log(mx), color = .group),
        linewidth = 0.5) +
      ggplot2::geom_ribbon(data = fc_comb,
        ggplot2::aes(x = year, ymin = log(mx_lo), ymax = log(mx_hi),
                     fill = .group), alpha = 0.12) +
      ggplot2::geom_line(data = fc_comb,
        ggplot2::aes(x = year, y = log(mx), color = .group),
        linewidth = 0.6, linetype = "dashed") +
      ggplot2::scale_color_brewer(palette = "RdYlGn", direction = -1, name = cc) +
      ggplot2::scale_fill_brewer(palette = "RdYlGn", direction = -1, name = cc)
  } else {
    p <- p +
      ggplot2::geom_point(data = obs_ag,
        ggplot2::aes(x = year, y = log(mx)),
        size = 0.5, alpha = 0.35) +
      ggplot2::geom_line(data = fit_ag,
        ggplot2::aes(x = year, y = log(mx)),
        linewidth = 0.5) +
      ggplot2::geom_ribbon(data = fc_comb,
        ggplot2::aes(x = year, ymin = log(mx_lo), ymax = log(mx_hi)),
        alpha = 0.12) +
      ggplot2::geom_line(data = fc_comb,
        ggplot2::aes(x = year, y = log(mx)),
        linewidth = 0.6, linetype = "dashed")
  }
  p <- p +
    ggplot2::geom_vline(xintercept = min(forecast_years) - 0.5,
                        linetype = "dotted", color = "grey50")
  if (!is.null(gs$facet_cols)) {
    facet_f <- stats::reformulate(gs$facet_cols, response = "age_group")
    p <- p + ggplot2::facet_grid(facet_f, scales = "free_y")
  } else {
    p <- p + ggplot2::facet_wrap(~age_group, scales = "free_y")
  }
  p <- p +
    ggplot2::labs(title = paste0(cause_label, " log(mx) by Age & Group"),
                  subtitle = if (has_ci) "Dots=observed, solid=fitted (95% Bayes CI), dashed=forecast (95% PI)"
                             else "Dots=observed, solid=fitted, dashed=forecast (95% PI)",
                  x = "Year", y = "log(mx)") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8),
                   legend.position = "bottom")
  ggplot2::ggsave(
    file.path(output_dir, paste0("plot_fc_trend_byage_", cause_short, ".png")),
    p, width = 12, height = 14, dpi = 150)
}

#' Forecast age profile for a single year
#' @param fc_dt data.table of forecasts.
#' @param output_dir,cause_label,cause_short Standard plot params.
#' @param age_groups Ordered character vector of age-group labels.
#' @param group_cols Character vector of group column names.
#' @export
plot_forecast_age <- function(fc_dt, output_dir, cause_label, cause_short,
                              age_groups, group_cols = character(0)) {
  gs <- .plot_group_setup(group_cols)
  cc <- gs$color_col
  last_yr <- max(fc_dt$year)
  fc_sub  <- fc_dt[year == last_yr]

  # If facet columns exist, filter to first level of first facet col
  if (!is.null(gs$facet_cols)) {
    first_facet <- gs$facet_cols[1]
    first_val <- fc_sub[[first_facet]][1]
    fc_sub <- fc_sub[fc_sub[[first_facet]] == first_val, ]
  }
  fc_sub[, age_group := factor(age_group, levels = age_groups)]

  if (!is.null(cc)) {
    fc_sub[, .group := factor(fc_sub[[cc]])]
    p <- ggplot2::ggplot(fc_sub,
          ggplot2::aes(x = age_group, group = .group)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = log(mx_lo95_adj), ymax = log(mx_hi95_adj),
                     fill = .group), alpha = 0.2) +
      ggplot2::geom_line(
        ggplot2::aes(y = log(mx_adjusted), color = .group), linewidth = 1) +
      ggplot2::scale_color_brewer(palette = "RdYlGn", direction = -1, name = cc) +
      ggplot2::scale_fill_brewer(palette = "RdYlGn", direction = -1, name = cc)
  } else {
    p <- ggplot2::ggplot(fc_sub,
          ggplot2::aes(x = age_group, group = 1)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = log(mx_lo95_adj), ymax = log(mx_hi95_adj)),
        alpha = 0.2) +
      ggplot2::geom_line(
        ggplot2::aes(y = log(mx_adjusted)), linewidth = 1)
  }
  p <- p +
    ggplot2::labs(title = paste0(cause_label, " Forecast log(mx) ", last_yr),
                  x = "Age Group", y = "log(mx)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(
    file.path(output_dir, paste0("plot_fc_age_", cause_short, ".png")),
    p, width = 9, height = 5, dpi = 150)
}

#' Plot disaggregated single-year mortality
#' @param mx_1yr_dt data.table with `age`, `mx_1yr`, `mx_1yr_lo95`, `mx_1yr_hi95`.
#' @param output_dir,cause_label,cause_short Standard plot params.
#' @param group_cols Character vector of group column names.
#' @export
plot_disaggregation <- function(mx_1yr_dt, output_dir, cause_label, cause_short,
                                group_cols = character(0)) {
  gs <- .plot_group_setup(group_cols)
  cc <- gs$color_col
  last_yr <- max(mx_1yr_dt$year)
  dsub <- mx_1yr_dt[year == last_yr]

  # Filter to first level of facet columns (if any)
  if (!is.null(gs$facet_cols)) {
    for (fc in gs$facet_cols) {
      first_val <- dsub[[fc]][1]
      dsub <- dsub[dsub[[fc]] == first_val, ]
    }
  }

  if (!is.null(cc)) {
    lvls <- sort(unique(dsub[[cc]]))
    if (length(lvls) > 2) lvls <- lvls[c(1, length(lvls))]
    dsub <- dsub[dsub[[cc]] %in% lvls, ]
    dsub[, .group := factor(dsub[[cc]])]

    p <- ggplot2::ggplot(dsub,
          ggplot2::aes(x = age, y = log(mx_1yr), color = .group)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = log(mx_1yr_lo95), ymax = log(mx_1yr_hi95),
                     fill = .group), alpha = 0.15) +
      ggplot2::scale_color_brewer(palette = "RdYlGn", direction = -1, name = cc) +
      ggplot2::scale_fill_brewer(palette = "RdYlGn", direction = -1, guide = "none")
  } else {
    p <- ggplot2::ggplot(dsub,
          ggplot2::aes(x = age, y = log(mx_1yr))) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = log(mx_1yr_lo95), ymax = log(mx_1yr_hi95)),
        alpha = 0.15)
  }
  p <- p +
    ggplot2::labs(title = paste0(cause_label, " Single-Year log(mx) (", last_yr, ")"),
                  x = "Age", y = "log(mx)") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(
    file.path(output_dir, paste0("plot_disagg_", cause_short, ".png")),
    p, width = 9, height = 5, dpi = 150)
}

#' Plot B-mode time components
#' @param fc_obj       Tucker3 forecast object (with `B` and `arima_fits`).
#' @param forecast_years Integer vector of future years.
#' @param sim_years    Integer vector of observed years.
#' @param output_dir,cause_label,cause_short Standard plot params.
#' @export
plot_B_components <- function(fc_obj, forecast_years, sim_years,
                              output_dir, cause_label, cause_short) {
  R_B <- ncol(fc_obj$B)
  B_list <- list()
  for (r in seq_len(R_B)) {
    B_list[[length(B_list) + 1]] <- data.table::data.table(
      year = sim_years, value = fc_obj$B[, r], comp = paste0("B", r),
      type = "observed", lo95 = NA_real_, hi95 = NA_real_)
    fc_r <- forecast::forecast(fc_obj$arima_fits[[r]],
                               h = length(forecast_years), level = 95)
    B_list[[length(B_list) + 1]] <- data.table::data.table(
      year = forecast_years, value = as.numeric(fc_r$mean),
      comp = paste0("B", r), type = "forecast",
      lo95 = as.numeric(fc_r$lower), hi95 = as.numeric(fc_r$upper))
  }
  B_dt <- data.table::rbindlist(B_list, fill = TRUE)

  p <- ggplot2::ggplot(B_dt, ggplot2::aes(x = year, y = value)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo95, ymax = hi95),
                         fill = "lightblue", alpha = 0.4, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(color = type), linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c(observed = "black", forecast = "red")) +
    ggplot2::facet_wrap(~comp, scales = "free_y") +
    ggplot2::labs(title = paste0("B-components: ", cause_label, " (R_B=", R_B, ")"),
                  x = "Year", y = "Score") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(
    file.path(output_dir, paste0("plot_B_components_", cause_short, ".png")),
    p, width = 10, height = 6, dpi = 150)
}

#' Scatterplot of raw vs smoothed log(mx)
#' @param mx_raw,mx_smooth 3D mortality arrays.
#' @param years  Integer vector.
#' @param ages   Integer vector.
#' @param output_dir,cause_short Standard plot params.
#' @param covid_years Years to exclude from comparison.
#' @param floor_threshold Cells below this are excluded.
#' @export
plot_smoothing_comparison <- function(mx_raw, mx_smooth, years, ages,
                                     output_dir, cause_short, covid_years,
                                     floor_threshold = 1e-7) {
  obs_idx <- which(!years %in% covid_years)
  raw_v <- as.vector(log(pmax(mx_raw[, obs_idx, ], 1e-10)))
  smo_v <- as.vector(log(pmax(mx_smooth[, obs_idx, ], 1e-10)))
  real_v <- as.vector(mx_raw[, obs_idx, ] >= floor_threshold)
  ok <- is.finite(raw_v) & is.finite(smo_v) & real_v
  p_dt <- data.table::data.table(raw = raw_v[ok], smooth = smo_v[ok])
  rho <- stats::cor(p_dt$raw, p_dt$smooth)

  p <- ggplot2::ggplot(p_dt, ggplot2::aes(x = raw, y = smooth)) +
    ggplot2::geom_point(alpha = 0.1, size = 0.5) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "red",
                         linetype = "dashed") +
    ggplot2::labs(title = paste0("Smoothing: log(mx) raw vs smoothed (r=",
                                 round(rho, 4), ")"),
                  x = "log(mx) raw", y = "log(mx) smoothed") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(
    file.path(output_dir, paste0("plot_smoothing_", cause_short, ".png")),
    p, width = 6, height = 6, dpi = 150)
}

#' Disaggregation method comparison plot
#' @param mx_grouped Grouped rates vector.
#' @param age_starts Age-group starts.
#' @param target_ages Single-year target ages.
#' @param output_dir,cause_short Standard plot params.
#' @param width Age-group width.
#' @export
plot_disagg_comparison <- function(mx_grouped, age_starts, target_ages,
                                   output_dir, cause_short, width = 5) {
  mx_hk <- disaggregate_hyman_kannisto(mx_grouped, age_starts, target_ages, width)
  mx_ng <- disaggregate_natural_gompertz(mx_grouped, age_starts, target_ages, width)
  dt <- data.table::rbindlist(list(
    data.table::data.table(age = target_ages, log_mx = log(mx_hk),
                           method = "Hyman+Kannisto"),
    data.table::data.table(age = target_ages, log_mx = log(mx_ng),
                           method = "Natural+Gompertz")))
  n_closed <- length(age_starts) - 1
  mids <- c(age_starts[1:n_closed] + width / 2, age_starts[length(age_starts)] + 5)
  pts <- data.table::data.table(age = mids, log_mx = log(mx_grouped), method = "Grouped")

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = age, y = log_mx, color = method)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(data = pts, size = 3, shape = 16) +
    ggplot2::labs(title = "Age Disaggregation Comparison", x = "Age", y = "log(mx)") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(
    file.path(output_dir, paste0("plot_disagg_methods_", cause_short, ".png")),
    p, width = 9, height = 5, dpi = 150)
}

#' Stacked area chart of cause shares (coherent mode)
#' @param shares_obs_4d  4D observed shares array.
#' @param shares_fc_4d   4D forecast shares array.
#' @param years_obs,years_fc Year vectors.
#' @param disease_order Character vector of disease names.
#' @param output_dir    Output directory.
#' @param group_index   Integer index into \code{groups} (default 1).
#' @param groups        data.frame of groups.
#' @param sex_filter,quintile_filter Deprecated strata filters (backward compat).
#' @export
plot_share_evolution <- function(shares_obs_4d, shares_fc_4d,
                                years_obs, years_fc,
                                disease_order, output_dir,
                                group_index = 1, groups,
                                sex_filter = NULL, quintile_filter = NULL) {
  K <- length(disease_order)
  cause_names <- c(disease_order, "remainder")

  # Backward compatibility: if sex_filter/quintile_filter given, find k_idx

  if (!is.null(sex_filter) && !is.null(quintile_filter) &&
      "sex" %in% names(groups) && "quintile" %in% names(groups)) {
    k_idx <- which(groups$sex == sex_filter &
                   groups$quintile == quintile_filter)
    if (length(k_idx) == 0) k_idx <- group_index
    grp_label <- paste0(sex_filter, ", Q", quintile_filter)
  } else {
    k_idx <- group_index
    grp_label <- paste(
      vapply(names(groups), function(g) groups[[g]][k_idx], character(1)),
      collapse = ", ")
  }

  obs_shares <- apply(shares_obs_4d[, , k_idx, , drop = FALSE], c(2, 4), mean)
  fc_shares  <- apply(shares_fc_4d[, , k_idx, , drop = FALSE], c(2, 4), mean)

  share_dt <- data.table::rbindlist(list(
    data.table::data.table(
      year  = rep(years_obs, each = K + 1),
      cause = rep(cause_names, length(years_obs)),
      share = as.vector(t(obs_shares)),
      type  = "Observed"),
    data.table::data.table(
      year  = rep(years_fc, each = K + 1),
      cause = rep(cause_names, length(years_fc)),
      share = as.vector(t(fc_shares)),
      type  = "Forecast")))
  share_dt[, cause := factor(cause, levels = rev(cause_names))]

  p <- ggplot2::ggplot(share_dt,
        ggplot2::aes(x = year, y = share, fill = cause)) +
    ggplot2::geom_area(alpha = 0.8) +
    ggplot2::geom_vline(xintercept = max(years_obs) + 0.5, linetype = "dotted") +
    ggplot2::scale_fill_brewer(palette = "Set3", name = "Cause") +
    ggplot2::labs(title = paste0("Cause-of-Death Share Evolution (",
                                  grp_label, ")"),
                  x = "Year", y = "Share of All-Cause Mortality") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right")
  ggplot2::ggsave(
    file.path(output_dir, "plot_share_evolution_coherent.png"),
    p, width = 12, height = 6, dpi = 150)
}

#' Consistency check boxplot (sum of diseases / all-cause)
#' @param fc_dt_combined Combined forecast data.table (all causes).
#' @param ac_fc_dt       All-cause forecast data.table.
#' @param output_dir     Output directory.
#' @param group_cols     Character vector of group column names.
#' @export
plot_consistency_check <- function(fc_dt_combined, ac_fc_dt, output_dir,
                                   group_cols = character(0)) {
  by_cols <- c("year", group_cols, "age_group")
  disease_sum <- fc_dt_combined[cause != "All Causes",
    .(mx_sum = sum(mx_point, na.rm = TRUE)),
    by = by_cols]
  ac_sel <- c("year", group_cols, "age_group", "mx_point")
  ac <- ac_fc_dt[, ac_sel, with = FALSE]
  data.table::setnames(ac, "mx_point", "mx_all")
  check <- merge(disease_sum, ac, by = c("year", group_cols, "age_group"))
  check[, ratio := mx_sum / mx_all]

  p <- ggplot2::ggplot(check, ggplot2::aes(x = factor(year), y = ratio)) +
    ggplot2::geom_boxplot(fill = "#66c2a5", alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::labs(title = "Coherence Check: Sum(Disease mx) / All-Cause mx",
                  x = "Year", y = "Ratio") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(
    file.path(output_dir, "plot_consistency_coherent.png"),
    p, width = 9, height = 5, dpi = 150)
}


#' Disaggregation profile diagnostic
#'
#' For a single group and forecast year, overlays grouped rates as
#' horizontal step bars, GAM raw, parametric, blended (pre-correction),
#' and final corrected curves.  A dashed vertical line marks the blend
#' centre.
#'
#' @param disagg      Return value from \code{\link{disaggregate_forecast_gam}}.
#' @param mx_grouped  Numeric vector of grouped rates for this group/year.
#' @param ages        Integer vector of age-group starts.
#' @param width       Age-group width (default 5).
#' @param fc_year_idx Integer index into forecast years (default 1).
#' @param group_idx   Integer group index (default 1).
#' @param open_mid    Blend centre age (default 95).
#' @param output_dir  Output directory.
#' @param cause_short Short cause label for filename.
#' @param group_label Character label for subtitle.
#' @export
plot_disagg_profile <- function(disagg, mx_grouped, ages, width = 5,
                                fc_year_idx = 1L, group_idx = 1L,
                                open_mid = 95, output_dir, cause_short,
                                group_label = "") {
  ta <- disagg$target_ages
  d  <- disagg$diag

  n_a <- length(ages)
  step_dt <- data.table::data.table(
    xmin = ages,
    xmax = c(ages[-1], max(ta) + 1),
    y    = log(pmax(mx_grouped, 1e-12)))

  lines_dt <- rbind(
    data.table::data.table(age = ta,
      log_mx = log(pmax(d$gam_raw[, fc_year_idx, group_idx], 1e-12)),
      stage = "1 GAM raw"),
    data.table::data.table(age = ta,
      log_mx = log(pmax(d$param_raw[, fc_year_idx, group_idx], 1e-12)),
      stage = "2 Parametric"),
    data.table::data.table(age = ta,
      log_mx = log(pmax(d$mx_pre_correction[, fc_year_idx, group_idx], 1e-12)),
      stage = "3 Blended"),
    data.table::data.table(age = ta,
      log_mx = log(pmax(disagg$mx_fc_1yr[, fc_year_idx, group_idx], 1e-12)),
      stage = "4 Final (corrected)"))

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = step_dt,
      ggplot2::aes(xmin = xmin, xmax = xmax,
                   ymin = y - 0.03, ymax = y + 0.03),
      fill = "grey70", alpha = 0.5) +
    ggplot2::geom_line(data = lines_dt,
      ggplot2::aes(x = age, y = log_mx, color = stage, linetype = stage),
      linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = open_mid, linetype = "dashed",
                        color = "grey40") +
    ggplot2::annotate("text", x = open_mid + 0.5,
      y = min(lines_dt$log_mx, na.rm = TRUE),
      label = paste0("blend centre = ", open_mid),
      hjust = 0, size = 3, color = "grey40") +
    ggplot2::scale_color_manual(
      values = c("1 GAM raw" = "#377EB8", "2 Parametric" = "#E41A1C",
                 "3 Blended" = "#FF7F00", "4 Final (corrected)" = "#4DAF4A"),
      name = "Stage") +
    ggplot2::scale_linetype_manual(
      values = c("1 GAM raw" = "solid", "2 Parametric" = "dotted",
                 "3 Blended" = "dashed", "4 Final (corrected)" = "solid"),
      name = "Stage") +
    ggplot2::labs(title = paste0("Disaggregation Profile: ", cause_short),
                  subtitle = group_label, x = "Age", y = "log(mx)") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(
    file.path(output_dir, paste0("plot_disagg_profile_", cause_short, ".png")),
    p, width = 10, height = 6, dpi = 150)
}


#' Disaggregation smoothness heatmap
#'
#' Displays \eqn{|d^2 \log(mx) / d(\text{age})^2|} across ages and
#' groups as a heatmap.  Bright spots indicate kinks or jumps.
#'
#' @param mx_1yr_dt  \code{data.table} with columns \code{age},
#'   \code{mx_1yr}, \code{sex}, \code{quintile}, \code{year}.
#' @param output_dir Output directory.
#' @param cause_short Short cause label for filename.
#' @export
plot_disagg_smoothness <- function(mx_1yr_dt, output_dir, cause_short) {
  sm_list <- list()
  for (yr in unique(mx_1yr_dt$year)) {
    for (s in unique(mx_1yr_dt$sex)) {
      for (q in unique(mx_1yr_dt$quintile)) {
        sub <- mx_1yr_dt[year == yr & sex == s & quintile == q][order(age)]
        if (nrow(sub) < 3) next
        d2 <- c(NA, abs(diff(log(pmax(sub$mx_1yr, 1e-12)),
                             differences = 2)), NA)
        sm_list[[length(sm_list) + 1]] <- data.table::data.table(
          age = sub$age, d2 = d2,
          grp = paste0(s, " Q", q, " ", yr))
      }
    }
  }
  if (length(sm_list) == 0) return(invisible(NULL))
  sm_dt <- data.table::rbindlist(sm_list)
  sm_dt <- sm_dt[!is.na(d2)]

  p <- ggplot2::ggplot(sm_dt, ggplot2::aes(x = age, y = grp,
                                            fill = pmin(d2, 1))) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "|d\u00b2log(mx)|",
                                  limits = c(0, 1), option = "inferno") +
    ggplot2::geom_vline(xintercept = 95, linetype = "dashed",
                        color = "white", linewidth = 0.3) +
    ggplot2::labs(title = paste0("Disaggregation Smoothness: ", cause_short),
                  subtitle = "Bright = kink/jump; dashed line = blend centre",
                  x = "Age", y = "") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
  n_grp <- length(unique(sm_dt$grp))
  ggplot2::ggsave(
    file.path(output_dir, paste0("plot_disagg_smoothness_", cause_short, ".png")),
    p, width = 12,
    height = min(49, max(4, n_grp * 0.25)), dpi = 150)
}


#' Temporal consistency diagnostic
#'
#' Shows mortality rates at selected ages across all years (observed +
#' forecast) for one population group.  Useful for verifying that the GAM
#' approach produces smooth temporal evolution at individual single-year
#' ages, including those near the blend boundary.
#'
#' @param disagg        Return value from
#'   \code{\link{disaggregate_forecast_gam}}.
#' @param years         Integer vector of observed years.
#' @param forecast_years Integer vector of forecast years.
#' @param group_idx     Integer group index (default 1).
#' @param output_dir    Output directory.
#' @param cause_short   Short cause label.
#' @param group_label   Character label for subtitle.
#' @param probe_ages    Ages to plot (default: 50, 70, 85, 90, 95, 100).
#' @export
plot_disagg_temporal <- function(disagg, years, forecast_years,
                                 group_idx = 1L, output_dir, cause_short,
                                 group_label = "",
                                 probe_ages = c(50, 70, 85, 90, 95, 100)) {
  ta <- disagg$target_ages
  probe_ages <- probe_ages[probe_ages %in% ta]
  n_t <- length(years); h <- length(forecast_years)

  lines <- list()
  for (a in probe_ages) {
    ai <- which(ta == a)
    obs_mx <- disagg$mx_obs_1yr[ai, , group_idx]
    fc_mx  <- disagg$mx_fc_1yr[ai, , group_idx]
    lines[[length(lines) + 1]] <- data.table::data.table(
      year = c(years, forecast_years),
      log_mx = log(pmax(c(obs_mx, fc_mx), 1e-12)),
      age = paste0("age ", a),
      type = c(rep("Observed", n_t), rep("Forecast", h)))
  }
  dt <- data.table::rbindlist(lines)

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = year, y = log_mx,
                                         color = age)) +
    ggplot2::geom_line(ggplot2::aes(linetype = type), linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = max(years) + 0.5,
                        linetype = "dotted", color = "grey50") +
    ggplot2::scale_linetype_manual(
      values = c(Observed = "solid", Forecast = "dashed")) +
    ggplot2::labs(title = paste0("Temporal Consistency: ", cause_short),
                  subtitle = group_label, x = "Year", y = "log(mx)") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(
    file.path(output_dir, paste0("plot_disagg_temporal_", cause_short, ".png")),
    p, width = 10, height = 6, dpi = 150)
}


#' Group-mean correction diagnostic
#'
#' Bar chart showing the ratio \code{mean(single-year) / grouped rate} per
#' age group, before and after correction.  Both bars should be close to
#' 1.0 after correction.
#'
#' @param disagg       Return value from
#'   \code{\link{disaggregate_forecast_gam}}.
#' @param mx_forecast  3D array \code{[n_a x h x n_g]} of grouped forecast
#'   rates.
#' @param ages         Integer vector of age-group starts.
#' @param fc_year_idx  Integer forecast year index (default 1).
#' @param group_idx    Integer group index (default 1).
#' @param output_dir   Output directory.
#' @param cause_short  Short cause label.
#' @param max_age      Maximum age (default 100).
#' @export
plot_disagg_correction <- function(disagg, mx_forecast, ages,
                                   fc_year_idx = 1L, group_idx = 1L,
                                   output_dir, cause_short, max_age = 100) {
  ta <- disagg$target_ages
  d  <- disagg$diag
  mx_grp  <- mx_forecast[, fc_year_idx, group_idx]
  mx_pre  <- d$mx_pre_correction[, fc_year_idx, group_idx]
  mx_post <- disagg$mx_fc_1yr[, fc_year_idx, group_idx]
  n_a <- length(ages)

  calc_factors <- function(mx_single) {
    f <- numeric(n_a)
    for (g in seq_along(ages)) {
      lower <- ages[g]
      upper <- if (g < n_a) ages[g + 1] - 1 else max_age
      idx <- which(ta >= lower & ta <= upper)
      if (length(idx) > 0 && mx_grp[g] > 0)
        f[g] <- mean(mx_single[idx]) / mx_grp[g]
      else f[g] <- 1
    }
    f
  }

  dt <- rbind(
    data.table::data.table(age_group = ages, ratio = calc_factors(mx_pre),
                           stage = "Before correction"),
    data.table::data.table(age_group = ages, ratio = calc_factors(mx_post),
                           stage = "After correction"))

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = factor(age_group), y = ratio,
                                         fill = stage)) +
    ggplot2::geom_col(position = "dodge", alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::scale_fill_manual(
      values = c("Before correction" = "#FC8D62",
                 "After correction" = "#66C2A5")) +
    ggplot2::labs(
      title = paste0("Group-Mean Correction: ", cause_short),
      subtitle = "Ratio = mean(single-year) / grouped rate; should be 1.0",
      x = "Age group", y = "Ratio") +
    ggplot2::theme_minimal(base_size = 12)
  ggplot2::ggsave(
    file.path(output_dir, paste0("plot_disagg_correction_", cause_short, ".png")),
    p, width = 9, height = 5, dpi = 150)
}
