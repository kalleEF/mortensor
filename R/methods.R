#' S3 methods for mort_data and mort_forecast objects
#'
#' @name methods
#' @importFrom data.table as.data.table
#' @importFrom stats median residuals
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon
#'   geom_tile facet_wrap facet_grid labs theme_minimal scale_color_brewer
#'   scale_fill_gradient2 ggsave
NULL

# ---- mort_data methods ----

#' Convert a mort_data or mort_forecast object to a long-format data.table
#'
#' Package-native alternative to `as.data.table()` that works after
#' `library(mortensor)` without needing `library(data.table)`. For users who
#' have data.table loaded, `as.data.table()` dispatches to the same methods.
#'
#' @param x A `mort_data` or `mort_forecast` object.
#' @param ... Passed to the underlying method.
#' @return A long-format `data.table`.
#' @export
to_data_table <- function(x, ...) UseMethod("to_data_table")

#' @export
to_data_table.mort_data <- function(x, ...) as.data.table.mort_data(x, ...)

#' @export
to_data_table.mort_forecast <- function(x, ...) as.data.table.mort_forecast(x, ...)

#' @describeIn to_data_table Convert a mort_data object
#' @param x A `mort_data` object.
#' @param ... Ignored.
#' @return A `data.table` with columns: year, age_group, age_start,
#'   group columns, mx, and optionally pop.
#' @export
as.data.table.mort_data <- function(x, ...) {
  rows <- vector("list", x$n_t * x$n_g)
  idx <- 0L
  for (j in seq_len(x$n_t)) {
    for (k in seq_len(x$n_g)) {
      idx <- idx + 1L
      row <- data.table::data.table(
        year      = x$years[j],
        age_group = x$age_labels,
        age_start = x$ages,
        mx        = x$mx[, j, k]
      )
      if (!is.null(x$pop))
        data.table::set(row, j = "pop", value = x$pop[, j, k])
      for (g in names(x$groups))
        data.table::set(row, j = g, value = x$groups[[g]][k])
      if (!is.null(x$cause))
        data.table::set(row, j = "cause", value = x$cause)
      rows[[idx]] <- row
    }
  }
  data.table::rbindlist(rows)
}

#' Plot a mort_data object
#'
#' @param x A `mort_data` object.
#' @param type One of `"heatmap"` (Lexis surface), `"age_profile"` (log-mx
#'   vs age for selected years), or `"trend"` (log-mx vs year for selected
#'   ages).
#' @param group Integer index of the group to plot (default 1).
#' @param years Integer vector of years to highlight (for `"age_profile"`).
#'   Defaults to 4 evenly-spaced years.
#' @param ages Integer vector of age starts to highlight (for `"trend"`).
#'   Defaults to every other age group.
#' @param ... Ignored.
#' @return A `ggplot` object (invisibly).
#' @export
plot.mort_data <- function(x, type = c("heatmap", "age_profile", "trend"),
                           group = 1L, years = NULL, ages = NULL, ...) {
  type <- match.arg(type)

  if (type == "heatmap") {
    dt <- expand.grid(age = x$ages, year = x$years)
    dt$log_mx <- as.vector(log(pmax(x$mx[, , group], 1e-12)))
    dt <- as.data.frame(dt)
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = year, y = age, fill = log_mx)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                     midpoint = median(dt$log_mx)) +
      ggplot2::labs(title = paste0("Mortality surface: ",
                                    x$cause %||% "all causes"),
                    subtitle = paste("Group", group),
                    x = "Year", y = "Age", fill = "log(mx)") +
      ggplot2::theme_minimal(base_size = 12)

  } else if (type == "age_profile") {
    if (is.null(years)) {
      idx <- round(seq(1, x$n_t, length.out = min(4, x$n_t)))
      years <- x$years[idx]
    }
    yrs_idx <- which(x$years %in% years)
    rows <- list()
    for (j in yrs_idx) {
      rows[[length(rows) + 1]] <- data.frame(
        age = x$ages, year = factor(x$years[j]),
        log_mx = log(pmax(x$mx[, j, group], 1e-12)))
    }
    dt <- do.call(rbind, rows)
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = age, y = log_mx,
                                           color = year, group = year)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::labs(title = paste0("Age profiles: ",
                                    x$cause %||% "all causes"),
                    x = "Age", y = "log(mx)", color = "Year") +
      ggplot2::theme_minimal(base_size = 12)

  } else {
    if (is.null(ages))
      ages <- x$ages[seq(1, x$n_a, by = max(1, x$n_a %/% 6))]
    ages_idx <- which(x$ages %in% ages)
    rows <- list()
    for (i in ages_idx) {
      rows[[length(rows) + 1]] <- data.frame(
        year = x$years,
        age_group = factor(x$age_labels[i]),
        log_mx = log(pmax(x$mx[i, , group], 1e-12)))
    }
    dt <- do.call(rbind, rows)
    p <- ggplot2::ggplot(dt, ggplot2::aes(x = year, y = log_mx,
                                           color = age_group)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::labs(title = paste0("Mortality trends: ",
                                    x$cause %||% "all causes"),
                    x = "Year", y = "log(mx)", color = "Age group") +
      ggplot2::theme_minimal(base_size = 12)
  }

  print(p)
  invisible(p)
}

# ---- export generic ----

#' Export mortality objects to files
#'
#' Writes a mort_forecast or mort_calibration to CSV/RDS files.
#'
#' @param x An object to export.
#' @param ... Additional arguments passed to methods.
#' @export
export <- function(x, ...) UseMethod("export")

#' @describeIn export Export a mort_forecast
#' @param dir Output directory (default `"."`).
#' @param prefix File name prefix (defaults to cause label or `"forecast"`).
#' @param format One of `"csv"` or `"rds"`.
#' @export
export.mort_forecast <- function(x, dir = ".", prefix = NULL,
                                  format = c("csv", "rds"), ...) {
  format <- match.arg(format)
  if (is.null(prefix)) prefix <- gsub("[^a-zA-Z0-9_]", "_",
                                       x$cause %||% "forecast")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  dt <- to_data_table(x)
  if (format == "csv") {
    out <- file.path(dir, paste0(prefix, "_forecast_grouped.csv"))
    data.table::fwrite(dt, out)
    message("Wrote ", out)
  } else {
    out <- file.path(dir, paste0(prefix, "_forecast_grouped.rds"))
    saveRDS(dt, out)
    message("Wrote ", out)
  }

  # Also export single-year if available
  mx_1yr <- x$extras$mx_1yr_dt
  if (!is.null(mx_1yr)) {
    if (format == "csv") {
      out2 <- file.path(dir, paste0(prefix, "_forecast_single_year.csv"))
      data.table::fwrite(mx_1yr, out2)
      message("Wrote ", out2)
    } else {
      out2 <- file.path(dir, paste0(prefix, "_forecast_single_year.rds"))
      saveRDS(mx_1yr, out2)
      message("Wrote ", out2)
    }
  }
  invisible(x)
}

#' @describeIn export Export a mort_calibration
#' @export
export.mort_calibration <- function(x, dir = ".", prefix = NULL,
                                     format = c("csv", "rds"), ...) {
  format <- match.arg(format)
  if (is.null(prefix)) prefix <- gsub("[^a-zA-Z0-9_]", "_",
                                       x$data$cause %||% "calibration")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  comp <- x$comparison
  if (!is.null(comp) && nrow(comp) > 0) {
    if (format == "csv") {
      out <- file.path(dir, paste0(prefix, "_model_comparison.csv"))
      data.table::fwrite(comp, out)
      message("Wrote ", out)
    } else {
      out <- file.path(dir, paste0(prefix, "_model_comparison.rds"))
      saveRDS(comp, out)
      message("Wrote ", out)
    }
  }
  invisible(x)
}

# ---- mort_forecast methods ----

#' Summary of a mort_forecast object
#'
#' @param object A `mort_forecast` object.
#' @param ... Ignored.
#' @return Invisibly, the object.
#' @export
summary.mort_forecast <- function(object, ...) {
  cat("<mort_forecast> summary\n")
  cat("  Model:    ", object$model, "\n")
  cat("  Cause:    ", object$cause %||% "(not set)", "\n")
  cat("  Training: ", min(object$train_years), "-",
      max(object$train_years), "\n")
  cat("  Forecast: ", min(object$forecast_years), "-",
      max(object$forecast_years), " (h =", object$h, ")\n")
  cat("  Ages:     ", object$n_a, "groups (",
      object$age_labels[1], " .. ", object$age_labels[object$n_a], ")\n")
  cat("  Groups:   ", object$n_g, "\n")

  # Point forecast range
  log_range <- range(log(pmax(object$mx_point, 1e-15)))
  cat("  log(mx) forecast range: [", round(log_range[1], 2), ",",
      round(log_range[2], 2), "]\n")

  has_pi <- !is.null(object$mx_lower) && !is.null(object$mx_upper)
  cat("  PI:       ", if (has_pi) "yes (95%)" else "no", "\n")
  cat("  Fitted:   ", if (!is.null(object$mx_fitted)) "yes" else "no", "\n")

  if (!is.null(object$fit_info$fit_pct))
    cat("  Fit pct:  ", round(object$fit_info$fit_pct, 1), "%\n")
  if (!is.null(object$fit_info$mase))
    cat("  CV MASE:  ", round(object$fit_info$mase, 4), "\n")
  if (!is.null(object$fit_info$model_comparison)) {
    cat("  Model comparison:\n")
    print(object$fit_info$model_comparison[
      order(object$fit_info$model_comparison$mase_combined), ])
  }
  invisible(object)
}

#' Extract fitted values from a mort_forecast
#'
#' @param object A `mort_forecast` object.
#' @param ... Ignored.
#' @return A 3D array `[n_ages x n_train_years x n_groups]` of fitted
#'   mortality rates, or `NULL` if not available.
#' @export
fitted.mort_forecast <- function(object, ...) {
  if (is.null(object$mx_fitted)) {
    message("No fitted values stored. Re-run forecast with return_fitted = TRUE.")
    return(invisible(NULL))
  }
  object$mx_fitted
}

#' Extract residuals from a mort_forecast
#'
#' @param object A `mort_forecast` object.
#' @param type One of `"response"` (raw scale: observed - fitted) or
#'   `"deviance"` (log scale: log(observed) - log(fitted)).
#' @param ... Ignored.
#' @return A 3D array of residuals, or `NULL` if fitted values are not
#'   available.
#' @export
residuals.mort_forecast <- function(object,
                                     type = c("response", "deviance"), ...) {
  type <- match.arg(type)
  if (is.null(object$mx_fitted)) {
    message("No fitted values stored; cannot compute residuals.")
    return(invisible(NULL))
  }
  raw <- object$mx_raw
  if (is.null(raw)) {
    message("No raw data stored; cannot compute residuals.")
    return(invisible(NULL))
  }
  if (type == "response") {
    raw - object$mx_fitted
  } else {
    log(pmax(raw, 1e-12)) - log(pmax(object$mx_fitted, 1e-12))
  }
}

#' Extract model coefficients from a mort_forecast
#'
#' For Tucker3-based models, returns the decomposition factors A (age),
#' B (time), C (group), and core tensor G. For other models, returns
#' any extras stored during fitting.
#'
#' @param object A `mort_forecast` object.
#' @param ... Ignored.
#' @return A named list of model components, or `NULL`.
#' @export
coef.mort_forecast <- function(object, ...) {
  ex <- object$extras
  if (is.null(ex) || length(ex) == 0) {
    message("No model coefficients stored.")
    return(invisible(NULL))
  }
  # For Tucker3, return the factor matrices and core tensor
  tucker_names <- c("A", "B", "C", "G_mat", "arima_fits")
  found <- intersect(names(ex), tucker_names)
  if (length(found) > 0) {
    return(ex[found])
  }
  ex
}

#' Plot a mort_calibration object
#'
#' @param x A `mort_calibration` object.
#' @param type One of `"comparison"` (bar chart of MASE by model) or
#'   `"ranks"` (Tucker3 rank selection heatmap, if available).
#' @param ... Ignored.
#' @return A `ggplot` object (invisibly).
#' @export
plot.mort_calibration <- function(x, type = c("comparison", "ranks"), ...) {
  type <- match.arg(type)
  if (type == "comparison") {
    .plot_calib_comparison(x)
  } else {
    .plot_calib_ranks(x)
  }
}

#' Plot a mort_forecast object
#'
#' @param x A `mort_forecast` object.
#' @param type One of `"trends"` (age-averaged trends with PIs),
#'   `"age"` (age profiles for forecast years), `"pi"` (PI fan chart),
#'   `"residuals"` (residual diagnostics), `"smoothing"` (observed vs
#'   smoothed), `"components"` (Tucker3 B-mode temporal components),
#'   `"model_comparison"` (bar chart from calibration), or `"disagg"`
#'   (single-year disaggregated profiles).
#' @param group Integer index or `NULL` for all groups.
#' @param ... Additional arguments passed to ggplot.
#' @return A `ggplot` object (invisibly).
#' @export
plot.mort_forecast <- function(x,
                                type = c("trends", "age", "pi", "residuals",
                                         "smoothing", "components",
                                         "model_comparison", "disagg"),
                                group = NULL, ...) {
  type <- match.arg(type)

  switch(type,
    trends           = .plot_fc_trends(x, group),
    age              = .plot_fc_age(x, group),
    pi               = .plot_fc_pi(x, group),
    residuals        = .plot_fc_residuals(x, group),
    smoothing        = .plot_fc_smoothing(x, group),
    components       = .plot_fc_components(x, group),
    model_comparison = .plot_fc_model_comparison(x),
    disagg           = .plot_fc_disagg(x, group)
  )
}

# ---- Internal plot helpers ----

#' @keywords internal
.plot_fc_trends <- function(fc, group = NULL) {
  groups_to_plot <- if (!is.null(group)) group else seq_len(fc$n_g)
  rows <- list()
  for (k in groups_to_plot) {
    mx_mean <- colMeans(fc$mx_point[, , k, drop = FALSE], dims = 1)
    row <- data.frame(
      year = fc$forecast_years, log_mx = log(pmax(mx_mean, 1e-12)),
      group = paste("Group", k), source = "forecast",
      log_lo = NA_real_, log_hi = NA_real_)
    if (!is.null(fc$mx_lower) && !is.null(fc$mx_upper)) {
      lo_mean <- colMeans(fc$mx_lower[, , k, drop = FALSE], dims = 1)
      hi_mean <- colMeans(fc$mx_upper[, , k, drop = FALSE], dims = 1)
      row$log_lo <- log(pmax(lo_mean, 1e-12))
      row$log_hi <- log(pmax(hi_mean, 1e-12))
    }
    rows[[length(rows) + 1]] <- row

    if (!is.null(fc$mx_fitted)) {
      fit_mean <- colMeans(fc$mx_fitted[, , k, drop = FALSE], dims = 1)
      rows[[length(rows) + 1]] <- data.frame(
        year = fc$train_years, log_mx = log(pmax(fit_mean, 1e-12)),
        group = paste("Group", k), source = "fitted",
        log_lo = NA_real_, log_hi = NA_real_)
    }
  }
  dt <- do.call(rbind, rows)

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = year, y = log_mx,
                                         color = source)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(title = paste0("Forecast trends: ", fc$model),
                  subtitle = fc$cause %||% "",
                  x = "Year", y = "log(mx)") +
    ggplot2::theme_minimal(base_size = 12)

  if (any(!is.na(dt$log_lo))) {
    fc_dt <- dt[dt$source == "forecast", ]
    p <- p + ggplot2::geom_ribbon(
      data = fc_dt,
      ggplot2::aes(ymin = log_lo, ymax = log_hi),
      alpha = 0.2, color = NA, fill = "steelblue")
  }

  if (length(groups_to_plot) > 1)
    p <- p + ggplot2::facet_wrap(~group)

  print(p)
  invisible(p)
}

#' @keywords internal
.plot_fc_age <- function(fc, group = NULL) {
  k <- if (!is.null(group)) group[1] else 1L
  yrs_idx <- c(1, max(1, fc$h %/% 2), fc$h)
  yrs_idx <- unique(pmin(yrs_idx, fc$h))

  rows <- list()
  for (j in yrs_idx) {
    rows[[length(rows) + 1]] <- data.frame(
      age = fc$ages, year = factor(fc$forecast_years[j]),
      log_mx = log(pmax(fc$mx_point[, j, k], 1e-12)))
  }
  dt <- do.call(rbind, rows)

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = age, y = log_mx,
                                         color = year)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::labs(title = paste0("Forecast age profiles: ", fc$model),
                  subtitle = paste("Group", k),
                  x = "Age", y = "log(mx)", color = "Year") +
    ggplot2::theme_minimal(base_size = 12)

  print(p)
  invisible(p)
}

#' @keywords internal
.plot_fc_pi <- function(fc, group = NULL) {
  if (is.null(fc$mx_lower) || is.null(fc$mx_upper)) {
    message("No prediction intervals available. Use add_pi() first.")
    return(invisible(NULL))
  }
  k <- if (!is.null(group)) group[1] else 1L
  age_idx <- c(1, max(1, fc$n_a %/% 2), fc$n_a)
  age_idx <- unique(pmin(age_idx, fc$n_a))

  rows <- list()
  for (i in age_idx) {
    rows[[length(rows) + 1]] <- data.frame(
      year = fc$forecast_years, age_group = fc$age_labels[i],
      log_mx = log(pmax(fc$mx_point[i, , k], 1e-12)),
      log_lo = log(pmax(fc$mx_lower[i, , k], 1e-12)),
      log_hi = log(pmax(fc$mx_upper[i, , k], 1e-12)))
  }
  dt <- do.call(rbind, rows)

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = year, y = log_mx)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = log_lo, ymax = log_hi),
                          alpha = 0.25, fill = "steelblue") +
    ggplot2::geom_line(linewidth = 0.8, color = "darkblue") +
    ggplot2::facet_wrap(~age_group, scales = "free_y") +
    ggplot2::labs(title = paste0("Prediction intervals: ", fc$model),
                  subtitle = paste("Group", k, "- 95% PI"),
                  x = "Year", y = "log(mx)") +
    ggplot2::theme_minimal(base_size = 12)

  print(p)
  invisible(p)
}

#' @keywords internal
.plot_fc_residuals <- function(fc, group = NULL) {
  resid <- residuals(fc, type = "deviance")
  if (is.null(resid)) return(invisible(NULL))

  k <- if (!is.null(group)) group[1] else 1L
  dt <- expand.grid(age = fc$ages, year = fc$train_years)
  dt$residual <- as.vector(resid[, , k])

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = year, y = age, fill = residual)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                   midpoint = 0) +
    ggplot2::labs(title = paste0("Deviance residuals: ", fc$model),
                  subtitle = paste("Group", k),
                  x = "Year", y = "Age", fill = "Residual\n(log scale)") +
    ggplot2::theme_minimal(base_size = 12)

  print(p)
  invisible(p)
}

# ---- mort_calibration plot helpers ----

#' @keywords internal
.plot_calib_comparison <- function(cal) {
  comp <- cal$comparison
  if (is.null(comp) || nrow(comp) == 0) {
    message("No model comparison data available.")
    return(invisible(NULL))
  }
  comp <- as.data.frame(comp)
  if (!"mase_combined" %in% names(comp)) {
    message("comparison table missing 'mase_combined' column.")
    return(invisible(NULL))
  }
  comp <- comp[order(comp$mase_combined), ]
  comp$model <- factor(comp$model, levels = comp$model)
  best <- comp$model[1]

  p <- ggplot2::ggplot(comp, ggplot2::aes(x = model, y = mase_combined)) +
    ggplot2::geom_col(fill = ifelse(comp$model == best, "steelblue", "grey70")) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = paste0("Model comparison: ",
                                  cal$data$cause %||% "all causes"),
                  subtitle = paste("Best:", best),
                  x = NULL, y = "MASE (combined)") +
    ggplot2::theme_minimal(base_size = 12)
  print(p)
  invisible(p)
}

#' @keywords internal
.plot_calib_ranks <- function(cal) {
  R_A <- R_B <- mase <- NULL  # NSE bindings
  calib <- cal$calib_coda %||% cal$calib_direct
  if (is.null(calib) || is.null(calib$rank_search)) {
    message("No rank search results available for rank plot.")
    return(invisible(NULL))
  }
  rs <- as.data.frame(calib$rank_search)
  if (!all(c("R_A", "R_B", "mase") %in% names(rs))) {
    message("rank_search missing expected columns.")
    return(invisible(NULL))
  }

  p <- ggplot2::ggplot(rs, ggplot2::aes(x = factor(R_A), y = factor(R_B),
                                          fill = mase)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "steelblue", high = "red") +
    ggplot2::labs(title = "Tucker3 rank selection",
                  x = "R_A (age rank)", y = "R_B (time rank)",
                  fill = "MASE") +
    ggplot2::theme_minimal(base_size = 12)
  print(p)
  invisible(p)
}

# ---- Additional mort_forecast plot helpers ----

#' @keywords internal
.plot_fc_smoothing <- function(fc, group = NULL) {
  mx_raw <- fc$mx_raw
  mx_smooth <- fc$extras$mx_smooth
  if (is.null(mx_raw) || is.null(mx_smooth)) {
    message("No raw/smoothed data available for smoothing plot. ",
            "Ensure mx_raw and extras$mx_smooth are stored.")
    return(invisible(NULL))
  }
  k <- if (!is.null(group)) group[1] else 1L
  raw_vec <- as.vector(log(pmax(mx_raw[, , k], 1e-12)))
  smo_vec <- as.vector(log(pmax(mx_smooth[, , k], 1e-12)))
  dt <- data.frame(raw = raw_vec, smooth = smo_vec)

  lims <- range(c(dt$raw, dt$smooth), na.rm = TRUE)
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = raw, y = smooth)) +
    ggplot2::geom_point(alpha = 0.3, size = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                          color = "red") +
    ggplot2::coord_equal(xlim = lims, ylim = lims) +
    ggplot2::labs(title = paste0("Observed vs smoothed: ", fc$cause %||% ""),
                  subtitle = paste("Group", k),
                  x = "log(mx) observed", y = "log(mx) smoothed") +
    ggplot2::theme_minimal(base_size = 12)
  print(p)
  invisible(p)
}

#' @keywords internal
.plot_fc_components <- function(fc, group = NULL) {
  B <- fc$extras$B
  arima_fits <- fc$extras$arima_fits
  if (is.null(B)) {
    message("No Tucker3 B-mode matrix available. ",
            "Ensure forecast used a Tucker3 model with return_full = TRUE.")
    return(invisible(NULL))
  }
  R_B <- ncol(B)
  train_years <- fc$train_years

  rows <- list()
  for (r in seq_len(R_B)) {
    rows[[length(rows) + 1]] <- data.frame(
      year = train_years, value = B[, r],
      component = paste("Component", r), type = "observed")
  }

  if (!is.null(arima_fits) && length(arima_fits) >= R_B) {
    for (r in seq_len(R_B)) {
      if (!is.null(arima_fits[[r]])) {
        fc_obj <- tryCatch(
          forecast::forecast(arima_fits[[r]], h = fc$h),
          error = function(e) NULL)
        if (!is.null(fc_obj)) {
          rows[[length(rows) + 1]] <- data.frame(
            year = fc$forecast_years, value = as.numeric(fc_obj$mean),
            component = paste("Component", r), type = "forecast")
        }
      }
    }
  }

  dt <- do.call(rbind, rows)
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = year, y = value,
                                          color = type, linetype = type)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~component, scales = "free_y") +
    ggplot2::labs(title = paste0("Tucker3 B-mode components: ", fc$model),
                  x = "Year", y = "Component value") +
    ggplot2::scale_linetype_manual(values = c(observed = "solid",
                                               forecast = "dashed")) +
    ggplot2::theme_minimal(base_size = 12)
  print(p)
  invisible(p)
}

#' @keywords internal
.plot_fc_model_comparison <- function(fc) {
  comp <- fc$fit_info$model_comparison
  if (is.null(comp) || nrow(comp) == 0) {
    message("No model comparison data in fit_info.")
    return(invisible(NULL))
  }
  comp <- as.data.frame(comp)
  comp <- comp[order(comp$mase_combined), ]
  comp$model <- factor(comp$model, levels = comp$model)

  p <- ggplot2::ggplot(comp, ggplot2::aes(x = model, y = mase_combined)) +
    ggplot2::geom_col(fill = ifelse(comp$model == comp$model[1],
                                     "steelblue", "grey70")) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = paste0("Model comparison: ", fc$cause %||% ""),
                  subtitle = paste("Selected:", fc$model),
                  x = NULL, y = "MASE (combined)") +
    ggplot2::theme_minimal(base_size = 12)
  print(p)
  invisible(p)
}

#' @keywords internal
.plot_fc_disagg <- function(fc, group = NULL) {
  mx_1yr_dt <- fc$extras$mx_1yr_dt
  if (is.null(mx_1yr_dt)) {
    message("No single-year disaggregated data available. ",
            "Use disaggregate() first or run the pipeline.")
    return(invisible(NULL))
  }
  k <- if (!is.null(group)) group[1] else 1L

  dt <- as.data.frame(mx_1yr_dt)
  grp_cols <- intersect(names(fc$groups), names(dt))
  if (length(grp_cols) > 0) {
    for (g in grp_cols)
      dt <- dt[dt[[g]] == fc$groups[[g]][k], , drop = FALSE]
  }

  all_yrs <- sort(unique(dt$year))
  show_yrs <- all_yrs[unique(c(1, max(1, length(all_yrs) %/% 2),
                                length(all_yrs)))]
  dt <- dt[dt$year %in% show_yrs, ]
  dt$year <- factor(dt$year)

  mx_col <- if ("mx_1yr" %in% names(dt)) "mx_1yr" else "mx_point"
  dt$log_mx <- log(pmax(dt[[mx_col]], 1e-12))
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = age, y = log_mx,
                                          color = year)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(title = paste0("Single-year age profiles: ",
                                  fc$cause %||% ""),
                  subtitle = paste("Group", k),
                  x = "Age", y = "log(mx)", color = "Year") +
    ggplot2::theme_minimal(base_size = 12)
  print(p)
  invisible(p)
}
