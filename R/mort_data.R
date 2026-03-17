# ---- S3 class: mort_data ----

#' Create a mortality data object
#'
#' Constructs a `mort_data` S3 object that encodes stratified mortality data
#' as a 3D array with metadata about ages, years, and groups (strata).
#'
#' @param mx 3D array `[n_ages x n_years x n_groups]` of mortality rates.
#' @param ages Integer vector of age-group start values (e.g. \code{c(30, 35, 40, 90)}).
#' @param years Integer vector of calendar years.
#' @param groups A data.frame with one row per group (slice of 3rd dimension).
#'   Must contain columns that define strata (e.g. `sex`, `deprivation`).
#' @param pop Optional 3D array of population exposure (same dims as `mx`).
#' @param cause Optional character label for the cause of death.
#' @param age_labels Optional character vector of age-group labels
#'   (e.g. \code{"30-34"}, \code{"35-39"}). Derived from \code{ages} if omitted.
#' @param smoothed Logical; has the data been smoothed?
#' @return A `mort_data` object (list with class attribute).
#' @export
mort_data <- function(mx, ages, years, groups,
                      pop = NULL, cause = NULL,
                      age_labels = NULL, smoothed = FALSE) {
  stopifnot(is.array(mx), length(dim(mx)) == 3)
  n_a <- dim(mx)[1]; n_t <- dim(mx)[2]; n_g <- dim(mx)[3]
  stopifnot(length(ages) == n_a)
  stopifnot(length(years) == n_t)
  if (!is.data.frame(groups)) groups <- as.data.frame(groups)
  stopifnot(nrow(groups) == n_g)
  if (!is.null(pop)) {
    stopifnot(is.array(pop), length(dim(pop)) == 3)
    stopifnot(identical(dim(pop), dim(mx)))
  }
  if (is.null(age_labels)) {
    age_labels <- paste0(ages, "-", ages + 4L)
    age_labels[length(age_labels)] <- paste0(ages[length(ages)], "+")
  }
  stopifnot(length(age_labels) == n_a)

  structure(
    list(
      mx         = mx,
      pop        = pop,
      ages       = as.integer(ages),
      years      = as.integer(years),
      groups     = groups,
      cause      = cause,
      age_labels = age_labels,
      smoothed   = smoothed,
      n_a        = n_a,
      n_t        = n_t,
      n_g        = n_g
    ),
    class = "mort_data"
  )
}

#' @export
print.mort_data <- function(x, ...) {
  cat("<mort_data>\n")
  cat("  Cause:  ", x$cause %||% "(not set)", "\n")
  cat("  Ages:   ", x$n_a, "groups (",
      x$age_labels[1], " .. ", x$age_labels[x$n_a], ")\n")
  cat("  Years:  ", min(x$years), "-", max(x$years),
      " (", x$n_t, " years)\n")
  cat("  Groups: ", x$n_g, " (", paste(names(x$groups), collapse = " x "), ")\n")
  cat("  Pop:    ", if (!is.null(x$pop)) "yes" else "no", "\n")
  cat("  Smoothed:", x$smoothed, "\n")
  invisible(x)
}

#' @export
summary.mort_data <- function(object, ...) {
  cat("<mort_data> summary\n")
  cat("  Cause:     ", object$cause %||% "(not set)", "\n")
  cat("  Dimensions:", object$n_a, "ages x", object$n_t, "years x",
      object$n_g, "groups\n")
  cat("  log(mx) range: [",
      round(min(log(pmax(object$mx, 1e-15))), 2), ", ",
      round(max(log(pmax(object$mx, 1e-15))), 2), "]\n")
  cat("  mx range:  [",
      formatC(min(object$mx, na.rm = TRUE), format = "e", digits = 2), ", ",
      formatC(max(object$mx, na.rm = TRUE), format = "e", digits = 2), "]\n")
  if (!is.null(object$pop)) {
    cat("  Pop range: [",
        formatC(min(object$pop, na.rm = TRUE), format = "g"),  ", ",
        formatC(max(object$pop, na.rm = TRUE), format = "g"), "]\n")
  }
  cat("  Strata:\n")
  for (nm in names(object$groups)) {
    cat("    ", nm, ": ", paste(unique(object$groups[[nm]]), collapse = ", "), "\n")
  }
  invisible(object)
}

#' Subset a mort_data object by group indices
#'
#' @param x A `mort_data` object.
#' @param i Age indices (integer or logical).
#' @param j Year indices (integer or logical).
#' @param k Group indices (integer or logical).
#' @param ... Ignored.
#' @return A new `mort_data` object.
#' @export
`[.mort_data` <- function(x, i = NULL, j = NULL, k = NULL, ...) {
  mx  <- x$mx
  pop <- x$pop
  ages <- x$ages; years <- x$years; groups <- x$groups
  age_labels <- x$age_labels

  if (!is.null(i)) {
    mx  <- mx[i, , , drop = FALSE]
    if (!is.null(pop)) pop <- pop[i, , , drop = FALSE]
    ages <- ages[i]; age_labels <- age_labels[i]
  }
  if (!is.null(j)) {
    mx  <- mx[, j, , drop = FALSE]
    if (!is.null(pop)) pop <- pop[, j, , drop = FALSE]
    years <- years[j]
  }
  if (!is.null(k)) {
    mx  <- mx[, , k, drop = FALSE]
    if (!is.null(pop)) pop <- pop[, , k, drop = FALSE]
    groups <- groups[k, , drop = FALSE]
  }

  mort_data(mx = mx, ages = ages, years = years, groups = groups,
            pop = pop, cause = x$cause, age_labels = age_labels,
            smoothed = x$smoothed)
}

# ---- S3 class: mort_forecast ----

#' Create a mortality forecast object
#'
#' @param mx_point 3D array `[n_ages x h x n_groups]` of point forecasts.
#' @param train A `mort_data` object used for training (provides metadata).
#' @param forecast_years Integer vector of forecast years.
#' @param model Character string naming the forecast model.
#' @param mx_median Optional 3D array of bootstrap median forecasts.
#' @param mx_lower Optional 3D array of lower PI bound.
#' @param mx_upper Optional 3D array of upper PI bound.
#' @param sims Optional 4D array `[n_ages x h x n_groups x n_sim]` of
#'   bootstrap simulations.
#' @param mx_fitted Optional 3D array `[n_ages x n_train_years x n_groups]`
#'   of in-sample fitted values.
#' @param mx_raw Optional 3D array of raw (unsmoothed) training mortality
#'   rates, used for residual computation.
#' @param fit_info Optional named list of model fit information, e.g.
#'   `fit_pct`, `mase`, `model_comparison`.
#' @param extras Named list of additional model-specific outputs (Tucker
#'   decomposition components, ARIMA fits, etc.).
#' @return A `mort_forecast` object.
#' @export
mort_forecast <- function(mx_point, train, forecast_years, model,
                          mx_median = NULL, mx_lower = NULL,
                          mx_upper = NULL, sims = NULL,
                          mx_fitted = NULL, mx_raw = NULL,
                          fit_info = list(), extras = list()) {
  stopifnot(inherits(train, "mort_data"))
  stopifnot(is.array(mx_point), length(dim(mx_point)) == 3)
  h <- dim(mx_point)[2]
  stopifnot(length(forecast_years) == h)

  structure(
    list(
      mx_point       = mx_point,
      mx_median      = mx_median,
      mx_lower       = mx_lower,
      mx_upper       = mx_upper,
      sims           = sims,
      mx_fitted      = mx_fitted,
      mx_raw         = mx_raw,
      model          = model,
      forecast_years = as.integer(forecast_years),
      ages           = train$ages,
      age_labels     = train$age_labels,
      groups         = train$groups,
      cause          = train$cause,
      train_years    = train$years,
      n_a            = train$n_a,
      h              = h,
      n_g            = train$n_g,
      fit_info       = fit_info,
      extras         = extras
    ),
    class = "mort_forecast"
  )
}

#' @export
print.mort_forecast <- function(x, ...) {
  cat("<mort_forecast>\n")
  cat("  Model:  ", x$model, "\n")
  cat("  Cause:  ", x$cause %||% "(not set)", "\n")
  cat("  Train:  ", min(x$train_years), "-", max(x$train_years),
      " (", length(x$train_years), " years)\n")
  cat("  Horizon:", x$h, "years (",
      min(x$forecast_years), "-", max(x$forecast_years), ")\n")
  cat("  Ages:   ", x$n_a, "groups (",
      x$age_labels[1], " .. ", x$age_labels[x$n_a], ")\n")
  cat("  Groups: ", x$n_g, " (",
      paste(names(x$groups), collapse = " x "), ")\n")

  has_pi <- !is.null(x$mx_lower) && !is.null(x$mx_upper)
  cat("  PI:     ", if (has_pi) "yes (95%)" else "no", "\n")
  cat("  Fitted: ", if (!is.null(x$mx_fitted)) "yes" else "no", "\n")
  if (!is.null(x$sims))
    cat("  Sims:   ", dim(x$sims)[4], "\n")
  if (!is.null(x$fit_info$fit_pct))
    cat("  Fit:    ", round(x$fit_info$fit_pct, 1), "%\n")
  if (!is.null(x$fit_info$mase))
    cat("  MASE:   ", round(x$fit_info$mase, 4), "\n")
  invisible(x)
}

#' Convert mort_forecast to a data.table
#'
#' @param x A `mort_forecast` object.
#' @param keep.rownames Ignored (required for generic compatibility).
#' @param keep_sims Logical; include simulation array? Default FALSE.
#' @param ... Ignored.
#' @return A `data.table` with columns: year, age_group, age_start, group
#'   columns, mx_point, and optionally mx_median, mx_lo95, mx_hi95.
#' @export
as.data.table.mort_forecast <- function(x, keep.rownames = FALSE, ..., keep_sims = FALSE) {
  rows <- list()
  for (j in seq_len(x$h)) {
    for (k in seq_len(x$n_g)) {
      row <- data.table::data.table(
        cause     = x$cause %||% NA_character_,
        year      = x$forecast_years[j],
        age_group = x$age_labels,
        age_start = x$ages,
        mx_point  = x$mx_point[, j, k]
      )
      # Append group columns
      for (g in names(x$groups))
        data.table::set(row, j = g, value = x$groups[[g]][k])

      if (!is.null(x$mx_median))
        data.table::set(row, j = "mx_median", value = x$mx_median[, j, k])
      if (!is.null(x$mx_lower))
        data.table::set(row, j = "mx_lo95", value = x$mx_lower[, j, k])
      if (!is.null(x$mx_upper))
        data.table::set(row, j = "mx_hi95", value = x$mx_upper[, j, k])
      rows[[length(rows) + 1L]] <- row
    }
  }
  data.table::rbindlist(rows)
}

# ---- S3 class: mort_config ----

#' Create a validated configuration object
#'
#' @param cause_cols Character vector of column names to sum for this cause.
#'   Set to `NULL` when providing data via [mort_data] objects directly
#'   (i.e., not using the ONS loading path).
#' @param cause_label Human-readable cause label.
#' @param cause_short Short tag for filenames (must be unique across configs).
#' @param is_cause_specific Logical; if TRUE, cap against all-causes is applied.
#' @param ... Additional named parameters that override defaults. Common
#'   overrides include `forecast_h`, `n_sim`, `smoothing_method`,
#'   `disagg_method`, `pi_method`, `covid_years`, `holdout_train`,
#'   `holdout_test`, `output_dir`.
#' @return A `mort_config` object.
#' @export
mort_config <- function(cause_cols = NULL, cause_label = "All causes",
                        cause_short = "all",
                        is_cause_specific = FALSE, ...) {
  defaults <- list(
    age_groups  = c("30-34","35-39","40-44","45-49","50-54","55-59",
                    "60-64","65-69","70-74","75-79","80-84","85-89","90+"),
    age_starts  = c(30L, 35L, 40L, 45L, 50L, 55L, 60L, 65L, 70L, 75L,
                    80L, 85L, 90L),
    age_width   = 5L,
    forecast_h  = 10L,
    n_sim       = 500L,
    default_R_A = 3L,
    default_R_B = 3L,
    default_R_C = 2L,
    default_ts  = "ets",
    sex_constraint       = FALSE,
    ensemble_mode        = TRUE,
    ensemble_top_k       = 3L,
    target_ages          = 30L:100L,
    seed                 = 42L,
    smoothing_method     = "gam",
    cv_n_origins         = 4L,
    calibrate_models     = TRUE,
    disagg_method        = "natural_gompertz",
    pi_method            = "bootstrap",
    output_max_age       = 100L,
    open_interval_midpoint = 95,
    extrapolation_method = "kannisto",
    blend_width          = 2,
    covid_years          = c(2020L, 2021L, 2022L),
    holdout_train        = NULL,
    holdout_test         = NULL
  )

  overrides <- list(...)
  cfg <- c(
    list(
      cause_cols        = cause_cols,
      cause_label       = cause_label,
      cause_short       = cause_short,
      is_cause_specific = is_cause_specific
    ),
    defaults
  )
  for (nm in names(overrides)) cfg[[nm]] <- overrides[[nm]]

  structure(cfg, class = "mort_config")
}

#' @export
print.mort_config <- function(x, ...) {
  cat("<mort_config>\n")
  cat("  Cause:   ", x$cause_label, " [", x$cause_short, "]\n")
  if (!is.null(x$cause_cols))
    cat("  Cols:    ", paste(x$cause_cols, collapse = " + "), "\n")
  cat("  Ages:    ", length(x$age_groups), "groups\n")
  cat("  h:       ", x$forecast_h, "  n_sim:", x$n_sim, "\n")
  cat("  Smooth:  ", x$smoothing_method, "\n")
  cat("  Disagg:  ", x$disagg_method, "\n")
  cat("  PI:      ", x$pi_method, "\n")
  invisible(x)
}

# ---- S3 class: mort_shared_config ----

#' Create a shared configuration for multi-disease pipelines
#'
#' Provides common parameters used by [run_reconciled_forecast()] and
#' [run_coherent_forecast()].
#'
#' @param output_dir Output directory for CSV and plot files.
#' @param forecast_h Forecast horizon in years.
#' @param n_sim Number of bootstrap simulations.
#' @param age_width Width of closed age groups.
#' @param seed Random seed.
#' @param target_ages Integer vector of single-year target ages for
#'   disaggregation.
#' @param sex_constraint Logical; enforce Female <= Male?
#' @param ... Additional named parameters (e.g. `holdout_train`,
#'   `holdout_test`, `cv_years`, `cv_h`, `cv_min_train`, `age_groups`).
#' @return A `mort_shared_config` object (named list with class).
#' @export
mort_shared_config <- function(output_dir = ".", forecast_h = 10L,
                                n_sim = 500L, age_width = 5L,
                                seed = 42L, target_ages = 30L:100L,
                                sex_constraint = FALSE, ...) {
  cfg <- c(
    list(
      output_dir     = output_dir,
      forecast_h     = as.integer(forecast_h),
      n_sim          = as.integer(n_sim),
      age_width      = as.integer(age_width),
      seed           = as.integer(seed),
      target_ages    = as.integer(target_ages),
      sex_constraint = sex_constraint
    ),
    list(...)
  )
  structure(cfg, class = "mort_shared_config")
}

#' @export
print.mort_shared_config <- function(x, ...) {
  cat("<mort_shared_config>\n")
  cat("  Output: ", x$output_dir, "\n")
  cat("  h:      ", x$forecast_h, "  n_sim:", x$n_sim, "\n")
  cat("  Ages:   ", min(x$target_ages), "-", max(x$target_ages), "\n")
  invisible(x)
}
