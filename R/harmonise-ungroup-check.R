#' Diagnostics of ungrouped data
#'
#' Checks of single-year (grid) data estimated from grouped observations by
#' [mt_ungroup()], for one method or several side by side.
#'
#' * `mt_plot_ungroup(type = "profile")`: the observed group rates (grey
#'   steps) and the ungrouped rates (lines) over age in one period. The
#'   ungrouped rates should pass through the steps and change smoothly.
#' * `type = "time"`: the ungrouped rates of a few ages over the periods
#'   (log scale), with the known rates of `reference` as points or, without
#'   `reference`, the rate of the observed group that contains each age
#'   (grey points). One-dimensional ungrouping fits every year separately,
#'   so single ages can jump from year to year; two-dimensional ungrouping
#'   smooths over time as well.
#' * `type = "error"`: where the single years are known (`reference`, for
#'   example data that were grouped for the check), the log ratio of the
#'   ungrouped to the known rates by age, one panel per method: the median
#'   over periods and strata (line) and the range of the middle 80%
#'   (coloured band). The known rates contain Poisson noise, so even perfect
#'   estimates scatter around zero; the grey band shows the middle 80% of
#'   that noise (plus or minus 1.28 standard errors of a log rate,
#'   `1 / sqrt(deaths)`). Errors within the grey band cannot be told apart
#'   from noise.
#' * `mt_ungroup_check()` summarises the fits: the share of fits whose
#'   smoothing parameter is at the boundary of its search range
#'   (`at_bound`), the largest relative gap between fitted and observed group
#'   totals before rescaling (`max_rel_gap`) and, with `reference`, by age
#'   band the median absolute log error (`mae_log`), the same for the
#'   Poisson noise of the known rates (`noise_log`, the error that perfect
#'   estimates would have) and the median signed log error (`bias_log`).
#'
#' @param data The grouped [mt_data] that were ungrouped.
#' @param ungrouped A result of [mt_ungroup()], or a named list of results
#'   (one per method or setting).
#' @param type `"profile"`, `"time"` or `"error"`.
#' @param period Period of the profile (default: the last).
#' @param strata Stratum labels to show (default: all).
#' @param ages Ages shown with `type = "time"` (default: four ages spread
#'   over the grid).
#' @param reference An [mt_data] object with the known deaths and exposures
#'   on the grid of `ungrouped` (for `type = "error"` and the errors of
#'   `mt_ungroup_check()`).
#' @param bands Lower bounds of the age bands of `mt_ungroup_check()`
#'   (default: 0, 15, 45, 65 and 85, as far as the grid covers them).
#' @return A ggplot object; `mt_ungroup_check()` a `data.table` with one
#'   row per method (and age band, with `reference`).
#' @seealso The vignette "Ungrouping age groups: methods and diagnostics".
#' @examples
#' x <- mt_example(years = 2010:2014)
#' u <- list(PCLM = mt_ungroup(x), "Hyman spline" = mt_ungroup(x, method = "hyman_kannisto"))
#' mt_ungroup_check(u)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   mt_plot_ungroup(x, u, strata = "female")
#'   mt_plot_ungroup(x, u, type = "time", strata = "female")
#' }
#' @name ungroup-diagnostics
NULL

#' @rdname ungroup-diagnostics
#' @export
mt_plot_ungroup <- function(data, ungrouped, period = NULL, strata = NULL, type = c("profile", "time", "error"),
                            ages = NULL, reference = NULL) {
  need_ggplot()
  type <- match.arg(type)
  check_class(data, "mt_data", "data")
  single <- inherits(ungrouped, "mt_data")
  us <- ungroup_list(ungrouped)
  U <- ungroup_long(us)
  if (!is.null(strata)) U <- U[label %in% strata]
  if (!nrow(U)) mt_abort("Nothing to show: check `strata`.", class = "argument")
  col_scale <- list(scale_colour_mortensor(type = "categorical"), scale_fill_mortensor(type = "categorical"))
  if (type == "error") {
    if (is.null(reference)) mt_abort("`type = \"error\"` needs the known rates in `reference`.", class = "argument")
    err <- ungroup_errors(U, reference)
    s <- err[, list(med = stats::median(lr), lo = stats::quantile(lr, 0.1), hi = stats::quantile(lr, 0.9),
                    v = stats::qnorm(0.9) * stats::median(noise)), by = c("method", "age_lo")]
    return(ggplot2::ggplot(s, ggplot2::aes(x = age_lo + 0.5)) +
             ggplot2::geom_ribbon(ggplot2::aes(ymin = -v, ymax = v), fill = "grey85") +
             ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
             ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = method), alpha = 0.3) +
             ggplot2::geom_line(ggplot2::aes(y = med, colour = method)) +
             ggplot2::facet_wrap(~method) +
             col_scale +
             ggplot2::labs(x = "Age", y = "log(ungrouped / known rate)", colour = NULL, fill = NULL,
                           subtitle = "Median over periods and strata; band: middle 80%; grey: Poisson noise") +
             theme_mortensor(legend = "none"))
  }
  periods <- if (type == "profile") period %||% max(U$period) else sort(unique(U$period))
  G <- group_rates(data, us[[1]], periods)
  if (!is.null(strata)) G <- G[label %in% strata]
  if (type == "profile") {
    U <- U[period %in% periods]
    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = G, ggplot2::aes(x = age_lo, xend = x_hi, y = rate, yend = rate), colour = "grey40")
    p <- if (single) {
      p + ggplot2::geom_line(data = U, ggplot2::aes(x = age_lo + 0.5, y = rate), colour = "firebrick")
    } else {
      p + ggplot2::geom_line(data = U, ggplot2::aes(x = age_lo + 0.5, y = rate, colour = method)) + col_scale
    }
    title <- if (single) sprintf("Grouped (grey) and ungrouped (red) rates, %s", periods) else
      sprintf("Grouped (grey) and ungrouped rates, %s", periods)
    return(p + ggplot2::facet_wrap(~label) + ggplot2::scale_y_log10() +
             ggplot2::labs(x = "Age", y = "Death rate (log scale)", title = title, colour = NULL) +
             theme_mortensor())
  }
  grid_ages <- sort(unique(U$age_lo))
  ages <- ages %||% grid_ages[unique(round(stats::quantile(seq_along(grid_ages), c(0.2, 0.45, 0.7, 0.95))))]
  if (!all(ages %in% grid_ages)) mt_abort("`ages` must be ages of the grid.", class = "argument")
  U <- U[age_lo %in% ages]
  U[, panel := factor(paste("age", age_lo), paste("age", ages))]
  if (is.null(reference)) {
    Pt <- rbindlist(lapply(ages, function(a) G[age_lo <= a & age_hi > a][, panel := paste("age", a)]))
    pt_col <- "grey50"
    what <- "grey points: rate of the observed age group"
  } else {
    check_class(reference, "mt_data", "reference")
    Pt <- grid_rates_long(reference)[age_lo %in% ages & label %in% U$label & period %in% U$period & rate > 0]
    Pt[, panel := paste("age", age_lo)]
    pt_col <- "black"
    what <- "points: known rates"
  }
  Pt[, panel := factor(panel, levels(U$panel))]
  ggplot2::ggplot(U, ggplot2::aes(x = period, y = rate)) +
    ggplot2::geom_point(data = Pt, colour = pt_col, size = 0.8) +
    ggplot2::geom_line(ggplot2::aes(colour = method)) +
    ggplot2::facet_grid(label ~ panel, scales = "free_y") +
    ggplot2::scale_y_log10() +
    col_scale +
    ggplot2::labs(x = NULL, y = "Death rate (log scale)", colour = NULL,
                  subtitle = sprintf("Lines: ungrouped rates; %s", what)) +
    theme_mortensor()
}

#' @rdname ungroup-diagnostics
#' @export
mt_ungroup_check <- function(ungrouped, reference = NULL, bands = NULL) {
  us <- ungroup_list(ungrouped)
  fits <- rbindlist(lapply(names(us), function(n) {
    a <- as.data.table(attr(us[[n]], "ungroup"))
    if ("measure" %in% names(a) && any(a$measure == "deaths")) a <- a[measure == "deaths"]
    data.table(method = n,
               at_bound = if ("at_bound" %in% names(a) && nrow(a)) mean(a$at_bound) else NA_real_,
               max_rel_gap = if ("max_rel_diff" %in% names(a) && nrow(a)) max(a$max_rel_diff) else NA_real_)
  }))
  fits[, method := factor(method, names(us))]
  if (is.null(reference)) return(fits[])
  err <- ungroup_errors(ungroup_long(us), reference)
  lows <- sort(unique(err$age_lo))
  bands <- bands %||% c(0, 15, 45, 65, 85)
  bands <- sort(unique(c(max(bands[bands <= min(lows)], min(lows)), bands[bands > min(lows) & bands <= max(lows)])))
  hi <- c(bands[-1], Inf)
  err[, band := factor(findInterval(age_lo, bands), seq_along(bands),
                       ifelse(is.finite(hi), sprintf("%g-%g", bands, hi - 1), sprintf("%g+", bands)))]
  out <- err[, list(mae_log = stats::median(abs(lr)), noise_log = stats::qnorm(0.75) * stats::median(noise),
                    bias_log = stats::median(lr)), by = c("method", "band")]
  out <- merge(out, fits, by = "method", sort = FALSE)
  setorderv(out, c("method", "band"))
  num <- c("mae_log", "noise_log", "bias_log", "max_rel_gap")
  out[, (num) := lapply(.SD, signif, 3), .SDcols = num]
  out[]
}

# A single ungrouped data set or a named list of them
ungroup_list <- function(ungrouped) {
  if (inherits(ungrouped, "mt_data")) return(list(ungrouped = ungrouped))
  ok <- is.list(ungrouped) && length(ungrouped) && all(vapply(ungrouped, inherits, TRUE, "mt_data"))
  if (!ok || is.null(names(ungrouped)) || any(!nzchar(names(ungrouped))) || anyDuplicated(names(ungrouped))) {
    mt_abort("`ungrouped` must be an mt_data object or a named list of them.", class = "argument")
  }
  ungrouped
}

# Rates on the grid of each ungrouped data set, with the method name
ungroup_long <- function(us) {
  rbindlist(lapply(names(us), function(n) {
    out <- grid_rates_long(us[[n]])
    out[, method := factor(n, names(us))]
  }))
}

grid_rates_long <- function(u) {
  D <- mt_obs(u, "deaths")
  E <- mt_obs(u, "exposure")
  out <- merge(D[, list(stratum, label, age_lo, age_hi, period = period_lo, deaths = value)],
               E[, list(stratum, age_lo, period = period_lo, exposure = value)], by = c("stratum", "age_lo", "period"))
  out[, rate := deaths / exposure]
  out[]
}

# Observed rates of the age groups of `data`, with the exposures of the
# ungrouped data summed over each group
group_rates <- function(data, u, periods) {
  D <- mt_obs(data, "deaths")[period_lo %in% periods]
  E <- mt_obs(u, "exposure")[period_lo %in% periods]
  top <- max(E$age_lo) + 1
  D[, exposure := vapply(seq_len(.N), function(i) {
    sum(E$value[E$stratum == D$stratum[i] & E$period_lo == D$period_lo[i] & E$age_lo >= D$age_lo[i] &
                  E$age_lo < D$age_hi[i]])
  }, 0)]
  D[, `:=`(rate = value / exposure, x_hi = pmin(age_hi, top), period = period_lo)]
  D[exposure > 0]
}

# Log ratios of ungrouped to known rates at the cells both have
ungroup_errors <- function(U, reference) {
  check_class(reference, "mt_data", "reference")
  R <- grid_rates_long(reference)[, list(label, age_lo, period, ref = rate, noise = 1 / sqrt(deaths))]
  err <- merge(U, R, by = c("label", "age_lo", "period"))
  err <- err[is.finite(log(rate)) & is.finite(log(ref))]
  if (!nrow(err)) mt_abort("`reference` shares no cells with the ungrouped data.", class = "argument")
  err[, lr := log(rate / ref)]
  err[]
}
