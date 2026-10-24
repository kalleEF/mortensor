# Disaggregation of forecast draws to finer ages (proportional allocation by
# a standard schedule), its age shapes and checks, and the representative
# ages of age cells used by the smoothers.

# Representative ages of age cells [lo, hi): the midpoint of a closed cell;
# for an open cell (hi = Inf) its lower bound plus half the width of the
# closed cell below it (one year when there is none), for example 92.5 for
# 90+ after 85-89. A convention of this package, used wherever a smoother or
# spline needs one age per cell (the PCLM instead expands an open interval
# into single years up to omega).
cell_ages <- function(lo, hi) {
  mid <- (lo + hi) / 2
  closed <- which(is.finite(hi))
  for (i in which(is.infinite(hi))) {
    below <- closed[hi[closed] <= lo[i] + 1e-9]
    w <- if (length(below)) {
      j <- below[which.max(hi[below])]
      hi[j] - lo[j]
    } else {
      1
    }
    mid[i] <- lo[i] + w / 2
  }
  mid
}

#' Disaggregate forecast draws to a finer age grid
#'
#' Maps draws of age groups to finer ages (for example single years), path
#' by path, when a model was forecast at the resolution of the observed
#' deaths.
#'
#' The method (`"shape"`) is proportional allocation by a standard
#' schedule. The rate of each finer age is the group rate of the draw times
#' a fixed age-specific factor taken from rates at the finer resolution
#' (`shape_rates`), the factors being scaled to an exposure-weighted mean of
#' 1 within each group, so that the exposure-weighted mean of the finer
#' rates equals the group rate of every draw. [mt_shape_rates()] estimates
#' the shape from the last observed years of data at the finer resolution;
#' with grouped training data, estimate those first with [mt_ungroup()]
#' (for example `mt_shape_rates(mt_ungroup(x), periods)`). A group whose
#' shape rates are all 0 gets a flat shape. Missing group rates (for
#' example of a stratum whose fit failed) give missing rates at the finer
#' ages of that group.
#'
#' @param draws An [mt_draws] object on age groups.
#' @param age The finer [mt_axis]; its cells must nest in the groups.
#' @param exposure Exposures on the finer ages, `[age, stratum]` or
#'   `[age, period, stratum]`, used as weights.
#' @param method `"shape"` (proportional allocation by a standard schedule).
#' @param shape_rates Rates on the finer ages, `[age, stratum]` or `[age,
#'   period, stratum]` for the forecast periods (see [mt_shape_rates()]).
#' @param omega Maximum age.
#' @return An `mt_draws` object on `age`.
#' @examples
#' x <- mt_example(years = 2000:2019, deaths = "single")
#' xg <- mt_regroup(x, age = c(seq(30, 85, 5), 90))
#' fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 3, n_draws = 20, seed = 1)
#' e <- mt_array(x, "exposure")[, 20, ]
#' sh <- mt_shape_rates(x, periods = 2020:2022)
#' fs <- mt_disaggregate(fc, mt_grid(x)$age, exposure = e, shape_rates = sh)
#' dim(fs$rates)
#' @export
mt_disaggregate <- function(draws, age, exposure, method = "shape", shape_rates = NULL, omega = 110) {
  check_class(draws, "mt_draws", "draws")
  check_class(age, "mt_axis", "age")
  match.arg(method)
  ga <- draws$age
  d <- dim(draws$rates)
  if (length(dim(exposure)) == 2) exposure <- array(exposure, c(length(age), 1, d[3]))
  if (dim(exposure)[1] != length(age) || dim(exposure)[3] != d[3]) {
    mt_abort("`exposure` must be [age, stratum] or [age, period, stratum] on the finer ages.", class = "argument")
  }
  if (dim(exposure)[2] == 1) exposure <- exposure[, rep(1, d[2]), , drop = FALSE]
  grp <- findInterval(age$lower, ga$lower)
  if (any(grp == 0) || any(pmin(age$upper, omega) > pmin(ga$upper[pmax(grp, 1)], omega) + 1e-9)) {
    mt_abort("The finer ages do not nest in the age groups of the draws.", class = "argument")
  }
  out <- array(NA_real_, c(length(age), d[2], d[3], d[4]))
  if (is.null(shape_rates)) {
    mt_abort(paste("`shape_rates` is needed: rates at the finer ages, for example mt_shape_rates() of data at the",
                   "finer resolution, or of grouped data ungrouped with mt_ungroup()."), class = "argument")
  }
  if (length(dim(shape_rates)) == 2) shape_rates <- array(shape_rates, c(nrow(shape_rates), 1, ncol(shape_rates)))
  if (dim(shape_rates)[1] != length(age) || dim(shape_rates)[3] != d[3] ||
      !dim(shape_rates)[2] %in% c(1, d[2])) {
    mt_abort("`shape_rates` must be [age, stratum] or [age, period, stratum] on the finer ages and forecast years.",
             class = "argument")
  }
  for (s in seq_len(d[3])) {
    for (t in seq_len(d[2])) {
      e <- exposure[, t, s]
      r <- shape_rates[, min(t, dim(shape_rates)[2]), s]
      norm <- tapply(e * r, grp, sum) / tapply(e, grp, sum)
      nk <- as.vector(norm[as.character(grp)])
      # groups whose shape rates (or exposures) are all 0 get a flat shape
      f <- ifelse(is.finite(nk) & nk > 0, r / nk, 1)
      out[, t, s, ] <- draws$rates[grp, t, s, ] * f
    }
  }
  res <- draws
  dimnames(out) <- list(age = age$labels, period = dimnames(draws$rates)$period,
                        stratum = dimnames(draws$rates)$stratum, draw = NULL)
  res$rates <- out
  res$age <- age
  res$observed <- NULL
  res$fitted <- NULL
  res$exposure_last <- exposure[, 1, , drop = FALSE]
  dimnames(res$exposure_last) <- list(age = age$labels, period = dimnames(exposure)[[2]][1],
                                      stratum = dimnames(draws$rates)$stratum)
  res
}

#' Age shapes within groups for disaggregating forecasts
#'
#' Rates on a fine age grid that give [mt_disaggregate()] the age pattern
#' within each age group of a forecast. The method (`"last"`) takes the
#' rates of the last `years` years (deaths and exposures summed), the same
#' for every forecast year: a standard schedule for proportional allocation
#' within groups.
#'
#' @param data An [mt_data] object at the fine age resolution (for example
#'   ungrouped training data).
#' @param periods Forecast years.
#' @param method `"last"` (the rates of the last `years` years).
#' @param years Number of last observed years used.
#' @return An array `[age, period, stratum]` of rates.
#' @examples
#' x <- mt_example(years = 2000:2019, deaths = "single")
#' sh <- mt_shape_rates(x, periods = 2020:2022)
#' xg <- mt_regroup(x, age = c(seq(30, 85, 5), 90))
#' fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 3, n_draws = 20, seed = 1)
#' fs <- mt_disaggregate(fc, mt_grid(x)$age, exposure = mt_array(x, "exposure")[, 20, ],
#'                       method = "shape", shape_rates = sh)
#' @export
mt_shape_rates <- function(data, periods, method = "last", years = 5) {
  check_class(data, "mt_data", "data")
  match.arg(method)
  check_number(years, "years", lower = 1, integerish = TRUE)
  g <- mt_grid(data)
  D <- mt_array(data, "deaths")
  E <- mt_array(data, "exposure")
  np <- dim(D)[2]
  use <- seq(max(1, np - years + 1), np)
  out <- array(NA_real_, c(dim(D)[1], length(periods), dim(D)[3]),
               dimnames = list(g$age$labels, as.character(periods), data$strata$label))
  r <- apply(D[, use, , drop = FALSE], c(1, 3), sum) / apply(E[, use, , drop = FALSE], c(1, 3), sum)
  out[] <- r[, rep(seq_len(dim(D)[3]), each = length(periods))]
  out
}

#' Check and plot disaggregated forecasts
#'
#' `mt_disaggregate_check()` compares fine-age draws with the grouped draws
#' they come from: the largest relative difference between the
#' exposure-weighted mean of the fine rates and the group rate (0 when the
#' group totals are kept) and the largest absolute second difference of
#' the median log rates over age (large values mean kinks or oscillations).
#' `mt_plot_disaggregation()` shows, for one forecast year, the group medians
#' as steps and the fine-age medians as lines (`type = "profile"`), or the
#' absolute second differences of the fine log rates by age as a heat map
#' (`type = "smoothness"`).
#'
#' @param grouped The [mt_draws] on age groups.
#' @param fine The result of [mt_disaggregate()] for `grouped`.
#' @param exposure The exposures used by [mt_disaggregate()].
#' @param period Forecast year shown (default: the first).
#' @param strata Stratum labels shown (default: all).
#' @param type `"profile"` or `"smoothness"`.
#' @return `mt_disaggregate_check()`: a `data.table` per stratum and period
#'   with `max_group_error` and `max_second_difference`.
#'   `mt_plot_disaggregation()`: a ggplot object.
#' @examples
#' x <- mt_example(years = 2000:2019, deaths = "single")
#' xg <- mt_regroup(x, age = c(seq(30, 85, 5), 90))
#' fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 3, n_draws = 20, seed = 1)
#' e <- mt_array(x, "exposure")[, 20, ]
#' sh <- mt_shape_rates(x, periods = 2020:2022)
#' fs <- mt_disaggregate(fc, mt_grid(x)$age, exposure = e, shape_rates = sh)
#' mt_disaggregate_check(fc, fs, e)
#' if (requireNamespace("ggplot2", quietly = TRUE)) mt_plot_disaggregation(fc, fs)
#' @export
mt_disaggregate_check <- function(grouped, fine, exposure) {
  check_class(grouped, "mt_draws", "grouped")
  check_class(fine, "mt_draws", "fine")
  d <- dim(fine$rates)
  if (length(dim(exposure)) == 2) exposure <- array(exposure, c(d[1], 1, d[3]))
  grp <- findInterval(fine$age$lower, grouped$age$lower)
  rows <- list()
  for (s in seq_len(d[3])) {
    for (t in seq_len(d[2])) {
      e <- exposure[, min(t, dim(exposure)[2]), s]
      Fm <- matrix(fine$rates[, t, s, ], d[1])
      wm <- rowsum(Fm * e, grp) / as.vector(rowsum(e, grp))
      G <- matrix(grouped$rates[sort(unique(grp)), t, s, ], length(unique(grp)))
      med <- log(apply(Fm, 1, stats::median))
      d2 <- if (d[1] > 2) max(abs(diff(med, differences = 2))) else NA_real_
      rows[[length(rows) + 1L]] <- data.table(stratum = s, period = fine$periods[t],
                                              max_group_error = max(abs(wm / G - 1)), max_second_difference = d2)
    }
  }
  out <- merge(fine$strata, rbindlist(rows), by = "stratum", sort = FALSE)
  out[, stratum := NULL]
  out[]
}

#' @rdname mt_disaggregate_check
#' @export
mt_plot_disaggregation <- function(grouped, fine, period = NULL, strata = NULL, type = c("profile", "smoothness")) {
  need_ggplot()
  check_class(grouped, "mt_draws", "grouped")
  check_class(fine, "mt_draws", "fine")
  type <- match.arg(type)
  period <- period %||% fine$periods[1]
  t <- match(period, fine$periods)
  if (is.na(t)) mt_abort("`period` is not a forecast year.", class = "argument")
  strata <- check_strata_labels(strata %||% fine$strata$label, fine$strata$label)
  si <- match(strata, fine$strata$label)
  fm <- apply(fine$rates[, t, si, , drop = FALSE], c(1, 3), stats::median)
  fa <- fine$age
  if (type == "smoothness") {
    sd2 <- apply(log(fm), 2, function(v) c(NA, abs(diff(v, differences = 2)), NA))
    tab <- data.table(age = rep(fa$lower, length(si)), label = rep(strata, each = length(fa)), value = as.vector(sd2))
    return(ggplot2::ggplot(tab, ggplot2::aes(x = age, y = label, fill = pmin(value, 1))) +
             ggplot2::geom_tile() +
             ggplot2::scale_fill_viridis_c(option = "inferno", limits = c(0, 1), name = "|second difference|") +
             ggplot2::labs(x = "Age", y = NULL, title = "Smoothness of the disaggregated rates",
                           subtitle = sprintf("Median log rates, %s; values near 0 are smooth", period)) +
             theme_mortensor(legend = "right"))
  }
  gm <- apply(grouped$rates[, t, si, , drop = FALSE], c(1, 3), stats::median)
  ga <- grouped$age
  hi <- pmin(ga$upper, max(fa$upper[is.finite(fa$upper)], ga$lower + 5))
  gt <- data.table(label = rep(strata, each = length(ga)), lo = rep(ga$lower, length(si)),
                   hi = rep(hi, length(si)), rate = as.vector(gm))
  ft <- data.table(label = rep(strata, each = length(fa)), age = rep(fa$lower + 0.5, length(si)), rate = as.vector(fm))
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = gt, ggplot2::aes(x = lo, xend = hi, y = rate, yend = rate), colour = "grey55",
                          linewidth = 0.8) +
    ggplot2::geom_line(data = ft, ggplot2::aes(x = age, y = rate), colour = "#D7191C", linewidth = 0.6) +
    ggplot2::facet_wrap(ggplot2::vars(label)) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Age", y = "Death rate (log scale)", title = "Disaggregated forecast",
                  subtitle = sprintf("%s, medians; grey: age groups; red: finer ages", period)) +
    theme_mortensor()
}
