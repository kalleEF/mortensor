#' Plots
#'
#' All plot functions return ggplot objects (package \pkg{ggplot2} must be
#' installed) and never print or write files.
#'
#' * `mt_plot_forecast()`: observed rates and forecast fan (median and
#'   central intervals) over time for selected ages, log scale;
#' * `mt_plot_age()`: age profiles of the forecast for selected years;
#' * `mt_plot_life_expectancy()`: life expectancy forecast fan;
#' * `mt_plot_rates()`: Lexis surface (age by period) of observed log rates;
#' * `mt_plot_ungroup()`: grouped rates as steps and ungrouped single-year
#'   rates for one period;
#' * `mt_plot_coverage()`: empirical interval coverage by horizon from a
#'   backtest, with the nominal levels.
#'
#' @param draws An [mt_draws] object.
#' @param ages Age labels to show (default: four ages across the range).
#' @param strata Stratum labels to show (default: all).
#' @param levels Central interval levels for the fans.
#' @param periods Forecast years (default: first and last).
#' @param at Age for life expectancy.
#' @param data An [mt_data] object.
#' @param ungrouped The [mt_ungroup()] result for `data`.
#' @param period Year shown.
#' @param backtest An [mt_backtest()] result.
#' @return A ggplot object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single")
#'   fc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 100, seed = 1)
#'   mt_plot_forecast(fc)
#' }
#' @name mt_plot
NULL

need_ggplot <- function() need_pkg("ggplot2", "for plots")

#' @rdname mt_plot
#' @export
mt_plot_forecast <- function(draws, ages = NULL, strata = NULL, levels = c(0.8, 0.95)) {
  need_ggplot()
  check_class(draws, "mt_draws", "draws")
  ages <- ages %||% draws$age$labels[unique(round(seq(1, length(draws$age), length.out = 4)))]
  strata <- strata %||% draws$strata$label
  probs <- sort(c((1 - levels) / 2, 0.5, 1 - (1 - levels) / 2))
  s <- mt_summarise(draws, probs)
  s <- s[age %in% ages & label %in% strata]
  obs <- observed_long(draws)[age %in% ages & label %in% strata]
  qn <- function(p) qnames(p)
  p <- ggplot2::ggplot(s, ggplot2::aes(x = period))
  for (lv in sort(levels, decreasing = TRUE)) {
    lo <- qn((1 - lv) / 2)
    hi <- qn(1 - (1 - lv) / 2)
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]], fill = label), alpha = 0.25)
  }
  p + ggplot2::geom_line(ggplot2::aes(y = q0.5, colour = label)) +
    ggplot2::geom_point(data = obs, ggplot2::aes(y = rate, colour = label), size = 0.6) +
    ggplot2::facet_wrap(~ factor(age, levels = ages), scales = "free_y") +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL, y = "Death rate (log scale)", colour = NULL, fill = NULL, title = draws$model,
                  subtitle = sprintf("Median and %s intervals", paste0(levels * 100, "%", collapse = ", ")))
}

observed_long <- function(draws) {
  o <- draws$observed
  if (is.null(o)) return(data.table(label = character(), age = character(), period = numeric(), rate = numeric()))
  d <- dim(o)
  data.table(label = rep(draws$strata$label, each = d[1] * d[2]), age = rep(draws$age$labels, d[2] * d[3]),
             period = rep(rep(draws$observed_periods, each = d[1]), d[3]), rate = as.vector(o))
}

#' @rdname mt_plot
#' @export
mt_plot_age <- function(draws, periods = NULL, strata = NULL, levels = 0.95) {
  need_ggplot()
  check_class(draws, "mt_draws", "draws")
  periods <- periods %||% range(draws$periods)
  strata <- strata %||% draws$strata$label
  probs <- c((1 - levels) / 2, 0.5, 1 - (1 - levels) / 2)
  s <- mt_summarise(draws, probs)[period %in% periods & label %in% strata]
  names(s)[names(s) %in% qnames(probs)] <- c("lo", "med", "hi")
  ggplot2::ggplot(s, ggplot2::aes(x = age_lo, group = factor(period))) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = factor(period)), alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = med, colour = factor(period))) +
    ggplot2::facet_wrap(~label) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Age", y = "Death rate (log scale)", colour = "Year", fill = "Year", title = draws$model)
}

#' @rdname mt_plot
#' @export
mt_plot_life_expectancy <- function(draws, at = NULL, levels = c(0.8, 0.95)) {
  need_ggplot()
  probs <- sort(c((1 - levels) / 2, 0.5, 1 - (1 - levels) / 2))
  e <- life_expectancy(draws, at = at, probs = probs)
  qn <- function(p) qnames(p)
  p <- ggplot2::ggplot(e, ggplot2::aes(x = period))
  for (lv in sort(levels, decreasing = TRUE)) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[qn((1 - lv) / 2)]], ymax = .data[[qn(1 - (1 - lv) / 2)]],
                                               fill = label), alpha = 0.25)
  }
  p + ggplot2::geom_line(ggplot2::aes(y = q0.5, colour = label)) +
    ggplot2::labs(x = NULL, y = sprintf("Life expectancy at age %s", e$at[1]), colour = NULL, fill = NULL,
                  title = draws$model)
}

#' @rdname mt_plot
#' @export
mt_plot_rates <- function(data, strata = NULL) {
  need_ggplot()
  check_class(data, "mt_data", "data")
  m <- mt_array(data, "deaths") / mt_array(data, "exposure")
  g <- mt_grid(data)
  d <- dim(m)
  long <- data.table(label = rep(data$strata$label, each = d[1] * d[2]), age = rep(g$age$lower, d[2] * d[3]),
                     period = rep(rep(g$period$lower, each = d[1]), d[3]), rate = as.vector(m))
  if (!is.null(strata)) long <- long[label %in% strata]
  ggplot2::ggplot(long, ggplot2::aes(x = period, y = age, fill = log(rate))) +
    ggplot2::geom_raster() +
    ggplot2::facet_wrap(~label) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(x = NULL, y = "Age", fill = "log rate")
}

#' @rdname mt_plot
#' @export
mt_plot_ungroup <- function(data, ungrouped, period = NULL, strata = NULL) {
  need_ggplot()
  g <- mt_grid(ungrouped)
  period <- period %||% max(g$period$lower)
  D <- mt_obs(data, "deaths")[period_lo == period]
  E <- mt_obs(ungrouped, "exposure")[period_lo == period]
  D[, exposure := vapply(seq_len(.N), function(i) {
    sum(E$value[E$stratum == D$stratum[i] & E$age_lo >= D$age_lo[i] & E$age_hi <= D$age_hi[i]])
  }, 0)]
  D[, `:=`(rate = value / exposure, x_hi = pmin(age_hi, max(g$age$lower) + 1))]
  U <- mt_obs(ungrouped, "deaths")[period_lo == period]
  U <- merge(U, E[, list(stratum, age_lo, exposure = value)], by = c("stratum", "age_lo"))
  U[, `:=`(rate = value / exposure, x = age_lo + 0.5)]
  if (!is.null(strata)) {
    D <- D[label %in% strata]
    U <- U[label %in% strata]
  }
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = D, ggplot2::aes(x = age_lo, xend = x_hi, y = rate, yend = rate), colour = "grey40") +
    ggplot2::geom_line(data = U, ggplot2::aes(x = x, y = rate), colour = "firebrick") +
    ggplot2::facet_wrap(~label) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Age", y = "Death rate (log scale)",
                  title = sprintf("Grouped (grey) and ungrouped (red) rates, %s", period))
}

#' @rdname mt_plot
#' @export
mt_plot_coverage <- function(backtest) {
  need_ggplot()
  s <- mt_score(backtest, by = c("model", "horizon"))
  long <- rbind(s[, list(model, horizon, level = 0.5, coverage = cov50)],
                s[, list(model, horizon, level = 0.8, coverage = cov80)],
                s[, list(model, horizon, level = 0.95, coverage = cov95)])
  ggplot2::ggplot(long, ggplot2::aes(x = horizon, y = coverage, colour = model)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = level), linetype = 2, colour = "grey50") +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ paste0(level * 100, "% interval")) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Horizon (years)", y = "Empirical coverage", colour = NULL)
}
