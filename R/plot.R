#' Plots
#'
#' All plot functions return ggplot objects (package \pkg{ggplot2} must be
#' installed) and never print or write files.
#'
#' * `mt_plot_forecast()`: observed rates and forecast fan (median and
#'   central intervals) over time for selected ages, log scale. The points
#'   are the training data after preprocessing (for example with
#'   interpolated shock years or smoothed rates); give the original data as
#'   `data` to show the rates as published;
#' * `mt_plot_age()`: age profiles of the forecast for selected years;
#' * `mt_plot_life_expectancy()`: life expectancy forecast fan with the
#'   observed history (from `data` if given, otherwise from the training
#'   data), optionally one panel per level of a stratum variable (`facet`);
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
#' @param data An [mt_data] object (for the forecast plots: the original
#'   data, whose observed rates are shown).
#' @param facet Stratum variable whose levels get one panel each (for
#'   example `"sex"`), or `NULL`.
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
mt_plot_forecast <- function(draws, ages = NULL, strata = NULL, levels = c(0.8, 0.95), data = NULL) {
  need_ggplot()
  check_class(draws, "mt_draws", "draws")
  ages <- ages %||% draws$age$labels[unique(round(seq(1, length(draws$age), length.out = 4)))]
  bad <- setdiff(ages, draws$age$labels)
  if (length(bad)) {
    mt_abort(sprintf("Unknown age label(s): %s. Available: %s.", paste(bad, collapse = ", "),
                     paste(draws$age$labels, collapse = ", ")), class = "argument")
  }
  strata <- check_strata_labels(strata %||% draws$strata$label, draws$strata$label)
  probs <- sort(c((1 - levels) / 2, 0.5, 1 - (1 - levels) / 2))
  s <- mt_summarise(draws, probs)
  s <- s[age %in% ages & label %in% strata]
  obs <- observed_long(draws, data)[age %in% ages & label %in% strata]
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

observed_long <- function(draws, data = NULL) {
  o <- observed_rates(draws, data)
  if (is.null(o)) return(data.table(label = character(), age = character(), period = numeric(), rate = numeric()))
  d <- dim(o)
  data.table(label = rep(dimnames(o)[[3]], each = d[1] * d[2]), age = rep(dimnames(o)[[1]], d[2] * d[3]),
             period = rep(rep(as.numeric(dimnames(o)[[2]]), each = d[1]), d[3]), rate = as.vector(o))
}

# Observed rates [age, period, stratum] on the grid of the draws: from the
# original data if given, otherwise the training data of the fit
observed_rates <- function(draws, data = NULL) {
  if (is.null(data)) {
    o <- draws$observed
    if (is.null(o)) return(NULL)
    dimnames(o) <- list(draws$age$labels, as.character(draws$observed_periods), draws$strata$label)
    return(o)
  }
  check_class(data, "mt_data", "data")
  y <- tryCatch(mt_set_grid(data, age = draws$age), error = function(e) {
    mt_abort("`data` must have the age groups of the forecast (or finer groups that add up to them).",
             class = "argument")
  })
  y <- mt_regroup(y, age = draws$age)
  o <- mt_array(y, "deaths") / mt_array(y, "exposure")
  o[, , intersect(dimnames(o)[[3]], draws$strata$label), drop = FALSE]
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
mt_plot_life_expectancy <- function(draws, at = NULL, levels = c(0.8, 0.95), strata = NULL, facet = NULL,
                                    data = NULL) {
  need_ggplot()
  check_class(draws, "mt_draws", "draws")
  strata <- check_strata_labels(strata %||% draws$strata$label, draws$strata$label)
  probs <- sort(c((1 - levels) / 2, 0.5, 1 - (1 - levels) / 2))
  e <- life_expectancy(draws, at = at, probs = probs)[label %in% strata]
  o <- observed_rates(draws, data)
  hist <- NULL
  if (!is.null(o)) {
    eo <- life_expectancy(o, draws$age, at = at)
    hist <- data.table(label = rep(colnames(eo), each = nrow(eo)), period = rep(as.numeric(rownames(eo)), ncol(eo)),
                       value = as.vector(eo))[label %in% strata]
  }
  if (!is.null(facet)) {
    if (!facet %in% names(draws$strata)) mt_abort("`facet` must be a stratum variable.", class = "argument")
    map <- stats::setNames(draws$strata[[facet]], draws$strata$label)
    e[, panel := map[label]]
    if (!is.null(hist)) hist[, panel := map[label]]
  }
  qn <- function(p) qnames(p)
  p <- ggplot2::ggplot(e, ggplot2::aes(x = period))
  for (lv in sort(levels, decreasing = TRUE)) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[qn((1 - lv) / 2)]], ymax = .data[[qn(1 - (1 - lv) / 2)]],
                                               fill = label), alpha = 0.25)
  }
  p <- p + ggplot2::geom_line(ggplot2::aes(y = q0.5, colour = label))
  if (!is.null(hist)) p <- p + ggplot2::geom_line(data = hist, ggplot2::aes(y = value, colour = label))
  if (!is.null(facet)) p <- p + ggplot2::facet_wrap(~panel)
  p + ggplot2::labs(x = NULL, y = sprintf("Life expectancy at age %s", e$at[1]), colour = NULL, fill = NULL,
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

#' Further plots: Tucker factors, cause shares and backtest scores
#'
#' * `mt_plot_tucker()`: the factor matrices of a [tucker3()] fit (age,
#'   period and stratum modes), one panel per mode and one line per factor;
#'   signs are arbitrary.
#' * `mt_plot_shares()`: shares of each cause in the deaths over the
#'   selected ages, observed and forecast (median and 95% interval), from a
#'   forecast with a cause variable among the strata (for example from
#'   [cause_shares()]).
#' * `mt_plot_scores()`: mean score of each model by horizon or age from a
#'   backtest.
#'
#' @param fit A fitted [tucker3()] model.
#' @param draws An [mt_draws] object whose strata include `var`.
#' @param var Stratum variable with the causes (default: the variable with
#'   shared exposures in the hierarchy of a reconciled forecast).
#' @param total_level Level of `var` for all causes (left out of the shares;
#'   default: its total in the hierarchy).
#' @param ages Age labels summed over (default all).
#' @param groups Panels to show: combinations of the other stratum variables
#'   joined by "." (for example `"female.Q1"`; default all).
#' @param causes Levels of `var` to show (default all; shares are always
#'   computed among all levels).
#' @param backtest An [mt_backtest()] result.
#' @param metric Score column, for example `"crps"`, `"wis"` or `"ae"`.
#' @param by `"horizon"` or `"age"`.
#' @return A ggplot object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 2000:2019)
#'   mt_plot_tucker(mt_fit(tucker3(), x))
#'   bt <- rbind(mt_backtest(rwd(), x, origins = c(2012, 2014), h = 4, n_draws = 50),
#'               mt_backtest(lee_carter(), x, origins = c(2012, 2014), h = 4, n_draws = 50))
#'   mt_plot_scores(bt)
#' }
#' @name mt_plot_more
NULL

#' @rdname mt_plot_more
#' @export
mt_plot_tucker <- function(fit) {
  need_ggplot()
  check_class(fit, "mt_fit", "fit")
  dec <- fit$fits$dec
  if (is.null(dec)) mt_abort("`fit` is not a Tucker3 fit.", class = "argument")
  labs <- list(age = axis_midpoints(fit$age, fit$omega), period = fit$periods, stratum = seq_len(nrow(fit$strata)))
  long <- rbindlist(lapply(1:3, function(n) {
    U <- dec$U[[n]]
    data.table(mode = names(labs)[n], index = rep(labs[[n]], ncol(U)),
               factor = factor(rep(seq_len(ncol(U)), each = nrow(U))), value = as.vector(U))
  }))
  long[, mode := factor(mode, levels = names(labs))]
  ggplot2::ggplot(long, ggplot2::aes(x = index, y = value, colour = factor)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = long[mode == "stratum"]) +
    ggplot2::facet_wrap(~ mode, scales = "free") +
    ggplot2::labs(x = NULL, y = "Loading", colour = "Factor", title = fit$spec$label,
                  caption = sprintf("Strata: %s", paste(seq_len(nrow(fit$strata)), fit$strata$label, sep = " = ",
                                                        collapse = ", ")))
}

#' @rdname mt_plot_more
#' @export
mt_plot_shares <- function(draws, var = NULL, total_level = NULL, ages = NULL, groups = NULL, causes = NULL) {
  need_ggplot()
  check_class(draws, "mt_draws", "draws")
  st <- draws$strata
  if ((is.null(var) || is.null(total_level)) && inherits(draws$spec, "mt_spec_cause_shares")) {
    var <- var %||% draws$spec$params$var
    total_level <- total_level %||% draws$spec$params$total_level
  }
  if (is.null(var) || is.null(total_level)) {
    h <- draws$hierarchy
    shared <- if (!is.null(h)) names(h$vars)[vapply(h$vars, function(v) identical(v$additive, "rates"), logical(1))]
    if (!length(shared)) {
      mt_abort(paste("Give `var` and `total_level` (for example var = \"cause\", total_level = \"all\"), or",
                     "plot a forecast of reconciled() with causes in its hierarchy."), class = "argument")
    }
    var <- var %||% shared[1]
    total_level <- total_level %||% h$vars[[var]]$total
  }
  if (!var %in% names(st)) mt_abort(sprintf("The strata have no variable `%s`.", var), class = "argument")
  ai <- if (is.null(ages)) seq_along(draws$age) else match(ages, draws$age$labels)
  parts <- which(st[[var]] != total_level)
  other <- setdiff(names(st), c("stratum", "label", var))
  grp <- if (length(other)) do.call(paste, c(as.list(st[, other, with = FALSE]), sep = ".")) else rep("all", nrow(st))
  share_tab <- function(Dsum, periods, kind) {
    # Dsum: [period, stratum] or [period, stratum, draw]
    rbindlist(lapply(unique(grp[parts]), function(g) {
      ii <- parts[grp[parts] == g]
      tot <- apply(Dsum[, ii, , drop = FALSE], c(1, 3), sum)
      rbindlist(lapply(ii, function(i) {
        sh <- Dsum[, i, , drop = FALSE][, 1, ] / tot
        sh <- matrix(sh, nrow = length(periods))
        q <- t(apply(sh, 1, stats::quantile, c(0.025, 0.5, 0.975), na.rm = TRUE))
        data.table(group = g, cause = st[[var]][i], period = periods, lo = q[, 1], med = q[, 2], hi = q[, 3],
                   kind = kind)
      }))
    }))
  }
  out <- list()
  if (!is.null(draws$train_D)) {
    Do <- apply(draws$train_D[ai, , , drop = FALSE], c(2, 3), sum)
    out$obs <- share_tab(array(Do, c(dim(Do), 1)), draws$observed_periods, "observed")
  }
  E <- matrix(draws$exposure_last[ai, 1, ], length(ai))
  Df <- apply(sweep(draws$rates[ai, , , , drop = FALSE], c(1, 3), E, "*"), c(2, 3, 4), sum)
  out$fc <- share_tab(Df, draws$periods, "forecast")
  tab <- rbindlist(out)
  if (!is.null(groups)) tab <- tab[group %in% groups]
  if (!is.null(causes)) tab <- tab[cause %in% causes]
  ggplot2::ggplot(tab, ggplot2::aes(x = period, y = med, colour = cause)) +
    ggplot2::geom_ribbon(data = tab[kind == "forecast"], ggplot2::aes(ymin = lo, ymax = hi, fill = cause),
                         alpha = 0.2, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(linetype = kind)) +
    ggplot2::facet_wrap(~ group) +
    ggplot2::labs(x = NULL, y = "Share of deaths", colour = NULL, fill = NULL, linetype = NULL, title = draws$model)
}

#' @rdname mt_plot_more
#' @export
mt_plot_scores <- function(backtest, metric = "crps", by = c("horizon", "age")) {
  need_ggplot()
  by <- match.arg(by)
  if (!metric %in% names(backtest)) mt_abort(sprintf("No column `%s` in the backtest.", metric), class = "argument")
  tab <- as.data.table(unclass_backtest(backtest))
  if (by == "age") {
    s <- tab[, list(score = mean(get(metric), na.rm = TRUE)), by = c("model", "age_lo")]
    setnames(s, "age_lo", "x")
  } else {
    s <- tab[, list(score = mean(get(metric), na.rm = TRUE)), by = c("model", "horizon")]
    setnames(s, "horizon", "x")
  }
  ggplot2::ggplot(s, ggplot2::aes(x = x, y = score, colour = model)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(x = if (by == "age") "Age" else "Horizon (years)", y = sprintf("Mean %s", toupper(metric)),
                  colour = NULL)
}
