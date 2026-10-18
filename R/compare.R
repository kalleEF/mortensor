#' Compare models and forecasts
#'
#' Helpers for the comparisons that most analyses make: several models in
#' one backtest, their scores side by side, and one quantity of several
#' forecasts (models, scenarios, jump-off choices) as a table or a figure.
#'
#' * `mt_compare()` runs [mt_backtest()] for each specification in a named
#'   list on the same data and origins and returns one backtest whose
#'   `model` column holds the names (in the order given). An element can be
#'   a model specification or `list(spec = , prep = )` with its own
#'   preprocessing. Without `origins`, the default origins of
#'   [mt_backtest()] are chosen once for all models, with enough training
#'   years for the model that needs the most.
#' * `mt_score_table()` gives the mean scores of a backtest for every
#'   target it holds side by side: log death rates without suffix, life
#'   expectancy with `_le` (years), the age-standardised rate with `_asr`,
#'   and differences between strata (`contrast` in [mt_backtest()]) with
#'   `_le_gap` and `_asr_gap`. The same model with different preprocessing
#'   gets one row per pipeline (column `prep`, as in [mt_score()]). With
#'   several models, only the cells that every model scored enter (the
#'   same origin, horizon, stratum, age group, year and target, as in
#'   [mt_dm_test()]), with a warning of class `mt_warning_backtest` when
#'   cells are left out: where a model failed at an origin or in a stratum,
#'   the scores of the other models for those cells would otherwise enter
#'   their means only, and the means would refer to different targets.
#'   `common = FALSE` scores every model on its own cells; a target of
#'   which no cell was scored by every model (for example with different
#'   age groups) is always scored so. The CRPS depends on the number of
#'   draws (see [mt_backtest()]); compare models
#'   backtested with the same `n_draws`, as `mt_compare()` does.
#' * `mt_compare_forecasts()` gives one quantity of several forecasts
#'   (`mt_draws`) or observed data (`mt_data`, with Poisson intervals) as
#'   "median (lower-upper)", one column per forecast (empty where a
#'   forecast has no value for the periods shown).
#' * `mt_plot_compare()` draws the same quantity of several forecasts in
#'   one figure, one colour per forecast: medians (dashed), intervals,
#'   fitted values (solid) and the observed values (points). With `periods`
#'   around the last observed year it shows the jump-off of each forecast.
#'
#' @param specs A named list of model specifications (or of
#'   `list(spec = , prep = )`).
#' @param data An [mt_data] object.
#' @param origins,h,n_draws,quantities Passed to [mt_backtest()].
#' @param ... Further arguments of [mt_backtest()].
#' @param backtest An [mt_backtest()] result.
#' @param metrics Score columns of [mt_score()], for example `"crps"`,
#'   `"rmse"`, `"cov80"`, `"cov95"`.
#' @param targets Targets to show (default: all in the backtest).
#' @param by Columns that define the rows (default `"model"`).
#' @param digits Significant digits (`mt_score_table()`) or decimals
#'   (`mt_compare_forecasts()`).
#' @param hierarchy Optional [mt_hierarchy()]: scores are then given
#'   separately for the bottom series, the aggregates and the grand total
#'   (column `level`).
#' @param common Score only the cells that every model scored (default
#'   `TRUE`); `FALSE` scores every model on its own cells.
#' @param forecasts A named list of [mt_draws] objects or [mt_data] objects
#'   (observed values).
#' @param quantity `"life_expectancy"`, `"asr"` (per 100,000, European
#'   Standard Population 2013, see [asr()]), `"deaths"`
#'   (summed over the ages, with the exposures of the last observed year
#'   unless `exposure` is given; forecasts only), `"rate"` (death rates of
#'   the ages in `ages`) or `"lifespan_variation"` (the index `measure` of
#'   [lifespan_variation()]).
#' @param measure Index of lifespan variation for `quantity =
#'   "lifespan_variation"` (see [lifespan_variation()]).
#' @param periods Years shown (for the plot a range `c(first, last)`, or
#'   one year: then medians and intervals per stratum side by side).
#' @param at Age for life expectancy and lifespan variation.
#' @param ages Age labels for `quantity = "rate"`.
#' @param strata Stratum labels shown (default all).
#' @param level Level of the central interval.
#' @param wide One column per forecast (`TRUE`) or a long table.
#' @param exposure Future exposures for `quantity = "deaths"` (see
#'   [mt_deaths()]).
#' @return `mt_compare()`: an `mt_backtest` table. `mt_score_table()`,
#'   `mt_compare_forecasts()`: a `data.table`. `mt_plot_compare()`: a ggplot
#'   object.
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2012)
#' models <- list("random walk" = rwd(), "Lee-Carter" = lee_carter())
#' bt <- mt_compare(models, x, origins = c(2005, 2007), h = 5, n_draws = 50,
#'                  quantities = "life_expectancy")
#' mt_score_table(bt)
#' fc <- list("random walk" = mt_forecast(mt_fit(rwd(), x), h = 10, n_draws = 100, seed = 1),
#'            "Lee-Carter" = mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 100, seed = 1))
#' mt_compare_forecasts(fc, periods = 2022)
#' if (requireNamespace("ggplot2", quietly = TRUE)) mt_plot_compare(fc, data = x)
#' @name mt_compare
NULL

#' @rdname mt_compare
#' @export
mt_compare <- function(specs, data, origins = NULL, h = 10, n_draws = 200, quantities = NULL, ...) {
  if (inherits(specs, "mt_spec")) specs <- list(specs)
  if (!is.list(specs) || !length(specs)) mt_abort("`specs` must be a list of model specifications.", class = "argument")
  items <- lapply(specs, spec_and_prep)
  nm <- names(specs) %||% rep("", length(specs))
  nm[!nzchar(nm)] <- vapply(items[!nzchar(nm)], function(it) it$spec$label, "")
  if (anyDuplicated(nm)) mt_abort("The names of `specs` must be unique.", class = "argument")
  check_class(data, "mt_data", "data")
  # the same origins for every model: enough training years for each
  if (is.null(origins)) {
    origins <- backtest_origins(data, h, list(...)[["n_origins"]] %||% 5,
                                max(vapply(items, function(it) it$spec$min_periods, 0)))
  }
  res <- lapply(seq_along(items), function(i) {
    bt <- mt_backtest(items[[i]]$spec, data, origins = origins, h = h, prep = items[[i]]$prep, n_draws = n_draws,
                      quantities = quantities, ...)
    bt[, model := nm[i]]
  })
  out <- rbindlist(res, fill = TRUE)
  out[, model := factor(model, levels = nm)]
  setattr(out, "class", c("mt_backtest", class(out)))
  out[]
}

# Rows of the cells that every model of a backtest scored (models as in
# mt_dm_test(): "model (prep)" when a model has several pipelines), with a
# warning when cells are left out
common_cells <- function(bt) {
  if (!all(c("model", "origin", "horizon") %in% names(bt))) return(bt)
  who <- dm_model_units(bt)
  n_units <- length(attr(who, "levels"))
  if (n_units < 2) return(bt)
  keys <- intersect(c("origin", "horizon", "stratum", "age_lo", "age_hi", "period", "target"), names(bt))
  cell <- do.call(paste, c(unname(as.list(bt[, keys, with = FALSE])), sep = "\r"))
  cnt <- data.table(cell = cell, unit = as.vector(who))[, list(n = data.table::uniqueN(unit)), by = "cell"]
  keep <- cnt$n[match(cell, cnt$cell)] == n_units
  if (all(keep)) return(bt)
  # a target without any cell scored by every model (for example different
  # age groups with step_regroup()) keeps the cells of each model
  tg <- if ("target" %in% names(bt)) bt[["target"]] else rep("rate", nrow(bt))
  none <- setdiff(unique(tg), unique(tg[keep]))
  keep[tg %in% none] <- TRUE
  msg <- c(if (!all(keep)) {
    sprintf(paste("The models did not score the same cells: %d of %d rows, of cells not scored by every model,",
                  "are left out so that the models are compared on the same cells (common = FALSE scores every",
                  "model on its own cells)."), sum(!keep), length(keep))
  }, if (length(none)) {
    sprintf(paste("No cell of target(s) %s was scored by every model (different age groups?); these targets are",
                  "scored on the cells of each model."), paste(none, collapse = ", "))
  })
  mt_warn(msg, class = "backtest", call = NULL)
  bt[keep]
}

# Level of each backtest row in a hierarchy: "total" where every variable is
# at its total, "aggregate" where some variable is at its total or a group,
# "bottom" otherwise
score_levels <- function(bt, hierarchy) {
  h <- as_hierarchy(hierarchy)
  vars <- names(h$vars)
  miss <- setdiff(vars, names(bt))
  if (length(miss)) {
    mt_abort(sprintf("The backtest has no stratum column %s of the hierarchy.", paste(miss, collapse = ", ")),
             class = "argument")
  }
  agg <- Reduce(`|`, lapply(vars, function(v) bt[[v]] %in% c(h$vars[[v]]$total, names(h$vars[[v]]$groups))))
  tot <- Reduce(`&`, lapply(vars, function(v) bt[[v]] == h$vars[[v]]$total))
  ifelse(tot, "total", ifelse(agg, "aggregate", "bottom"))
}

# Specification and preprocessing of a list element
spec_and_prep <- function(x) {
  if (inherits(x, "mt_spec")) return(list(spec = x, prep = NULL))
  if (is.list(x) && inherits(x$spec, "mt_spec")) return(list(spec = x$spec, prep = x$prep))
  mt_abort("Every element must be a model specification or list(spec = , prep = ).", class = "argument")
}

#' @rdname mt_compare
#' @export
mt_score_table <- function(backtest, metrics = c("crps", "cov95"), targets = NULL, by = "model", digits = 3,
                           hierarchy = NULL, common = TRUE) {
  check_flag(common, "common")
  bt <- as.data.table(backtest)
  if (common) bt <- common_cells(bt)
  if (!is.null(hierarchy)) {
    bt <- copy(bt)
    bt[, level := factor(score_levels(bt, hierarchy), c("bottom", "aggregate", "total"))]
    if (!"level" %in% by) by <- c(by, "level")
  }
  backtest <- bt
  have <- if ("target" %in% names(bt)) unique(bt$target) else "rate"
  targets <- targets %||% intersect(c("rate", "life_expectancy", "asr", "life_expectancy_gap", "asr_gap"), have)
  # the rows of mt_score(), the same for every target (with `prep` when a
  # model was backtested with several pipelines)
  by <- score_by(bt, by)
  suffix <- c(rate = "", life_expectancy = "_le", asr = "_asr", life_expectancy_gap = "_le_gap", asr_gap = "_asr_gap")
  tabs <- lapply(targets, function(tg) {
    s <- mt_score(backtest, by = by, target = tg)
    miss <- setdiff(metrics, names(s))
    if (length(miss)) mt_abort(sprintf("Unknown metric(s): %s.", paste(miss, collapse = ", ")), class = "argument")
    s <- s[, c(by, metrics), with = FALSE]
    setnames(s, metrics, paste0(metrics, suffix[[tg]]))
    s
  })
  out <- Reduce(function(a, b) merge(a, b, by = by, sort = FALSE, all = TRUE), tabs)
  num <- setdiff(names(out), by)
  out[, (num) := lapply(.SD, signif, digits), .SDcols = num]
  out[]
}

#' @rdname mt_compare
#' @export
mt_compare_forecasts <- function(forecasts, quantity = c("life_expectancy", "asr", "deaths", "rate",
                                                         "lifespan_variation"),
                                 periods = NULL, at = NULL, ages = NULL, strata = NULL, level = 0.95, digits = 1,
                                 wide = TRUE, exposure = NULL, measure = "edagger") {
  quantity <- match.arg(quantity)
  measure <- check_compare_measure(quantity, measure)
  tab <- compare_long(forecasts, quantity, at = at, ages = ages, level = level, exposure = exposure,
                      measure = measure)
  if (!is.null(periods)) tab <- tab[period %in% periods]
  if (!is.null(strata)) tab <- tab[label %in% strata]
  if (!nrow(tab)) mt_abort("Nothing left to compare: check `periods` and `strata`.", class = "argument")
  fmt <- function(v) formatC(v, digits = digits, format = "f")
  tab[, value := sprintf("%s (%s-%s)", fmt(med), fmt(lo), fmt(hi))]
  keys <- c("label", if ("age" %in% names(tab)) "age", "period")
  if (!wide) return(tab[, c("forecast", keys, "med", "lo", "hi", "value"), with = FALSE][])
  out <- dcast(tab, stats::as.formula(paste(paste(keys, collapse = " + "), "~ forecast")), value.var = "value")
  # a forecast without values for the periods shown keeps an empty column
  for (f in setdiff(levels(tab$forecast), names(out))) set(out, j = f, value = NA_character_)
  setcolorder(out, c(keys, levels(tab$forecast)))
  out[]
}

# One quantity of several forecasts: columns forecast, label, (age),
# period, med, lo, hi
compare_long <- function(forecasts, quantity, at = NULL, ages = NULL, level = 0.95, exposure = NULL,
                         measure = "edagger") {
  if (inherits(forecasts, c("mt_draws", "mt_data"))) forecasts <- list(forecasts)
  if (!is.list(forecasts) || !length(forecasts)) {
    mt_abort("`forecasts` must be a list of mt_draws or mt_data objects.", class = "argument")
  }
  nm <- names(forecasts) %||% rep("", length(forecasts))
  nm[!nzchar(nm)] <- vapply(forecasts[!nzchar(nm)], function(f) f$model %||% "observed", "")
  if (anyDuplicated(nm)) mt_abort("The names of `forecasts` must be unique.", class = "argument")
  probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
  qn <- qnames(probs)
  one <- function(f) {
    if (quantity == "deaths") {
      if (!inherits(f, "mt_draws")) mt_abort("Deaths are compared for forecasts (mt_draws) only.", class = "argument")
      D <- mt_deaths(f, exposure = exposure, noise = FALSE)
      tot <- apply(D, c(2, 3, 4), sum)
      M <- matrix(tot, ncol = dim(tot)[3])
      Q <- row_quantiles(M, probs)
      return(data.table(label = rep(f$strata$label, each = dim(tot)[1]), period = rep(f$periods, dim(tot)[2]),
                        lo = Q[, 1], med = Q[, 2], hi = Q[, 3]))
    }
    s <- switch(quantity,
                life_expectancy = life_expectancy(f, at = at, probs = probs),
                asr = asr(f, probs = probs),
                lifespan_variation = lifespan_variation(f, measure = measure, at = at, probs = probs),
                rate = {
                  d <- as_draws(f)
                  r <- mt_summarise(d, probs = probs)
                  if (!is.null(ages)) r <- r[age %in% ages]
                  r
                })
    out <- data.table(label = s$label, period = s$period, lo = s[[qn[1]]], med = s[[qn[2]]], hi = s[[qn[3]]])
    if (quantity == "rate") out[, age := factor(s$age, levels = unique(s$age))]
    out
  }
  tab <- rbindlist(lapply(seq_along(forecasts), function(i) one(forecasts[[i]])[, forecast := nm[i]]))
  tab[, forecast := factor(forecast, levels = nm)]
  tab[]
}

#' @rdname mt_compare
#' @param data Observed data shown as points (`mt_plot_compare()`; default:
#'   the training data of the first forecast).
#' @export
mt_plot_compare <- function(forecasts, quantity = c("life_expectancy", "rate", "asr", "lifespan_variation"),
                            at = NULL, ages = NULL, strata = NULL, periods = NULL, level = 0.95, data = NULL,
                            measure = "edagger") {
  need_ggplot()
  quantity <- match.arg(quantity)
  measure <- check_compare_measure(quantity, measure)
  if (inherits(forecasts, "mt_draws")) forecasts <- list(forecasts)
  first <- forecasts[[1]]
  check_class(first, "mt_draws", "forecasts[[1]]")
  if (quantity == "rate" && is.null(ages)) ages <- first$age$labels[unique(round(stats::quantile(
    seq_along(first$age), c(0.25, 0.75))))]
  tab <- compare_long(forecasts, quantity, at = at, ages = ages, level = level, measure = measure)
  hist <- compare_history(first, quantity, at = at, ages = ages, data = data, measure = measure)
  fit <- rbindlist(lapply(seq_along(forecasts), function(i) {
    f <- forecasts[[i]]
    if (!inherits(f, "mt_draws") || is.null(f$fitted)) return(NULL)
    fr <- f$fitted
    dimnames(fr) <- list(f$age$labels, as.character(f$observed_periods), f$strata$label)
    E <- if (quantity == "asr" && identical(dim(f$train_E), dim(fr))) f$train_E
    h <- compare_values(fr, f$age, quantity, at, ages, measure, E)
    if (is.null(h)) NULL else h[, forecast := names(forecasts)[i] %||% f$model]
  }), fill = TRUE)
  if (nrow(fit)) fit[, forecast := factor(forecast, levels = levels(tab$forecast))]
  keep <- function(d) {
    if (is.null(d) || !nrow(d)) return(d)
    if (!is.null(strata)) d <- d[label %in% strata]
    if (!is.null(periods)) d <- d[period >= periods[1] & period <= periods[length(periods)]]
    d
  }
  tab <- keep(tab)
  hist <- keep(hist)
  fit <- keep(fit)
  if (quantity == "rate") {
    lab <- function(d) if (!is.null(d) && nrow(d)) d[, panel := paste0(label, ", age ", age)]
    lab(tab)
    lab(hist)
    lab(fit)
  }
  ylab <- switch(quantity, life_expectancy = sprintf("Life expectancy at age %s", at %||% first$age$lower[1]),
                 rate = "Death rate (log scale)", asr = "Age-standardised rate per 100,000",
                 lifespan_variation = sprintf("%s from age %s", lifespan_labels[[measure]],
                                              at %||% first$age$lower[1]))
  if (length(periods) == 1) {
    # one year: median and interval per stratum (and age), one colour per forecast
    lev <- unique(if (quantity == "rate") tab$panel else tab$label)
    if (!is.null(strata)) {
      first <- if (quantity == "rate") tab$label[match(lev, tab$panel)] else lev
      lev <- lev[order(match(first, strata))]
    }
    tab[, row := factor(if (quantity == "rate") panel else label, levels = rev(lev))]
    return(ggplot2::ggplot(tab, ggplot2::aes(x = med, y = row, colour = forecast)) +
             ggplot2::geom_pointrange(ggplot2::aes(xmin = lo, xmax = hi), size = 0.25,
                                      position = ggplot2::position_dodge(width = 0.7)) +
             (if (quantity == "rate") ggplot2::scale_x_log10()) +
             scale_colour_mortensor(type = "categorical") +
             ggplot2::labs(x = sprintf("%s in %s", ylab, periods), y = NULL, colour = NULL,
                           subtitle = sprintf("Median and %g%% interval", 100 * level)) +
             theme_mortensor(legend = "right"))
  }
  panel <- if (quantity == "rate") ggplot2::facet_wrap(ggplot2::vars(panel), scales = "free_y") else
    if (length(unique(tab$label)) > 1) ggplot2::facet_wrap(ggplot2::vars(label))
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = period, colour = forecast, fill = forecast)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(y = med), linetype = "dashed", linewidth = 0.6)
  if (!is.null(fit) && nrow(fit)) p <- p + ggplot2::geom_line(data = fit, ggplot2::aes(y = value), linewidth = 0.5)
  if (!is.null(hist) && nrow(hist)) {
    p <- p + ggplot2::geom_point(data = hist, ggplot2::aes(x = period, y = value), inherit.aes = FALSE,
                                 colour = "grey30", size = 0.9, alpha = 0.7)
  }
  p <- p + ggplot2::geom_vline(xintercept = first$origin + 0.5, linetype = "dotted", colour = "grey50") + panel
  if (quantity == "rate") p <- p + ggplot2::scale_y_log10()
  p + scale_colour_mortensor(type = "categorical") + scale_fill_mortensor(type = "categorical") +
    ggplot2::labs(x = NULL, y = ylab, colour = NULL, fill = NULL,
                  subtitle = sprintf("Points: observed; solid: fitted; dashed: forecast median with %g%% interval",
                                     100 * level)) +
    theme_mortensor()
}

# Observed values of a quantity, from `data` or the training data of a forecast
compare_history <- function(draws, quantity, at = NULL, ages = NULL, data = NULL, measure = "edagger") {
  o <- observed_rates(draws, data)
  if (is.null(o)) return(NULL)
  E <- if (quantity == "asr") observed_exposures(draws, data, o)
  compare_values(o, draws$age, quantity, at, ages, measure, E)
}

# Exposures [age, period, stratum] that go with observed_rates(): of `data`
# or of the training data of the draws; NULL if there are none
observed_exposures <- function(draws, data, o) {
  if (is.null(data)) {
    E <- draws$train_E
    if (is.null(E) || !identical(dim(E), dim(o))) return(NULL)
    dimnames(E) <- dimnames(o)
    return(E)
  }
  y <- mt_regroup(mt_set_grid(data, age = draws$age), age = draws$age)
  E <- mt_array(y, "exposure")
  E[, dimnames(o)[[2]], dimnames(o)[[3]], drop = FALSE]
}

# One index of lifespan variation for quantity = "lifespan_variation"
check_compare_measure <- function(quantity, measure) {
  if (quantity != "lifespan_variation") return(measure)
  measure <- check_lifespan_measure(measure)
  if (length(measure) != 1) mt_abort("Give one `measure` of lifespan variation.", class = "argument")
  measure
}

# A quantity of a rate array [age, period, stratum] (point values); E: the
# exposures of the same cells (age-standardised rates)
compare_values <- function(m, age, quantity, at = NULL, ages = NULL, measure = "edagger", E = NULL) {
  dn <- dimnames(m)
  if (quantity == "rate") {
    d <- dim(m)
    out <- data.table(label = rep(dn[[3]], each = d[1] * d[2]), age = rep(dn[[1]], d[2] * d[3]),
                      period = rep(rep(as.numeric(dn[[2]]), each = d[1]), d[3]), value = as.vector(m))
    out <- out[age %in% (ages %||% dn[[1]])]
    out[, age := factor(age, levels = unique(age))]
    return(out)
  }
  v <- if (quantity == "life_expectancy") {
    life_expectancy(m, age, at = at)
  } else if (quantity == "lifespan_variation") {
    lifespan_variation(m, measure, age, at = at)
  } else {
    M <- matrix(m, nrow = dim(m)[1])
    array(asr_columns(M, age, standard_population("esp2013"), E, 1e5), dim(m)[-1])
  }
  data.table(label = rep(dn[[3]], each = length(dn[[2]])), period = rep(as.numeric(dn[[2]]), length(dn[[3]])),
             value = as.vector(v))
}

#' Didactic plots of backtests and tuning
#'
#' * `mt_plot_backtest_design()` shows the training and test years of each
#'   origin of a rolling-origin backtest: every model is fitted to the
#'   years up to the origin and scored on the years after it.
#' * `mt_plot_pit()` shows the histogram of the probability integral
#'   transform (PIT: the share of forecast draws below the observed value)
#'   per model. Calibrated forecasts give a flat histogram; a U shape means
#'   intervals that are too narrow, a hump in the middle intervals that are
#'   too wide, a slope a bias.
#' * `mt_plot_tuning()` shows the scores of the settings of [mt_tune()] with
#'   one standard error, the chosen setting highlighted.
#'
#' @param data An [mt_data] object (its years).
#' @param origins Forecast origins (last training years).
#' @param h Horizon.
#' @param tune The result of [mt_tune()].
#' @param backtest An [mt_backtest()] result.
#' @param target `"rate"`, `"life_expectancy"`, `"asr"`,
#'   `"life_expectancy_gap"` or `"asr_gap"` (see `quantities` and `contrast`
#'   in [mt_backtest()]).
#' @param bins Number of histogram bins.
#' @return A ggplot object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 1990:2012)
#'   mt_plot_backtest_design(x, origins = 2003:2007, h = 5)
#' }
#' @name mt_plot_design
NULL

#' @rdname mt_plot_design
#' @export
mt_plot_backtest_design <- function(data, origins, h) {
  need_ggplot()
  check_class(data, "mt_data", "data")
  yrs <- mt_grid(data)$period$lower
  first <- min(yrs)
  last <- max(yrs)
  d <- rbindlist(lapply(origins, function(o) {
    data.table(origin = factor(o, levels = rev(origins)),
               part = factor(c("training", "test"), levels = c("training", "test")),
               from = c(first, o + 1) - 0.5, to = c(o, min(o + h, last)) + 0.5)
  }))
  ggplot2::ggplot(d, ggplot2::aes(y = origin, colour = part)) +
    ggplot2::geom_segment(ggplot2::aes(x = from, xend = to, yend = origin), linewidth = 4) +
    ggplot2::scale_colour_manual(values = c(training = "#9ECAE1", test = "#D95F02")) +
    ggplot2::labs(x = NULL, y = "Origin", colour = NULL,
                  subtitle = "Each model is fitted to the training years and scored on the test years") +
    theme_mortensor()
}

#' @rdname mt_plot_design
#' @export
mt_plot_tuning <- function(tune) {
  need_ggplot()
  s <- copy(tune$scores)
  if (is.null(s$selected)) mt_abort("`tune` must be a result of mt_tune().", class = "argument")
  metric <- intersect(c("crps", "wis", "mae", "rmse", "mase"), names(s))[1]
  set_cols <- setdiff(names(s), c("group", metric, "se", "selected"))
  s[, setting := do.call(paste, c(lapply(set_cols, function(cl) paste(cl, "=", s[[cl]])), sep = ", "))]
  s[, setting := factor(setting, levels = unique(setting))]
  s[, score := s[[metric]]]
  p <- ggplot2::ggplot(s, ggplot2::aes(x = setting, y = score, colour = selected)) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = score - se, ymax = score + se)) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "#D95F02"),
                                 labels = c(`FALSE` = "other", `TRUE` = "chosen")) +
    ggplot2::coord_flip()
  if ("group" %in% names(s)) p <- p + ggplot2::facet_wrap(ggplot2::vars(group))
  p + ggplot2::labs(x = NULL, y = sprintf("%s in the inner backtest (+/- one standard error)", toupper(metric)),
                    colour = NULL) +
    theme_mortensor()
}

#' @rdname mt_plot_design
#' @export
mt_plot_pit <- function(backtest, target = c("rate", "life_expectancy", "asr", "life_expectancy_gap", "asr_gap"),
                        bins = 10) {
  need_ggplot()
  target <- match.arg(target)
  bt <- as.data.table(backtest)
  if ("target" %in% names(bt)) bt <- bt[bt[["target"]] == target]
  if (!nrow(bt)) mt_abort(sprintf("The backtest holds no scores for target \"%s\".", target), class = "argument")
  ggplot2::ggplot(bt, ggplot2::aes(x = pit)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), breaks = seq(0, 1, length.out = bins + 1),
                            fill = "#7570B3", colour = "white") +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    ggplot2::facet_wrap(ggplot2::vars(model)) +
    ggplot2::labs(x = "PIT (share of forecast draws below the observed value)", y = "Density",
                  subtitle = "Calibrated forecasts give a flat histogram (dashed line)") +
    theme_mortensor()
}

#' Plot gaps between strata over time
#'
#' Draws the result of [mt_contrast()] over time: forecast medians (dashed)
#' with intervals and, if given, the observed gaps (points). Several
#' forecasts (models, scenarios) can be compared as a named list, one
#' colour each; otherwise the levels of the contrast get the colours. Other
#' stratum variables (for example sex) give the panels.
#'
#' @param contrast The result of [mt_contrast()] for a forecast, or a named
#'   list of such results.
#' @param observed The result of [mt_contrast()] for observed data, or
#'   `NULL`. Without `contrast`, the observed gaps are drawn with their
#'   Poisson intervals.
#' @param levels Levels of the compared variable to show (default all).
#' @return A ggplot object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 1995:2012)
#'   fc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 100, seed = 1)
#'   mt_plot_contrast(mt_contrast(fc, var = "sex", ref = "male", at = 30),
#'                    observed = mt_contrast(x, var = "sex", ref = "male", at = 30))
#' }
#' @export
mt_plot_contrast <- function(contrast = NULL, observed = NULL, levels = NULL) {
  need_ggplot()
  if (is.null(contrast)) {
    if (is.null(observed)) mt_abort("Give `contrast`, `observed` or both.", class = "argument")
    return(plot_contrast_observed(as.data.table(observed), levels))
  }
  multi <- is.list(contrast) && !is.data.frame(contrast)
  tabs <- if (multi) contrast else list(contrast)
  if (multi && (is.null(names(tabs)) || any(!nzchar(names(tabs))))) {
    mt_abort("Give several contrasts as a named list.", class = "argument")
  }
  tab <- rbindlist(lapply(seq_along(tabs), function(i) {
    as.data.table(tabs[[i]])[, forecast := if (multi) names(tabs)[i] else ""]
  }), fill = TRUE)
  need <- c("level", "ref", "period", "q0.5", "quantity", "type")
  if (!all(need %in% names(tab))) mt_abort("`contrast` must be a result of mt_contrast().", class = "argument")
  qb <- quantile_bounds(names(tab))
  lo <- qb[1]
  hi <- qb[2]
  if (!is.null(levels)) tab <- tab[level %in% levels]
  obs <- if (!is.null(observed)) as.data.table(observed)[if (is.null(levels)) TRUE else level %in% levels]
  pvars <- setdiff(names(tab), c(need, "mean", lo, hi, grep("^q0", names(tab), value = TRUE), "p_positive",
                                 "forecast", "at", "from", "to", "ages", "standard", "per", "age"))
  tab[, col := if (multi) factor(forecast, levels = names(tabs)) else factor(level, levels = unique(level))]
  if (!is.null(obs)) obs[, col := if (multi) NA else factor(level, levels = levels(tab$col))]
  what <- gsub("_", " ", tab$quantity[1])
  ylab <- sprintf("%s in %s", if (tab$type[1] == "ratio") "Ratio" else "Difference", what)
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = period, group = interaction(col, level))) +
    ggplot2::geom_hline(yintercept = if (tab$type[1] == "ratio") 1 else 0, colour = "grey70") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]], fill = col), alpha = 0.15) +
    ggplot2::geom_line(ggplot2::aes(y = q0.5, colour = col), linetype = "dashed", linewidth = 0.6)
  if (!is.null(obs) && nrow(obs)) {
    p <- p + if (multi) {
      ggplot2::geom_point(data = obs, ggplot2::aes(x = period, y = q0.5), inherit.aes = FALSE, colour = "grey30",
                          size = 0.9, alpha = 0.7)
    } else {
      ggplot2::geom_point(data = obs, ggplot2::aes(x = period, y = q0.5, colour = col), inherit.aes = FALSE,
                          size = 0.9, alpha = 0.7)
    }
  }
  if (length(pvars)) {
    p <- p + ggplot2::facet_wrap(pvars)
  }
  pal <- if (multi) "categorical" else if (looks_ordered(levels(tab$col))) "ordered" else "categorical"
  p + scale_colour_mortensor(type = pal) + scale_fill_mortensor(type = pal) +
    ggplot2::labs(x = NULL, y = ylab, colour = NULL, fill = NULL,
                  subtitle = sprintf("Compared with %s; %sdashed: forecast median with interval", tab$ref[1],
                                     if (!is.null(obs)) "points: observed; " else "")) +
    theme_mortensor()
}

# Observed gaps with their (Poisson) intervals
plot_contrast_observed <- function(obs, levels = NULL) {
  if (!all(c("level", "period", "q0.5") %in% names(obs))) {
    mt_abort("`observed` must be a result of mt_contrast().", class = "argument")
  }
  if (!is.null(levels)) obs <- obs[level %in% levels]
  qs <- grep("^q0", names(obs), value = TRUE)
  qb <- quantile_bounds(qs)
  lo <- qb[1]
  hi <- qb[2]
  pvars <- setdiff(names(obs), c("level", "ref", "period", "mean", "p_positive", "quantity", "type", "at", "from",
                                 "to", "ages", "standard", "per", "age", qs))
  obs[, col := factor(level, levels = unique(level))]
  pal <- if (looks_ordered(levels(obs$col))) "ordered" else "categorical"
  p <- ggplot2::ggplot(obs, ggplot2::aes(x = period, y = q0.5, colour = col)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = .data[[lo]], ymax = .data[[hi]]), alpha = 0.6) +
    ggplot2::geom_point(size = 1)
  if (length(pvars)) p <- p + ggplot2::facet_wrap(pvars)
  p + scale_colour_mortensor(type = pal) +
    ggplot2::labs(x = NULL, y = sprintf("%s in %s", if (obs$type[1] == "ratio") "Ratio" else "Difference",
                                        gsub("_", " ", obs$quantity[1])), colour = NULL,
                  subtitle = sprintf("Compared with %s; observed, with Poisson intervals", obs$ref[1])) +
    theme_mortensor()
}

# The lowest and the highest quantile columns (q0.025, q0.1, ...) of a table
quantile_bounds <- function(nms) {
  qs <- grep("^q[0-9.]+$", nms, value = TRUE)
  p <- suppressWarnings(as.numeric(substring(qs, 2)))
  qs <- qs[!is.na(p)]
  p <- p[!is.na(p)]
  c(qs[which.min(p)], qs[which.max(p)])
}
