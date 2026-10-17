#' Plots
#'
#' All plot functions return ggplot objects (package \pkg{ggplot2} must be
#' installed) and never print or write files.
#'
#' The forecast plots share one layout (see [mt_plot_style]): observed
#' values as points, fitted values as solid lines, forecast medians as dashed
#' lines with interval ribbons, and a dotted line before the first forecast
#' year. A step between the solid and the dashed line shows the jump from the
#' fitted rates to the forecast. One stratum variable gives the colours (by
#' default the one with the most levels among the shown strata, or
#' `colour`), the others give the panels.
#'
#' * `mt_plot_forecast()`: rates over time, one row of panels per age and one
#'   column per combination of the other stratum variables. The points are
#'   the training data after preprocessing (for example with interpolated
#'   shock years or smoothed rates); give the original data as `data` to show
#'   the rates as published;
#' * `mt_plot_age()`: age profiles of the forecast for selected years
#'   (default: the last);
#' * `mt_plot_life_expectancy()`: life expectancy, observed, fitted and
#'   forecast (observed from `data` if given, otherwise from the training
#'   data), with one panel per level of `facet` if given;
#' * `mt_plot_rates()`: Lexis surface (age by period) of observed log rates;
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
#' @param colour Stratum variable that gives the colours (default: the one
#'   with the most levels among the shown strata).
#' @param log_values Show natural logarithms of the rates on a linear axis
#'   (`log(death rate)`) instead of rates on a logarithmic axis.
#' @param facet Stratum variable whose levels get one panel each (for
#'   example `"sex"`), or `NULL`.
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
mt_plot_forecast <- function(draws, ages = NULL, strata = NULL, levels = 0.95, data = NULL, colour = NULL,
                             log_values = FALSE) {
  need_ggplot()
  check_class(draws, "mt_draws", "draws")
  ages <- ages %||% draws$age$labels[unique(round(seq(1, length(draws$age), length.out = 4)))]
  bad <- setdiff(ages, draws$age$labels)
  if (length(bad)) {
    mt_abort(sprintf("Unknown age label(s): %s. Available: %s.", paste(bad, collapse = ", "),
                     paste(draws$age$labels, collapse = ", ")), class = "argument")
  }
  strata <- check_strata_labels(strata %||% draws$strata$label, draws$strata$label)
  lay <- plot_layout(draws$strata, strata, colour)
  probs <- sort(c((1 - levels) / 2, 0.5, 1 - (1 - levels) / 2))
  s <- mt_summarise(draws, probs)[age %in% ages & label %in% strata]
  obs <- observed_long(draws, data)[age %in% ages & label %in% strata]
  fit <- fitted_long(draws)[age %in% ages & label %in% strata]
  tabs <- lapply(list(s = s, obs = obs, fit = fit), function(d) {
    d <- merge(d, lay$map, by = "label", sort = FALSE)
    d[, age_f := factor(age, levels = ages)]
  })
  tr <- if (log_values) log else identity
  qn <- function(p) qnames(p)
  p <- ggplot2::ggplot(tabs$s, ggplot2::aes(x = period, group = label))
  for (lv in sort(levels, decreasing = TRUE)) {
    lo <- qn((1 - lv) / 2)
    hi <- qn(1 - (1 - lv) / 2)
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = tr(.data[[lo]]), ymax = tr(.data[[hi]]), fill = col),
                                  alpha = 0.3 / length(levels))
  }
  p <- p + ggplot2::geom_line(ggplot2::aes(y = tr(q0.5), colour = col), linetype = "dashed", linewidth = 0.6)
  if (nrow(tabs$fit)) {
    p <- p + ggplot2::geom_line(data = tabs$fit, ggplot2::aes(y = tr(rate), colour = col), linewidth = 0.5)
  }
  p <- p + ggplot2::geom_point(data = tabs$obs, ggplot2::aes(y = tr(rate), colour = col), size = 0.6, alpha = 0.4) +
    ggplot2::geom_vline(xintercept = draws$origin + 0.5, linetype = "dotted", colour = "grey50")
  p <- p + if (lay$has_panel) {
    ggplot2::facet_grid(ggplot2::vars(age_f), ggplot2::vars(panel), scales = "free_y")
  } else {
    ggplot2::facet_wrap(ggplot2::vars(age_f), scales = "free_y")
  }
  if (!log_values) p <- p + ggplot2::scale_y_log10()
  p + plot_colours(lay) +
    ggplot2::labs(x = NULL, y = if (log_values) "log(death rate)" else "Death rate (log scale)",
                  colour = lay$name, fill = lay$name, title = draws$model,
                  subtitle = forecast_subtitle(levels, nrow(tabs$fit) > 0)) +
    theme_mortensor(base_size = if (length(ages) > 2 || lay$has_panel) 9 else 11)
}

forecast_subtitle <- function(levels, fitted = TRUE) {
  sprintf("Points: observed; %sdashed: forecast median with %s interval%s", if (fitted) "solid: fitted; " else "",
          paste0(levels * 100, "%", collapse = " and "), if (length(levels) > 1) "s" else "")
}

observed_long <- function(draws, data = NULL) {
  array_long(observed_rates(draws, data))
}

fitted_long <- function(draws) {
  f <- draws$fitted
  if (!is.null(f)) {
    dimnames(f) <- list(draws$age$labels, as.character(draws$observed_periods), draws$strata$label)
  }
  array_long(f)
}

array_long <- function(o) {
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

# Colour and panels of the strata in a plot: the stratum variable with the
# most levels among the shown strata gives the colour (unless `colour` is
# given), the other varying variables the panels
plot_layout <- function(strata, labels, colour = NULL, panel = NULL) {
  st <- strata[match(labels, strata$label)]
  vars <- setdiff(names(st), c("stratum", "label"))
  nlev <- vapply(vars, function(v) length(unique(st[[v]])), 0)
  if (!is.null(colour) && !colour %in% c(vars, "label")) {
    mt_abort(sprintf("`colour` must be a stratum variable (%s).", paste(vars, collapse = ", ")), class = "argument")
  }
  if (is.null(colour)) {
    cand <- setdiff(vars[nlev > 1], panel)
    colour <- if (length(cand)) cand[max(which(nlev[cand] == max(nlev[cand])))] else "label"
  }
  col <- if (colour == "label") st$label else as.character(st[[colour]])
  pvars <- panel %||% setdiff(vars[nlev > 1], colour)
  pan <- if (length(pvars)) do.call(paste, c(as.list(st[, pvars, with = FALSE]), sep = ".")) else rep("", nrow(st))
  lev <- unique(col)
  map <- data.table(label = st$label, col = factor(col, levels = lev), panel = factor(pan, levels = unique(pan)))
  list(map = map, name = if (colour == "label") NULL else colour, levels = lev, has_panel = length(pvars) > 0,
       ordered = looks_ordered(lev))
}

# Levels such as Q1..Q5, 1..10 or D1..D10 are treated as ordered
looks_ordered <- function(lev) length(lev) > 1 && all(grepl("^[A-Za-z]{0,3}[0-9]+$", lev))

plot_colours <- function(lay) {
  vals <- stats::setNames(mt_palette(length(lay$levels), if (lay$ordered) "ordered" else "categorical"), lay$levels)
  list(ggplot2::scale_colour_manual(values = vals, breaks = lay$levels, aesthetics = c("colour", "fill")))
}

#' @rdname mt_plot
#' @export
mt_plot_age <- function(draws, periods = NULL, strata = NULL, levels = 0.95, colour = NULL, log_values = FALSE) {
  need_ggplot()
  check_class(draws, "mt_draws", "draws")
  periods <- periods %||% max(draws$periods)
  strata <- check_strata_labels(strata %||% draws$strata$label, draws$strata$label)
  lay <- plot_layout(draws$strata, strata, colour)
  probs <- c((1 - levels) / 2, 0.5, 1 - (1 - levels) / 2)
  s <- mt_summarise(draws, probs)[period %in% periods & label %in% strata]
  names(s)[names(s) %in% qnames(probs)] <- c("lo", "med", "hi")
  s <- merge(s, lay$map, by = "label", sort = FALSE)
  s[, panel2 := if (lay$has_panel) paste(panel, period) else as.character(period)]
  tr <- if (log_values) log else identity
  p <- ggplot2::ggplot(s, ggplot2::aes(x = age_lo, group = label)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = tr(lo), ymax = tr(hi), fill = col), alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = tr(med), colour = col), linewidth = 0.8)
  if (length(unique(s$panel2)) > 1) p <- p + ggplot2::facet_wrap(ggplot2::vars(panel2))
  if (!log_values) p <- p + ggplot2::scale_y_log10()
  p + plot_colours(lay) +
    ggplot2::labs(x = "Age", y = if (log_values) "log(death rate)" else "Death rate (log scale)",
                  colour = lay$name, fill = lay$name, title = draws$model,
                  subtitle = sprintf("Forecast median and %s interval, %s", paste0(levels * 100, "%"),
                                     paste(periods, collapse = ", "))) +
    theme_mortensor(legend = "right")
}

#' @rdname mt_plot
#' @export
mt_plot_life_expectancy <- function(draws, at = NULL, levels = 0.95, strata = NULL, facet = NULL,
                                    data = NULL, colour = NULL) {
  need_ggplot()
  check_class(draws, "mt_draws", "draws")
  strata <- check_strata_labels(strata %||% draws$strata$label, draws$strata$label)
  if (!is.null(facet) && !facet %in% names(draws$strata)) {
    mt_abort("`facet` must be a stratum variable.", class = "argument")
  }
  lay <- plot_layout(draws$strata, strata, colour, panel = facet)
  probs <- sort(c((1 - levels) / 2, 0.5, 1 - (1 - levels) / 2))
  e <- life_expectancy(draws, at = at, probs = probs)[label %in% strata]
  e0_long <- function(m) {
    if (is.null(m)) return(NULL)
    eo <- life_expectancy(m, draws$age, at = at)
    data.table(label = rep(colnames(eo), each = nrow(eo)), period = rep(as.numeric(rownames(eo)), ncol(eo)),
               value = as.vector(eo))[label %in% strata]
  }
  hist <- e0_long(observed_rates(draws, data))
  fr <- draws$fitted
  if (!is.null(fr)) dimnames(fr) <- list(draws$age$labels, as.character(draws$observed_periods), draws$strata$label)
  fit <- e0_long(fr)
  add <- function(d) if (is.null(d)) NULL else merge(d, lay$map, by = "label", sort = FALSE)
  e <- add(e)
  hist <- add(hist)
  fit <- add(fit)
  qn <- function(p) qnames(p)
  p <- ggplot2::ggplot(e, ggplot2::aes(x = period, group = label))
  for (lv in sort(levels, decreasing = TRUE)) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[qn((1 - lv) / 2)]], ymax = .data[[qn(1 - (1 - lv) / 2)]],
                                               fill = col), alpha = 0.3 / length(levels))
  }
  p <- p + ggplot2::geom_line(ggplot2::aes(y = q0.5, colour = col), linetype = "dashed", linewidth = 0.6)
  if (!is.null(fit) && nrow(fit)) p <- p + ggplot2::geom_line(data = fit, ggplot2::aes(y = value, colour = col),
                                                              linewidth = 0.5)
  if (!is.null(hist)) p <- p + ggplot2::geom_point(data = hist, ggplot2::aes(y = value, colour = col), size = 0.8,
                                                   alpha = 0.5)
  p <- p + ggplot2::geom_vline(xintercept = draws$origin + 0.5, linetype = "dotted", colour = "grey50")
  if (lay$has_panel) p <- p + ggplot2::facet_wrap(ggplot2::vars(panel))
  p + plot_colours(lay) +
    ggplot2::labs(x = NULL, y = sprintf("Life expectancy at age %s", e$at[1]), colour = lay$name, fill = lay$name,
                  title = draws$model, subtitle = forecast_subtitle(levels, !is.null(fit) && nrow(fit) > 0)) +
    theme_mortensor()
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
    ggplot2::labs(x = NULL, y = "Age", fill = "log rate") +
    theme_mortensor(legend = "right")
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
    scale_colour_mortensor(type = "categorical") +
    ggplot2::labs(x = "Horizon (years)", y = "Empirical coverage", colour = NULL,
                  subtitle = "Dashed: nominal level") +
    theme_mortensor()
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
#'   backtest, or one bar per model (`by = "model"`); `target` selects
#'   rates or a derived quantity scored by the backtest.
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
#' @param by `"horizon"`, `"age"` or `"model"`.
#' @param target `"rate"`, `"life_expectancy"` or `"asr"` (see
#'   [mt_backtest()]).
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
    scale_colour_mortensor(type = "categorical") +
    theme_mortensor(legend = "right") +
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
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_line(ggplot2::aes(linetype = kind), linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = draws$origin + 0.5, linetype = "dotted", colour = "grey50") +
    ggplot2::facet_wrap(~ group) +
    ggplot2::scale_linetype_manual(values = c(observed = "solid", forecast = "dashed"),
                                   guide = ggplot2::guide_legend(order = 2)) +
    scale_colour_mortensor(type = "categorical", guide = ggplot2::guide_legend(order = 1)) +
    scale_fill_mortensor(type = "categorical", guide = ggplot2::guide_legend(order = 1)) +
    ggplot2::labs(x = NULL, y = "Share of deaths", colour = NULL, fill = NULL, linetype = NULL, title = draws$model,
                  subtitle = "Solid: observed; dashed: forecast median with 95% interval") +
    theme_mortensor()
}

#' @rdname mt_plot_more
#' @export
mt_plot_scores <- function(backtest, metric = "crps", by = c("horizon", "age", "model"),
                           target = c("rate", "life_expectancy", "asr", "life_expectancy_gap", "asr_gap")) {
  need_ggplot()
  by <- match.arg(by)
  target <- match.arg(target)
  if (!metric %in% names(backtest)) mt_abort(sprintf("No column `%s` in the backtest.", metric), class = "argument")
  tab <- as.data.table(unclass_backtest(backtest))
  if ("target" %in% names(tab)) {
    keep <- tab[["target"]] == target
    tab <- tab[keep]
  }
  if (!nrow(tab)) mt_abort(sprintf("The backtest holds no scores for target \"%s\".", target), class = "argument")
  if (by == "model") {
    s <- tab[, list(score = mean(get(metric), na.rm = TRUE)), by = "model"]
    s[, model := stats::reorder(model, score)]
    return(ggplot2::ggplot(s, ggplot2::aes(x = model, y = score, fill = model)) +
             ggplot2::geom_col(show.legend = FALSE, alpha = 0.9) +
             ggplot2::coord_flip() +
             scale_fill_mortensor(type = "categorical") +
             ggplot2::labs(x = NULL, y = sprintf("Mean %s", toupper(metric)),
                           title = sprintf("Model comparison (%s)", toupper(metric)),
                           subtitle = sprintf("Backtest of %s; lower is better", gsub("_", " ", target))) +
             theme_mortensor())
  }
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
    scale_colour_mortensor(type = "categorical") +
    ggplot2::labs(x = if (by == "age") "Age" else "Horizon (years)", y = sprintf("Mean %s", toupper(metric)),
                  colour = NULL) +
    theme_mortensor()
}
