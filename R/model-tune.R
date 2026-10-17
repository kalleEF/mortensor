#' Tune model and preprocessing settings by forecast accuracy
#'
#' `mt_tune()` evaluates every combination of settings in `grid` by how well
#' it forecasts the last years of `data` (an inner backtest) and returns the
#' scores and the best specification. `tuned()` wraps the same procedure
#' into a model specification: each time it is fitted (for example at every
#' origin of an outer [mt_backtest()]), the tuning uses only the training
#' data it receives, so the selection never sees the years being forecast.
#' The chosen model is then fitted to all training data.
#'
#' **Design.** Each candidate is refitted at `inner_origins` consecutive
#' origins and forecasts `inner_h` years from each; the last origin is
#' `inner_h` years before the end of the data. `inner_origins = 1` is a
#' single hold-out: fit to the years up to `T - inner_h` and forecast the
#' last `inner_h` years (with `inner_h = 5` or `10`, the usual choice of
#' smoothing and model parameters for `demography::fdm()`). Several origins
#' average over more periods and make the choice depend less on the events
#' of one period (Tashman 2000).
#'
#' **What can be tuned.** `spec_fun` receives one row of `grid` as named
#' arguments and returns a model specification, or `list(spec = , prep = )`
#' to tune preprocessing together with the model: for example the smoothing
#' parameter of [step_smooth()], the fitting window of [step_window()] and
#' the number of components of [hyndman_ullah()]. Preprocessing is given
#' either by `spec_fun` or as `prep`.
#'
#' **Criterion.** `metric` is computed for `target`: the log death rates of
#' the observed cells (`"rate"`), life expectancy at birth
#' (`"life_expectancy"`, years) or the log age-standardised rate (`"asr"`),
#' per test year and stratum. `"crps"` and `"wis"` score the whole forecast
#' distribution; `"mae"` and `"rmse"` the error of the median (point
#' accuracy); `"mase"` scales the absolute error by the mean year-to-year
#' change of the training data (rates only).
#'
#' **Selection.** `select = "min"` takes the lowest score. `select =
#' "one_se"` takes the first row of `grid` whose score is within one
#' standard error of the lowest (the one-standard-error rule of Hastie,
#' Tibshirani and Friedman 2009, Section 7.10), so order `grid` from the
#' simplest to the most complex setting. The standard error comes from the
#' mean scores per origin and stratum (per test year for a hold-out of a
#' single stratum); these blocks are correlated, so it is approximate and
#' tends to be too small. With few
#' test years, differences between settings are often within this error:
#' the selection is then uncertain, and `scores` shows it.
#'
#' **Groups.** With `by`, the strata are split by these stratum variables
#' and each group is tuned and fitted on its own (for example women and men
#' with their own number of components). A joint model such as [li_lee()]
#' then pools only within each group.
#'
#' @param spec_fun A function that returns a model specification (or
#'   `list(spec = , prep = )`) from named arguments, for example
#'   `function(order) hyndman_ullah(order = order)`.
#' @param grid A data frame of settings (columns match the arguments of
#'   `spec_fun`), ordered from the simplest to the most complex setting for
#'   `select = "one_se"`.
#' @param data An [mt_data] object.
#' @param prep Preprocessing pipeline used by every candidate (when
#'   `spec_fun` returns specifications only).
#' @param inner_origins Number of inner origins (`1` for a single hold-out).
#' @param inner_h Horizon of the inner backtest (the number of years held
#'   out from each origin).
#' @param metric `"crps"`, `"wis"`, `"mae"`, `"rmse"` or `"mase"`.
#' @param n_draws Draws per inner forecast.
#' @param strategy `"full"` evaluates every row of `grid`; `"two_stage"`
#'   first evaluates a coarse grid (every second value of each column,
#'   including the first and the last) and then the rows of `grid` within one
#'   step of the best coarse combination in every column. For large grids,
#'   such as Tucker ranks, this needs far fewer fits.
#' @param target `"rate"`, `"life_expectancy"` or `"asr"`.
#' @param select `"min"` or `"one_se"`.
#' @param by Stratum variables: tune each group of strata separately.
#' @return `mt_tune()`: a list with `scores` (one row per setting and group:
#'   the settings, the score, its standard error `se` and whether it was
#'   `selected`), `by_horizon` (scores by forecast horizon), `best` (the
#'   chosen specification; with `by`, a list per group), `best_prep` (the
#'   preprocessing of the chosen setting), `chosen` and `origins`.
#'   `tuned()`: a model specification; its fit carries `tuning` (the scores)
#'   and `tuning_choice`.
#' @references
#' Tashman LJ (2000). Out-of-sample tests of forecasting accuracy: an
#' analysis and review. *International Journal of Forecasting* 16(4):
#' 437-450.
#'
#' Hastie T, Tibshirani R, Friedman J (2009). *The Elements of Statistical
#' Learning*, 2nd edition. Springer, New York.
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2012)
#' tn <- mt_tune(function(drift_window) rwd(drift_window = drift_window),
#'               data.frame(drift_window = c(5, 10, 20)), x, inner_origins = 2, inner_h = 3)
#' tn$scores
#' # a single hold-out of the last 5 years, tuning the smoothing too
#' smoothed_lc <- function(lambda) {
#'   list(spec = lee_carter(), prep = mt_prep(step_smooth(lambda = lambda)))
#' }
#' tn2 <- mt_tune(smoothed_lc, data.frame(lambda = c(1, 100)), x, inner_origins = 1, inner_h = 5,
#'                metric = "rmse", n_draws = 50)
#' tn2$scores
#' @export
mt_tune <- function(spec_fun, grid, data, prep = NULL, inner_origins = 3, inner_h = 5, metric = "crps",
                    n_draws = 100, strategy = c("full", "two_stage"), target = c("rate", "life_expectancy", "asr"),
                    select = c("min", "one_se"), by = NULL) {
  if (!is.function(spec_fun)) mt_abort("`spec_fun` must be a function.", class = "argument")
  check_class(data, "mt_data", "data")
  check_not_ungrouped(data)
  strategy <- match.arg(strategy)
  target <- match.arg(target)
  select <- match.arg(select)
  check_choice(metric, c("crps", "wis", "mae", "rmse", "mase"), "metric")
  if (metric == "mase" && target != "rate") {
    mt_abort("`metric = \"mase\"` is defined for rates only.", class = "argument")
  }
  check_number(inner_origins, "inner_origins", lower = 1, integerish = TRUE)
  check_number(inner_h, "inner_h", lower = 1, integerish = TRUE)
  grid <- as.data.frame(grid)
  if (!nrow(grid)) mt_abort("`grid` has no rows.", class = "argument")
  opts <- list(spec_fun = spec_fun, grid = grid, prep = prep, inner_origins = inner_origins, inner_h = inner_h,
               metric = metric, target = target, select = select, n_draws = n_draws, strategy = strategy)
  if (is.null(by)) return(tune_one(data, opts))
  groups <- strata_groups(data$strata, by)
  res <- lapply(groups, function(idx) tune_one(strata_subset(data, idx), opts))
  with_group <- function(field) {
    rbindlist(lapply(names(res), function(g) cbind(data.table(group = g), res[[g]][[field]])))
  }
  list(scores = with_group("scores"), by_horizon = with_group("by_horizon"), chosen = with_group("chosen"),
       best = lapply(res, `[[`, "best"), best_prep = lapply(res, function(r) r$best_prep),
       origins = res[[1]]$origins, by = by, groups = groups)
}

# Indices of the strata in each group of the stratum variables `by`
strata_groups <- function(strata, by) {
  miss <- setdiff(by, names(strata))
  if (length(miss)) {
    mt_abort(sprintf("`by` names unknown stratum variables: %s.", paste(miss, collapse = ", ")), class = "argument")
  }
  key <- do.call(paste, c(as.list(strata[, by, with = FALSE]), sep = "."))
  split(seq_len(nrow(strata)), factor(key, levels = unique(key)))
}

# Specification and preprocessing of one candidate
tune_candidate <- function(spec_fun, args, prep) {
  out <- do.call(spec_fun, args)
  if (inherits(out, "mt_spec")) return(list(spec = out, prep = prep))
  if (is.list(out) && inherits(out$spec, "mt_spec")) {
    if (!is.null(out$prep) && !is.null(prep)) {
      mt_abort("Give the preprocessing either as `prep` or from `spec_fun`.", class = "argument")
    }
    return(list(spec = out$spec, prep = out$prep %||% prep))
  }
  mt_abort("`spec_fun` must return a model specification or list(spec = , prep = ).", class = "argument")
}

# Score of one inner backtest: mean loss, its standard error over the
# blocks origin x stratum, and the score by horizon
tune_score <- function(bt, tgt, metric) {
  none <- list(value = NA_real_, se = NA_real_, by_h = NULL)
  if (is.null(bt)) return(none)
  b <- bt[target == tgt]
  if (!nrow(b)) return(none)
  lv <- switch(metric, crps = b$crps, wis = b$wis, mae = b$ae, rmse = b$ae^2, mase = b$ae / b$scale)
  ok <- is.finite(lv)
  if (!any(ok)) return(none)
  b <- b[ok]
  lv <- lv[ok]
  agg <- function(v) if (metric == "rmse") sqrt(mean(v)) else mean(v)
  val <- agg(lv)
  blocks <- tapply(lv, paste(b$origin, b$stratum), mean)
  # one origin and one stratum (a hold-out of one population): blocks of
  # test years, whose errors are correlated, so the error is a rough guide
  if (length(blocks) < 2) blocks <- tapply(lv, paste(b$origin, b$stratum, b$period), mean)
  se <- if (length(blocks) > 1) stats::sd(blocks) / sqrt(length(blocks)) else NA_real_
  if (metric == "rmse") se <- se / (2 * val)
  by_h <- data.table(horizon = b$horizon, loss = lv)[, list(value = agg(loss)), keyby = horizon]
  list(value = val, se = se, by_h = by_h)
}

# Tuning of one set of strata
tune_one <- function(data, o) {
  last <- max(mt_grid(data)$period$lower)
  origins <- seq(last - o$inner_origins - o$inner_h + 1, last - o$inner_h)
  quantities <- if (o$target == "rate") NULL else o$target
  acc <- new.env(parent = emptyenv())
  acc$hz <- list()
  score_rows <- function(rows) {
    rbindlist(lapply(rows, function(i) {
      args <- as.list(o$grid[i, , drop = FALSE])
      cand <- tune_candidate(o$spec_fun, args, o$prep)
      bt <- suppressWarnings(tryCatch(mt_backtest(cand$spec, data, origins = origins, h = o$inner_h,
                                                  prep = cand$prep, n_draws = o$n_draws, quantities = quantities),
                                      error = function(e) NULL))
      sc <- tune_score(bt, o$target, o$metric)
      if (!is.null(sc$by_h)) acc$hz[[length(acc$hz) + 1L]] <- cbind(data.table(row = i), sc$by_h)
      cbind(data.table(row = i), as.data.table(o$grid[i, , drop = FALSE]), value = sc$value, se = sc$se)
    }))
  }
  if (o$strategy == "full") {
    scores <- score_rows(seq_len(nrow(o$grid)))
  } else {
    pos <- lapply(o$grid, function(v) match(v, sort(unique(v))))
    coarse <- Reduce(`&`, lapply(pos, function(k) {
      n <- max(k)
      k %% 2 == 1 | k == n
    }))
    scores <- score_rows(which(coarse))
    if (all(is.na(scores$value))) mt_abort("Every setting failed in the inner backtest.", class = "tune")
    b <- scores$row[which.min(scores$value)]
    near <- Reduce(`&`, lapply(pos, function(k) abs(k - k[b]) <= 1))
    scores <- rbind(scores, score_rows(setdiff(which(near), scores$row)))
  }
  if (all(is.na(scores$value))) mt_abort("Every setting failed in the inner backtest.", class = "tune")
  setorder(scores, row)
  best_i <- tune_select(scores, o$select)
  cand <- tune_candidate(o$spec_fun, as.list(o$grid[best_i, , drop = FALSE]), o$prep)
  scores[, selected := row == best_i]
  by_h <- merge(rbindlist(acc$hz), scores[, !c("value", "se", "selected")], by = "row")
  setcolorder(by_h, c(setdiff(names(by_h), c("horizon", "value")), "horizon", "value"))
  chosen <- scores[row == best_i]
  tidy <- function(tab) {
    setnames(tab, "value", o$metric)
    tab[, row := NULL]
  }
  tidy(scores)
  tidy(by_h)
  tidy(chosen)
  chosen[, selected := NULL]
  list(scores = scores[], by_horizon = by_h[], chosen = chosen[], best = cand$spec, best_prep = cand$prep,
       origins = origins)
}

# Row of the chosen setting: the lowest score, or the first row within one
# standard error of it
tune_select <- function(scores, select) {
  i <- which.min(scores$value)
  if (select == "min") return(scores$row[i])
  thr <- scores$value[i] + (if (is.finite(scores$se[i])) scores$se[i] else 0)
  scores$row[which(!is.na(scores$value) & scores$value <= thr)[1]]
}

#' @rdname mt_tune
#' @param label Display name.
#' @export
tuned <- function(spec_fun, grid, inner_origins = 3, inner_h = 5, metric = "crps", n_draws = 100,
                  label = "Tuned model", strategy = c("full", "two_stage"),
                  target = c("rate", "life_expectancy", "asr"), select = c("min", "one_se"), by = NULL) {
  if (!is.function(spec_fun)) mt_abort("`spec_fun` must be a function.", class = "argument")
  check_choice(metric, c("crps", "wis", "mae", "rmse", "mase"), "metric")
  grid <- as.data.frame(grid)
  new_mt_spec("tuned", list(spec_fun = spec_fun, grid = grid, inner_origins = inner_origins,
                            inner_h = inner_h, metric = metric, n_draws = n_draws, strategy = match.arg(strategy),
                            target = match.arg(target), select = match.arg(select), by = by),
              min_periods = inner_origins + inner_h + 5, label = label, smooth = tuned_smooth(spec_fun, grid))
}

# Whether every candidate of a tuned() spec describes smooth rates (TRUE, the
# default of new_mt_spec(), when the candidates cannot be created here; the
# fit carries the chosen spec, from which predictive draws take their noise)
tuned_smooth <- function(spec_fun, grid) {
  tryCatch(all(vapply(seq_len(nrow(grid)), function(i) {
    isTRUE(tune_candidate(spec_fun, as.list(grid[i, , drop = FALSE]), NULL)$spec$smooth)
  }, logical(1))), error = function(e) TRUE)
}

#' @export
mt_fit.mt_spec_tuned <- function(spec, data, prep = NULL, cores = 1, ...) {
  p <- spec$params
  tn <- mt_tune(p$spec_fun, p$grid, data, prep = prep, inner_origins = p$inner_origins, inner_h = p$inner_h,
                metric = p$metric, n_draws = p$n_draws, strategy = p$strategy %||% "full",
                target = p$target %||% "rate", select = p$select %||% "min", by = p$by)
  if (is.null(p$by)) {
    fit <- mt_fit(tn$best, data, prep = tn$best_prep, cores = cores)
    fit$tuning <- tn$scores
    fit$tuning_choice <- tn$chosen
    fit$spec$label <- sprintf("%s: %s (%s)", spec$label, tn$best$label, setting_label(tn$chosen, names(p$grid)))
    return(fit)
  }
  fits <- lapply(names(tn$groups), function(g) {
    mt_fit(tn$best[[g]], strata_subset(data, tn$groups[[g]]), prep = tn$best_prep[[g]], cores = cores)
  })
  spec$label <- sprintf("%s by %s", spec$label, paste(p$by, collapse = " x "))
  fit <- combine_group_fits(fits, tn$groups, data, spec)
  fit$tuning <- tn$scores
  fit$tuning_choice <- tn$chosen
  fit
}

# "name = value" of the chosen settings
setting_label <- function(chosen, cols) {
  paste(vapply(cols, function(cl) paste(cl, "=", format(chosen[[cl]][1])), ""), collapse = ", ")
}

# One fit object from fits of groups of strata (in the order of data$strata)
combine_group_fits <- function(fits, groups, data, spec) {
  st <- data$strata
  ns <- nrow(st)
  f1 <- fits[[1]]
  same <- vapply(fits, function(f) identical(f$age$labels, f1$age$labels) && identical(f$periods, f1$periods),
                 logical(1))
  if (!all(same)) mt_abort("The groups were fitted on different age or period grids.", class = "model")
  bind <- function(field, fallback = field) {
    a1 <- f1[[field]] %||% f1[[fallback]]
    if (is.null(a1)) return(NULL)
    dn <- dimnames(a1)
    if (length(dim(a1)) == 2) {
      out <- matrix(NA_real_, nrow(a1), ns, dimnames = list(dn[[1]], st$label))
      for (g in seq_along(fits)) out[, groups[[g]]] <- fits[[g]][[field]] %||% fits[[g]][[fallback]]
    } else {
      out <- array(NA_real_, c(dim(a1)[1:2], ns), dimnames = c(dn[1:2], list(st$label)))
      for (g in seq_along(fits)) out[, , groups[[g]]] <- fits[[g]][[field]] %||% fits[[g]][[fallback]]
    }
    out
  }
  any_raw <- any(vapply(fits, function(f) !is.null(f$raw_D), logical(1)))
  # diagnostics of all groups in the order of the strata (an ensemble has one
  # row per component and stratum, a bootstrap one more for all strata)
  diag <- rbindlist(lapply(fits, `[[`, "diagnostics"), fill = TRUE)
  diag <- diag[order(match(diag$stratum, st$label))]
  structure(list(spec = spec, prep = NULL, fits = NULL, group_fits = fits, groups = groups,
                 diagnostics = diag, age = f1$age, periods = f1$periods,
                 strata = st, omega = data$omega, observed = bind("observed"), train_D = bind("train_D"),
                 raw_D = if (any_raw) bind("raw_D", "train_D"), train_E = bind("train_E"),
                 exposure_last = bind("exposure_last"), level_totals = data$level_totals,
                 noise_sd = bind("noise_sd")),
            class = c("mt_fit_split", "mt_fit"))
}

#' @export
mt_forecast.mt_fit_split <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                                     exposure = NULL, ...) {
  last <- max(fit$periods)
  if (!is.null(to)) {
    if (!missing(h)) mt_warn("Both `h` and `to` are given; `to` is used.", class = "argument")
    h <- to - last
  }
  check_number(h, "h", lower = 1, integerish = TRUE)
  ns <- nrow(fit$strata)
  arr <- array(NA_real_, c(length(fit$age), h, ns, n_draws))
  with_seed(seed, {
    for (g in seq_along(fit$group_fits)) {
      idx <- fit$groups[[g]]
      fc <- suppressWarnings(mt_forecast(fit$group_fits[[g]], h = h, n_draws = n_draws, predictive = predictive,
                                         exposure = exposure_strata(exposure, idx, fit$strata$label)))
      arr[, , idx, ] <- fc$rates
    }
  })
  check_draws_plausible(arr, fit$age, fit$strata, fit$spec$label)
  out <- new_mt_draws(arr, fit$age, last + seq_len(h), fit$strata, model = fit$spec$label, origin = last,
                      fit = fit, seed = seed)
  out$predictive <- predictive
  out
}

# Exposures of some strata ([age, stratum] or [age, period, stratum], by
# label when the array has stratum names)
exposure_strata <- function(exposure, idx, labels) {
  if (is.null(exposure)) return(NULL)
  E <- if (inherits(exposure, "mt_data")) mt_array(exposure, "exposure") else exposure
  k <- length(dim(E))
  if (!k %in% 2:3) return(E)
  dn <- dimnames(E)[[k]]
  sel <- if (!is.null(dn)) match(labels[idx], dn) else idx
  if (anyNA(sel)) {
    mt_abort(sprintf("`exposure` lacks strata of the forecast: %s.", paste(labels[idx][is.na(sel)], collapse = ", ")),
             class = "argument")
  }
  if (k == 2) E[, sel, drop = FALSE] else E[, , sel, drop = FALSE]
}

#' @export
fitted.mt_fit_split <- function(object, ...) {
  out <- array(NA_real_, dim(object$train_D), dimnames(object$train_D))
  for (g in seq_along(object$group_fits)) {
    f <- tryCatch(stats::fitted(object$group_fits[[g]]), error = function(e) NULL)
    if (!is.null(f)) out[, , object$groups[[g]]] <- f
  }
  out
}

#' @export
print.mt_fit_split <- function(x, ...) {
  NextMethod()
  ch <- x$tuning_choice
  if (!is.null(ch)) {
    cols <- setdiff(names(ch), c("group", "se", names(ch)[ncol(ch) - 1L]))
    lab <- vapply(seq_along(x$group_fits), function(g) {
      sprintf("%s (%s)", x$group_fits[[g]]$spec$label, setting_label(ch[g], cols))
    }, "")
    cat(paste0("  ", names(x$groups), ": ", lab, "\n"), sep = "")
  }
  invisible(x)
}
