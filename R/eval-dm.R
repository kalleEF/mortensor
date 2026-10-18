#' Test differences in forecast accuracy between models (Diebold-Mariano)
#'
#' Tests whether models of a backtest differ in their expected score with
#' the Diebold-Mariano test (Diebold and Mariano 1995) and the small-sample
#' correction of Harvey, Leybourne and Newbold (1997). It complements
#' [mt_score_table()], which shows the mean scores: the test says whether a
#' difference is larger than the variation between forecast origins.
#'
#' The test uses one loss per forecast origin, the time index of a
#' rolling-origin backtest. Scores of the cells of one origin and horizon
#' (ages, strata) concern the same test year and are strongly dependent, so
#' they are averaged first:
#'
#' * per horizon `h`: the mean score over the cells of each origin at that
#'   horizon gives one loss per origin and model. The differences between
#'   the losses of a model and of `ref` form a series over the origins,
#'   tested with the autocovariances up to lag `h - 1` (forecasts `h` years
#'   ahead from consecutive origins share test years and errors).
#' * overall (`horizon = "all"`): the mean of these per-horizon losses over
#'   all horizons gives one loss per origin, tested with the largest horizon
#'   `H` (lags up to `H - 1`). Only origins scored at every horizon enter.
#'
#' Only cells scored for both models enter (the same origin, horizon,
#' stratum, age group, year and target). The statistic, its variance (sum of
#' autocovariances, or Bartlett weights with `varestimator = "bartlett"`),
#' the correction factor and the t distribution with `T - 1` degrees of
#' freedom (`T` origins) are those of `forecast::dm.test(e1, e2, h = h,
#' power = 1)` with `e1` and `e2` the losses per origin (scores are not
#' negative, so the loss is the score itself). With origins spaced `s`
#' years apart (for example `origins = seq(1981, 2001, by = 10)`), test years
#' overlap only up to lag `ceiling(h / s) - 1`, which is used instead; with
#' unequal spacing the smallest gap counts.
#'
#' With few origins the test has little power: a backtest with 6 origins
#' finds only large and consistent differences. When the number of origins
#' `T` is not larger than the number of lags needed (`T <= h` for
#' consecutive origins), or the variance estimate is not positive, the test
#' falls back to lag 0 (column `lag`) with a warning of class
#' `mt_warning_dm`. Loss differences of overlapping test years are mostly
#' positively correlated, so these p-values are usually too small. For long
#' horizons, use more origins or origins spaced `h` years apart. The p-values
#' are not adjusted for the number of tests (models, targets, horizons); see
#' [stats::p.adjust()].
#'
#' @param backtest An [mt_backtest()] or [mt_compare()] result with at least
#'   two models. Models are identified by `model`, and by
#'   `"model (prep)"` when a model was backtested with several
#'   preprocessing pipelines.
#' @param ref Reference model (default: the first model). Every other model
#'   is compared with it.
#' @param metrics Scores to compare: `"crps"`, `"wis"`, `"ae"` (absolute
#'   error of the median) and `"se"` (squared error of the median). CRPS and
#'   WIS compare whole predictive distributions; `"ae"` and `"se"` the
#'   medians only.
#' @param targets Targets to test (default: all in the backtest; see
#'   `quantities` and `contrast` in [mt_backtest()]).
#' @param by Optional stratum columns (for example `"sex"`): the losses per
#'   origin and the tests are then computed separately for each of their
#'   values.
#' @param alternative `"two.sided"`, `"less"` (the model has a lower expected
#'   score than `ref`, that is, it is more accurate) or `"greater"`.
#' @param varestimator `"acf"` (sum of the autocovariances, default) or
#'   `"bartlett"` (Bartlett weights), as in `forecast::dm.test()`.
#' @return A `data.table` with one row per model (other than `ref`), value
#'   of `by`, target, metric and horizon: `model`, `ref`, the `by` columns,
#'   `target`, `metric`, `horizon` (a factor: the horizons, then `"all"`),
#'   `n_origins` (number of losses tested), `lag` (largest autocovariance
#'   lag used), `score` and `score_ref` (mean loss over the origins),
#'   `mean_diff` (`score - score_ref`; negative when the model is more
#'   accurate), `statistic` (Diebold-Mariano statistic with the correction
#'   of Harvey et al.) and `p_value`.
#' @references
#' Diebold FX, Mariano RS (1995). Comparing predictive accuracy. *Journal of
#' Business & Economic Statistics* 13(3): 253-263.
#'
#' Harvey D, Leybourne S, Newbold P (1997). Testing the equality of
#' prediction mean squared errors. *International Journal of Forecasting*
#' 13(2): 281-291.
#' @seealso [mt_compare()], [mt_score_table()], [mt_score()].
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2012)
#' bt <- mt_compare(list(rw = rwd(), lc = lee_carter()), x, origins = 2001:2007, h = 3, n_draws = 50)
#' mt_dm_test(bt)
#' @export
mt_dm_test <- function(backtest, ref = NULL, metrics = "crps", targets = NULL, by = NULL,
                       alternative = c("two.sided", "less", "greater"), varestimator = c("acf", "bartlett")) {
  alternative <- match.arg(alternative)
  varestimator <- match.arg(varestimator)
  if (!is.data.frame(backtest)) mt_abort("`backtest` must be a backtest table.", class = "argument")
  bt <- as.data.table(backtest)
  miss <- setdiff(c("model", "origin", "horizon"), names(bt))
  if (length(miss)) {
    mt_abort(sprintf("`backtest` lacks column(s) %s; give a result of mt_backtest().", paste(miss, collapse = ", ")),
             class = "argument")
  }
  ok_metrics <- c("crps", "wis", "ae", "se")
  if (!is.character(metrics) || !length(metrics) || !all(metrics %in% ok_metrics)) {
    mt_abort(sprintf("`metrics` can be %s.", paste0("\"", ok_metrics, "\"", collapse = ", ")), class = "argument")
  }
  have <- setdiff(names(bt), "se")
  miss <- setdiff(ifelse(metrics == "se", "ae", metrics), have)
  if (length(miss)) {
    mt_abort(sprintf("The backtest has no score column(s) %s.", paste(unique(miss), collapse = ", ")),
             class = "argument")
  }
  reserved <- c("model", "prep", "origin", "horizon", "target", "period", "age_lo", "age_hi")
  if (length(by) && (!is.character(by) || !all(by %in% names(bt)) || any(by %in% reserved))) {
    mt_abort("`by` must name stratum columns of the backtest.", class = "argument")
  }
  if (!"target" %in% names(bt)) bt[, target := "rate"]
  if (!is.null(targets)) {
    bad <- setdiff(targets, unique(bt$target))
    if (length(bad)) {
      mt_abort(sprintf("The backtest holds no scores for target(s) %s.", paste(bad, collapse = ", ")),
               class = "argument")
    }
    bt <- bt[target %in% targets]
  }
  unit <- dm_model_units(bt)
  units <- attr(unit, "levels")
  if (length(units) < 2) mt_abort("The backtest must hold at least two models.", class = "argument")
  ref <- ref %||% units[1]
  if (!is.character(ref) || length(ref) != 1 || !ref %in% units) {
    mt_abort(sprintf("`ref` must be one of the models: %s.", paste0("\"", units, "\"", collapse = ", ")),
             class = "argument")
  }
  keys <- intersect(c("origin", "horizon", "stratum", "age_lo", "age_hi", "period", "target"), names(bt))
  cell <- do.call(paste, c(unname(as.list(bt[, keys, with = FALSE])), sep = "\r"))
  for (u in units) {
    if (anyDuplicated(cell[unit == u])) {
      mt_abort(sprintf(paste("The backtest holds cells of model \"%s\" more than once; give each backtest its own",
                             "model name (mt_compare())."), u), class = "argument")
    }
  }
  ir <- which(unit == ref)
  groups <- c(by, "target")
  rows <- list()
  for (u in setdiff(units, ref)) {
    iu <- which(unit == u)
    jr <- match(cell[iu], cell[ir])
    iu <- iu[!is.na(jr)]
    jr <- ir[jr[!is.na(jr)]]
    if (!length(iu)) {
      mt_warn(sprintf("Models \"%s\" and \"%s\" have no scored cells in common.", u, ref), class = "dm")
      next
    }
    for (m in metrics) {
      d <- bt[iu, c(by, "target", "origin", "horizon"), with = FALSE]
      set(d, j = "e1", value = dm_loss(bt, iu, m))
      set(d, j = "e2", value = dm_loss(bt, jr, m))
      res <- dm_tests(d, groups, alternative, varestimator)
      set(res, j = "model", value = u)
      set(res, j = "ref", value = ref)
      set(res, j = "metric", value = m)
      rows[[length(rows) + 1L]] <- res
    }
  }
  if (!length(rows)) mt_abort("No pair of models has scored cells in common.", class = "dm")
  out <- rbindlist(rows)
  dm_warn_fallbacks(out)
  # order: models, `by`, targets, metrics as given, horizons
  hs <- sort(unique(suppressWarnings(as.numeric(setdiff(out$horizon, "all")))))
  set(out, j = "horizon", value = factor(out$horizon, levels = c(as.character(hs), "all")))
  tg <- c("rate", "life_expectancy", "asr", "life_expectancy_gap", "asr_gap")
  tg <- c(intersect(tg, out$target), setdiff(unique(out$target), tg))
  set(out, j = ".mo", value = match(out$model, units))
  set(out, j = ".tg", value = match(out$target, tg))
  set(out, j = ".me", value = match(out$metric, metrics))
  setorderv(out, c(".mo", by, ".tg", ".me", "horizon"))
  cols <- c("model", "ref", by, "target", "metric", "horizon", "n_origins", "lag", "score", "score_ref", "mean_diff",
            "statistic", "p_value")
  out[, cols, with = FALSE][]
}

# Loss of the rows `i` of a backtest: a score column, or the squared error
# of the median ("se")
dm_loss <- function(bt, i, metric) {
  if (metric == "se") bt$ae[i]^2 else bt[[metric]][i]
}

# Model of each row of a backtest: `model`, or "model (prep)" when a model
# was backtested with several preprocessing pipelines; the levels in the
# order of the models
dm_model_units <- function(bt) {
  mod <- as.character(bt$model)
  lev <- if (is.factor(bt$model)) intersect(levels(bt$model), mod) else unique(mod)
  unit <- mod
  if ("prep" %in% names(bt)) {
    n_prep <- tapply(as.character(bt$prep), mod, function(p) length(unique(p)))
    several <- mod %in% names(n_prep)[n_prep > 1]
    unit[several] <- sprintf("%s (%s)", mod[several], as.character(bt$prep)[several])
  }
  ord <- order(match(mod, lev), seq_along(mod))
  structure(unit, levels = unique(unit[ord]))
}

# Tests per horizon and overall for the losses `e1` (model) and `e2` (ref)
# of the cells in `d` (columns origin, horizon and `groups`)
dm_tests <- function(d, groups, alternative, varestimator) {
  po <- d[, list(e1 = mean(e1), e2 = mean(e2)), keyby = c(groups, "horizon", "origin")]
  per_h <- po[, dm_series(e1, e2, origin, .BY$horizon, alternative, varestimator), by = c(groups, "horizon")]
  per_h[, horizon := as.character(horizon)]
  # overall: per origin, the mean of the per-horizon losses; origins scored
  # at every horizon of the group
  po[, `:=`(n_h = length(unique(horizon)), h_max = max(horizon)), by = groups]
  pa <- po[, list(e1 = mean(e1), e2 = mean(e2), full = .N == n_h[1], h_max = h_max[1]), by = c(groups, "origin")]
  pa <- pa[full == TRUE]
  if (!nrow(pa)) return(per_h)
  all_h <- pa[, dm_series(e1, e2, origin, h_max[1], alternative, varestimator), by = groups]
  all_h[, horizon := "all"]
  rbind(per_h, all_h, use.names = TRUE)
}

# Diebold-Mariano test of one series of losses over origins. `h` is the
# horizon in years; with origins `s` years apart, ceiling(h / s) steps
dm_series <- function(e1, e2, origin, h, alternative, varestimator) {
  o <- order(origin)
  e1 <- e1[o]
  e2 <- e2[o]
  step <- if (length(origin) > 1) min(diff(origin[o])) else 1
  r <- dm_stat(e1 - e2, max(1L, as.integer(ceiling(h / step - 1e-9))), alternative, varestimator)
  c(list(n_origins = length(e1), score = mean(e1), score_ref = mean(e2), mean_diff = mean(e1 - e2)), r)
}

# Diebold-Mariano statistic with the correction of Harvey, Leybourne and
# Newbold (1997) for the loss differences `d` and `h` steps (lags 0 to
# h - 1), as forecast::dm.test(). Falls back to lag 0 when there are too few
# losses for the lags (`fallback = "origins"`) or the variance estimate is
# not positive (`"variance"`); NA without a variance (`"none"`)
dm_stat <- function(d, h, alternative = "two.sided", varestimator = "acf") {
  n <- length(d)
  na <- list(lag = NA_integer_, statistic = NA_real_, p_value = NA_real_, fallback = "none")
  if (n < 2 || anyNA(d)) return(na)
  fallback <- ""
  if (h >= n) {
    fallback <- "origins"
    h <- 1L
  }
  v <- dm_variance(d, h, varestimator)
  # a constant difference leaves only rounding error in the variance
  tiny <- (1e-8 * max(abs(d)))^2 / n
  if (!(v > tiny) && h > 1L) {
    fallback <- "variance"
    h <- 1L
    v <- dm_variance(d, h, varestimator)
  }
  if (!(v > tiny)) return(na)
  k <- sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)
  stat <- mean(d) / sqrt(v) * k
  p <- switch(alternative,
              two.sided = 2 * stats::pt(-abs(stat), df = n - 1),
              less = stats::pt(stat, df = n - 1),
              greater = stats::pt(stat, df = n - 1, lower.tail = FALSE))
  list(lag = h - 1L, statistic = stat, p_value = p, fallback = fallback)
}

# Variance of the mean loss difference from the autocovariances up to lag
# h - 1 (denominator n, as stats::acf())
dm_variance <- function(d, h, varestimator) {
  n <- length(d)
  dc <- d - mean(d)
  gam <- vapply(seq_len(h) - 1L, function(k) sum(dc[seq_len(n - k)] * dc[seq_len(n - k) + k]) / n, 0)
  w <- if (varestimator == "acf" || h == 1L) rep(1, h - 1L) else 1 - seq_len(h - 1L) / h
  (gam[1] + 2 * sum(w * gam[-1])) / n
}

# One warning for the tests that fell back to lag 0 or have no statistic
dm_warn_fallbacks <- function(out) {
  horizons <- function(f) {
    h <- unique(as.character(out$horizon[out$fallback == f]))
    if (!length(h)) return(NULL)
    num <- suppressWarnings(as.numeric(h))
    paste(c(sort(num[!is.na(num)]), h[is.na(num)]), collapse = ", ")
  }
  h_orig <- horizons("origins")
  h_var <- horizons("variance")
  h_none <- horizons("none")
  msg <- c(
    if (!is.null(h_orig)) {
      sprintf(paste("Too few origins for the autocovariances up to lag h - 1 at horizon(s) %s: lag 0 was used",
                    "(column `lag`). Losses of overlapping test years are correlated, so these p-values are too",
                    "small; use more origins or origins spaced h years apart."), h_orig)
    },
    if (!is.null(h_var)) {
      sprintf(paste("The variance estimate with lags up to h - 1 was not positive at horizon(s) %s: lag 0 was",
                    "used, as by forecast::dm.test(). Try varestimator = \"bartlett\"."), h_var)
    },
    if (!is.null(h_none)) {
      sprintf("No test at horizon(s) %s: fewer than 2 origins or no variation of the loss difference.", h_none)
    })
  if (length(msg)) mt_warn(msg, class = "dm", call = NULL)
  invisible(NULL)
}
