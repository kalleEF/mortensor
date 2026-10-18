#' Test differences in forecast accuracy between models (Diebold-Mariano)
#'
#' Tests whether models of a backtest differ in their expected score with
#' the Diebold-Mariano test (Diebold and Mariano 1995). By default the
#' long-run variance is the equal-weighted cosine (EWC) estimate with
#' fixed-smoothing critical values from a t distribution, as Lazarus, Lewis,
#' Stock and Watson (2018) recommend for tests of a single restriction; the
#' Diebold-Mariano statistic is such a test, of the mean loss difference.
#' The fixed-smoothing test of Coroneo and Iacone (2020) for the
#' Diebold-Mariano test, with a Daniell-kernel periodogram estimate, and
#' the version of `forecast::dm.test()`, with the small-sample correction of
#' Harvey, Leybourne and Newbold (1997), are options. It complements
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
#'   the losses of a model and of `ref` form a series over the origins. The
#'   errors of forecasts `h` years ahead from consecutive origins build up
#'   over overlapping years, so the differences are autocorrelated.
#' * overall (`horizon = "all"`): the mean of these per-horizon losses over
#'   all horizons gives one loss per origin. Only origins scored at every
#'   horizon enter.
#'
#' Only cells scored for both models enter (the same origin, horizon,
#' stratum, age group, year and target). Scores are not negative, so the
#' loss is the score itself.
#'
#' @section Variance estimators:
#' With `d_t` the loss differences of the `T` origins (`t = 1, ..., T`, in
#' the order of the origins) and `dbar` their mean:
#'
#' * `"ewc"` (default): `Lambda_j = sqrt(2 / T) * sum_t (d_t - dbar) *
#'   cos(pi * j * (t - 1/2) / T)` for `j = 1, ..., nu`, the long-run
#'   variance `Omega = mean(Lambda_j^2)`, the statistic `sqrt(T) * dbar /
#'   sqrt(Omega)` and p-values from the t distribution with `nu` degrees of
#'   freedom, `nu = 0.4 * T^(2/3)` cosine terms as Lazarus et al. (2018)
#'   recommend for the EWC test, rounded down and at least 1 (1 for up to 11
#'   origins, 2 for 12 to 20, 3 for 21 to 31). There is no further
#'   small-sample correction. The estimate does not use the horizon, so
#'   there is a test for every horizon with at least 2 origins; `lag` is
#'   `NA`.
#' * `"daniell"`: the weighted periodogram estimate with the Daniell kernel
#'   of Coroneo and Iacone (2020, Eq. 8), `Omega = (1 / m) * sum_j I_j` for
#'   the first `m` Fourier frequencies `lambda_j = 2 * pi * j / T`, with
#'   `I_j = |sum_t d_t exp(-i lambda_j t)|^2 / T`, and p-values from the t
#'   distribution with `2 m` degrees of freedom (their Eq. 9). `m =
#'   floor(T^(1/3))`, at least 1, as in their simulations and application
#'   (Section 5: the best combination of size and power in a preliminary
#'   study): 1 for up to 7 origins, 2 for 8 to 26. `lag` is `NA`. Coroneo
#'   and Iacone (2020) also propose the Bartlett-kernel estimate with
#'   bandwidth `floor(T^(1/2))` and fixed-b critical values given as a
#'   polynomial for two quantiles, which give no p-values; it is not
#'   implemented. They do not use the equal-weighted cosine estimate.
#' * `"acf"`: the test of `forecast::dm.test(e1, e2, h = h, power = 1)` with
#'   `e1` and `e2` the losses per origin: the variance is the sum of the
#'   autocovariances up to lag `h - 1`, the statistic has the correction of
#'   Harvey et al. (1997), and p-values come from the t distribution with
#'   `T - 1` degrees of freedom. With origins spaced `s` years apart (for
#'   example `origins = seq(1981, 2001, by = 10)`), test years overlap only
#'   up to lag `ceiling(h / s) - 1`, which is used instead; with unequal
#'   spacing the smallest gap counts. When `T` is not larger than the
#'   number of lags needed (`T <= h` for consecutive origins), or when the
#'   variance estimate is not positive, there is no test (`statistic` and
#'   `p_value` are `NA`, with a warning of class `mt_warning_dm`);
#'   `forecast::dm.test()` falls back to lag 0 in the second case, which
#'   rejects a true null hypothesis about half the time.
#' * `"bartlett"`: as `"acf"` with Bartlett weights `1 - k / h` for the
#'   autocovariances at lags `k = 1, ..., h - 1`, as in `forecast::dm.test()`.
#'
#' @section Size and power with few origins:
#' Backtests have few origins, and with consecutive origins the test years
#' of horizons above 1 overlap. In a simulation under the null hypothesis
#' (20,000 series of loss differences that are sums of `h` consecutive
#' independent normal errors, a moving average of order `h - 1` as for
#' forecasts `h` years ahead from consecutive origins), the `"acf"` test at
#' the nominal 5% level rejected in 10% (`h = 2`), 16% (`h = 3`) and 19%
#' (`h = 5`) of the samples with 7 origins, and in 16% to 17% with 15 or
#' 16 origins at `h = 5`; `"bartlett"` in 10% to 15% with 7 origins. Their
#' p-values are too small with few overlapping origins. The fallback to lag
#' 0 of `forecast::dm.test()`, needed in 4% of the samples with 7 origins
#' at `h = 3` and in 22% at `h = 5`, rejected in 39% to 55% of them, so
#' that the test with the fallback rejected in 27% at `h = 5` with 7
#' origins. The EWC test rejected in 4.9% to 7.2% of the samples with 7
#' origins and in 4.8% to 7.0% with 15 or 16 origins (`h` from 1 to 5). Its
#' size is nearly correct, but with few origins it has few degrees of
#' freedom and low power: with 7 origins (`nu = 1`) it finds only large
#' and consistent differences. The `"daniell"` test of Coroneo and Iacone
#' (2020) has more degrees of freedom (2 with 7 origins, 4 with 8 to 26)
#' and rejected in 4.9%, 6.0%, 7.6%, 10.5% and 15.3% of the samples with 7
#' origins at `h` from 1 to 5, and in 4.8% to 13.8% with 15 or 16 origins:
#' its size is close to 5% for `h` up to 2, where test years overlap little,
#' and too large beyond. Coroneo and Iacone (2020) studied one-step
#' forecasts with 40 or 120 origins. More origins, or origins spaced `h`
#' years apart, which avoid overlapping test years, give more informative
#' tests.
#' The p-values are not adjusted for the number of tests (models, targets,
#' horizons); see [stats::p.adjust()].
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
#' @param varestimator `"ewc"` (equal-weighted cosine estimate with
#'   fixed-smoothing t critical values, default), `"daniell"` (Daniell-kernel
#'   periodogram estimate with fixed-smoothing t critical values, Coroneo and
#'   Iacone 2020), `"acf"` (sum of the autocovariances, as
#'   `forecast::dm.test()`) or `"bartlett"` (Bartlett weights, as
#'   `forecast::dm.test()`). See Variance estimators.
#' @return A `data.table` with one row per model (other than `ref`), value
#'   of `by`, target, metric and horizon: `model`, `ref`, the `by` columns,
#'   `target`, `metric`, `horizon` (a factor: the horizons, then `"all"`),
#'   `n_origins` (number of losses tested), `lag` (largest autocovariance
#'   lag used; `NA` for `"ewc"` and `"daniell"`), `df` (degrees of freedom of the t
#'   distribution), `score` and `score_ref` (mean loss over the origins),
#'   `mean_diff` (`score - score_ref`; negative when the model is more
#'   accurate), `statistic` and `p_value`.
#' @references
#' Coroneo L, Iacone F (2020). Comparing predictive accuracy in small
#' samples using fixed-smoothing asymptotics. *Journal of Applied
#' Econometrics* 35(4): 391-409. \doi{10.1002/jae.2756}
#'
#' Diebold FX, Mariano RS (1995). Comparing predictive accuracy. *Journal of
#' Business & Economic Statistics* 13(3): 253-263.
#'
#' Harvey D, Leybourne S, Newbold P (1997). Testing the equality of
#' prediction mean squared errors. *International Journal of Forecasting*
#' 13(2): 281-291.
#'
#' Lazarus E, Lewis DJ, Stock JH, Watson MW (2018). HAR inference:
#' recommendations for practice. *Journal of Business & Economic
#' Statistics* 36(4): 541-559.
#' @seealso [mt_compare()], [mt_score_table()], [mt_score()].
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2012)
#' bt <- mt_compare(list(rw = rwd(), lc = lee_carter()), x, origins = 2001:2007, h = 3, n_draws = 30)
#' mt_dm_test(bt)
#' # the version of forecast::dm.test()
#' mt_dm_test(bt, varestimator = "acf")
#' @export
mt_dm_test <- function(backtest, ref = NULL, metrics = "crps", targets = NULL, by = NULL,
                       alternative = c("two.sided", "less", "greater"),
                       varestimator = c("ewc", "daniell", "acf", "bartlett")) {
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
  cols <- c("model", "ref", by, "target", "metric", "horizon", "n_origins", "lag", "df", "score", "score_ref",
            "mean_diff", "statistic", "p_value")
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
# (autocovariance estimators only)
dm_series <- function(e1, e2, origin, h, alternative, varestimator) {
  o <- order(origin)
  e1 <- e1[o]
  e2 <- e2[o]
  step <- if (length(origin) > 1) min(diff(origin[o])) else 1
  r <- dm_stat(e1 - e2, max(1L, as.integer(ceiling(h / step - 1e-9))), alternative, varestimator)
  c(list(n_origins = length(e1), score = mean(e1), score_ref = mean(e2), mean_diff = mean(e1 - e2)), r)
}

# Diebold-Mariano test for the loss differences `d` (in the order of the
# origins). "ewc": equal-weighted cosine long-run variance with nu =
# max(1, floor(0.4 T^(2/3))) cosines and the t distribution with nu degrees
# of freedom (Lazarus et al. 2018). "daniell": Daniell-kernel periodogram
# estimate with m = max(1, floor(T^(1/3))) Fourier frequencies and the t
# distribution with 2 m degrees of freedom (Coroneo and Iacone 2020, Eqs.
# 8-9). "acf", "bartlett": `h` steps (lags 0 to
# h - 1) and the correction of Harvey, Leybourne and Newbold (1997), as
# forecast::dm.test(), but NA instead of lag 0 when the variance estimate is
# not positive. `fallback` records why there is no test: "origins" (too few
# losses for the lags), "variance" (variance estimate not positive) or
# "none" (fewer than 2 losses or no variation)
dm_stat <- function(d, h, alternative = "two.sided", varestimator = "ewc") {
  n <- length(d)
  na <- list(lag = NA_integer_, df = NA_real_, statistic = NA_real_, p_value = NA_real_, fallback = "none")
  if (n < 2 || anyNA(d)) return(na)
  # a constant difference leaves only rounding error in the variance
  tiny <- (1e-8 * max(abs(d)))^2
  if (varestimator == "ewc") {
    nu <- dm_ewc_df(n)
    omega <- dm_ewc_variance(d, nu)
    if (!(omega > tiny)) return(na)
    stat <- sqrt(n) * mean(d) / sqrt(omega)
    return(list(lag = NA_integer_, df = nu, statistic = stat, p_value = dm_p_value(stat, nu, alternative),
                fallback = ""))
  }
  if (varestimator == "daniell") {
    dn <- dm_daniell(d)
    if (!(dn$omega > tiny)) return(na)
    stat <- sqrt(n) * mean(d) / sqrt(dn$omega)
    return(list(lag = NA_integer_, df = dn$df, statistic = stat, p_value = dm_p_value(stat, dn$df, alternative),
                fallback = ""))
  }
  # too few origins for the autocovariances of overlapping test years: no
  # test (lag 0 would understate the variance and give too small p-values)
  if (h >= n) return(utils::modifyList(na, list(fallback = "origins")))
  v <- dm_variance(d, h, varestimator)
  if (!(v > tiny / n)) {
    no_var <- h == 1L || !(dm_variance(d, 1L, varestimator) > tiny / n)
    return(utils::modifyList(na, list(fallback = if (no_var) "none" else "variance")))
  }
  k <- sqrt((n + 1 - 2 * h + h * (h - 1) / n) / n)
  stat <- mean(d) / sqrt(v) * k
  list(lag = h - 1L, df = n - 1, statistic = stat, p_value = dm_p_value(stat, n - 1, alternative), fallback = "")
}

dm_p_value <- function(stat, df, alternative) {
  switch(alternative,
         two.sided = 2 * stats::pt(-abs(stat), df = df),
         less = stats::pt(stat, df = df),
         greater = stats::pt(stat, df = df, lower.tail = FALSE))
}

# Number of cosine terms (degrees of freedom) of the EWC test for T losses,
# nu = 0.4 T^(2/3) (Lazarus et al. 2018), at least 1
dm_ewc_df <- function(n) max(1, floor(0.4 * n^(2 / 3)))

# Equal-weighted cosine estimate of the long-run variance of `d` with `nu`
# cosine terms: the mean of Lambda_j^2, Lambda_j = sqrt(2 / T) sum_t
# (d_t - dbar) cos(pi j (t - 1/2) / T)
dm_ewc_variance <- function(d, nu) {
  n <- length(d)
  dc <- d - mean(d)
  lam <- vapply(seq_len(nu), function(j) sqrt(2 / n) * sum(dc * cos(pi * j * (seq_len(n) - 0.5) / n)), 0)
  mean(lam^2)
}

# Weighted periodogram estimate of the long-run variance of `d` with the
# Daniell kernel (Coroneo and Iacone 2020, Eq. 8): 2 pi times the mean of
# the periodogram I(lambda_j) = |sum_t d_t exp(-i lambda_j t)|^2 / (2 pi T)
# at the Fourier frequencies lambda_j = 2 pi j / T, j = 1, ..., m, with m =
# floor(T^(1/3)) (their Section 5), at least 1. Each frequency below pi
# carries 2 degrees of freedom (Eq. 9: t with 2 m); the frequency pi (T = 2)
# carries one
dm_daniell <- function(d) {
  n <- length(d)
  m <- max(1, floor(n^(1 / 3) + 1e-9))
  tt <- seq_len(n)
  dc <- d - mean(d)
  pg <- vapply(seq_len(m), function(j) {
    lam <- 2 * pi * j / n
    (sum(dc * cos(lam * tt))^2 + sum(dc * sin(lam * tt))^2) / n
  }, 0)
  list(omega = mean(pg), df = 2 * m - (2 * m >= n), m = m)
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

# One warning for the tests without a statistic
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
      sprintf(paste("No test at horizon(s) %s: too few origins for the autocovariances up to lag h - 1 of",
                    "overlapping test years. Use varestimator = \"ewc\", more origins or origins spaced h years",
                    "apart."), h_orig)
    },
    if (!is.null(h_var)) {
      sprintf(paste("No test at horizon(s) %s: the variance estimate with lags up to h - 1 was not positive",
                    "(forecast::dm.test() would use lag 0, which gives too small p-values). Use",
                    "varestimator = \"ewc\"."), h_var)
    },
    if (!is.null(h_none)) {
      sprintf("No test at horizon(s) %s: fewer than 2 origins or no variation of the loss difference.", h_none)
    })
  if (length(msg)) mt_warn(msg, class = "dm", call = NULL)
  invisible(NULL)
}
