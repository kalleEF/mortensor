#' Rolling-origin backtest
#'
#' For each forecast origin, the data are cut at the origin, the
#' preprocessing pipeline and the model are fitted on the
#' training window only, and predictive draws for the next `h` years are
#' compared with the observations (rolling-origin evaluation; Tashman
#' 2000). Nothing after the origin enters the fit, including ungrouping and
#' smoothing (they are part of `prep`).
#'
#' Evaluation happens at the resolution of the observed deaths: grid draws of
#' deaths (Poisson given the exposures of the test years, which are treated
#' as known) are summed to each observed cell, so models on single years of
#' age can be scored against data in 5-year groups. Exposures of the test
#' years are ungrouped if needed (no deaths of the test years are used for
#' that).
#'
#' @param spec A model specification.
#' @param data An [mt_data] object.
#' @param origins Last training years. Default: the last `n_origins` years
#'   whose full horizon `h` is observed, so that every horizon is scored
#'   equally often.
#' @param h Forecast horizon.
#' @param prep Preprocessing pipeline ([mt_prep()]).
#' @param n_draws Number of draws per origin.
#' @param n_origins Number of origins when `origins` is `NULL`.
#' @param seed Random seed (forecast draws and the randomised PIT); `NULL`
#'   for none.
#' @param score_estimated Score cells whose deaths were estimated (for
#'   example suppressed cells filled by [mt_fill_suppressed()]); by default
#'   only observed cells are scored.
#' @param cores Cores used over origins (forked processes; 1 on Windows).
#' @param quantities Derived quantities scored in addition to the rates:
#'   `"life_expectancy"` (at the first age of the grid, in years) and/or
#'   `"asr"` (age-standardised rate, European Standard Population 2013, on
#'   the log scale; direct standardisation on the age groups of the
#'   standard with the exposures of the test year, see [asr()]), computed
#'   from the predictive draws of deaths aggregated to the observed age
#'   groups of each stratum and test year, like the observed values. Years
#'   whose observed cells do not cover the full age range without gaps
#'   (missing, suppressed or estimated cells) are skipped. A zero ASR (no
#'   deaths at any age) is replaced by the smallest value that half a death
#'   gives. They are added as rows with `target` `"life_expectancy"` or
#'   `"asr"`; summarise them with `mt_score(..., target = )`.
#' @param contrast Optional `c(var = , ref = )`: a stratum variable and its
#'   reference level. With `quantities`, the differences between every other
#'   level and the reference level (with the same values of the other
#'   stratum variables) are scored too, path by path: targets
#'   `"life_expectancy_gap"` (difference in years) and `"asr_gap"`
#'   (difference of log ASR, the log rate ratio). This checks forecasts of
#'   inequality, whose uncertainty depends on how the paths of the strata
#'   move together; the rows carry the stratum of the compared level.
#' @details
#' Data whose deaths were ungrouped beforehand with [mt_ungroup()], and
#' data smoothed with [mt_smooth()], are refused: the forecasts would be scored against
#' estimates instead of observations, and an ungrouping or smoothing over
#' time has seen the test years. Ungroup or smooth inside each fit instead
#' (`prep = mt_prep(step_ungroup())`, `step_smooth()`). With
#' `step_regroup()`, the test data are summed to the age groups of the fit.
#'
#' The log rate of an observed cell is `log(max(D, 0.5) / E)`: a cell
#' without deaths counts as half a death, so that its log rate is finite (a
#' convention of this package; reference unknown). The predictive draws of
#' deaths are first summed from the grid to the observed cell and then
#' floored in the same way, so that draws and observations are treated
#' alike; a floor per grid cell would put the draws of a grouped cell with
#' few deaths above the observation. The column `floored` marks the cells
#' whose observed deaths are below half a death, and [mt_score()] counts
#' them (`n_floored`). The scores of floored cells can dominate the means
#' of sparse data. The derived quantities (`quantities = "asr"` or
#' `"life_expectancy"`), which aggregate the deaths over age before
#' scoring, or coarser age groups (`step_regroup()`) avoid the floor.
#'
#' Strata without forecast draws at an origin (the model failed in them or
#' gave missing rates) are not scored, with a warning of class
#' `mt_warning_backtest`. Compare models on the cells all of them scored
#' ([mt_score_table()] does so for [mt_compare()] results).
#'
#' Scores (see [mt_score()] for their summaries): the CRPS is the
#' empirical-CDF estimator from the draws (Grimit et al. 2006), the mean
#' absolute difference between the draws and the observation minus half the
#' mean absolute difference between pairs of draws. It is the CRPS of the
#' empirical distribution of the draws, which exceeds that of the
#' predictive distribution by about `E|X - X'| / (2 n_draws)` on average,
#' so compare models backtested with the same `n_draws`. The weighted
#' interval score uses the 23 quantiles of the COVID-19 forecast hubs
#' (Bracher et al. 2021); the PIT is randomised where draws equal the
#' observation (Czado, Gneiting and Held 2009), as with zero deaths.
#' @references
#' Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
#' forecasts in an interval format. *PLoS Computational Biology* 17(2):
#' e1008618. \doi{10.1371/journal.pcbi.1008618}
#'
#' Czado C, Gneiting T, Held L (2009). Predictive model assessment for count
#' data. *Biometrics* 65(4): 1254-1261.
#'
#' Gneiting T, Raftery AE (2007). Strictly proper scoring rules, prediction,
#' and estimation. *Journal of the American Statistical Association*
#' 102(477): 359-378. \doi{10.1198/016214506000001437}
#'
#' Grimit EP, Gneiting T, Berrocal VJ, Johnson NA (2006). The continuous
#' ranked probability score for circular variables and its application to
#' mesoscale forecast ensemble verification. *Quarterly Journal of the Royal
#' Meteorological Society* 132: 2925-2942.
#'
#' Tashman LJ (2000). Out-of-sample tests of forecasting accuracy: an
#' analysis and review. *International Journal of Forecasting* 16(4):
#' 437-450.
#' @return An object of class `mt_backtest`: a `data.table` with one row per
#'   observed cell, origin and horizon: `model`, `prep` (the preprocessing
#'   steps, `"none"` without), `origin`, `horizon`,
#'   stratum columns, `age_lo`, `age_hi`, `period`, `observed` (log rate),
#'   `floored` (`TRUE` where the observed deaths are below 0.5, so that the
#'   log rate uses the floor of half a death; `NA` on rows of derived
#'   quantities), the MASE scale `scale` (mean absolute year-to-year change in the
#'   training data), `target` (`"rate"`, or the derived quantity), the
#'   forecast `median` and `mean`, and the scores `crps`, `wis`, `ae`
#'   (absolute error of the median), `cov50`, `cov80`, `cov95` and `pit`
#'   (randomised for ties), and for the opt-in summaries of [mt_score()]
#'   (`measures`) the standard deviation of the draws `sd`, `cov90` and the
#'   half-widths `hw80`, `hw90` and `hw95` of the central intervals. Rows of
#'   derived quantities have no `age_lo`, `age_hi` and `scale`. Summarise
#'   with [mt_score()].
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2012)
#' bt <- mt_backtest(rwd(), x, origins = c(2005, 2007), h = 5, n_draws = 100)
#' mt_score(bt)
#' @export
mt_backtest <- function(spec, data, origins = NULL, h = 10, prep = NULL, n_draws = 200, n_origins = 5,
                        seed = 1, cores = 1, score_estimated = FALSE, quantities = NULL, contrast = NULL) {
  backtest_impl(spec, data, origins, h, prep, n_draws, n_origins, seed, cores, score_estimated, quantities,
                contrast, call = sys.call())
}

# The backtest of mt_backtest(). With `keep_draws = TRUE` (inner backtests of
# ensemble()), the rows of rates also carry the draws of the scored log rates
# (list column `draws`) and the rows of derived quantities the quantity of the
# mean forecast rates (column `plugin`)
backtest_impl <- function(spec, data, origins = NULL, h = 10, prep = NULL, n_draws = 200, n_origins = 5,
                          seed = 1, cores = 1, score_estimated = FALSE, quantities = NULL, contrast = NULL,
                          keep_draws = FALSE, call = sys.call(-1)) {
  if (length(quantities)) {
    bad <- setdiff(quantities, c("life_expectancy", "asr"))
    if (length(bad)) {
      mt_abort("`quantities` can be \"life_expectancy\" and \"asr\".", class = "argument", call = call)
    }
  }
  check_class(data, "mt_data", "data")
  contrast <- check_backtest_contrast(contrast, data, quantities)
  if (!is.null(attr(data, "smooth"))) {
    mt_abort(paste("`data` were smoothed with mt_smooth(): the forecasts would be scored against smoothed values,",
                   "and the smoother has seen the test years. Pass the raw data and smooth inside each fit:",
                   "prep = mt_prep(step_smooth(...))."), class = "smoothed_data", call = call)
  }
  check_not_ungrouped(data, call = call)
  last <- max(mt_grid(data)$period$lower)
  origins <- origins %||% backtest_origins(data, h, n_origins, spec$min_periods, call = call)
  origins <- origins[origins < last]
  if (!length(origins)) mt_abort("No origin leaves test data.", class = "backtest", call = call)
  one <- function(o) {
    res <- tryCatch(backtest_origin(spec, data, o, h, prep, n_draws, if (is.null(seed)) NULL else seed + o,
                                    score_estimated, quantities, contrast, keep_draws), error = function(e) e)
    if (inherits(res, "error")) {
      mt_warn(sprintf("Origin %s failed: %s", o, conditionMessage(res)), class = "backtest_origin")
      return(NULL)
    }
    res
  }
  # one random number stream per origin: fits that use random numbers give the
  # same results with any number of cores
  res <- par_lapply(origins, one, cores, streams = TRUE)
  backtest_warn_skipped(res, origins)
  out <- rbindlist(res, fill = TRUE)
  if (!nrow(out)) mt_abort("All origins failed.", class = "backtest", call = call)
  setattr(out, "class", c("mt_backtest", class(out)))
  out
}

# One warning for the strata that were not scored at some origins because
# their draws were missing (attribute "skipped" of the origin results; in
# the main process, so that it is not lost with forked processes)
backtest_warn_skipped <- function(res, origins) {
  sk <- lapply(res, function(r) if (is.null(r)) character() else attr(r, "skipped"))
  hit <- lengths(sk) > 0
  if (!any(hit)) return(invisible(NULL))
  labs <- unique(unlist(sk))
  mt_warn(sprintf(paste("Strata without forecast draws were not scored (the model failed in them or gave missing",
                        "rates): %s at origin(s) %s. Compare models on the cells all of them scored."),
                  paste(utils::head(labs, 5), collapse = ", "), paste(origins[hit], collapse = ", ")),
          class = "backtest", call = NULL)
  invisible(NULL)
}

# Default origins: the last `n_origins` years whose full horizon `h` is
# observed and that leave `min_periods` training years
backtest_origins <- function(data, h, n_origins, min_periods, call = sys.call(-1)) {
  cl <- call
  yrs <- mt_grid(data)$period$lower
  last <- max(yrs)
  origins <- seq(last - h - n_origins + 1, last - h)
  first_ok <- min(yrs) + min_periods - 1
  if (any(origins < first_ok)) {
    origins <- origins[origins >= first_ok]
    if (!length(origins)) {
      mt_abort(sprintf(paste("The data are too short for %d origins with a full horizon of %d years; give",
                             "`origins` or a shorter `h`."), n_origins, h), class = "backtest", call = cl)
    }
    mt_warn(sprintf("Only %d origin(s) have a full horizon of %d years and enough training years.",
                    length(origins), h), class = "backtest_origin", call = cl)
  }
  origins
}

# Backtests score observed cells: refuse data whose deaths were ungrouped
# beforehand. Inner backtests of a fit (tuned(), ensemble()) see the data the
# fit was given
check_not_ungrouped <- function(data, call = sys.call(-1)) {
  if (!isTRUE(.mt_state$fitting) && ungrouped_deaths(data)) {
    mt_abort(paste("The deaths in `data` were ungrouped with mt_ungroup(): the forecasts would be scored against",
                   "estimated cells instead of the observed age groups, and an ungrouping over time has seen the",
                   "test years. Pass the grouped data and ungroup inside each fit: prep = mt_prep(step_ungroup())."),
             class = "ungrouped_data", call = call)
  }
}

# Were the deaths of `data` ungrouped with mt_ungroup()?
ungrouped_deaths <- function(data) {
  ug <- attr(data, "ungroup")
  !is.null(ug) && "measure" %in% names(ug) && any(ug$measure == "deaths")
}

backtest_origin <- function(spec, data, origin, h, prep, n_draws, seed, score_estimated = FALSE, quantities = NULL,
                            contrast = NULL, keep_draws = FALSE) {
  train <- mt_window(data, to = origin)
  fit <- mt_fit(spec, train, prep = prep)
  last_fit <- max(fit$periods)
  to <- min(origin + h, max(mt_grid(data)$period$lower))
  test <- mt_window(data, from = last_fit + 1, to = to)
  test <- test_on_age_axis(test, fit$age)
  E_grid <- test_exposures(test, fit)
  fc <- mt_forecast(fit, to = to, n_draws = n_draws, seed = seed, predictive = TRUE, exposure = E_grid)
  # the strata of the data, in their order (a reconciled fit of bottom-level
  # data also forecasts the aggregates)
  sf <- match(test$strata$label, fc$strata$label)
  if (anyNA(sf)) mt_abort("The forecast lacks strata of the data.", class = "backtest")
  D_draws <- fc$rates[, , sf, , drop = FALSE] * as.vector(E_grid)
  obs <- test$obs$deaths
  obs <- obs[period_lo > origin]
  if (!score_estimated) obs <- obs[provenance != "estimated"]
  gt <- mt_grid(test)
  idx <- obs_grid_index(obs, gt, "deaths", strict = TRUE)
  ok <- !is.na(idx$a1)
  obs <- obs[ok]
  idx <- idx[ok]
  na <- length(gt$age)
  np <- length(gt$period)
  ages <- c(gt$age$lower[1], gt$age$upper[na])
  probs <- wis_probs()
  rows <- list()
  qdraws <- list()
  # the randomised PIT draws uniforms: seeded like the forecast
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  if (!is.null(seed)) seed_default(seed + 1L)
  on.exit(if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  skipped <- character()
  for (s in unique(obs$stratum)) {
    if (anyNA(D_draws[, , s, ])) {
      skipped <- c(skipped, test$strata$label[s])
      next
    }
    sel <- obs$stratum == s
    C <- build_operator(idx[sel], na, np)
    Ds <- matrix(D_draws[, , s, ], ncol = n_draws)
    Ecell <- as.numeric(C %*% as.vector(E_grid[, , s]))
    # deaths summed to the observed cells, then floored like the observed deaths
    Dcell <- as.matrix(C %*% Ds)
    X <- log(pmax(Dcell, 0.5) / Ecell)
    o <- obs[sel]
    y <- log(pmax(o$value, 0.5) / Ecell)
    sc <- score_samples(X, y, probs)
    scale <- mase_scale(fit, match(test$strata$label[s], fit$strata$label), o)
    r <- cbind(
      data.table(model = spec$label, prep = prep_label(prep), origin = origin, horizon = o$period_lo - origin,
                 stratum = s, age_lo = o$age_lo, age_hi = o$age_hi, period = o$period_lo, observed = y,
                 floored = o$value < 0.5, scale = scale, target = "rate"),
      sc)
    if (keep_draws) r[, draws := lapply(seq_len(nrow(X)), function(i) X[i, ])]
    rows[[length(rows) + 1L]] <- r
    if (length(quantities)) {
      qr <- quantity_scores(Dcell, o, Ecell, quantities, probs, ages, plugin = keep_draws)
      if (!is.null(qr)) {
        rows[[length(rows) + 1L]] <- quantity_rows(qr, spec, prep, origin, s)
        qdraws[[as.character(s)]] <- qr
      }
    }
  }
  if (!is.null(contrast) && length(qdraws)) {
    for (qr in contrast_quantity_scores(qdraws, data$strata, contrast, probs)) {
      rows[[length(rows) + 1L]] <- quantity_rows(qr, spec, prep, origin, qr$stratum)
    }
  }
  if (!length(rows)) mt_abort("The model failed in every stratum.", class = "backtest")
  # rows of rates and of derived quantities differ in the columns kept for
  # inner backtests (draws, plugin)
  out <- rbindlist(rows, fill = keep_draws)
  out <- merge(data$strata, out, by = "stratum", sort = FALSE)
  if (length(skipped)) setattr(out, "skipped", skipped)
  out[]
}

# Rows of a backtest table for scored derived quantities
quantity_rows <- function(qr, spec, prep, origin, s) {
  out <- cbind(data.table(model = spec$label, prep = prep_label(prep), origin = origin, horizon = qr$period - origin,
                          stratum = s, age_lo = NA_real_, age_hi = NA_real_, period = qr$period,
                          observed = qr$observed, floored = NA, scale = NA_real_, target = qr$target),
               qr$scores)
  if (!is.null(qr$plugin)) out[, plugin := qr$plugin]
  out
}

# `contrast` of mt_backtest(): c(var = , ref = ) with a stratum variable and
# one of its levels
check_backtest_contrast <- function(contrast, data, quantities) {
  if (is.null(contrast)) return(NULL)
  if (!length(quantities)) mt_abort("`contrast` needs `quantities` to score.", class = "argument")
  if (!is.character(contrast) || !all(c("var", "ref") %in% names(contrast))) {
    mt_abort("`contrast` must be c(var = , ref = ).", class = "argument")
  }
  st <- data$strata
  if (!contrast[["var"]] %in% setdiff(names(st), c("stratum", "label"))) {
    mt_abort(sprintf("`contrast`: \"%s\" is not a stratum variable.", contrast[["var"]]), class = "argument")
  }
  if (!contrast[["ref"]] %in% st[[contrast[["var"]]]]) {
    mt_abort(sprintf("`contrast`: \"%s\" is not a level of %s.", contrast[["ref"]], contrast[["var"]]),
             class = "argument")
  }
  contrast[c("var", "ref")]
}

# Scores of the differences of derived quantities between each level of a
# stratum variable and its reference level, path by path
contrast_quantity_scores <- function(qdraws, strata, contrast, probs) {
  st <- as.data.table(strata)
  var <- contrast[["var"]]
  others <- setdiff(names(st), c("stratum", "label", var))
  key <- if (length(others)) do.call(paste, c(unname(as.list(st[, others, with = FALSE])), sep = "\r")) else
    rep("", nrow(st))
  out <- list()
  for (s in names(qdraws)) {
    i <- match(as.integer(s), st$stratum)
    if (st[[var]][i] == contrast[["ref"]]) next
    r <- st$stratum[st[[var]] == contrast[["ref"]] & key == key[i]]
    qr <- qdraws[[as.character(r[1])]]
    if (length(r) != 1 || is.null(qr)) next
    qs <- qdraws[[s]]
    # the same year, quantity and first age on both sides
    j <- match(paste(qs$period, qs$target, qs$at), paste(qr$period, qr$target, qr$at))
    k <- which(!is.na(j))
    if (!length(k)) next
    gaps <- lapply(k, function(m) {
      g <- qs$draws[[m]] - qr$draws[[j[m]]]
      go <- qs$observed[m] - qr$observed[j[m]]
      list(observed = go, scores = score_samples(matrix(g, 1), go, probs))
    })
    out[[length(out) + 1L]] <- list(period = qs$period[k], target = paste0(qs$target[k], "_gap"),
                                    observed = vapply(gaps, `[[`, 0, "observed"),
                                    scores = rbindlist(lapply(gaps, `[[`, "scores")), stratum = as.integer(s))
  }
  out
}

# Scores of derived quantities per test year from the draws of deaths Dcell
# at the observed cells o, treated like the observed deaths (no floor per
# cell). Only years whose cells cover the full age range `ages` (lowest
# age, upper end) without gaps are scored. With `plugin = TRUE`, also the
# quantity of the mean of the draws of the rates in each cell (`plugin`)
quantity_scores <- function(Dcell, o, Ecell, quantities, probs, ages, plugin = FALSE) {
  same <- function(a, b) a == b || abs(a - b) < 1e-9
  out <- list()
  for (pp in sort(unique(o$period_lo))) {
    k <- which(o$period_lo == pp)
    k <- k[order(o$age_lo[k])]
    lo <- o$age_lo[k]
    hi <- o$age_hi[k]
    if (length(k) < 2 || !same(lo[1], ages[1]) || !same(hi[length(hi)], ages[2]) ||
        any(abs(hi[-length(hi)] - lo[-1]) > 1e-9)) next
    ax <- mt_axis(lo, hi, kind = "age")
    R <- Dcell[k, , drop = FALSE] / Ecell[k]
    r_obs <- o$value[k] / Ecell[k]
    r_mean <- rowMeans(R)
    qp <- NA_real_
    for (q in quantities) {
      if (q == "life_expectancy") {
        qd <- as.vector(life_expectancy(array(R, c(length(k), 1, 1, ncol(R))), ax))
        qo <- as.vector(life_expectancy(array(r_obs, c(length(k), 1, 1)), ax))
        if (plugin) qp <- as.vector(life_expectancy(array(r_mean, c(length(k), 1, 1)), ax))
      } else {
        # direct standardisation on the age groups of the ESP 2013 with the
        # exposures of the test year
        esp <- standard_population("esp2013")
        g <- asr_groups(ax, esp)
        ig <- which(!is.na(g$group))
        Eg <- as.vector(rowsum(Ecell[k][ig], g$group[ig], reorder = TRUE))
        # no deaths at any age: the smallest ASR that half a death gives
        pos <- g$weight > 0 & Eg > 0
        floor_asr <- min(0.5 / Eg[pos] * g$weight[pos]) / sum(g$weight)
        qd <- log(pmax(asr_columns(R, ax, esp, Ecell[k], per = 1), floor_asr))
        qo <- log(max(asr_columns(matrix(r_obs), ax, esp, Ecell[k], per = 1), floor_asr))
        if (plugin) qp <- log(max(asr_columns(matrix(r_mean), ax, esp, Ecell[k], per = 1), floor_asr))
      }
      out[[length(out) + 1L]] <- list(period = pp, target = q, at = lo[1], observed = qo, draws = qd,
                                      scores = score_samples(matrix(qd, 1), qo, probs), plugin = qp)
    }
  }
  if (!length(out)) return(NULL)
  list(period = vapply(out, `[[`, 0, "period"), target = vapply(out, `[[`, "", "target"),
       at = vapply(out, `[[`, 0, "at"), observed = vapply(out, `[[`, 0, "observed"),
       draws = lapply(out, `[[`, "draws"), scores = rbindlist(lapply(out, `[[`, "scores")),
       plugin = if (plugin) vapply(out, `[[`, 0, "plugin"))
}

# The test data on the age axis of the fit. Cells finer than the age groups
# of the fit (preprocessing with step_regroup()) are summed to them; a group
# that contains an estimated cell counts as estimated
test_on_age_axis <- function(test, age) {
  res <- tryCatch(mt_set_grid(test, age = age), mt_error_grid = function(e) NULL)
  if (!is.null(res)) return(res)
  d <- test$obs$deaths
  est <- if (is.null(d)) d else d[provenance == "estimated"]
  test <- mt_regroup(test, age = age)
  if (!is.null(est) && nrow(est)) {
    ia <- axis_locate(age, est$age_lo, est$age_hi)
    hit <- paste(est$stratum, age$lower[ia], est$period_lo)
    test$obs$deaths[paste(stratum, age_lo, period_lo) %in% hit, provenance := "estimated"]
  }
  mt_set_grid(test, age = age)
}

# Exposures of the forecast periods on the grid (ungrouped if needed; deaths not used)
test_exposures <- function(test, fit) {
  t2 <- test
  t2$obs$deaths <- NULL
  g <- mt_grid(t2)
  # the exposures convert forecast rates to deaths: keep their group totals
  if (needs_ungroup(t2$obs$exposure, g)) t2 <- mt_ungroup(t2, preserve_totals = TRUE)
  E <- mt_array(t2, "exposure")
  if (anyNA(E)) mt_abort("Exposures of the test periods are incomplete.", class = "backtest")
  E
}

# Mean absolute year-to-year change of the observed log rates in the training
# data (after preprocessing), aggregated to each scored age cell
mase_scale <- function(fit, s, o) {
  cells <- unique(data.table(age_lo = o$age_lo, age_hi = o$age_hi))
  a1 <- match_bound(cells$age_lo, fit$age$lower)
  a2 <- match_bound(cells$age_hi, fit$age$upper)
  D <- (fit$raw_D %||% fit$train_D)[, , s]
  E <- fit$train_E[, , s]
  cells[, scale := vapply(seq_len(.N), function(i) {
    if (is.na(a1[i]) || is.na(a2[i])) return(NA_real_)
    r <- a1[i]:a2[i]
    lr <- log(pmax(colSums(D[r, , drop = FALSE]), 0.5) / colSums(E[r, , drop = FALSE]))
    mean(abs(diff(lr)))
  }, 0)]
  cells$scale[match(paste(o$age_lo, o$age_hi), paste(cells$age_lo, cells$age_hi))]
}

wis_probs <- function() {
  alpha <- c(0.02, 0.05, seq(0.1, 0.9, 0.1))
  sort(unique(c(alpha / 2, 0.5, 1 - alpha / 2)))
}

#' Proper scores from forecast samples
#'
#' CRPS (the empirical-CDF estimator, Grimit et al. 2006: the CRPS of the
#' empirical distribution of the samples, which depends on their number),
#' weighted interval score (Bracher et al. 2021, with the 23 quantiles of
#' the COVID-19 forecast hubs: 11 central intervals plus the median),
#' interval coverage, the PIT (randomised where samples equal the
#' observation, Czado, Gneiting and Held 2009), the absolute error of the
#' median, the standard deviation of the samples and the half-widths of the
#' central 80%, 90% and 95% intervals.
#'
#' @param X Matrix of samples (one row per target).
#' @param y Observed values.
#' @param probs Quantile levels (must contain 0.5 and symmetric pairs,
#'   among them 0.025, 0.05, 0.1, 0.25 and their complements).
#' @return A `data.table` with `median`, `mean`, `crps`, `wis`, `ae`,
#'   `cov50`, `cov80`, `cov95`, `pit`, `sd`, `cov90`, `hw80`, `hw90`,
#'   `hw95`.
#' @references
#' Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
#' forecasts in an interval format. *PLoS Computational Biology* 17(2):
#' e1008618. \doi{10.1371/journal.pcbi.1008618}
#'
#' Czado C, Gneiting T, Held L (2009). Predictive model assessment for count
#' data. *Biometrics* 65(4): 1254-1261.
#'
#' Gneiting T, Raftery AE (2007). Strictly proper scoring rules, prediction,
#' and estimation. *Journal of the American Statistical Association*
#' 102(477): 359-378. \doi{10.1198/016214506000001437}
#'
#' Grimit EP, Gneiting T, Berrocal VJ, Johnson NA (2006). The continuous
#' ranked probability score for circular variables and its application to
#' mesoscale forecast ensemble verification. *Quarterly Journal of the Royal
#' Meteorological Society* 132: 2925-2942.
#' @keywords internal
score_samples <- function(X, y, probs = wis_probs()) {
  n <- ncol(X)
  S <- t(apply(X, 1, sort))
  if (n == 1L) S <- t(S)
  Q <- row_quantiles(X, probs)
  colnames(Q) <- paste0("q", probs)
  w <- (2 * seq_len(n) - n - 1)
  crps <- rowMeans(abs(X - y)) - as.numeric(S %*% w) / n^2
  med <- Q[, which.min(abs(probs - 0.5))]
  wis <- wis_from_quantiles(Q, y, probs)
  inside <- function(lo, hi) as.numeric(y >= Q[, lo] & y <= Q[, hi])
  pq <- function(p) which.min(abs(probs - p))
  less <- rowSums(X < y)
  eq <- rowSums(X == y)
  pit <- (less + stats::runif(length(y)) * eq) / n
  half <- function(lo, hi) (Q[, hi] - Q[, lo]) / 2
  sdev <- if (n > 1) sqrt(rowSums((X - rowMeans(X))^2) / (n - 1)) else rep(NA_real_, length(y))
  out <- data.table(median = med, mean = rowMeans(X), crps = crps, wis = wis, ae = abs(med - y),
                    cov50 = inside(pq(0.25), pq(0.75)), cov80 = inside(pq(0.1), pq(0.9)),
                    cov95 = inside(pq(0.025), pq(0.975)), pit = pit,
                    sd = sdev, cov90 = inside(pq(0.05), pq(0.95)), hw80 = half(pq(0.1), pq(0.9)),
                    hw90 = half(pq(0.05), pq(0.95)), hw95 = half(pq(0.025), pq(0.975)))
  out[]
}

wis_from_quantiles <- function(Q, y, probs) {
  lower <- probs[probs < 0.5]
  K <- length(lower)
  med <- Q[, which.min(abs(probs - 0.5))]
  total <- 0.5 * abs(y - med)
  for (p in lower) {
    alpha <- 2 * p
    l <- Q[, which.min(abs(probs - p))]
    u <- Q[, which.min(abs(probs - (1 - p)))]
    is <- (u - l) + (2 / alpha) * (l - y) * (y < l) + (2 / alpha) * (y - u) * (y > u)
    total <- total + (alpha / 2) * is
  }
  total / (K + 0.5)
}

#' Summarise backtest scores
#'
#' @param backtest An [mt_backtest()] result (or several combined with
#'   `rbind()`).
#' @param by Grouping columns, for example `c("model", "horizon")`.
#' @param target `"rate"` (log death rates per cell, default),
#'   `"life_expectancy"` or `"asr"` (when the backtest scored these
#'   quantities, see `quantities` in [mt_backtest()]), or
#'   `"life_expectancy_gap"` or `"asr_gap"` (differences between strata, see
#'   `contrast` in [mt_backtest()]).
#' @param measures Optional further summaries (none by default):
#'   * `"sape"`: the mean standardised absolute predictive error of Raftery
#'     et al. (2013), `sqrt(pi / 2) |y - median| / sd` with the standard
#'     deviation `sd` of the predictive draws. It is scaled so that its
#'     expected value is 1 when the predictive distribution is correct and
#'     normal; a mean above 1 indicates intervals that are too narrow, below
#'     1 too wide. Raftery et al. define it for life expectancy, so it is
#'     available for `target = "life_expectancy"` only. Rows whose draws
#'     have no spread (`sd = 0`) are left out;
#'   * `"cov90"`: the empirical coverage of the central 90% intervals.
#'     Raftery et al. (2013) report the coverage of the 80%, 90% and 95%
#'     intervals with their mean half-widths; with `by = c("model",
#'     "horizon")` they are given by horizon;
#'   * `"halfwidth"`: the mean half-widths `hw80`, `hw90`, `hw95` of these
#'     intervals (Raftery et al. 2013, Table 1);
#'   * `"revision"`: the mean absolute revision, the absolute change of the
#'     forecast median of a cell (stratum, age group, year) between
#'     consecutive origins of the backtest, for every cell forecast from both.
#'     Stoeldraijer et al. (2013) name robustness as a purpose of a forecast:
#'     small changes in the observations when the data are updated should
#'     result in only modest changes in the forecasts. Stoeldraijer et al.
#'     (2018) measure it by the standard deviation of the changes of the
#'     forecast for one target year between successive fitting periods. The
#'     mean absolute revision over all cells forecast from consecutive
#'     origins, each of which adds a year of data with the default origins,
#'     is a package-specific summary of these changes (reference unknown).
#'     A revision belongs to the later
#'     origin and its horizon; the first origin has none. Origins further
#'     apart than one year (for example `origins = c(2005, 2007)`) give
#'     revisions over several years of data.
#' @return A `data.table` with the mean CRPS, WIS and absolute error, the
#'   root mean squared error of the median (`rmse`), MASE
#'   (mean of absolute error over the scale, the mean absolute year-to-year
#'   change of the training data, as the mean absolute scaled error of
#'   Hyndman and Koehler 2006; rates only), empirical coverage
#'   of the 50%, 80% and 95% intervals, mean PIT, the number of scored
#'   cells (`n`), the number of them whose observed deaths were floored at
#'   half a death (`n_floored`, see [mt_backtest()]; for `target = "rate"`,
#'   `NA` for other targets and for backtest tables without the column
#'   `floored`) and the `unit` of CRPS, WIS and absolute error (log rates,
#'   or years for life expectancy). The scores of floored cells can
#'   dominate the means of sparse data; the derived quantities and coarser
#'   age groups avoid the floor. Backtests of the same model with different
#'   preprocessing (column `prep`) are scored separately. With `measures`,
#'   the columns `sape`, `cov90`, `hw80`, `hw90`, `hw95` and `revision`
#'   follow.
#' @references
#' Hyndman RJ, Koehler AB (2006). Another look at measures of forecast
#' accuracy. *International Journal of Forecasting* 22(4): 679-688.
#'
#' Raftery AE, Chunn JL, Gerland P, Sevcikova H (2013). Bayesian
#' probabilistic projections of life expectancy for all countries.
#' *Demography* 50(3): 777-801. \doi{10.1007/s13524-012-0193-x}
#'
#' Stoeldraijer L, van Duin C, van Wissen L, Janssen F (2013). Impact of
#' different mortality forecasting methods and explicit assumptions on
#' projected future life expectancy: the case of the Netherlands.
#' *Demographic Research* 29(13): 323-354.
#' \doi{10.4054/DemRes.2013.29.13}
#'
#' Stoeldraijer L, van Duin C, van Wissen L, Janssen F (2018). Comparing
#' strategies for matching mortality forecasts to the most recently observed
#' data: exploring the trade-off between accuracy and robustness. *Genus*
#' 74: 16. \doi{10.1186/s41118-018-0041-y}
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2012)
#' bt <- mt_backtest(rwd(), x, origins = 2007, h = 5, n_draws = 100)
#' mt_score(bt, by = c("model", "horizon"))
#' # life expectancy: SAPE and 90% coverage by horizon (Raftery et al. 2013),
#' # revisions between consecutive origins
#' bl <- mt_backtest(rwd(), x, origins = 2004:2007, h = 5, n_draws = 100,
#'                   quantities = "life_expectancy")
#' mt_score(bl, by = c("model", "horizon"), target = "life_expectancy",
#'          measures = c("sape", "cov90", "revision"))
#' @export
mt_score <- function(backtest, by = "model",
                     target = c("rate", "life_expectancy", "asr", "life_expectancy_gap", "asr_gap"),
                     measures = NULL) {
  if (!is.data.frame(backtest)) mt_abort("`backtest` must be a backtest table.", class = "argument")
  target <- match.arg(target)
  measures <- check_score_measures(measures, target)
  bt <- as.data.table(backtest)
  if ("target" %in% names(bt)) {
    keep <- bt[["target"]] == target
    bt <- bt[keep]
  } else if (target != "rate") {
    mt_abort("The backtest has no derived quantities; use `quantities` in mt_backtest().", class = "argument")
  }
  if (!nrow(bt)) mt_abort(sprintf("The backtest holds no scores for target \"%s\".", target), class = "argument")
  rates_only <- identical(target, "rate")
  # cells whose observed deaths were floored at half a death (rates only)
  count_floored <- rates_only && "floored" %in% names(bt)
  by <- score_by(bt, by)
  out <- bt[, list(crps = mean(crps), wis = mean(wis), mae = mean(ae), rmse = sqrt(mean(ae^2)),
                   mase = if (rates_only) mean(ae / scale, na.rm = TRUE) else NA_real_,
                   cov50 = mean(cov50), cov80 = mean(cov80), cov95 = mean(cov95), pit_mean = mean(pit), n = .N,
                   n_floored = if (count_floored) sum(floored %in% TRUE) else NA_integer_),
            by = by]
  out[, unit := switch(target, rate = "log rate", life_expectancy = "years", asr = "log rate (ASR)",
                       life_expectancy_gap = "years (difference)", asr_gap = "log rate ratio (ASR)")]
  if (length(measures)) {
    ex <- score_extra(bt, by, measures)
    out <- cbind(out, ex[, setdiff(names(ex), by), with = FALSE])
  }
  out[]
}

score_measures <- c("sape", "cov90", "halfwidth", "revision")

check_score_measures <- function(measures, target) {
  if (is.null(measures)) return(character())
  if (!is.character(measures) || !all(measures %in% score_measures)) {
    mt_abort(sprintf("`measures` can be %s.", paste0("\"", score_measures, "\"", collapse = ", ")),
             class = "argument")
  }
  if ("sape" %in% measures && target != "life_expectancy") {
    mt_abort(paste("SAPE is defined for life expectancy (Raftery et al. 2013): use target = \"life_expectancy\"",
                   "(and quantities = \"life_expectancy\" in mt_backtest())."), class = "argument")
  }
  score_measures[score_measures %in% measures]
}

# The opt-in summaries of mt_score() for the rows of one target, grouped like
# the main scores (the same `by` on the same table gives the same group order)
score_extra <- function(bt, by, measures) {
  need <- c(sape = "sd", cov90 = "cov90", halfwidth = "hw80")[intersect(measures, c("sape", "cov90", "halfwidth"))]
  if (length(need) && !all(need %in% names(bt))) {
    mt_abort(sprintf("The backtest lacks the column(s) %s of `measures`; run it again with mt_backtest().",
                     paste(setdiff(need, names(bt)), collapse = ", ")), class = "argument")
  }
  ex <- data.table(.row = seq_len(nrow(bt)))
  if ("sape" %in% measures) {
    s <- bt[["sd"]]
    ex[, sape_ := ifelse(is.finite(s) & s > 0, sqrt(pi / 2) * bt[["ae"]] / s, NA_real_)]
  }
  if ("revision" %in% measures) ex[, revision_ := backtest_revisions(bt)]
  for (v in c(by, intersect(c("cov90", "hw80", "hw90", "hw95"), names(bt)))) set(ex, j = v, value = bt[[v]])
  na_mean <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  ex[, c(if ("sape" %in% measures) list(sape = na_mean(sape_)),
         if ("cov90" %in% measures) list(cov90 = mean(cov90)),
         if ("halfwidth" %in% measures) list(hw80 = mean(hw80), hw90 = mean(hw90), hw95 = mean(hw95)),
         if ("revision" %in% measures) list(revision = na_mean(abs(revision_)))),
     by = by]
}

# Revision of each row: its forecast median minus the median for the same
# cell (model, pipeline, stratum, age group, year, target) at the previous
# origin of the backtest of that model; NA at the first origin and for
# cells that the previous origin did not score
backtest_revisions <- function(bt) {
  key <- function(cols) {
    cols <- intersect(cols, names(bt))
    if (!length(cols)) return(rep("", nrow(bt)))
    do.call(paste, c(unname(as.list(bt[, cols, with = FALSE])), sep = "\r"))
  }
  unit <- key(c("model", "prep"))
  cell <- key(c("model", "prep", "stratum", "age_lo", "age_hi", "period", "target"))
  o_rank <- integer(nrow(bt))
  for (u in unique(unit)) {
    i <- which(unit == u)
    o_rank[i] <- match(bt$origin[i], sort(unique(bt$origin[i])))
  }
  j <- match(paste(cell, o_rank - 1L, sep = "\r"), paste(cell, o_rank, sep = "\r"))
  bt$median - bt$median[j]
}

# Grouping columns of the scores: the same model with different
# preprocessing is scored separately
score_by <- function(bt, by) {
  if ("model" %in% by && "prep" %in% names(bt) && !"prep" %in% by &&
      bt[, data.table::uniqueN(prep), by = model][, any(V1 > 1)]) {
    by <- c(by, "prep")
  }
  by
}

# Short description of a preprocessing pipeline for backtest tables
prep_label <- function(prep) {
  if (is.null(prep) || !length(prep$steps)) return("none")
  paste(vapply(prep$steps, function(st) {
    a <- st$args
    det <- c(a$along, a$method, if (!is.null(a$dims) && a$dims == 2) "2D")
    if (length(det)) sprintf("%s(%s)", st$type, paste(det, collapse = ", ")) else st$type
  }, ""), collapse = " + ")
}

#' @export
print.mt_backtest <- function(x, ...) {
  if (!all(c("model", "origin", "crps", "cov95") %in% names(x))) {
    # a derived table (for example x[, .N, by = horizon]) keeps the class
    print(unclass_backtest(x), ...)
    return(invisible(x))
  }
  cat(sprintf("<mt_backtest> %s; origins %s; %d scored cells\n", paste(unique(x$model), collapse = ", "),
              paste(unique(x$origin), collapse = ", "), nrow(x)))
  targets <- if ("target" %in% names(x)) unique(x$target) else "rate"
  print(mt_score(x, target = if ("rate" %in% targets) "rate" else targets[1]))
  if (length(targets) > 1) {
    cat(sprintf("Other targets: %s (mt_score(target = ))\n", paste(setdiff(targets, "rate"), collapse = ", ")))
  }
  invisible(x)
}

unclass_backtest <- function(x) {
  y <- copy(x)
  setattr(y, "class", setdiff(class(y), "mt_backtest"))
  y
}
