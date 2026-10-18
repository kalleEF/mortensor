#' Ensembles of models
#'
#' A mixture (linear pool, Gneiting and Ranjan 2013) of the predictive
#' distributions of several models: each draw of the ensemble is a draw of
#' one component, chosen with the weights as probabilities. Mixing draws
#' (rather than averaging quantiles) keeps the paths coherent, so derived
#' quantities such as life expectancy can be computed from the ensemble
#' draws.
#'
#' Weights:
#'
#' * `"equal"`;
#' * a numeric vector (normalised to sum to 1);
#' * `"crps"`: proportional to the inverse mean CRPS of each component in an
#'   inner rolling-origin backtest on the last years of the training window
#'   (`inner_origins`, `inner_h`), so the weights never use data after the
#'   forecast origin. This is a package heuristic (reference unknown),
#'   analogous to the inverse mean squared error weights of Bates and
#'   Granger (1969) with the mean CRPS in place of the mean squared error;
#' * `"stack"`: weights that minimise the squared error of the weighted mean
#'   of the components' median forecasts (log rates) in the same inner
#'   backtest, constrained to be non-negative and to sum to 1. This is
#'   stacking of point forecasts (Breiman 1996; for mortality, Kessy,
#'   Sherris, Villegas and Ziveyi 2022). Using these weights as the
#'   probabilities of the mixture is a deviation: they are optimal for the
#'   linear combination of the medians, not for the mixture distribution,
#'   whose median is not the weighted mean of the medians.
#'
#' Each draw of the ensemble comes from one component for all strata, so
#' paths keep the dependence between strata and the coherence of reconciled
#' components (`reconciled()`). A component that fails in a stratum is
#' replaced there by the other components in proportion to their weights;
#' this is recorded in the diagnostics. Any model specification can be a
#' component, including joint models, Bayesian models, tuned, bootstrapped
#' and reconciled models, and other ensembles. The components must forecast
#' the same strata: with [reconciled()] components, add the aggregates of the
#' hierarchy to the data first ([mt_aggregate_strata()]). Weights from the
#' inner backtest (`"crps"`, `"stack"`, `top_k`) compare the components on
#' the cells that all of them forecast.
#'
#' A mixture is wider than its components when they disagree, which often
#' corrects the too narrow intervals of single models; a linear pool of
#' different calibrated components is overdispersed (Gneiting and Ranjan
#' 2013). Its CRPS is never larger than the weighted mean of the CRPS of its
#' components (the CRPS is convex in the forecast distribution). The
#' vignette "Choosing and combining models" compares ensembles and single
#' models in a backtest. In populations with few deaths per cell, use
#' `rwd(observation_noise = TRUE)` as the random-walk component: the plain
#' random walk accumulates the Poisson noise of observed rates and its
#' long-horizon paths become implausible.
#'
#' @param ... Model specifications; names are used as labels.
#' @param weights `"equal"`, `"crps"`, `"stack"` or a numeric vector.
#' @param inner_origins,inner_h Inner backtest for `weights = "crps"` and
#'   `weights = "stack"`, and for `top_k`: the number of origins and the
#'   horizon of each inner forecast.
#' @param top_k Keep only the `top_k` components with the lowest CRPS in the
#'   inner backtest on the training years (the others get weight 0); `NULL`
#'   keeps all.
#' @param inner_n_draws Draws per inner forecast of the inner backtest.
#' @return A model specification.
#' @references
#' Gneiting T, Ranjan R (2013). Combining predictive distributions.
#' *Electronic Journal of Statistics* 7: 1747-1782.
#'
#' Bates JM, Granger CWJ (1969). The combination of forecasts. *Operational
#' Research Quarterly* 20(4): 451-468.
#'
#' Breiman L (1996). Stacked regressions. *Machine Learning* 24(1): 49-64.
#'
#' Kessy SR, Sherris M, Villegas AM, Ziveyi J (2022). Mortality forecasting
#' using stacked regression ensembles. *Scandinavian Actuarial Journal*
#' 2022(7): 591-626.
#' @examples
#' x <- mt_example(deaths = "single")
#' ens <- ensemble(lc = lee_carter(), rw = rwd())
#' fc <- mt_forecast(mt_fit(ens, x), h = 5, n_draws = 100, seed = 1)
#' @export
ensemble <- function(..., weights = "equal", inner_origins = 3, inner_h = 5, top_k = NULL, inner_n_draws = 100) {
  comps <- list(...)
  if (!length(comps) || !all(vapply(comps, inherits, logical(1), what = "mt_spec"))) {
    mt_abort("ensemble() needs model specifications.", class = "argument")
  }
  nm <- names(comps) %||% rep("", length(comps))
  nm[!nzchar(nm)] <- vapply(comps[!nzchar(nm)], `[[`, "", "label")
  names(comps) <- make.unique(nm)
  if (is.numeric(weights)) {
    if (length(weights) != length(comps) || any(weights < 0)) mt_abort("Invalid weights.", class = "argument")
  } else {
    check_choice(weights, c("equal", "crps", "stack"), "weights")
  }
  if (!is.null(top_k)) check_number(top_k, "top_k", lower = 1, upper = length(comps), integerish = TRUE)
  check_number(inner_origins, "inner_origins", lower = 1, integerish = TRUE)
  check_number(inner_h, "inner_h", lower = 1, integerish = TRUE)
  check_number(inner_n_draws, "inner_n_draws", lower = 10, integerish = TRUE)
  # draws describe smooth rates only when those of every component do (the
  # noise of predictive draws follows the component of each draw)
  sp <- new_mt_spec("ensemble", list(components = comps, weights = weights, inner_origins = inner_origins,
                                     inner_h = inner_h, top_k = top_k, inner_n_draws = inner_n_draws),
                    min_periods = max(vapply(comps, `[[`, 0, "min_periods")),
                    smooth = all(vapply(comps, function(s) isTRUE(s$smooth), TRUE)),
                    label = sprintf("Ensemble (%s)", paste(names(comps), collapse = ", ")))
  sp
}

#' @export
mt_fit.mt_spec_ensemble <- function(spec, data, prep = NULL, cores = 1, ...) {
  p <- spec$params
  fits <- lapply(p$components, function(cs) {
    tryCatch(mt_fit(cs, data, prep = prep, cores = cores), error = function(e) e)
  })
  ok <- !vapply(fits, inherits, logical(1), what = "error")
  if (!any(ok)) mt_abort("Every ensemble component failed.", class = "model_fit")
  ref <- fits[[which(ok)[1]]]
  same <- vapply(fits[ok], function(f) identical(f$strata$label, ref$strata$label), logical(1))
  if (!all(same)) {
    mt_abort(paste("The ensemble components forecast different strata. reconciled() adds the aggregates of its",
                   "hierarchy to data with bottom strata only; add them to the data for every component with",
                   "mt_aggregate_strata()."), class = "model_fit")
  }
  inner <- NULL
  if (!is.numeric(p$weights) && p$weights != "equal" || !is.null(p$top_k)) {
    inner <- ensemble_inner_backtests(p, data, prep)
  }
  w <- if (is.numeric(p$weights)) {
    p$weights
  } else if (p$weights == "equal") {
    rep(1, length(fits))
  } else {
    ensemble_inner_weights(p, inner)
  }
  if (!is.null(p$top_k) && !is.null(inner)) {
    # components compared on the cells that all of them forecast
    common <- ensemble_inner_common(inner)
    crps <- vapply(seq_along(inner), function(i) {
      if (is.null(common) || is.null(inner[[i]])) Inf else mean(common[[paste0("crps", i)]])
    }, 0)
    w[rank(crps, ties.method = "first") > p$top_k] <- 0
  }
  w[!ok] <- 0
  if (sum(w) == 0) w[ok] <- 1
  w <- w / sum(w)
  diag <- rbindlist(lapply(names(fits), function(n) {
    f <- fits[[n]]
    if (inherits(f, "error")) {
      data.table(component = n, stratum = ref$strata$label, status = "failed", message = conditionMessage(f),
                 seconds = NA_real_)
    } else {
      cbind(data.table(component = n), f$diagnostics)
    }
  }), fill = TRUE)
  # the data of the first component, without fields of wrappers around it
  # (calibration, anchor, tuning, ...): they do not describe the ensemble
  fields <- c("age", "periods", "strata", "omega", "observed", "train_D", "raw_D", "train_E", "exposure_last",
              "level_totals", "noise_sd")
  out <- c(list(spec = spec, prep = prep), ref[intersect(fields, names(ref))],
           list(components = fits[ok], weights = stats::setNames(w[ok], names(fits)[ok]), diagnostics = diag,
                reconcile = common_reconcile(fits[ok])))
  class(out) <- c("mt_fit_ensemble", "mt_fit")
  out
}

# The hierarchy of the components when all of them are reconciled over the
# same one: their draws, and so those of the ensemble, are then coherent, and
# wrappers around the ensemble (anchored(), recalibrated()) keep them so
common_reconcile <- function(fits) {
  rc <- lapply(fits, `[[`, "reconcile")
  if (any(vapply(rc, is.null, logical(1)))) return(NULL)
  same <- vapply(rc, function(r) {
    identical(r$hierarchy, rc[[1]]$hierarchy) && identical(r$structure$S, rc[[1]]$structure$S)
  }, logical(1))
  if (all(same)) rc[[1]] else NULL
}

ensemble_inner_backtests <- function(p, data, prep) {
  last <- max(mt_grid(data)$period$lower)
  origins <- seq(last - p$inner_origins - p$inner_h + 1, last - p$inner_h)
  lapply(p$components, function(cs) {
    suppressWarnings(tryCatch(mt_backtest(cs, data, origins = origins, h = p$inner_h, prep = prep,
                                          n_draws = p$inner_n_draws %||% 100), error = function(e) NULL))
  })
}

# Scores of the inner backtests on the cells (origin, stratum, age, year) that
# every successful one forecast, matched by these keys: columns crps<i> and
# median<i> for component i, and the observed value (NULL if there are none)
ensemble_inner_common <- function(bts) {
  keys <- c("origin", "stratum", "age_lo", "age_hi", "period", "target")
  ok <- which(!vapply(bts, is.null, logical(1)))
  if (!length(ok)) return(NULL)
  tabs <- lapply(ok, function(i) {
    b <- as.data.table(bts[[i]])[, c(keys, "observed", "crps", "median"), with = FALSE]
    setnames(b, c("crps", "median"), paste0(c("crps", "median"), i))
    if (i != ok[1]) b[, observed := NULL]
    b
  })
  common <- Reduce(function(a, b) merge(a, b, by = keys), tabs)
  if (!nrow(common)) NULL else common
}

ensemble_inner_weights <- function(p, bts) {
  ok <- !vapply(bts, is.null, logical(1))
  common <- ensemble_inner_common(bts)
  if (is.null(common)) return(rep(1, length(bts)))
  w <- numeric(length(bts))
  if (p$weights == "crps") {
    w[ok] <- vapply(which(ok), function(i) 1 / mean(common[[paste0("crps", i)]]), 0)
    return(if (all(w == 0)) rep(1, length(w)) else w)
  }
  M <- as.matrix(common[, paste0("median", which(ok)), with = FALSE])
  w[ok] <- simplex_least_squares(M, common$observed)
  w
}

# Least squares with weights on the simplex (softmax parametrisation): the
# stacking weights of the median forecasts
simplex_least_squares <- function(M, y) {
  k <- ncol(M)
  if (k == 1) return(1)
  sm <- function(z) exp(z - max(z)) / sum(exp(z - max(z)))
  f <- function(z) mean((M %*% sm(z) - y)^2)
  sm(stats::optim(rep(0, k), f, method = "BFGS")$par)
}

#' @export
mt_forecast.mt_fit_ensemble <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                                        exposure = NULL, ...) {
  last <- max(fit$periods)
  if (!is.null(to)) h <- to - last
  w <- fit$weights
  out <- with_seed(seed, {
    comp_draws <- lapply(names(fit$components), function(n) {
      mt_forecast(fit$components[[n]], h = h, n_draws = n_draws, predictive = predictive, exposure = exposure)$rates
    })
    d <- dim(comp_draws[[1]])
    arr <- array(NA_real_, d)
    # one component per draw for all strata, so that each path keeps the
    # dependence between strata (and the coherence of reconciled components)
    avail <- vapply(comp_draws, function(a) apply(a[, , , 1, drop = FALSE], 3, function(z) !anyNA(z)),
                    logical(d[3]))
    avail <- matrix(avail, d[3])
    pick <- sample(seq_along(w), n_draws, replace = TRUE, prob = w / sum(w))
    for (k in unique(pick)) {
      j <- which(pick == k)
      arr[, , , j] <- comp_draws[[k]][, , , j, drop = FALSE]
    }
    # a component that failed in a stratum is replaced there by the others
    for (s in which(rowSums(!avail) > 0)) {
      ws <- w * avail[s, ]
      if (sum(ws) == 0) next
      j <- which(!avail[s, pick])
      if (!length(j)) next
      alt <- if (sum(ws > 0) == 1) rep(which(ws > 0), length(j)) else
        sample(seq_along(ws), length(j), replace = TRUE, prob = ws / sum(ws))
      for (k in unique(alt)) arr[, , s, j[alt == k]] <- comp_draws[[k]][, , s, j[alt == k]]
    }
    list(arr = arr, pick = pick)
  })
  res <- new_mt_draws(out$arr, fit$age, last + seq_len(h), fit$strata, model = fit$spec$label, origin = last,
                      fit = fit, seed = seed)
  res$predictive <- predictive
  res$weights <- w
  # component of each draw (strata where it failed use another component)
  res$component <- out$pick
  res
}
