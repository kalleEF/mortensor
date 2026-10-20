#' Ensembles of models
#'
#' A mixture (linear pool, Gneiting and Ranjan 2013) of the predictive
#' distributions of several models: each draw of the ensemble is a draw of
#' one component, chosen with the weights as probabilities. Mixing draws
#' keeps the paths coherent, so derived
#' quantities such as life expectancy can be computed from the ensemble
#' draws.
#'
#' Weights:
#'
#' * `"equal"` (default);
#' * a numeric vector (normalised to sum to 1);
#' * `"crps_stack"`: stacking of predictive distributions (Yao, Vehtari,
#'   Simpson and Gelman 2018, Section 2.1, eq. 2.2): the weights on the
#'   simplex that optimise a proper scoring rule of the mixture, averaged
#'   over held-out predictions. The scoring rule is the CRPS. Yao et al.
#'   define stacking for any proper scoring rule (their Section 2.1 lists
#'   the CRPS among them) and use the logarithmic score, which needs
#'   predictive densities; the CRPS, which draws give directly, is a choice
#'   of this package. The CRPS is proper, and strictly proper for
#'   distributions with a finite mean (Gneiting and Raftery 2007, Section
#'   4.2). The held-out
#'   predictions are the forecasts of the inner rolling-origin backtest on
#'   the training years (leave-future-out validation, as in Barigou et al.
#'   2023) in place of the leave-one-out predictions of Yao et al.; they are
#'   scored on the observed log death rates, as in [mt_backtest()]. With
#'   draws `X_k` of component `k`, the CRPS of the mixture at an
#'   observation `y` is `sum_k w_k E|X_k - y| - 1/2 sum_k sum_l w_k w_l
#'   E|X_k - X_l'|` with `X_k` and `X_l'` independent (the kernel form of
#'   Gneiting and Raftery 2007, eq. 21, applied to the mixture). The
#'   expectations are means over all draws and all pairs of draws, so this
#'   is the CRPS of the pooled draws with weight `w_k / n` for each of the
#'   `n` draws of component `k`, the estimator of the backtest. Its mean
#'   over the scored cells is a quadratic function of the weights, convex
#'   on the simplex (the CRPS is the integral of a squared difference of
#'   distribution functions, Gneiting and Raftery 2007, eq. 20). The exact
#'   minimum of this quadratic program is found by an active-set method.
#'   Weights can be 0, and the minimum need not be unique (for example with
#'   identical components). Yao et al. (2018, Section 2.3) note that the
#'   weights can be unstable when the held-out sample is small; the cells
#'   of an inner backtest with few origins are correlated;
#' * `"bias"`: the weights of Kontis et al. (2017, Methods, "Bayesian model
#'   averaging", steps 1 to 3). They withhold the last 13 years of data, fit
#'   each model to the earlier years and project the withheld years. The
#'   projection bias of a model is the observed minus the projected life
#'   expectancy, averaged over the withheld years. A model gets the weight
#'   `exp(-|bias|)`, and the weights are normalised to sum to 1. The final
#'   projections pool draws of the models fitted to all data in proportion
#'   to the weights. Here the withheld years are those of the inner
#'   backtest with one origin 13 years before the last training year
#'   (`inner_origins = 1`, `inner_h = 13`, the defaults with `weights =
#'   "bias"`). The projected life expectancy is that of the mean forecast
#'   rates in each age group and year, here the means of the predictive
#'   draws of the inner backtest (the code published with the paper uses
#'   the life expectancy of the posterior mean rates). The observed life
#'   expectancy is that of the observed rates. Both are computed by
#'   [life_expectancy()] in the observed age groups at the first age of the
#'   data, so at birth only when the data start at age 0 (the code published
#'   with the paper uses life expectancy at birth). Withheld years whose
#'   observed age groups do
#'   not cover the age range are skipped, as in [mt_backtest()].
#'   Package-specific (reference unknown): with several strata, one weight
#'   per component from the mean over the strata of `|bias|` (the paper
#'   weights each country and sex separately; fit one ensemble per
#'   population with [stratified()] to do so); with `inner_origins > 1`, the
#'   mean of `|bias|` over the origins. The draws of the ensemble come from
#'   the components with the weights as probabilities; the paper takes
#'   numbers of draws proportional to the weights;
#' * `"crps"`: proportional to the inverse mean CRPS of each component in
#'   the inner backtest. Shang (2012, Section 3.2.1) weights point
#'   forecasts of age-specific life expectancies in inverse proportion to
#'   their mean absolute forecast errors in a validation period, following Bates and
#'   Granger (1969), and normalises the weights to sum to 1. He combines
#'   all models, or the two best ones, this way and compares these
#'   combinations with equal and Bayesian model averaging weights (Tables
#'   18 and 19). The CRPS of a point forecast is its absolute error
#'   (Gneiting and Raftery 2007, Section 4.2), so the inverse mean CRPS
#'   extends these weights to predictive distributions. This extension, the
#'   log rates as the scored quantity and the use of the weights as mixture
#'   probabilities are package-specific. `top_k = 2` keeps the two
#'   components with the lowest inner CRPS, as Shang's combination of the
#'   two best models;
#' * `"stack"` (deprecated; use `"crps_stack"`): weights that minimise the
#'   squared error of the weighted mean of the components' median forecasts
#'   (log rates) in the same inner backtest, constrained to be non-negative
#'   and to sum to 1. This is stacking of point forecasts (Breiman 1996;
#'   for mortality, Kessy, Sherris, Villegas and Ziveyi 2022). Using these
#'   weights as the probabilities of the mixture is a deviation: they are
#'   optimal for the linear combination of the medians, not for the mixture
#'   distribution, whose median is not the weighted mean of the medians. It
#'   still works and warns (class `mt_warning_deprecated`).
#'
#' Each draw of the ensemble comes from one component for all strata, so
#' paths keep the dependence between strata and the coherence of reconciled
#' components (`reconciled()`). A component that fails in a stratum is
#' replaced there by the other components in proportion to their weights;
#' this is recorded in the diagnostics. Any model specification can be a
#' component, including joint models, Bayesian models, tuned, bootstrapped
#' and reconciled models, and other ensembles. The components must forecast
#' the same strata: with [reconciled()] components, add the aggregates of the
#' hierarchy to the data first ([mt_aggregate_strata()]).
#'
#' Weights from the inner backtest (`"crps_stack"`, `"bias"`, `"crps"`,
#' `"stack"`, `top_k`) use only the training data, so they never use data
#' after the forecast origin, also inside an outer backtest. They compare
#' the components on the cells that all of them forecast. A component whose
#' inner backtest fails gets weight 0; when every inner backtest fails, the
#' weights are equal, with a warning of class `mt_warning_ensemble_weights`.
#' With `top_k`, the stacking weights of `"crps_stack"` are computed among
#' the kept components; with the other weights, the dropped components get
#' weight 0.
#'
#' A mixture is wider than its components when they disagree, which often
#' corrects the too narrow intervals of single models; a linear pool of
#' different calibrated components is overdispersed (Gneiting and Ranjan
#' 2013). Its CRPS is never larger than the weighted mean of the CRPS of its
#' components (the CRPS is convex in the forecast distribution). The
#' vignette "Choosing and combining models" compares ensembles and single
#' models in a backtest. In populations with few deaths per cell, the plain
#' random walk accumulates the Poisson noise of observed rates and its
#' long-horizon paths become implausible; the local level model with drift,
#' the random walk observed with noise, is `arima_series(order = c(0, 1, 1),
#' drift = TRUE)` (its reduced form, Harvey 1989).
#'
#' @param ... Model specifications; names are used as labels.
#' @param weights `"equal"`, `"crps_stack"`, `"bias"`, `"crps"`, the
#'   deprecated `"stack"`, or a numeric vector.
#' @param inner_origins,inner_h Inner backtest for the weights other than
#'   `"equal"` and numeric ones, and for `top_k`: the number of origins and
#'   the horizon of each inner forecast. The last origin is `inner_h` years
#'   before the last training year. Defaults: 3 origins and 5 years; with
#'   `weights = "bias"`, one origin and 13 years (Kontis et al. 2017).
#' @param top_k Keep only the `top_k` components with the lowest CRPS in the
#'   inner backtest on the training years (the others get weight 0); `NULL`
#'   keeps all.
#' @param inner_n_draws Draws per inner forecast of the inner backtest.
#' @return A model specification.
#' @references
#' Barigou K, Goffard P-O, Loisel S, Salhi Y (2023). Bayesian model
#' averaging for mortality forecasting using leave-future-out validation.
#' *International Journal of Forecasting* 39(2): 674-690.
#' \doi{10.1016/j.ijforecast.2022.01.011}
#'
#' Bates JM, Granger CWJ (1969). The combination of forecasts. *Operational
#' Research Quarterly* 20(4): 451-468.
#'
#' Breiman L (1996). Stacked regressions. *Machine Learning* 24(1): 49-64.
#'
#' Gneiting T, Raftery AE (2007). Strictly proper scoring rules, prediction,
#' and estimation. *Journal of the American Statistical Association*
#' 102(477): 359-378. \doi{10.1198/016214506000001437}
#'
#' Gneiting T, Ranjan R (2013). Combining predictive distributions.
#' *Electronic Journal of Statistics* 7: 1747-1782.
#'
#' Harvey AC (1989). *Forecasting, Structural Time Series Models and the
#' Kalman Filter*. Cambridge University Press.
#'
#' Kessy SR, Sherris M, Villegas AM, Ziveyi J (2022). Mortality forecasting
#' using stacked regression ensembles. *Scandinavian Actuarial Journal*
#' 2022(7): 591-626.
#'
#' Kontis V, Bennett JE, Mathers CD, Li G, Foreman K, Ezzati M (2017).
#' Future life expectancy in 35 industrialised countries: projections with a
#' Bayesian model ensemble. *The Lancet* 389(10076): 1323-1335.
#' \doi{10.1016/S0140-6736(16)32381-9}
#'
#' Shang HL (2012). Point and interval forecasts of age-specific life
#' expectancies: a model averaging approach. *Demographic Research* 27(21):
#' 593-644. \doi{10.4054/DemRes.2012.27.21}
#'
#' Yao Y, Vehtari A, Simpson D, Gelman A (2018). Using stacking to average
#' Bayesian predictive distributions (with discussion). *Bayesian Analysis*
#' 13(3): 917-1007. \doi{10.1214/17-BA1091}
#' @examples
#' x <- mt_example(deaths = "single")
#' ens <- ensemble(lc = lee_carter(), rw = rwd())
#' fc <- mt_forecast(mt_fit(ens, x), h = 5, n_draws = 100, seed = 1)
#' \donttest{
#' # stacking weights from an inner backtest on the training years
#' st <- ensemble(lc = lee_carter(), rw = rwd(), weights = "crps_stack",
#'                inner_origins = 2, inner_h = 3)
#' mt_fit(st, x)$weights
#' }
#' @export
ensemble <- function(..., weights = "equal", inner_origins = NULL, inner_h = NULL, top_k = NULL,
                     inner_n_draws = 100) {
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
    check_choice(weights, c("equal", "crps_stack", "bias", "crps", "stack"), "weights")
    if (weights == "stack") {
      mt_warn(paste("ensemble(weights = \"stack\") is deprecated: its weights optimise the weighted mean of the",
                    "median forecasts and are used as mixture probabilities. Use weights = \"crps_stack\" (stacking",
                    "of the predictive distributions)."), class = "deprecated", call = NULL)
    }
  }
  # Kontis et al. (2017) withhold the last 13 years
  bias <- identical(weights, "bias")
  inner_origins <- inner_origins %||% if (bias) 1 else 3
  inner_h <- inner_h %||% if (bias) 13 else 5
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
    ensemble_inner_check(inner, p, names(p$components))
  }
  # top_k: the components with the lowest mean CRPS on the cells that all of
  # them forecast
  keep <- rep(TRUE, length(fits))
  if (!is.null(p$top_k) && !is.null(inner)) {
    common <- ensemble_inner_common(inner)
    crps <- vapply(seq_along(inner), function(i) {
      if (is.null(common) || is.null(inner[[i]])) Inf else mean(common[[paste0("crps", i)]])
    }, 0)
    keep <- rank(crps, ties.method = "first") <= p$top_k
  }
  w <- if (is.numeric(p$weights)) {
    p$weights
  } else if (p$weights == "equal") {
    rep(1, length(fits))
  } else {
    ensemble_inner_weights(p, inner, keep)
  }
  w[!keep] <- 0
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
              "level_totals")
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

# Inner rolling-origin backtests of the components on the training data:
# `inner_origins` origins whose horizon `inner_h` ends in the last training
# year. Stacking keeps the draws of the scored log rates; the weights of
# Kontis et al. (2017) score life expectancy and keep its value at the mean
# forecast rates
ensemble_inner_backtests <- function(p, data, prep) {
  last <- max(mt_grid(data)$period$lower)
  origins <- seq(last - p$inner_origins - p$inner_h + 1, last - p$inner_h)
  n_draws <- p$inner_n_draws %||% 100
  keep <- identical(p$weights, "crps_stack") || identical(p$weights, "bias")
  quantities <- if (identical(p$weights, "bias")) "life_expectancy"
  lapply(p$components, function(cs) {
    suppressWarnings(tryCatch({
      if (keep) {
        backtest_impl(cs, data, origins = origins, h = p$inner_h, prep = prep, n_draws = n_draws,
                      quantities = quantities, keep_draws = TRUE)
      } else {
        mt_backtest(cs, data, origins = origins, h = p$inner_h, prep = prep, n_draws = n_draws)
      }
    }, error = function(e) NULL))
  })
}

# Warns when the inner backtest of a component failed or scored fewer
# origins than requested (too few training years at the first inner origin,
# or failed fits), so that the weights rest on fewer forecasts than set
ensemble_inner_check <- function(inner, p, labels) {
  n_req <- p$inner_origins
  got <- vapply(inner, function(b) if (is.null(b)) 0L else data.table::uniqueN(as.data.table(b)$origin), 0L)
  short <- which(got < n_req)
  if (length(short)) {
    mt_warn(sprintf(paste("The inner backtests for the ensemble weights scored fewer origins than",
                          "`inner_origins` = %d for %s; the weights rest on fewer forecasts. Use more training",
                          "years or smaller `inner_origins` and `inner_h`."), n_req,
                    paste(sprintf("%s (%d)", labels[short], got[short]), collapse = ", ")),
            class = "ensemble_weights")
  }
  invisible(got)
}

# Rows of the inner backtests on the cells (origin, stratum, age, year) that
# every successful one forecast, matched by these keys, for the target `tgt`
# ("rate", or a derived quantity matched without ages): columns crps<i>,
# median<i> and, when kept, draws<i> and plugin<i> for component i, and the
# observed value (NULL if there are none)
ensemble_inner_common <- function(bts, tgt = "rate") {
  keys <- c("origin", "stratum", if (tgt == "rate") c("age_lo", "age_hi"), "period", "target")
  ok <- which(!vapply(bts, is.null, logical(1)))
  if (!length(ok)) return(NULL)
  tabs <- lapply(ok, function(i) {
    b <- as.data.table(bts[[i]])
    b <- b[target == tgt]
    cols <- c("crps", "median", intersect(c("draws", "plugin"), names(b)))
    b <- b[, c(keys, "observed", cols), with = FALSE]
    setnames(b, cols, paste0(cols, i))
    if (i != ok[1]) b[, observed := NULL]
    b
  })
  common <- Reduce(function(a, b) merge(a, b, by = keys), tabs)
  if (!nrow(common)) NULL else common
}

# Weights from the inner backtests (0 for components whose inner backtest
# failed); `keep` are the components kept by top_k, among which the stacking
# weights are computed
ensemble_inner_weights <- function(p, bts, keep = NULL) {
  if (identical(p$weights, "crps_stack") && !is.null(keep)) bts[!keep] <- list(NULL)
  ok <- !vapply(bts, is.null, logical(1))
  common <- ensemble_inner_common(bts, if (identical(p$weights, "bias")) "life_expectancy" else "rate")
  if (is.null(common)) {
    mt_warn(sprintf(paste("The inner backtests for weights = \"%s\" scored no cells that every component",
                          "forecast (too few training years for inner_origins = %d and inner_h = %d, or",
                          "failed fits); the weights are equal."), p$weights, p$inner_origins %||% NA_integer_,
                    p$inner_h %||% NA_integer_), class = "ensemble_weights", call = NULL)
    return(rep(1, length(bts)))
  }
  w <- numeric(length(bts))
  idx <- which(ok)
  if (p$weights == "crps") {
    w[ok] <- vapply(idx, function(i) 1 / mean(common[[paste0("crps", i)]]), 0)
    return(if (all(w == 0)) rep(1, length(w)) else w)
  }
  if (p$weights == "crps_stack") {
    X <- lapply(idx, function(i) do.call(rbind, common[[paste0("draws", i)]]))
    tm <- crps_mixture_terms(X, common$observed)
    w[ok] <- simplex_qp(-tm$B, tm$a)
    return(w)
  }
  if (p$weights == "bias") {
    w[ok] <- vapply(idx, function(i) exp(-mean_abs_bias(common, paste0("plugin", i))), 0)
    return(w)
  }
  M <- as.matrix(common[, paste0("median", idx), with = FALSE])
  w[ok] <- simplex_least_squares(M, common$observed)
  w
}

# Kontis et al. (2017), step 1: the projection bias, observed minus projected
# life expectancy averaged over the withheld years, of each origin and
# stratum; the mean of its absolute value over origins and strata
# (package-specific when there are several). Inf when there is none
mean_abs_bias <- function(common, col) {
  b <- tapply(common$observed - common[[col]], list(common$origin, common$stratum), mean)
  b <- abs(b[is.finite(b)])
  if (!length(b)) Inf else mean(b)
}

# Terms of the CRPS of a mixture of the empirical distributions of draws
# (kernel form, Gneiting and Raftery 2007, eq. 21): for weights w on the
# simplex, the mean over cells of the CRPS of the mixture is
# sum(a * w) - t(w) %*% B %*% w / 2, with a[k] the mean of |X_k - y| and
# B[k, l] the mean of |X_k - X_l'| over all draws and pairs of draws. `X`:
# list of matrices of draws (cells x draws), `y`: observations
crps_mixture_terms <- function(X, y) {
  K <- length(X)
  a <- vapply(X, function(x) mean(abs(x - y)), 0)
  self <- lapply(X, pair_abs_sums)
  B <- matrix(0, K, K)
  for (k in seq_len(K)) {
    B[k, k] <- mean(self[[k]]) / ncol(X[[k]])^2
    for (l in seq_len(k - 1)) {
      # pairs within the pooled draws minus pairs within each component
      pooled <- pair_abs_sums(cbind(X[[k]], X[[l]]))
      B[k, l] <- B[l, k] <- mean(pooled - self[[k]] - self[[l]]) / (2 * ncol(X[[k]]) * ncol(X[[l]]))
    }
  }
  list(a = a, B = B)
}

# Sums of |x_i - x_j| over all ordered pairs of the values in each row, from
# the order statistics: 2 sum_i (2 i - n - 1) x_(i)
pair_abs_sums <- function(x) {
  n <- ncol(x)
  s <- if (n == 1L) x else t(apply(x, 1, sort))
  2 * as.vector(s %*% (2 * seq_len(n) - n - 1))
}

# Minimum of t(w) %*% Q %*% w / 2 + sum(q * w) over the simplex (w >= 0,
# sum(w) = 1), for Q positive semidefinite on the directions that sum to 0,
# so that the problem is convex (Q itself need not be): a primal active-set
# method over the faces of the simplex, starting from the best vertex. On
# each face, the step minimises the objective in an orthonormal basis of
# the face; directions of zero curvature that decrease the objective are
# followed to the boundary, so that singular problems (identical
# components) end at one of their minima
simplex_qp <- function(Q, q, max_iter = 100L * length(q)) {
  K <- length(q)
  if (K == 1L) return(1)
  Q <- (Q + t(Q)) / 2
  tol <- 1e-10 * max(abs(Q), abs(q), .Machine$double.xmin)
  w <- numeric(K)
  free <- which.min(diag(Q) / 2 + q)
  w[free] <- 1
  for (it in seq_len(max_iter)) {
    g <- as.vector(Q %*% w) + q
    step <- simplex_face_step(Q, g, free, tol)
    if (is.null(step)) {
      # stationary on the face: multipliers of the bounds w_j >= 0
      mu <- g - mean(g[free])
      mu[free] <- Inf
      if (min(mu) >= -tol) return(w)
      free <- sort(c(free, which.min(mu)))
      next
    }
    d <- numeric(K)
    d[free] <- step$d
    neg <- which(d < 0)
    ratio <- -w[neg] / d[neg]
    r <- if (length(neg)) min(ratio) else Inf
    alpha <- if (step$newton) min(1, r) else r
    if (!is.finite(alpha)) break
    w <- w + alpha * d
    if (alpha == r) {
      hit <- neg[which.min(ratio)]
      w[hit] <- 0
      free <- setdiff(free, hit)
    }
    w <- pmax(w, 0)
    w <- w / sum(w)
  }
  mt_warn("The quadratic program of the stacking weights did not converge.", class = "ensemble_weights",
          call = NULL)
  w
}

# Step on the face of the simplex with the variables `free`: the minimum of
# the quadratic on the face (`newton`), or a descent direction of zero
# curvature; NULL when the point is stationary on the face
simplex_face_step <- function(Q, g, free, tol) {
  nf <- length(free)
  if (nf < 2L) return(NULL)
  Z <- qr.Q(qr(matrix(1, nf, 1)), complete = TRUE)[, -1, drop = FALSE]
  H <- crossprod(Z, Q[free, free, drop = FALSE] %*% Z)
  r <- as.vector(crossprod(Z, g[free]))
  e <- eigen((H + t(H)) / 2, symmetric = TRUE)
  pos <- e$values > tol
  V <- e$vectors
  if (any(!pos)) {
    V0 <- V[, !pos, drop = FALSE]
    r0 <- as.vector(V0 %*% crossprod(V0, r))
    if (sqrt(sum(r0^2)) > tol) return(list(d = -as.vector(Z %*% r0), newton = FALSE))
  }
  if (!any(pos)) return(NULL)
  Vp <- V[, pos, drop = FALSE]
  d <- -as.vector(Z %*% (Vp %*% (crossprod(Vp, r) / e$values[pos])))
  if (max(abs(d)) <= 1e-12) return(NULL)
  list(d = d, newton = TRUE)
}

# Least squares with weights on the simplex (softmax parametrisation): the
# stacking weights of the median forecasts (deprecated weights = "stack")
simplex_least_squares <- function(M, y) {
  k <- ncol(M)
  if (k == 1) return(1)
  sm <- function(z) exp(z - max(z)) / sum(exp(z - max(z)))
  f <- function(z) mean((M %*% sm(z) - y)^2)
  sm(stats::optim(rep(0, k), f, method = "BFGS")$par)
}

# The weights of an ensemble in one sentence (model_assumptions())
ensemble_weights_text <- function(p) {
  w <- p$weights
  inner <- sprintf("an inner backtest on the training years (%d %s, horizon %d years)", p$inner_origins,
                   if (isTRUE(p$inner_origins == 1)) "origin" else "origins", p$inner_h)
  txt <- if (is.numeric(w)) {
    "Weights: given."
  } else {
    switch(w,
           equal = "Weights: equal.",
           crps_stack = paste("Weights: stacking of the predictive distributions (Yao et al. 2018) with the CRPS:",
                              "the weights minimise the mean CRPS of the mixture in", paste0(inner, ".")),
           bias = paste("Weights: exp(-|bias|) of the life expectancy at the first age projected for the",
                        "withheld years of", paste0(inner, " (Kontis et al. 2017).")),
           crps = paste("Weights: inverse mean CRPS in", inner, "(inverse past errors, Shang 2012; for",
                        "predictive distributions package-specific)."),
           stack = paste("Weights: stacking of the median forecasts in", inner, "(deprecated)."))
  }
  if (!is.null(p$top_k)) txt <- paste(txt, sprintf("Only the %d components with the lowest inner CRPS are kept.",
                                                   p$top_k))
  txt
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
