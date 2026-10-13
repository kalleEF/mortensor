#' Ensembles of models
#'
#' A mixture of the predictive distributions of several models: each draw of
#' the ensemble is a draw of one component, chosen in proportion to the
#' weights. Mixing draws (rather than averaging quantiles) keeps the paths
#' coherent, so derived quantities such as life expectancy can be computed
#' from the ensemble draws.
#'
#' Weights:
#'
#' * `"equal"`;
#' * a numeric vector (normalised to sum to 1);
#' * `"crps"`: inverse mean CRPS of an inner rolling-origin backtest on the
#'   last `inner_origins` years of the training window, so the weights never
#'   use data after the forecast origin;
#' * `"stack"`: weights on the simplex that minimise the squared error of the
#'   weighted median forecasts (log scale) in the same inner backtest (the
#'   non-negative stacking of mortensor 0.1, fitted on training data only).
#'
#' A component that fails in a stratum is dropped for that stratum and the
#' remaining weights are renormalised; this is recorded in the diagnostics.
#'
#' A mixture is wider than its components when they disagree. In backtests
#' on HMD data (England & Wales, France) and on the ONS deprivation data,
#' Lee-Carter alone was too narrow (95% coverage 0.70 to 0.84), while an
#' equal-weight mixture of `lee_carter()` and `rwd()` had 95% coverage of
#' 0.94 to 0.98 and the best or second-best CRPS (`dev/calibration.md` in
#' the source repository).
#'
#' @param ... Model specifications; names are used as labels.
#' @param weights `"equal"`, `"crps"` or a numeric vector.
#' @param inner_origins,inner_h Inner backtest for `weights = "crps"`.
#' @return A model specification.
#' @examples
#' x <- mt_example(deaths = "single")
#' ens <- mt_ensemble(lc = lee_carter(), rw = rwd())
#' fc <- mt_forecast(mt_fit(ens, x), h = 5, n_draws = 100, seed = 1)
#' @export
mt_ensemble <- function(..., weights = "equal", inner_origins = 3, inner_h = 5) {
  comps <- list(...)
  if (!length(comps) || !all(vapply(comps, inherits, logical(1), what = "mt_spec"))) {
    mt_abort("mt_ensemble() needs model specifications.", class = "argument")
  }
  nm <- names(comps) %||% rep("", length(comps))
  nm[!nzchar(nm)] <- vapply(comps[!nzchar(nm)], `[[`, "", "label")
  names(comps) <- make.unique(nm)
  if (is.numeric(weights)) {
    if (length(weights) != length(comps) || any(weights < 0)) mt_abort("Invalid weights.", class = "argument")
  } else {
    check_choice(weights, c("equal", "crps", "stack"), "weights")
  }
  sp <- new_mt_spec("ensemble", list(components = comps, weights = weights, inner_origins = inner_origins,
                                     inner_h = inner_h),
                    min_periods = max(vapply(comps, `[[`, 0, "min_periods")),
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
  w <- if (is.numeric(p$weights)) {
    p$weights
  } else if (p$weights == "equal") {
    rep(1, length(fits))
  } else {
    ensemble_inner_weights(p, data, prep)
  }
  w[!ok] <- 0
  w <- w / sum(w)
  ref <- fits[[which(ok)[1]]]
  diag <- rbindlist(lapply(names(fits), function(n) {
    f <- fits[[n]]
    if (inherits(f, "error")) {
      data.table(component = n, stratum = ref$strata$label, status = "failed", message = conditionMessage(f),
                 seconds = NA_real_)
    } else {
      cbind(data.table(component = n), f$diagnostics)
    }
  }), fill = TRUE)
  out <- ref
  out$spec <- spec
  out$prep <- prep
  out$components <- fits[ok]
  out$weights <- stats::setNames(w[ok], names(fits)[ok])
  out$diagnostics <- diag
  class(out) <- c("mt_fit_ensemble", "mt_fit")
  out
}

ensemble_inner_weights <- function(p, data, prep) {
  last <- max(mt_grid(data)$period$lower)
  origins <- seq(last - p$inner_origins - p$inner_h + 1, last - p$inner_h)
  bts <- lapply(p$components, function(cs) {
    suppressWarnings(tryCatch(mt_backtest(cs, data, origins = origins, h = p$inner_h, prep = prep, n_draws = 100),
                              error = function(e) NULL))
  })
  ok <- !vapply(bts, is.null, logical(1))
  if (!any(ok)) return(rep(1, length(bts)))
  if (p$weights == "crps") {
    w <- vapply(bts, function(b) if (is.null(b)) 0 else 1 / mean(b$crps), 0)
    return(if (all(w == 0)) rep(1, length(w)) else w)
  }
  M <- do.call(cbind, lapply(bts[ok], `[[`, "median"))
  y <- bts[[which(ok)[1]]]$observed
  w <- numeric(length(bts))
  w[ok] <- simplex_least_squares(M, y)
  w
}

# Least squares with weights on the simplex (softmax parametrisation)
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
    for (s in seq_len(d[3])) {
      avail <- vapply(comp_draws, function(a) !anyNA(a[, , s, 1]), logical(1))
      ws <- w * avail
      if (sum(ws) == 0) next
      pick <- sample(seq_along(ws), n_draws, replace = TRUE, prob = ws / sum(ws))
      for (k in unique(pick)) {
        j <- which(pick == k)
        arr[, , s, j] <- comp_draws[[k]][, , s, j]
      }
    }
    arr
  })
  res <- new_mt_draws(out, fit$age, last + seq_len(h), fit$strata, model = fit$spec$label, origin = last,
                      fit = fit, seed = seed)
  res$predictive <- predictive
  res$weights <- w
  res
}
