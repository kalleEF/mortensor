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
#'   use data after the forecast origin.
#'
#' A component that fails in a stratum is dropped for that stratum and the
#' remaining weights are renormalised; this is recorded in the diagnostics.
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
    check_choice(weights, c("equal", "crps"), "weights")
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
    ensemble_crps_weights(p, data, prep)
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

ensemble_crps_weights <- function(p, data, prep) {
  last <- max(mt_grid(data)$period$lower)
  origins <- seq(last - p$inner_origins - p$inner_h + 1, last - p$inner_h)
  crps <- vapply(p$components, function(cs) {
    bt <- tryCatch(mt_backtest(cs, data, origins = origins, h = p$inner_h, prep = prep, n_draws = 100),
                   error = function(e) NULL)
    if (is.null(bt)) NA_real_ else mean(bt$crps)
  }, 0)
  w <- 1 / crps
  w[!is.finite(w)] <- 0
  if (all(w == 0)) w[] <- 1
  w
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
