#' Per-age ARIMA and ETS models
#'
#' Benchmarks that model the log death rate of every age as a separate time
#' series with \pkg{forecast}: `arima_series()` selects an ARIMA model with
#' `forecast::auto.arima()` (drift allowed) or fits a given order with
#' `forecast::Arima()` (`order`), `ets_series()` an exponential smoothing
#' state space model without seasonality with `forecast::ets()` (damped
#' trends allowed). Draws are simulated future paths of the fitted models
#' (`simulate(..., future = TRUE)`), so uncertainty grows with the horizon.
#' Ages are simulated independently; parameter uncertainty is not included.
#' The models describe observed rates, so predictive draws add nothing.
#'
#' `arima_series(order = c(0, 1, 1), drift = TRUE)` is the local level model
#' with drift (a random walk with drift observed with noise; Harvey 1989) in
#' its reduced form: the year-to-year changes of a random walk with drift
#' plus independent noise are an MA(1) process with mean equal to the drift.
#' The parameters are estimated by maximum likelihood (`method = "ML"` of
#' `forecast::Arima()`, unless `...` gives another `method`). Where deaths
#' are few, it separates the noise of the observed rates from the trend,
#' which the plain random walk ([rwd()]) accumulates.
#'
#' @param ... Arguments passed to `forecast::auto.arima()`,
#'   `forecast::Arima()` (with `order`) or `forecast::ets()`.
#' @param order ARIMA order `c(p, d, q)` fitted to every age with
#'   `forecast::Arima()`; `NULL` (default) selects the order with
#'   `forecast::auto.arima()`.
#' @param drift With `order`: include a drift (a constant in the differenced
#'   series, `include.drift` of `forecast::Arima()`). Without `order`:
#'   whether `forecast::auto.arima()` may select a drift (`allowdrift`;
#'   `NULL` keeps its default).
#' @param bootstrap Simulate with resampled residuals instead of normal
#'   errors.
#' @return A model specification.
#' @references
#' Hyndman RJ, Khandakar Y (2008). Automatic time series forecasting: the
#' forecast package for R. *Journal of Statistical Software* 27(3): 1-22.
#' \doi{10.18637/jss.v027.i03}
#'
#' Harvey AC (1989). *Forecasting, Structural Time Series Models and the
#' Kalman Filter*. Cambridge University Press.
#' @examples
#' if (requireNamespace("forecast", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 2000:2019)
#'   fc <- mt_forecast(mt_fit(ets_series(), x), h = 3, n_draws = 20, seed = 1)
#'   # local level model with drift (random walk with drift plus noise)
#'   fc2 <- mt_forecast(mt_fit(arima_series(order = c(0, 1, 1), drift = TRUE), x), h = 3, n_draws = 20,
#'                      seed = 1)
#' }
#' @export
arima_series <- function(..., order = NULL, drift = NULL, bootstrap = FALSE) {
  need_pkg("forecast", "for arima_series()")
  if (!is.null(order) && (!is.numeric(order) || length(order) != 3 || anyNA(order) || any(order < 0) ||
                            any(order != round(order)))) {
    mt_abort("`order` must be an ARIMA order c(p, d, q) of non-negative whole numbers.", class = "argument")
  }
  if (!is.null(drift)) check_flag(drift, "drift")
  label <- if (is.null(order)) "ARIMA per age" else
    sprintf("ARIMA(%s)%s per age", paste(order, collapse = ","), if (isTRUE(drift)) " with drift" else "")
  new_mt_spec("series", list(method = "arima", args = list(...), order = order, drift = drift,
                             bootstrap = bootstrap),
              min_periods = 8, label = label, smooth = FALSE)
}

#' @rdname arima_series
#' @export
ets_series <- function(..., bootstrap = FALSE) {
  need_pkg("forecast", "for ets_series()")
  new_mt_spec("series", list(method = "ets", args = list(...), bootstrap = bootstrap), min_periods = 8,
              label = "ETS per age", smooth = FALSE)
}

#' @export
fit_model.mt_spec_series <- function(spec, D, E, ages, periods, ...) {
  y <- log_rates(D, E)
  p <- spec$params
  models <- lapply(seq_len(nrow(y)), function(i) {
    yi <- stats::ts(y[i, ], start = periods[1])
    if (p$method == "arima" && !is.null(p$order)) {
      do.call(forecast::Arima, c(list(yi), utils::modifyList(list(order = p$order, include.drift = isTRUE(p$drift),
                                                                  method = "ML"), p$args)))
    } else if (p$method == "arima") {
      do.call(forecast::auto.arima, c(list(yi), if (!is.null(p$drift)) list(allowdrift = p$drift), p$args))
    } else {
      do.call(forecast::ets, c(list(yi), utils::modifyList(list(model = "ZZN"), p$args)))
    }
  })
  list(models = models, bootstrap = p$bootstrap)
}

#' @export
simulate_model.mt_spec_series <- function(spec, fit, h, n_draws, ...) {
  na <- length(fit$models)
  out <- array(NA_real_, c(na, h, n_draws))
  for (i in seq_len(na)) {
    m <- fit$models[[i]]
    for (j in seq_len(n_draws)) {
      out[i, , j] <- as.numeric(stats::simulate(m, nsim = h, future = TRUE, bootstrap = fit$bootstrap))
    }
  }
  exp(out)
}
