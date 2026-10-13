#' Per-age ARIMA and ETS models
#'
#' Benchmarks that model the log death rate of every age as a separate time
#' series with \pkg{forecast}: `arima_series()` selects an ARIMA model with
#' `forecast::auto.arima()` (drift allowed), `ets_series()` an exponential
#' smoothing state space model without seasonality with `forecast::ets()`
#' (damped trends allowed). Draws are simulated future paths of the fitted
#' models (`simulate(..., future = TRUE)`), so uncertainty grows with the
#' horizon. Ages are simulated independently; parameter uncertainty is not
#' included.
#'
#' @param ... Arguments passed to `forecast::auto.arima()` or
#'   `forecast::ets()`.
#' @param bootstrap Simulate with resampled residuals instead of normal
#'   errors.
#' @return A model specification.
#' @references
#' Hyndman RJ, Khandakar Y (2008). Automatic time series forecasting: the
#' forecast package for R. *Journal of Statistical Software* 27(3): 1-22.
#' \doi{10.18637/jss.v027.i03}
#' @examples
#' if (requireNamespace("forecast", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 2000:2019)
#'   fc <- mt_forecast(mt_fit(ets_series(), x), h = 3, n_draws = 20, seed = 1)
#' }
#' @export
arima_series <- function(..., bootstrap = FALSE) {
  need_pkg("forecast", "for arima_series()")
  new_mt_spec("series", list(method = "arima", args = list(...), bootstrap = bootstrap), min_periods = 8,
              label = "ARIMA per age", smooth = FALSE)
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
    if (p$method == "arima") {
      do.call(forecast::auto.arima, c(list(yi), p$args))
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
