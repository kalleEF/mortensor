#' Random Walk with Drift forecast (per-series)
#'
#' Applies a random walk with drift model to each (age, group) series of
#' log(mx) independently.
#'
#' @param mx_array 3D array `[n_a x n_t x n_g]` of mortality rates.
#' @param years Integer vector of calendar years.
#' @param h Integer forecast horizon.
#' @param drift_window Integer; if > 0, use only the last `drift_window` years
#'   to estimate the drift (default 0 = all years).
#' @param return_fitted Logical; return fitted values?
#' @return A list with `mx_forecast` `[n_a x h x n_g]` and optionally
#'   `mx_fitted`.
#' @export
rwd_forecast <- function(mx_array, years, h, drift_window = 0L,
                         return_fitted = FALSE) {
  n_a <- dim(mx_array)[1]; n_t <- dim(mx_array)[2]; n_g <- dim(mx_array)[3]
  rwd_arr <- array(NA_real_, dim = c(n_a, h, n_g))
  rwd_fit <- if (return_fitted) array(NA_real_, dim = c(n_a, n_t, n_g)) else NULL

  for (i in seq_len(n_a)) {
    for (k in seq_len(n_g)) {
      ls <- log(pmax(mx_array[i, , k], 1e-12))
      win <- if (drift_window > 0 && length(ls) > drift_window)
        utils::tail(ls, drift_window) else ls
      drift  <- mean(diff(win))
      last_v <- utils::tail(ls, 1)
      for (j in seq_len(h))
        rwd_arr[i, j, k] <- exp(last_v + drift * j)
      if (return_fitted) rwd_fit[i, , k] <- mx_array[i, , k]
    }
  }
  list(mx_forecast = rwd_arr, mx_fitted = rwd_fit)
}

#' Per-series ARIMA forecast
#'
#' Fits an auto.arima model to each (age, group) series of log(mx).
#'
#' @param mx_array 3D array of mortality rates.
#' @param years Integer vector of calendar years.
#' @param h Integer forecast horizon.
#' @param max.p,max.q,max.d Maximum ARIMA orders.
#' @param return_fitted Logical; return fitted values?
#' @return A list with `mx_forecast` and optionally `mx_fitted`.
#' @export
arima_forecast <- function(mx_array, years, h,
                           max.p = 2L, max.q = 1L, max.d = 1L,
                           return_fitted = FALSE) {
  n_a <- dim(mx_array)[1]; n_t <- dim(mx_array)[2]; n_g <- dim(mx_array)[3]
  as_arr <- array(NA_real_, dim = c(n_a, h, n_g))
  as_fit <- if (return_fitted) array(NA_real_, dim = c(n_a, n_t, n_g)) else NULL

  for (i in seq_len(n_a)) {
    for (k in seq_len(n_g)) {
      ls <- log(pmax(mx_array[i, , k], 1e-12))
      ts_r <- stats::ts(ls, start = min(years), frequency = 1)
      fit_r <- safe_auto_arima(ts_r, max.p = max.p, max.q = max.q, max.d = max.d)
      fc_r <- forecast::forecast(fit_r, h = h)
      as_arr[i, , k] <- exp(as.numeric(fc_r$mean))
      if (return_fitted) as_fit[i, , k] <- exp(as.numeric(stats::fitted(fit_r)))
    }
  }
  list(mx_forecast = as_arr, mx_fitted = as_fit)
}

#' Damped ETS forecast (per-series)
#'
#' Fits an ETS (exponential smoothing) model to each (age, group) series of
#' log(mx) with optional damping.
#'
#' @param mx_array 3D array of mortality rates.
#' @param years Integer vector of calendar years.
#' @param h Integer forecast horizon.
#' @param ets_model Character ETS model specification (default `"AAN"`).
#' @param damped Logical; apply trend damping?
#' @param return_fitted Logical; return fitted values?
#' @return A list with `mx_forecast` and optionally `mx_fitted`.
#' @export
ets_forecast <- function(mx_array, years, h,
                         ets_model = "AAN", damped = TRUE,
                         return_fitted = FALSE) {
  n_a <- dim(mx_array)[1]; n_t <- dim(mx_array)[2]; n_g <- dim(mx_array)[3]
  ed_arr <- array(NA_real_, dim = c(n_a, h, n_g))
  ed_fit <- if (return_fitted) array(NA_real_, dim = c(n_a, n_t, n_g)) else NULL

  for (i in seq_len(n_a)) {
    for (k in seq_len(n_g)) {
      ls <- log(pmax(mx_array[i, , k], 1e-12))
      ts_r <- stats::ts(ls, start = min(years), frequency = 1)
      fit_r <- tryCatch(
        forecast::ets(ts_r, model = ets_model, damped = damped),
        error = function(e) safe_auto_arima(ts_r))
      fc_r <- forecast::forecast(fit_r, h = h)
      ed_arr[i, , k] <- exp(as.numeric(fc_r$mean))
      if (return_fitted) ed_fit[i, , k] <- exp(as.numeric(stats::fitted(fit_r)))
    }
  }
  list(mx_forecast = ed_arr, mx_fitted = ed_fit)
}


#' Safe ARIMA fitting with multiple fallbacks
#'
#' Tries `auto.arima` first, then falls back to `Arima(0,1,0)`, then to a
#' trivial constant model `Arima(0,0,0)`.
#'
#' @param ts_r A `ts` object.
#' @param stepwise Logical; stepwise search?
#' @param approximation Logical; use approximations?
#' @param max.p,max.q,max.d Maximum ARIMA orders.
#' @return An ARIMA model object.
#' @keywords internal
safe_auto_arima <- function(ts_r, stepwise = TRUE, approximation = TRUE,
                            max.p = 2, max.q = 1, max.d = 1) {
  tryCatch(
    forecast::auto.arima(ts_r, stepwise = stepwise,
                         approximation = approximation,
                         max.p = max.p, max.q = max.q, max.d = max.d),
    error = function(e1) tryCatch(
      forecast::Arima(ts_r, order = c(0, 1, 0)),
      error = function(e2) {
        forecast::Arima(ts_r, order = c(0, 0, 0), include.mean = TRUE)
      }
    )
  )
}
