#' Tucker3-CoDa forecast
#'
#' Forecasts mortality by applying Tucker3 tensor decomposition to the
#' Centered Log-Ratio (CLR) transform of the death distribution.
#' Pipeline: mx -> dx -> CLR -> Tucker3 -> forecast B-mode -> inv-CLR -> mx.
#'
#' @param mx_array 3D array `[n_a x n_t x n_g]` of mortality rates.
#' @param years Integer vector of calendar years.
#' @param ages Integer vector of age-group start values.
#' @param h Integer forecast horizon.
#' @param R_A,R_B,R_C Integer Tucker3 ranks for age, time, group modes.
#' @param ts_method Time-series method for B-mode: `"ets"`, `"arima"`, or
#'   `"rwd"` (random walk with drift).
#' @param fast Logical; use faster but approximate fitting?
#' @param return_full Logical; if TRUE, return decomposition details.
#' @param width Integer age-group width (default 5).
#' @return A list with `mx_forecast` `[n_a x h x n_g]` and `fit_pct`.
#'   If `return_full = TRUE`, also includes decomposition factors, residuals, etc.
#' @export
tucker3_coda_forecast <- function(mx_array, years, ages, h,
                                  R_A = 3, R_B = 3, R_C = 2,
                                  ts_method = "ets", fast = TRUE,
                                  return_full = FALSE, width = 5) {
  n_a <- dim(mx_array)[1]
  n_t <- dim(mx_array)[2]
  n_g <- dim(mx_array)[3]
  R_A <- min(R_A, n_a - 1); R_B <- min(R_B, n_t - 1); R_C <- min(R_C, n_g - 1)

  # mx -> dx -> CLR
  X_clr <- array(NA_real_, dim = c(n_a, n_t, n_g))
  for (j in seq_len(n_t))
    for (k in seq_len(n_g)) {
      dx <- compute_dx(mx_array[, j, k], ages, width)
      X_clr[, j, k] <- as.numeric(compositions::clr(dx))
    }

  # Tucker3 decomposition
  decomp <- rTensor::tucker(rTensor::as.tensor(X_clr),
                            ranks = c(R_A, R_B, R_C),
                            max_iter = ifelse(fast, 50, 200),
                            tol = ifelse(fast, 1e-6, 1e-10))

  A     <- decomp$U[[1]]
  B     <- decomp$U[[2]]
  C     <- decomp$U[[3]]
  G_mat <- rTensor::k_unfold(decomp$Z, m = 1)@data

  # Forecast B-mode columns
  bm <- forecast_B_mode(B, years, h, ts_method, fast)
  B_fc <- bm$B_fc

  B_full     <- rbind(B, B_fc)
  X_clr_full <- reconstruct_array(A, B_full, C, G_mat)

  # inv-CLR -> dx -> mx
  mx_fc <- array(NA_real_, dim = c(n_a, h, n_g))
  for (j in seq_len(h))
    for (k in seq_len(n_g)) {
      dx_fc <- inv_clr_to_dx(X_clr_full[, n_t + j, k])
      mx_fc[, j, k] <- dx_to_mx(dx_fc, ages, width)
    }

  result <- list(mx_forecast = mx_fc, fit_pct = decomp$norm_percent)
  if (return_full) {
    result$decomp      <- decomp
    result$A <- A; result$B <- B; result$C <- C; result$G_mat <- G_mat
    result$B_full      <- B_full
    result$X_clr       <- X_clr
    result$X_clr_full  <- X_clr_full
    result$resid_array <- X_clr - reconstruct_array(A, B, C, G_mat)
    result$n_t         <- n_t
    result$arima_fits  <- bm$arima_fits
  }
  result
}

#' Tucker3-Direct forecast
#'
#' Forecasts mortality by applying Tucker3 decomposition directly to log(mx).
#' Pipeline: log(mx) -> Tucker3 -> forecast B-mode -> exp -> mx.
#'
#' @inheritParams tucker3_coda_forecast
#' @return A list with `mx_forecast` and `fit_pct`.
#' @export
tucker3_direct_forecast <- function(mx_array, years, h,
                                    R_A = 3, R_B = 3, R_C = 2,
                                    ts_method = "ets", fast = TRUE,
                                    return_full = FALSE) {
  n_a <- dim(mx_array)[1]
  n_t <- dim(mx_array)[2]
  n_g <- dim(mx_array)[3]
  R_A <- min(R_A, n_a - 1); R_B <- min(R_B, n_t - 1); R_C <- min(R_C, n_g - 1)

  X_log  <- log(pmax(mx_array, 1e-12))
  decomp <- rTensor::tucker(rTensor::as.tensor(X_log),
                            ranks = c(R_A, R_B, R_C),
                            max_iter = ifelse(fast, 50, 200),
                            tol = ifelse(fast, 1e-6, 1e-10))

  A     <- decomp$U[[1]]
  B     <- decomp$U[[2]]
  C     <- decomp$U[[3]]
  G_mat <- rTensor::k_unfold(decomp$Z, m = 1)@data

  bm <- forecast_B_mode(B, years, h, ts_method, fast)
  B_fc <- bm$B_fc

  B_full     <- rbind(B, B_fc)
  X_log_full <- reconstruct_array(A, B_full, C, G_mat)

  mx_fc <- array(NA_real_, dim = c(n_a, h, n_g))
  for (j in seq_len(h))
    for (k in seq_len(n_g))
      mx_fc[, j, k] <- exp(X_log_full[, n_t + j, k])

  result <- list(mx_forecast = mx_fc, fit_pct = decomp$norm_percent)
  if (return_full) {
    result$decomp      <- decomp
    result$A <- A; result$B <- B; result$C <- C; result$G_mat <- G_mat
    result$B_full      <- B_full
    result$X_log       <- X_log
    result$X_log_full  <- X_log_full
    result$resid_array <- X_log - reconstruct_array(A, B, C, G_mat)
    result$n_t         <- n_t
    result$arima_fits  <- bm$arima_fits
  }
  result
}


# ---- Internal: B-mode time-series forecasting ----

#' Forecast all columns of the B factor matrix
#' @keywords internal
forecast_B_mode <- function(B, years, h, ts_method = "ets", fast = TRUE) {
  R_B <- ncol(B)
  arima_fits <- vector("list", R_B)
  B_fc <- matrix(NA_real_, h, R_B)
  for (r in seq_len(R_B)) {
    ts_r <- stats::ts(B[, r], start = min(years), frequency = 1)
    fit_r <- tryCatch({
      if (ts_method == "arima")
        forecast::auto.arima(ts_r, stepwise = fast, approximation = fast)
      else if (ts_method == "ets")
        forecast::ets(ts_r, model = "ZZN")
      else
        forecast::Arima(ts_r, order = c(0, 1, 0), include.drift = TRUE)
    }, error = function(e) tryCatch(
      forecast::Arima(ts_r, order = c(0, 1, 0)),
      error = function(e2) forecast::Arima(ts_r, order = c(0, 0, 0),
                                           include.mean = TRUE)))
    arima_fits[[r]] <- fit_r
    B_fc[, r] <- as.numeric(forecast::forecast(fit_r, h = h)$mean)
  }
  list(B_fc = B_fc, arima_fits = arima_fits)
}
