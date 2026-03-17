#' RMSE on the log scale
#'
#' @param actual Numeric vector of observed mortality rates.
#' @param predicted Numeric vector of predicted mortality rates.
#' @return Scalar RMSE of log(actual) vs log(predicted).
#' @export
compute_rmse_log <- function(actual, predicted) {
  ok <- actual > 0 & predicted > 0 & is.finite(actual) & is.finite(predicted)
  if (sum(ok) == 0) return(NA_real_)
  sqrt(mean((log(actual[ok]) - log(predicted[ok]))^2))
}

#' Mean Absolute Error
#'
#' @param actual Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values.
#' @return Scalar MAE.
#' @export
compute_mae <- function(actual, predicted) {
  ok <- is.finite(actual) & is.finite(predicted)
  if (sum(ok) == 0) return(NA_real_)
  mean(abs(actual[ok] - predicted[ok]))
}

#' Mean Absolute Scaled Error (MASE) on log(mx)
#'
#' Compares forecast accuracy to a naive random-walk baseline.
#'
#' @param actual_test Numeric vector of observed test-period rates.
#' @param predicted_test Numeric vector of predicted test-period rates.
#' @param train_series Numeric vector of training-period rates.
#' @return Scalar MASE.
#' @export
compute_mase <- function(actual_test, predicted_test, train_series) {
  ok_fc <- actual_test > 0 & predicted_test > 0 &
           is.finite(actual_test) & is.finite(predicted_test)
  if (sum(ok_fc) == 0) return(NA_real_)
  mae_fc <- mean(abs(log(actual_test[ok_fc]) - log(predicted_test[ok_fc])))
  ok_tr  <- train_series > 0 & is.finite(train_series)
  tr_log <- log(train_series[ok_tr])
  if (length(tr_log) < 2) return(NA_real_)
  naive_mae <- mean(abs(diff(tr_log)))
  if (naive_mae < 1e-12) return(NA_real_)
  mae_fc / naive_mae
}

#' Cell-level MASE across a 3D array
#'
#' Computes [compute_mase()] for each (age, group) series and returns the mean.
#'
#' @param arr_test 3D array of observed test-period rates.
#' @param fc_arr 3D array of forecast rates (same dims as test).
#' @param arr_train 3D array of training-period rates (same ages/groups).
#' @return Scalar mean MASE.
#' @export
compute_mase_array <- function(arr_test, fc_arr, arr_train) {
  mase_vals <- c()
  for (i in seq_len(dim(arr_test)[1])) {
    for (k in seq_len(dim(arr_test)[3])) {
      m <- compute_mase(arr_test[i, , k], fc_arr[i, , k], arr_train[i, , k])
      if (!is.na(m) && is.finite(m)) mase_vals <- c(mase_vals, m)
    }
  }
  if (length(mase_vals) == 0) NA_real_ else mean(mase_vals)
}

#' Group-trend MASE
#'
#' For each group, averages mx across ages first, then computes MASE on the
#' group-level trend. Measures whether the model captures group-specific
#' temporal dynamics.
#'
#' @param arr_test 3D array of observed test-period rates.
#' @param fc_arr 3D array of forecast rates.
#' @param arr_train 3D array of training-period rates.
#' @return Scalar mean group-trend MASE.
#' @export
compute_mase_group_trend <- function(arr_test, fc_arr, arr_train) {
  n_g <- dim(arr_test)[3]
  mase_vals <- c()
  for (k in seq_len(n_g)) {
    obs_avg   <- colMeans(arr_test[, , k, drop = FALSE], na.rm = TRUE)
    fc_avg    <- colMeans(fc_arr[, , k, drop = FALSE],  na.rm = TRUE)
    train_avg <- colMeans(arr_train[, , k, drop = FALSE], na.rm = TRUE)
    m <- compute_mase(obs_avg, fc_avg, train_avg)
    if (!is.na(m) && is.finite(m)) mase_vals <- c(mase_vals, m)
  }
  if (length(mase_vals) == 0) NA_real_ else mean(mase_vals)
}

#' Combined MASE (cell + trend)
#'
#' Balances cell-level and group-trend accuracy with weights 0.7 / 0.3.
#'
#' @param arr_test 3D array of observed test-period rates.
#' @param fc_arr 3D array of forecast rates.
#' @param arr_train 3D array of training-period rates.
#' @return Scalar combined MASE.
#' @export
compute_mase_combined <- function(arr_test, fc_arr, arr_train) {
  cell_mase  <- compute_mase_array(arr_test, fc_arr, arr_train)
  trend_mase <- compute_mase_group_trend(arr_test, fc_arr, arr_train)
  if (is.na(cell_mase))  return(trend_mase)
  if (is.na(trend_mase)) return(cell_mase)
  0.7 * cell_mase + 0.3 * trend_mase
}

#' CRPS for a Gaussian predictive distribution (closed-form)
#'
#' Operates on the log scale: `obs` and `mu` should be log(mx), `sigma` is
#' the log-scale standard deviation.
#'
#' @param obs Numeric vector of observed log(mx).
#' @param mu Numeric vector of predictive mean log(mx).
#' @param sigma Numeric vector of predictive SD (log scale).
#' @return Scalar mean CRPS.
#' @export
compute_crps_gaussian <- function(obs, mu, sigma) {
  ok <- is.finite(obs) & is.finite(mu) & sigma > 0 & is.finite(sigma)
  if (sum(ok) == 0) return(NA_real_)
  z <- (obs[ok] - mu[ok]) / sigma[ok]
  crps_i <- sigma[ok] * (z * (2 * stats::pnorm(z) - 1) +
                          2 * stats::dnorm(z) - 1 / sqrt(pi))
  mean(crps_i)
}

#' Negative log-score for Gaussian predictive distribution
#'
#' @param obs Numeric vector of observations.
#' @param mu Numeric vector of predictive means.
#' @param sigma Numeric vector of predictive standard deviations.
#' @return Scalar mean negative log-score.
#' @export
compute_log_score_gaussian <- function(obs, mu, sigma) {
  ok <- is.finite(obs) & is.finite(mu) & sigma > 0 & is.finite(sigma)
  if (sum(ok) == 0) return(NA_real_)
  nll <- 0.5 * log(2 * pi) + log(sigma[ok]) +
    0.5 * ((obs[ok] - mu[ok]) / sigma[ok])^2
  mean(nll)
}

#' Interval Score
#'
#' Computes the interval score (IS) for prediction intervals.
#' The IS penalises interval width and observations falling outside the
#' interval. Operates on the log scale for mortality rates.
#'
#' @param obs Numeric vector of observed values (log scale recommended).
#' @param lower Numeric vector of lower prediction interval bounds.
#' @param upper Numeric vector of upper prediction interval bounds.
#' @param alpha Numeric; nominal miscoverage rate. Default 0.05 for a 95
#'   percent PI.
#' @return Scalar mean interval score (lower is better).
#' @references Gneiting, T. and Raftery, A.E. (2007). Strictly proper scoring
#'   rules, prediction, and estimation. JASA 102(477), 359-378.
#' @export
compute_interval_score <- function(obs, lower, upper, alpha = 0.05) {
  ok <- is.finite(obs) & is.finite(lower) & is.finite(upper)
  if (sum(ok) == 0) return(NA_real_)
  obs <- obs[ok]; lower <- lower[ok]; upper <- upper[ok]

  width   <- upper - lower
  penalty <- (2 / alpha) * pmax(lower - obs, 0) +
             (2 / alpha) * pmax(obs - upper, 0)
  mean(width + penalty)
}

#' Weighted Interval Score (WIS)
#'
#' Computes the weighted interval score across multiple quantile levels.
#' A generalisation of the interval score that evaluates the full predictive
#' distribution via nested intervals.
#'
#' @param obs Numeric vector of observed values.
#' @param quantile_matrix Numeric matrix \code{[n x K]} of quantile forecasts,
#'   where rows correspond to observations and columns to quantile levels.
#' @param probs Numeric vector of length K with the quantile probabilities
#'   (e.g., \code{c(0.025, 0.25, 0.5, 0.75, 0.975)}).
#' @return Scalar mean WIS (lower is better).
#' @export
compute_wis <- function(obs, quantile_matrix, probs) {
  stopifnot(length(probs) == ncol(quantile_matrix),
            length(obs) == nrow(quantile_matrix))
  n <- length(obs)
  K <- length(probs)

  # Find symmetric pairs
  pairs <- list()
  for (i in seq_len(K %/% 2)) {
    p_lo <- probs[i]
    p_hi <- probs[K - i + 1]
    if (abs((p_lo + p_hi) - 1) < 1e-8) {
      alpha_k <- 2 * p_lo
      pairs[[length(pairs) + 1]] <- list(
        lo_col = i, hi_col = K - i + 1, alpha = alpha_k
      )
    }
  }

  if (length(pairs) == 0) return(NA_real_)

  # Compute IS for each pair, weighted by 0.5 * alpha
  total_score <- numeric(n)
  total_weight <- 0
  for (p in pairs) {
    w <- p$alpha / 2
    is_k <- (quantile_matrix[, p$hi_col] - quantile_matrix[, p$lo_col]) +
      (2 / p$alpha) * pmax(quantile_matrix[, p$lo_col] - obs, 0) +
      (2 / p$alpha) * pmax(obs - quantile_matrix[, p$hi_col], 0)
    total_score <- total_score + w * is_k
    total_weight <- total_weight + w
  }

  # Add absolute error for median if present
  med_idx <- which(abs(probs - 0.5) < 1e-8)
  if (length(med_idx) == 1) {
    total_score <- total_score + 0.5 * abs(obs - quantile_matrix[, med_idx])
    total_weight <- total_weight + 0.5
  }

  mean(total_score / total_weight)
}

#' Interval Score for 3D mortality arrays
#'
#' Convenience wrapper that computes the mean interval score across all
#' cells of 3D mortality rate arrays on the log scale.
#'
#' @param obs_arr 3D array \code{[n_a x n_t x n_g]} of observed rates.
#' @param lower_arr 3D array of lower PI bounds (same dims).
#' @param upper_arr 3D array of upper PI bounds (same dims).
#' @param alpha Nominal miscoverage rate. Default 0.05.
#' @return Scalar mean interval score.
#' @export
compute_interval_score_array <- function(obs_arr, lower_arr, upper_arr,
                                          alpha = 0.05) {
  obs_log   <- log(pmax(obs_arr, 1e-15))
  lower_log <- log(pmax(lower_arr, 1e-15))
  upper_log <- log(pmax(upper_arr, 1e-15))
  compute_interval_score(as.vector(obs_log), as.vector(lower_log),
                          as.vector(upper_log), alpha = alpha)
}

#' Empirical CRPS from bootstrap simulations (log scale)
#'
#' Computes CRPS using the energy-score formulation on the empirical CDF
#' of bootstrap samples.
#'
#' @param obs_arr 3D array `[n_a x n_t x n_g]` of observed rates.
#' @param sims_arr 4D array `[n_a x n_t x n_g x n_sim]` of simulated rates.
#' @return Scalar mean CRPS.
#' @export
compute_crps_empirical <- function(obs_arr, sims_arr) {
  crps_vals <- c()
  n_a <- dim(obs_arr)[1]; n_t <- dim(obs_arr)[2]; n_g <- dim(obs_arr)[3]
  for (i in seq_len(n_a)) {
    for (j in seq_len(n_t)) {
      for (k in seq_len(n_g)) {
        y <- log(pmax(obs_arr[i, j, k], 1e-15))
        x <- log(pmax(sims_arr[i, j, k, ], 1e-15))
        if (!is.finite(y) || all(!is.finite(x))) next
        x <- sort(x[is.finite(x)])
        m <- length(x)
        term1 <- mean(abs(x - y))
        term2 <- sum(abs(outer(x, x, "-"))) / (m * m)
        crps_vals <- c(crps_vals, term1 - 0.5 * term2)
      }
    }
  }
  if (length(crps_vals) == 0) NA_real_ else mean(crps_vals)
}
