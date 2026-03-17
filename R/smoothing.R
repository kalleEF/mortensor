#' Smooth mortality data
#'
#' Applies either a 2D tensor-product P-spline (via [mgcv::gam()]) or a 1D
#' Whittaker smoother to a `mort_data` object, smoothing log(mx) over age and
#' year (or year only for Whittaker).
#'
#' @param x A `mort_data` object, or a 3D array `[n_a x n_t x n_g]`.
#' @param method One of `"gam"` (2D P-spline, recommended) or `"whittaker"`
#'   (1D time-only fallback).
#' @param ages Integer vector of ages (required if `x` is a raw array).
#' @param years Integer vector of years (required if `x` is a raw array).
#' @param covid_years Integer vector of years to exclude (weight = 0).
#' @param lambda Optional smoothing parameter for Whittaker (auto-selected
#'   if NULL).
#' @param floor_threshold Mortality rates below this are treated as zero-count.
#' @return If `x` is `mort_data`, a new `mort_data` with `smoothed = TRUE`.
#'   If `x` is a raw array, a list with `arr` (smoothed array) and `method`.
#' @export
smooth_mortality <- function(x,
                             method = c("gam", "whittaker"),
                             ages = NULL, years = NULL,
                             covid_years = integer(0),
                             lambda = NULL,
                             floor_threshold = 1e-7) {
  method <- match.arg(method)
  if (inherits(x, "mort_data")) {
    result <- smooth_mortality(x$mx, method = method,
                               ages = x$ages, years = x$years,
                               covid_years = covid_years,
                               lambda = lambda,
                               floor_threshold = floor_threshold)
    x$mx <- result$arr
    x$smoothed <- TRUE
    x$gam_fits <- result$gam_fits  # store GAM fits for gam_bayes PIs (NULL if whittaker)
    return(x)
  }

  stopifnot(is.array(x), length(dim(x)) == 3)
  stopifnot(!is.null(ages), !is.null(years))

  if (method == "gam") {
    smooth_2d_array(x, ages, years, covid_years, floor_threshold)
  } else {
    smooth_time_array(x, years, covid_years, lambda, floor_threshold)
  }
}

#' Check smoothing quality
#'
#' Diagnostic checks comparing raw and smoothed arrays: MAD, max deviation,
#' trend preservation, and correlation.
#'
#' @param mx_raw 3D array of raw mortality rates.
#' @param mx_smooth 3D array of smoothed mortality rates.
#' @param years Integer vector of calendar years.
#' @param covid_years Years to exclude from comparison.
#' @param mad_threshold MAD threshold for over-smoothing warning.
#' @param floor_threshold Cells below this are excluded.
#' @return A list with `mad`, `max_dev`, `trend_ratio`, `correlation`, `pass`.
#' @export
check_smoothing <- function(mx_raw, mx_smooth, years,
                            covid_years = c(2020, 2021),
                            mad_threshold = 0.2,
                            floor_threshold = 1e-7) {
  obs_idx <- which(!years %in% covid_years)
  raw_obs    <- log(pmax(mx_raw[, obs_idx, ], 1e-10))
  smooth_obs <- log(pmax(mx_smooth[, obs_idx, ], 1e-10))

  real_mask <- mx_raw[, obs_idx, ] >= floor_threshold
  if (sum(real_mask) < 10) {
    return(list(mad = NA, max_dev = NA, trend_ratio = NA,
                correlation = NA, pass = TRUE))
  }

  deviations <- abs(raw_obs - smooth_obs)
  mad_val <- stats::median(deviations[real_mask])
  max_dev <- max(deviations[real_mask])

  n_a <- dim(mx_raw)[1]; n_g <- dim(mx_raw)[3]
  ratios <- c()
  for (i in seq_len(n_a)) {
    for (k in seq_len(n_g)) {
      if (any(mx_raw[i, obs_idx, k] < floor_threshold)) next
      r <- raw_obs[i, , k]; s <- smooth_obs[i, , k]
      if (all(is.finite(r)) && all(is.finite(s))) {
        v_raw <- stats::var(diff(r))
        if (v_raw > 1e-20) ratios <- c(ratios, stats::var(diff(s)) / v_raw)
      }
    }
  }
  mean_trend_ratio <- if (length(ratios) > 0) mean(ratios) else NA
  corr_raw_smooth  <- stats::cor(raw_obs[real_mask], smooth_obs[real_mask])

  list(mad = mad_val, max_dev = max_dev,
       trend_ratio = mean_trend_ratio, correlation = corr_raw_smooth,
       pass = mad_val <= mad_threshold &&
              !is.na(mean_trend_ratio) && mean_trend_ratio >= 0.15)
}


# ---- Internal smoothing functions ----

#' @keywords internal
whittaker_smooth <- function(y, lambda, w = rep(1, length(y)), d = 2L) {
  n <- length(y)
  E <- diag(n)
  D <- diff(E, differences = d)
  W <- diag(w)
  z <- solve(W + lambda * crossprod(D), W %*% y)
  as.numeric(z)
}

#' @keywords internal
select_lambda_robust <- function(log_mx_array, w_vec, d = 2L,
                                 lambdas = 10^seq(-1, 4, by = 0.25),
                                 min_trend_ratio = 0.15,
                                 floor_threshold = log(1e-7)) {
  n <- length(w_vec)
  obs_idx <- which(w_vec > 0)
  E <- diag(n)
  D_pen <- diff(E, differences = d)
  W <- diag(w_vec)

  n_lam        <- length(lambdas)
  loo_scores   <- rep(Inf, n_lam)
  trend_ratios <- rep(0,   n_lam)

  for (il in seq_len(n_lam)) {
    lam <- lambdas[il]
    H   <- solve(W + lam * crossprod(D_pen)) %*% W
    h_diag <- diag(H)
    denom_loo <- 1 - h_diag[obs_idx]
    if (any(abs(denom_loo) < 1e-10)) next
    total_loo <- 0; total_trend <- 0; n_series <- 0
    for (i in seq_len(dim(log_mx_array)[1])) {
      for (k in seq_len(dim(log_mx_array)[3])) {
        y <- log_mx_array[i, , k]
        if (!all(is.finite(y))) next
        if (any(y < floor_threshold)) next
        z <- as.numeric(H %*% y)
        loo_resid <- (y[obs_idx] - z[obs_idx]) / denom_loo
        total_loo <- total_loo + mean(loo_resid^2)
        d1_raw <- diff(y[obs_idx]); d1_smo <- diff(z[obs_idx])
        v_raw  <- stats::var(d1_raw)
        total_trend <- total_trend + if (v_raw > 1e-20) stats::var(d1_smo) / v_raw else 1
        n_series <- n_series + 1
      }
    }
    if (n_series > 0) {
      loo_scores[il]   <- total_loo   / n_series
      trend_ratios[il] <- total_trend / n_series
    }
  }
  loo_idx <- which.min(loo_scores)
  trend_ok <- which(trend_ratios >= min_trend_ratio)
  trend_idx <- if (length(trend_ok) > 0) max(trend_ok) else 1L
  best_idx <- min(loo_idx, trend_idx)
  chosen_by <- if (loo_idx <= trend_idx) "LOO-CV" else "trend-cap"

  list(lambda = lambdas[best_idx], chosen_by = chosen_by,
       loo_optimal = lambdas[loo_idx], trend_cap = lambdas[trend_idx],
       trend_ratio = trend_ratios[best_idx])
}

#' @keywords internal
smooth_time_array <- function(mx_array, years,
                              covid_years = c(2020L, 2021L),
                              lambda = NULL,
                              floor_threshold = 1e-7) {
  n_a <- dim(mx_array)[1]; n_t <- dim(mx_array)[2]; n_g <- dim(mx_array)[3]
  log_arr <- log(pmax(mx_array, 1e-10))
  w_base <- rep(1, n_t); w_base[years %in% covid_years] <- 0

  if (is.null(lambda)) {
    sel <- select_lambda_robust(log_arr, w_base,
                                floor_threshold = log(floor_threshold))
    lambda <- sel$lambda
  }

  smooth_arr <- mx_array
  for (i in seq_len(n_a)) {
    for (k in seq_len(n_g)) {
      y <- log_arr[i, , k]
      if (!all(is.finite(y))) next
      w_cell <- w_base
      is_floor <- mx_array[i, , k] < floor_threshold
      if (any(is_floor)) w_cell[is_floor] <- 0
      if (sum(w_cell > 0) < 3L) next
      z <- whittaker_smooth(y, lambda, w_cell)
      smooth_arr[i, , k] <- exp(z)
    }
  }
  list(arr = smooth_arr, lambda = lambda, method = "whittaker")
}

#' @keywords internal
smooth_2d_array <- function(mx_array, ages, years,
                            covid_years = integer(0),
                            floor_threshold = 1e-7) {
  n_a <- dim(mx_array)[1]; n_t <- dim(mx_array)[2]; n_g <- dim(mx_array)[3]
  smooth_arr <- mx_array
  gam_fits <- vector("list", n_g)

  grid <- expand.grid(age = ages, year = years)
  k_age  <- min(10, n_a - 1)
  k_year <- min(10, n_t - 1)

  for (k in seq_len(n_g)) {
    log_mx <- as.vector(log(pmax(mx_array[, , k], 1e-10)))
    w_vec <- rep(1, n_a * n_t)
    for (j in seq_len(n_t)) {
      if (years[j] %in% covid_years)
        w_vec[((j - 1) * n_a + 1):(j * n_a)] <- 0
    }
    is_floor <- as.vector(mx_array[, , k]) < floor_threshold
    w_vec[is_floor] <- 0

    df <- data.frame(age = grid$age, year = grid$year,
                     log_mx = log_mx, w = w_vec)

    gam_fit <- tryCatch({
      mgcv::gam(log_mx ~ te(age, year, bs = "ps", k = c(k_age, k_year)),
                data = df, weights = w, method = "REML")
    }, error = function(e) NULL)

    if (!is.null(gam_fit)) {
      pred <- stats::predict(gam_fit, newdata = df)
      smooth_arr[, , k] <- matrix(exp(pred), nrow = n_a, ncol = n_t)
      gam_fits[[k]] <- gam_fit
    }
  }
  list(arr = smooth_arr, method = "gam", gam_fits = gam_fits)
}

#' Compute pointwise Bayesian 95% credible intervals from GAM posteriors
#'
#' Uses `predict(gam_fit, se.fit = TRUE)` for efficient closed-form CIs.
#'
#' @param gam_fits List of [mgcv::gam] objects (one per group).
#' @param ages     Numeric vector of ages.
#' @param years    Integer vector of years.
#' @return List with `mx_lo95` and `mx_hi95` (3D arrays \[n_a x n_t x n_g\]).
#' @export
compute_gam_ci_arrays <- function(gam_fits, ages, years) {
  n_a <- length(ages); n_t <- length(years); n_g <- length(gam_fits)
  mx_lo <- array(NA_real_, dim = c(n_a, n_t, n_g))
  mx_hi <- array(NA_real_, dim = c(n_a, n_t, n_g))
  grid_df <- expand.grid(age = ages, year = years)
  for (k in seq_len(n_g)) {
    if (!is.null(gam_fits[[k]])) {
      pred <- stats::predict(gam_fits[[k]], newdata = grid_df, se.fit = TRUE)
      log_lo <- pred$fit - 1.96 * pred$se.fit
      log_hi <- pred$fit + 1.96 * pred$se.fit
      mx_lo[, , k] <- matrix(exp(log_lo), nrow = n_a, ncol = n_t)
      mx_hi[, , k] <- matrix(exp(log_hi), nrow = n_a, ncol = n_t)
    }
  }
  list(mx_lo95 = mx_lo, mx_hi95 = mx_hi)
}
