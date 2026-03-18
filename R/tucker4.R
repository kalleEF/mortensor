#' Tucker4-ILR coherent multi-disease forecasting
#'
#' Joint forecast of disease-specific mortality as compositional shares
#' of all-cause mortality.  ILR transform on the K+1-part cause
#' composition, then Tucker4 decomposition on the 4D tensor
#' \[ages x years x groups x ILR-coords\].
#'
#' @name tucker4
#' @importFrom compositions ilr ilrInv acomp
NULL

# ---- Cause-share computation ----

#' Compute cause shares from disease-specific and all-cause mortality
#'
#' @param mx_disease_list Named list of 3D arrays, one per disease.
#' @param mx_all          3D array of all-cause mortality.
#' @param disease_order   Character vector of disease names.
#' @return 4D array \[n_a x n_t x n_g x (K+1)\] of shares summing to 1.
#' @export
compute_cause_shares <- function(mx_disease_list, mx_all, disease_order) {
  n_a <- dim(mx_all)[1]; n_t <- dim(mx_all)[2]; n_g <- dim(mx_all)[3]
  K <- length(disease_order)
  shares <- array(NA_real_, dim = c(n_a, n_t, n_g, K + 1))

  for (d in seq_along(disease_order)) {
    nm <- disease_order[d]
    shares[, , , d] <- mx_disease_list[[nm]] / pmax(mx_all, 1e-12)
  }

  share_sum <- apply(shares[, , , 1:K, drop = FALSE], c(1, 2, 3), sum, na.rm = TRUE)
  remainder <- pmax(1 - share_sum, 1e-6)
  shares[, , , K + 1] <- remainder

  # Renormalise
  total <- apply(shares, c(1, 2, 3), sum, na.rm = TRUE)
  for (d in seq_len(K + 1))
    shares[, , , d] <- shares[, , , d] / total

  # Multiplicative zero replacement (Aitchison-style)
  delta <- 0.65 / (n_a * n_t * n_g)^2
  for (i in seq_len(n_a))
    for (j in seq_len(n_t))
      for (k in seq_len(n_g)) {
        sv <- shares[i, j, k, ]
        zeros <- sv < delta
        if (any(zeros)) {
          n_z <- sum(zeros)
          sv[zeros] <- delta
          sv[!zeros] <- sv[!zeros] * (1 - n_z * delta) / sum(sv[!zeros])
          shares[i, j, k, ] <- sv
        }
      }
  shares
}

# ---- ILR transform wrappers ----

#' ILR transform along the 4th mode of a shares array
#'
#' @param shares_4d 4D array \[n_a x n_t x n_g x K+1\] of shares.
#' @return 4D array \[n_a x n_t x n_g x K\] of ILR coordinates.
#' @export
ilr_cause_transform <- function(shares_4d) {
  dims <- dim(shares_4d)
  n_a <- dims[1]; n_t <- dims[2]; n_g <- dims[3]; Kp1 <- dims[4]
  K <- Kp1 - 1
  X_ilr <- array(NA_real_, dim = c(n_a, n_t, n_g, K))
  for (i in seq_len(n_a))
    for (j in seq_len(n_t))
      for (k in seq_len(n_g))
        X_ilr[i, j, k, ] <- as.numeric(
          compositions::ilr(compositions::acomp(shares_4d[i, j, k, ])))
  X_ilr
}

#' Inverse ILR: unconstrained coords -> compositional shares
#'
#' @param ilr_vec Numeric vector of length K.
#' @return Numeric vector of length K+1 summing to 1.
#' @keywords internal
ilr_cause_inverse <- function(ilr_vec) {
  as.numeric(compositions::ilrInv(ilr_vec))
}

# ---- Tucker4 reconstruction ----

#' 4D Tucker reconstruction
#'
#' Computes X_(1) = A G_(1) (D kron C kron B)^T.
#'
#' @param A     Mode-1 factor matrix (age).
#' @param B_mat Mode-2 factor matrix (time).
#' @param C     Mode-3 factor matrix (group).
#' @param D_mat Mode-4 factor matrix (cause).
#' @param G_mat Mode-1 unfolded core tensor.
#' @return 4D array \[n_a x n_t x n_g x n_d\].
#' @export
reconstruct_array4 <- function(A, B_mat, C, D_mat, G_mat) {
  X_mat <- A %*% G_mat %*% t(kronecker(D_mat, kronecker(C, B_mat)))
  array(X_mat, dim = c(nrow(A), nrow(B_mat), nrow(C), nrow(D_mat)))
}

# ---- Tucker4 coherent forecast ----

#' Tucker4 forecast on ILR cause-shares
#'
#' Decomposes the 4D ILR tensor, forecasts the B (time) mode,
#' reconstructs, and inverts ILR to produce share forecasts.
#'
#' @param X_ilr   4D array \[n_a x n_t x n_g x K\] of ILR coordinates.
#' @param years   Integer vector of observed years.
#' @param h       Forecast horizon.
#' @param R_A,R_B,R_C,R_D Tucker4 ranks.
#' @param ts_method `"ets"`, `"arima"`, or `"rwd"`.
#' @param fast    Use fast convergence settings?
#' @param return_full Return decomposition components for bootstrap?
#' @return List with `shares_forecast` \[n_a x h x n_g x K+1\] and `fit_pct`.
#' @export
tucker4_coherent_forecast <- function(X_ilr, years, h,
                                      R_A = 3, R_B = 3, R_C = 2, R_D = 3,
                                      ts_method = "ets", fast = TRUE,
                                      return_full = FALSE) {
  n_a <- dim(X_ilr)[1]; n_t <- dim(X_ilr)[2]
  n_g <- dim(X_ilr)[3]; K   <- dim(X_ilr)[4]
  R_A <- min(R_A, n_a - 1); R_B <- min(R_B, n_t - 1)
  R_C <- min(R_C, n_g - 1); R_D <- min(R_D, K)

  decomp <- rTensor::tucker(
    rTensor::as.tensor(X_ilr),
    ranks = c(R_A, R_B, R_C, R_D),
    max_iter = ifelse(fast, 50, 200),
    tol = ifelse(fast, 1e-6, 1e-10))

  A     <- decomp$U[[1]]
  B     <- decomp$U[[2]]
  C     <- decomp$U[[3]]
  D_mat <- decomp$U[[4]]
  G_mat <- rTensor::k_unfold(decomp$Z, m = 1)@data

  # Forecast B-columns
  arima_fits <- vector("list", R_B)
  B_fc <- matrix(NA_real_, h, R_B)
  for (r in seq_len(R_B)) {
    ts_r <- stats::ts(B[, r], start = min(years), frequency = 1)
    arima_fits[[r]] <- tryCatch({
      if (ts_method == "arima")
        forecast::auto.arima(ts_r, stepwise = fast, approximation = fast)
      else if (ts_method == "ets")
        forecast::ets(ts_r, model = "ZZN")
      else
        forecast::Arima(ts_r, order = c(0, 1, 0), include.drift = TRUE)
    }, error = function(e) tryCatch(
      forecast::Arima(ts_r, order = c(0, 1, 0)),
      error = function(e2)
        forecast::Arima(ts_r, order = c(0, 0, 0), include.mean = TRUE)))
    B_fc[, r] <- as.numeric(forecast::forecast(arima_fits[[r]], h = h)$mean)
  }
  B_full <- rbind(B, B_fc)

  X_ilr_full <- reconstruct_array4(A, B_full, C, D_mat, G_mat)

  # Invert ILR -> shares
  shares_fc <- array(NA_real_, dim = c(n_a, h, n_g, K + 1))
  for (i in seq_len(n_a))
    for (j in seq_len(h))
      for (k in seq_len(n_g))
        shares_fc[i, j, k, ] <- ilr_cause_inverse(X_ilr_full[i, n_t + j, k, ])

  result <- list(shares_forecast = shares_fc, fit_pct = decomp$norm_percent)
  if (return_full) {
    result$decomp     <- decomp
    result$arima_fits <- arima_fits
    result$A <- A; result$B <- B; result$C <- C; result$D_mat <- D_mat
    result$G_mat      <- G_mat
    result$B_full     <- B_full
    result$X_ilr      <- X_ilr
    result$X_ilr_full <- X_ilr_full
    result$resid_array <- X_ilr - reconstruct_array4(A, B, C, D_mat, G_mat)
    result$n_t        <- n_t
    result$K          <- K
  }
  result
}

# ---- Tucker4 cross-validation & calibration ----

#' Rolling-origin CV for Tucker4 coherent forecast
#'
#' @param X_ilr           4D ILR array.
#' @param mx_all          3D all-cause array.
#' @param mx_disease_list Named list of 3D disease arrays.
#' @param disease_order   Character vector.
#' @param years           Integer vector.
#' @param R_A,R_B,R_C,R_D Tucker4 ranks.
#' @param ts_method       Time-series method.
#' @param h_val           Validation horizon.
#' @param min_train       Minimum training years.
#' @param metric_fn Scoring function with signature
#'   `function(arr_test, fc_arr, arr_train)`. Default: [compute_mase_array()].
#' @return Scalar mean CV score.
#' @export
cv_rolling_origin_tucker4 <- function(X_ilr, mx_all, mx_disease_list,
                                      disease_order, years,
                                      R_A, R_B, R_C, R_D, ts_method,
                                      h_val = 3, min_train = 12,
                                      metric_fn = compute_mase_array) {
  n_yr <- length(years)
  last_origin  <- n_yr - h_val
  first_origin <- max(min_train, last_origin - 5)
  if (first_origin > last_origin) return(Inf)

  origins   <- first_origin:last_origin
  mase_vals <- numeric(length(origins))
  K <- length(disease_order)

  for (o in seq_along(origins)) {
    t_end    <- origins[o]
    h_actual <- min(h_val, n_yr - t_end)
    if (h_actual < 1) { mase_vals[o] <- Inf; next }

    ilr_train <- X_ilr[, 1:t_end, , , drop = FALSE]
    yrs_train <- years[1:t_end]

    result <- tryCatch(
      tucker4_coherent_forecast(ilr_train, yrs_train, h_actual,
                                R_A, R_B, R_C, R_D, ts_method, fast = TRUE),
      error = function(e) NULL)
    if (is.null(result)) { mase_vals[o] <- Inf; next }

    test_idx <- (t_end + 1):(t_end + h_actual)
    mx_all_test <- mx_all[, test_idx, , drop = FALSE]
    per_disease_mase <- c()

    for (d in seq_len(K)) {
      mx_fc_d <- result$shares_forecast[, , , d] * mx_all_test
      mx_obs_d <- mx_disease_list[[disease_order[d]]][, test_idx, , drop = FALSE]
      mx_train_d <- mx_disease_list[[disease_order[d]]][, 1:t_end, , drop = FALSE]
      m <- metric_fn(mx_obs_d, mx_fc_d, mx_train_d)
      if (!is.na(m) && is.finite(m)) per_disease_mase <- c(per_disease_mase, m)
    }
    mase_vals[o] <- if (length(per_disease_mase) == 0) Inf else mean(per_disease_mase)
  }
  finite <- mase_vals[is.finite(mase_vals)]
  if (length(finite) == 0) Inf else mean(finite)
}

#' Calibrate Tucker4 via two-stage grid search
#'
#' Stage 1: coarse grid. Stage 2: fine grid around best from stage 1.
#'
#' @param X_ilr           4D ILR array.
#' @param mx_all          3D all-cause array.
#' @param mx_disease_list Named list of disease arrays.
#' @param disease_order   Disease name vector.
#' @param years           Integer vector.
#' @param h_val           Validation horizon.
#' @param min_train       Minimum training years.
#' @param cv_metric Character; scoring metric for cross-validation.
#'   \code{"mase"} (default), \code{"mase_combined"}, or \code{"mase_trend"}.
#'   See [resolve_cv_metric()].
#' @return Named list of best ranks and `cv_score`.
#' @export
calibrate_tucker4 <- function(X_ilr, mx_all, mx_disease_list, disease_order,
                              years, h_val = 3, min_train = 12,
                              cv_metric = c("mase", "mase_combined", "mase_trend")) {
  metric_fn <- resolve_cv_metric(cv_metric)
  n_a <- dim(X_ilr)[1]; n_t <- dim(X_ilr)[2]
  n_g <- dim(X_ilr)[3]; K   <- dim(X_ilr)[4]

  ra_coarse <- seq(2, min(8, n_a - 1), by = 2)
  rb_coarse <- seq(2, min(6, n_t - 1), by = 2)
  rc_coarse <- seq(2, min(5, n_g - 1), by = 1)
  rd_coarse <- seq(2, min(7, K),       by = 2)
  ts_grid   <- c("ets", "arima", "rwd")

  best_mase   <- Inf
  best_params <- list(R_A = 3, R_B = 3, R_C = 2, R_D = 3, ts_method = "ets")

  eval_one <- function(ra, rb, rc, rd, ts_m) {
    tryCatch(
      cv_rolling_origin_tucker4(X_ilr, mx_all, mx_disease_list,
                                disease_order, years,
                                ra, rb, rc, rd, ts_m,
                                h_val = h_val, min_train = min_train,
                                metric_fn = metric_fn),
      error = function(e) Inf)
  }

  # Stage 1: coarse
  for (ra in ra_coarse)
    for (rb in rb_coarse)
      for (rc in rc_coarse)
        for (rd in rd_coarse)
          for (ts_m in ts_grid) {
            cv_m <- eval_one(ra, rb, rc, rd, ts_m)
            if (is.finite(cv_m) && cv_m < best_mase) {
              best_mase   <- cv_m
              best_params <- list(R_A = ra, R_B = rb, R_C = rc,
                                  R_D = rd, ts_method = ts_m)
            }
          }

  # Stage 2: fine (+/-1 around coarse best)
  refine <- function(val, lo, hi)
    unique(pmin(pmax(seq(val - 1, val + 1), lo), hi))
  ra_fine <- refine(best_params$R_A, 2, min(10, n_a - 1))
  rb_fine <- refine(best_params$R_B, 2, min(6,  n_t - 1))
  rc_fine <- refine(best_params$R_C, 2, min(5,  n_g - 1))
  rd_fine <- refine(best_params$R_D, 2, min(8,  K))

  for (ra in ra_fine)
    for (rb in rb_fine)
      for (rc in rc_fine)
        for (rd in rd_fine)
          for (ts_m in ts_grid) {
            cv_m <- eval_one(ra, rb, rc, rd, ts_m)
            if (is.finite(cv_m) && cv_m < best_mase) {
              best_mase   <- cv_m
              best_params <- list(R_A = ra, R_B = rb, R_C = rc,
                                  R_D = rd, ts_method = ts_m)
            }
          }
  best_params$cv_score <- best_mase
  best_params
}

# ---- Coherent holdout comparison ----

#' Compare coherent Tucker4 vs independent per-disease Tucker3
#'
#' @param X_ilr                 4D ILR array.
#' @param mx_all                3D all-cause array.
#' @param mx_disease_list       Named list of disease arrays (raw).
#' @param mx_disease_smooth_list Named list of smoothed disease arrays.
#' @param pop_all               3D population array.
#' @param disease_order         Character vector.
#' @param disease_configs       Named list of mort_config objects.
#' @param years,ages,groups     Study dimensions.
#' @param train_idx,test_idx    Integer index vectors.
#' @param calib_tucker4         Output of [calibrate_tucker4()].
#' @param width                 Age-group width.
#' @return data.table with per-disease coherent vs independent MASE.
#' @export
compare_coherent_vs_independent <- function(
    X_ilr, mx_all, mx_disease_list, mx_disease_smooth_list,
    pop_all, disease_order, disease_configs,
    years, ages, groups, train_idx, test_idx,
    calib_tucker4, width = 5) {

  h_total  <- max(test_idx) - max(train_idx)
  fc_steps <- test_idx - max(train_idx)
  K <- length(disease_order)

  results_list <- list()

  # Coherent Tucker4-ILR
  ilr_train <- X_ilr[, train_idx, , , drop = FALSE]
  yrs_train <- years[train_idx]
  coh_fc <- tryCatch(
    tucker4_coherent_forecast(ilr_train, yrs_train, h_total,
                              calib_tucker4$R_A, calib_tucker4$R_B,
                              calib_tucker4$R_C, calib_tucker4$R_D,
                              calib_tucker4$ts_method, fast = TRUE),
    error = function(e) NULL)

  mx_all_test <- mx_all[, test_idx, , drop = FALSE]

  for (d in seq_len(K)) {
    nm <- disease_order[d]
    mx_obs_test  <- mx_disease_list[[nm]][, test_idx, , drop = FALSE]
    mx_train_d   <- mx_disease_list[[nm]][, train_idx, , drop = FALSE]

    # Coherent MASE
    coh_mase <- NA_real_
    if (!is.null(coh_fc)) {
      shares_fc_d <- coh_fc$shares_forecast[, fc_steps, , d, drop = FALSE]
      dim(shares_fc_d) <- dim(mx_obs_test)
      mx_fc_coh <- shares_fc_d * mx_all_test
      coh_mase <- compute_mase_array(mx_obs_test, mx_fc_coh, mx_train_d)
    }

    # Independent Tucker3-Direct MASE
    indep_mase <- NA_real_
    tryCatch({
      mx_smooth_d <- mx_disease_smooth_list[[nm]]
      calib_d <- calibrate_tucker3(mx_smooth_d[, train_idx, , drop = FALSE],
                                   yrs_train, ages,
                                   method = "direct", width = width,
                                   h_val = length(test_idx), min_train = 12)
      fc_d <- tucker3_direct_forecast(mx_smooth_d[, train_idx, , drop = FALSE],
                                      yrs_train, h_total,
                                      calib_d$R_A, calib_d$R_B, calib_d$R_C,
                                      calib_d$ts_method, fast = TRUE)
      mx_fc_ind <- fc_d$mx_forecast[, fc_steps, , drop = FALSE]
      indep_mase <- compute_mase_array(mx_obs_test, mx_fc_ind, mx_train_d)
    }, error = function(e) NULL)

    results_list[[length(results_list) + 1]] <- data.table::data.table(
      disease = disease_configs[[nm]]$cause_label,
      disease_short = nm,
      coherent_mase = coh_mase,
      independent_mase = indep_mase)
  }
  data.table::rbindlist(results_list)
}

# ---- Coherent bootstrap PIs ----

#' Bootstrap PIs for Tucker4 coherent forecast
#'
#' @param fc_obj    Tucker4 forecast object (return_full = TRUE).
#' @param mx_all_fc 3D array \[n_a x h x n_g\] of all-cause point forecast.
#' @param n_sim     Number of simulations.
#' @param h         Forecast horizon.
#' @return List with 4D arrays: `mx_lower`, `mx_upper`, `mx_median`
#'   \[n_a x h x n_g x K\].
#' @export
bootstrap_coherent_forecast <- function(fc_obj, mx_all_fc, n_sim, h) {
  A <- fc_obj$A; B <- fc_obj$B; C <- fc_obj$C; D_mat <- fc_obj$D_mat
  G_mat <- fc_obj$G_mat
  n_a <- nrow(A); n_t <- fc_obj$n_t; R_B <- ncol(B)
  K <- fc_obj$K; n_g <- nrow(C)
  resid <- fc_obj$resid_array

  mx_sims <- array(NA_real_, dim = c(n_a, h, n_g, K, n_sim))
  for (s in seq_len(n_sim)) {
    B_sim <- matrix(NA_real_, h, R_B)
    for (r in seq_len(R_B)) {
      B_sim[, r] <- tryCatch(
        simulate(fc_obj$arima_fits[[r]], nsim = h,
                 future = TRUE, bootstrap = TRUE),
        error = function(e) {
          as.numeric(forecast::forecast(fc_obj$arima_fits[[r]], h = h)$mean) +
            stats::rnorm(h, 0, 0.1)
        })
    }
    B_full_sim <- rbind(B, B_sim)
    X_sim <- reconstruct_array4(A, B_full_sim, C, D_mat, G_mat)

    for (j in seq_len(h)) {
      ry <- sample(1:n_t, 1)
      for (k in seq_len(n_g)) {
        ilr_noisy <- X_sim[, n_t + j, k, ] + resid[, ry, k, ]
        shares_mat <- matrix(NA_real_, n_a, K + 1)
        for (i in seq_len(n_a))
          shares_mat[i, ] <- ilr_cause_inverse(ilr_noisy[i, ])
        for (d in seq_len(K))
          mx_sims[, j, k, d, s] <- pmax(shares_mat[, d] * mx_all_fc[, j, k], 1e-12)
      }
    }
  }

  list(
    mx_lower  = apply(mx_sims, c(1, 2, 3, 4), stats::quantile, 0.025, na.rm = TRUE),
    mx_upper  = apply(mx_sims, c(1, 2, 3, 4), stats::quantile, 0.975, na.rm = TRUE),
    mx_median = apply(mx_sims, c(1, 2, 3, 4), stats::quantile, 0.50,  na.rm = TRUE)
  )
}
