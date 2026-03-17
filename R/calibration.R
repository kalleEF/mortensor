#' Cross-validate Tucker3 with rolling origin
#'
#' @param mx_array 3D array of mortality rates.
#' @param years Integer vector of years.
#' @param ages Integer vector of age-group starts.
#' @param R_A,R_B,R_C Tucker3 ranks.
#' @param ts_method Time-series method.
#' @param method `"coda"` or `"direct"`.
#' @param h_val Validation horizon.
#' @param min_train Minimum training years.
#' @param width Age-group width.
#' @return Scalar mean MASE.
#' @export
cv_rolling_origin <- function(mx_array, years, ages,
                              R_A, R_B, R_C, ts_method,
                              method = "coda",
                              h_val = 3, min_train = 12, width = 5) {
  n_yr <- length(years)
  last_origin  <- n_yr - h_val
  first_origin <- max(min_train, last_origin - 5)
  if (first_origin > last_origin) return(Inf)

  origins   <- first_origin:last_origin
  mase_vals <- numeric(length(origins))

  for (i in seq_along(origins)) {
    t_end    <- origins[i]
    h_actual <- min(h_val, n_yr - t_end)
    if (h_actual < 1) { mase_vals[i] <- Inf; next }

    arr_train <- mx_array[, 1:t_end, , drop = FALSE]
    yrs_train <- years[1:t_end]
    arr_test  <- mx_array[, (t_end + 1):(t_end + h_actual), , drop = FALSE]

    result <- tryCatch({
      if (method == "coda")
        tucker3_coda_forecast(arr_train, yrs_train, ages, h_actual,
                              R_A, R_B, R_C, ts_method, fast = TRUE, width = width)
      else
        tucker3_direct_forecast(arr_train, yrs_train, h_actual,
                                R_A, R_B, R_C, ts_method, fast = TRUE)
    }, error = function(e) NULL)

    if (is.null(result)) { mase_vals[i] <- Inf; next }
    m <- compute_mase_array(arr_test, result$mx_forecast, arr_train)
    mase_vals[i] <- if (is.na(m) || !is.finite(m)) Inf else m
  }
  finite <- mase_vals[is.finite(mase_vals)]
  if (length(finite) == 0) Inf else mean(finite)
}

#' Information criterion for Tucker3 rank selection
#'
#' Computes BIC or AIC for a Tucker3 decomposition with given ranks.
#' Much faster than cross-validation and useful for initial rank screening.
#'
#' @param mx_array 3D array of mortality rates.
#' @param R_A,R_B,R_C Integer Tucker3 ranks.
#' @param method `"coda"` or `"direct"`.
#' @param criterion `"bic"` (default) or `"aic"`.
#' @param ages Integer vector of age-group starts (required for `"coda"`).
#' @param width Age-group width (default 5).
#' @return A list with `criterion` (BIC/AIC value), `fit_pct`, `n_params`, and
#'   `ranks`.
#' @export
tucker3_info_criterion <- function(mx_array, R_A, R_B, R_C,
                                    method = c("coda", "direct"),
                                    criterion = c("bic", "aic"),
                                    ages = NULL, width = 5) {
  method    <- match.arg(method)
  criterion <- match.arg(criterion)

  n_a <- dim(mx_array)[1]
  n_t <- dim(mx_array)[2]
  n_g <- dim(mx_array)[3]
  R_A <- min(R_A, n_a - 1)
  R_B <- min(R_B, n_t - 1)
  R_C <- min(R_C, n_g - 1)

  # Build input array
  if (method == "coda") {
    stopifnot(!is.null(ages))
    X <- array(NA_real_, dim = c(n_a, n_t, n_g))
    for (j in seq_len(n_t))
      for (k in seq_len(n_g)) {
        dx <- compute_dx(mx_array[, j, k], ages, width)
        X[, j, k] <- as.numeric(compositions::clr(dx))
      }
  } else {
    X <- log(pmax(mx_array, 1e-15))
  }

  # Decompose
  decomp <- rTensor::tucker(rTensor::as.tensor(X),
                             ranks = c(R_A, R_B, R_C),
                             max_iter = 50, tol = 1e-6)

  # Reconstruct and compute RSS
  A <- decomp$U[[1]]
  B <- decomp$U[[2]]
  C <- decomp$U[[3]]
  G_mat <- rTensor::k_unfold(decomp$Z, m = 1)@data
  X_hat <- reconstruct_array(A, B, C, G_mat)
  rss <- sum((X - X_hat)^2)

  # Parameter count: factor matrices + core tensor
  n_params <- R_A * n_a + R_B * n_t + R_C * n_g + R_A * R_B * R_C
  n_obs    <- prod(dim(X))

  # Information criterion
  ic <- if (criterion == "bic") {
    n_obs * log(rss / n_obs) + log(n_obs) * n_params
  } else {
    n_obs * log(rss / n_obs) + 2 * n_params
  }

  list(
    criterion = ic,
    fit_pct   = decomp$norm_percent,
    n_params  = n_params,
    ranks     = c(R_A = R_A, R_B = R_B, R_C = R_C)
  )
}

#' Calibrate Tucker3 via grid search
#'
#' Searches over rank combinations and time-series methods to find the best
#' Tucker3 configuration. Supports two selection strategies: rolling-origin
#' cross-validation (\code{"cv"}, slower but more robust) or information
#' criteria (\code{"bic"}/\code{"aic"}, fast screening).
#'
#' @param mx_array 3D array of mortality rates.
#' @param years Integer vector of years.
#' @param ages Integer vector of age-group starts.
#' @param method `"coda"` or `"direct"`.
#' @param width Age-group width.
#' @param h_val Validation horizon (used only for \code{rank_selection = "cv"}).
#' @param min_train Minimum training years (used only for \code{rank_selection = "cv"}).
#' @param rank_selection Character; \code{"cv"} (default) for rolling-origin
#'   cross-validation, \code{"bic"} or \code{"aic"} for information criterion
#'   based selection. When using IC, the ts_method is selected via a smaller CV
#'   using the best ranks.
#' @return Named list with `R_A`, `R_B`, `R_C`, `ts_method`, `cv_mase`
#'   (or `ic_value` for BIC/AIC selection).
#' @export
calibrate_tucker3 <- function(mx_array, years, ages,
                              method = "coda", width = 5,
                              h_val = 3, min_train = 12,
                              rank_selection = c("cv", "bic", "aic")) {
  rank_selection <- match.arg(rank_selection)
  n_a <- dim(mx_array)[1]; n_t <- dim(mx_array)[2]; n_g <- dim(mx_array)[3]

  ra_grid <- 2:min(6, n_a - 1)
  rb_grid <- 2:min(8, n_t - 1)
  rc_grid <- 2:min(4, n_g - 1)
  ts_grid <- c("ets", "arima", "rwd")

  if (rank_selection %in% c("bic", "aic")) {
    # Phase 1: Select ranks via information criterion
    best_ic     <- Inf
    best_ranks  <- c(R_A = 3, R_B = 3, R_C = 2)
    for (ra in ra_grid) {
      for (rb in rb_grid) {
        for (rc in rc_grid) {
          res <- tryCatch(
            tucker3_info_criterion(mx_array, ra, rb, rc,
                                    method = method, criterion = rank_selection,
                                    ages = ages, width = width),
            error = function(e) NULL)
          if (!is.null(res) && is.finite(res$criterion) && res$criterion < best_ic) {
            best_ic    <- res$criterion
            best_ranks <- c(R_A = ra, R_B = rb, R_C = rc)
          }
        }
      }
    }

    # Phase 2: Select ts_method via CV with the best ranks
    best_mase   <- Inf
    best_ts     <- "ets"
    for (ts_m in ts_grid) {
      cv_m <- tryCatch(
        cv_rolling_origin(mx_array, years, ages,
                          best_ranks["R_A"], best_ranks["R_B"], best_ranks["R_C"],
                          ts_m, method = method, h_val = h_val,
                          min_train = min_train, width = width),
        error = function(e) Inf)
      if (is.finite(cv_m) && cv_m < best_mase) {
        best_mase <- cv_m
        best_ts   <- ts_m
      }
    }

    return(list(R_A = best_ranks["R_A"], R_B = best_ranks["R_B"],
                R_C = best_ranks["R_C"], ts_method = best_ts,
                cv_mase = best_mase, ic_value = best_ic,
                rank_selection = rank_selection))
  }

  # Standard CV-based selection
  best_mase   <- Inf
  best_params <- list(R_A = 3, R_B = 3, R_C = 2, ts_method = "ets")

  for (ra in ra_grid) {
    for (rb in rb_grid) {
      for (rc in rc_grid) {
        for (ts_m in ts_grid) {
          cv_m <- tryCatch(
            cv_rolling_origin(mx_array, years, ages, ra, rb, rc, ts_m,
                              method = method, h_val = h_val,
                              min_train = min_train, width = width),
            error = function(e) Inf)
          if (is.finite(cv_m) && cv_m < best_mase) {
            best_mase   <- cv_m
            best_params <- list(R_A = ra, R_B = rb, R_C = rc, ts_method = ts_m)
          }
        }
      }
    }
  }
  best_params$cv_mase <- best_mase
  best_params
}

#' Expanding-window cross-validation for a generic model
#'
#' @param mx_array 3D array of mortality rates.
#' @param pop_array 3D array of population exposure.
#' @param years Integer vector.
#' @param ages Integer vector.
#' @param groups data.frame of group strata.
#' @param fit_fn A function with signature
#'   \code{function(arr_train, pop_train, yrs_train, ages, groups, h, ...)}
#'   returning a list with element \code{mx_forecast}.
#' @param h_val Validation horizon.
#' @param n_origins Number of expanding origins.
#' @param ... Passed to `fit_fn`.
#' @return Scalar mean MASE.
#' @export
cv_expanding_model <- function(mx_array, pop_array, years, ages, groups,
                               fit_fn, h_val = 4, n_origins = 4, ...) {
  n_yr <- length(years)
  last_origin  <- n_yr - h_val
  first_origin <- max(h_val + 5, last_origin - n_origins + 1)
  if (first_origin > last_origin) return(Inf)

  origins <- first_origin:last_origin
  mase_vals <- numeric(length(origins))
  for (i in seq_along(origins)) {
    t_end <- origins[i]
    h_actual <- min(h_val, n_yr - t_end)
    if (h_actual < 1) { mase_vals[i] <- Inf; next }
    arr_train <- mx_array[, 1:t_end, , drop = FALSE]
    arr_test  <- mx_array[, (t_end + 1):(t_end + h_actual), , drop = FALSE]
    pop_train <- pop_array[, 1:t_end, , drop = FALSE]
    yrs_train <- years[1:t_end]
    result <- tryCatch(
      fit_fn(arr_train, pop_train, yrs_train, ages, groups, h_actual, ...),
      error = function(e) NULL)
    if (is.null(result) || is.null(result$mx_forecast)) {
      mase_vals[i] <- Inf; next
    }
    m <- compute_mase_array(arr_test, result$mx_forecast, arr_train)
    mase_vals[i] <- if (is.na(m) || !is.finite(m)) Inf else m
  }
  finite <- mase_vals[is.finite(mase_vals)]
  if (length(finite) == 0) Inf else mean(finite)
}

#' Calibrate all non-Tucker3 models
#'
#' Grid-searches over hyperparameters for Lee-Carter, FDM, StMoMo-APC,
#' RWD, ARIMA, and ETS using expanding-window CV.
#'
#' @param mx_array 3D array of mortality rates.
#' @param pop_array 3D array of population exposure.
#' @param years Integer vector.
#' @param ages Integer vector.
#' @param groups data.frame of group strata.
#' @param h_val Validation horizon.
#' @param n_origins Number of CV origins.
#' @param width Age-group width.
#' @return Named list of best parameters per model.
#' @export
calibrate_all_models <- function(mx_array, pop_array, years, ages, groups,
                                 h_val = 4, n_origins = 4, width = 5) {
  calib <- list()
  n_a <- dim(mx_array)[1]

  # Lee-Carter
  best_lc <- list(adjust = "dt", chooseperiod = FALSE, cv_mase = Inf)
  if (requireNamespace("demography", quietly = TRUE)) {
    for (adj in c("dt", "dxt", "e0", "none")) {
      for (cp in c(FALSE, TRUE)) {
        lc_fn <- function(arr, pop, yrs, ag, gr, h_a, ...)
          lca_forecast(arr, pop, yrs, ag, gr, h_a, adjust = adj, chooseperiod = cp)
        m <- tryCatch(cv_expanding_model(mx_array, pop_array, years, ages, groups,
                                         lc_fn, h_val = h_val, n_origins = n_origins),
                      error = function(e) Inf)
        if (is.finite(m) && m < best_lc$cv_mase)
          best_lc <- list(adjust = adj, chooseperiod = cp, cv_mase = m)
      }
    }
  }
  calib$lc <- best_lc

  # FDM
  best_fdm <- list(order = 6, method = "arima", cv_mase = Inf)
  if (requireNamespace("demography", quietly = TRUE)) {
    for (ord in c(4, 6)) {
      for (meth in c("arima", "ets")) {
        fdm_fn <- function(arr, pop, yrs, ag, gr, h_a, ...)
          fdm_forecast(arr, pop, yrs, ag, gr, h_a, order = ord, method = meth)
        m <- tryCatch(cv_expanding_model(mx_array, pop_array, years, ages, groups,
                                         fdm_fn, h_val = h_val, n_origins = n_origins),
                      error = function(e) Inf)
        if (is.finite(m) && m < best_fdm$cv_mase)
          best_fdm <- list(order = ord, method = meth, cv_mase = m)
      }
    }
  }
  calib$fdm <- best_fdm

  # StMoMo
  best_stmomo <- list(model_type = "APC", cv_mase = Inf)
  if (requireNamespace("StMoMo", quietly = TRUE)) {
    for (mt in c("APC", "CBD")) {
      stm_fn <- function(arr, pop, yrs, ag, gr, h_a, ...)
        stmomo_forecast(arr, pop, yrs, ag, gr, h_a, model_type = mt)
      m <- tryCatch(cv_expanding_model(mx_array, pop_array, years, ages, groups,
                                       stm_fn, h_val = h_val, n_origins = n_origins),
                    error = function(e) Inf)
      if (is.finite(m) && m < best_stmomo$cv_mase)
        best_stmomo <- list(model_type = mt, cv_mase = m)
    }
  }
  calib$stmomo <- best_stmomo

  # RWD
  best_rwd <- list(drift_window = 0, cv_mase = Inf)
  for (dw in c(0, length(years), max(1, length(years) - 5))) {
    rwd_fn <- function(arr, pop, yrs, ag, gr, h_a, ...)
      rwd_forecast(arr, yrs, h_a, drift_window = dw)
    m <- tryCatch(cv_expanding_model(mx_array, pop_array, years, ages, groups,
                                     rwd_fn, h_val = h_val, n_origins = n_origins),
                  error = function(e) Inf)
    if (is.finite(m) && m < best_rwd$cv_mase)
      best_rwd <- list(drift_window = dw, cv_mase = m)
  }
  calib$rwd <- best_rwd

  # ARIMA
  best_arima <- list(max.p = 2, max.q = 1, max.d = 1, cv_mase = Inf)
  for (mp in c(1, 2)) {
    for (mq in c(0, 1)) {
      for (md in c(1, 2)) {
        arima_fn <- function(arr, pop, yrs, ag, gr, h_a, ...)
          arima_forecast(arr, yrs, h_a, max.p = mp, max.q = mq, max.d = md)
        m <- tryCatch(cv_expanding_model(mx_array, pop_array, years, ages, groups,
                                         arima_fn, h_val = h_val, n_origins = n_origins),
                      error = function(e) Inf)
        if (is.finite(m) && m < best_arima$cv_mase)
          best_arima <- list(max.p = mp, max.q = mq, max.d = md, cv_mase = m)
      }
    }
  }
  calib$arima <- best_arima

  # ETS
  best_ets <- list(ets_model = "AAN", damped = TRUE, cv_mase = Inf)
  for (em in c("AAN", "AAdN", "ZZN", "MMN", "ANN")) {
    ets_fn <- function(arr, pop, yrs, ag, gr, h_a, ...)
      ets_forecast(arr, yrs, h_a, ets_model = em, damped = TRUE)
    m <- tryCatch(cv_expanding_model(mx_array, pop_array, years, ages, groups,
                                     ets_fn, h_val = h_val, n_origins = n_origins),
                  error = function(e) Inf)
    if (is.finite(m) && m < best_ets$cv_mase)
      best_ets <- list(ets_model = em, damped = TRUE, cv_mase = m)
  }
  calib$ets <- best_ets

  calib
}
