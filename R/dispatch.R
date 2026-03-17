#' Dispatch forecast or bootstrap for any supported model
#'
#' @name dispatch
NULL

#' Dispatch a full-horizon forecast for any model
#'
#' Routes to the appropriate forecasting function based on `model_name`
#' and returns a standardised result list.
#'
#' @param model_name   One of `"Tucker3-CoDa"`, `"Tucker3-Direct"`,
#'   `"Lee-Carter"`, `"FDM"`, `"StMoMo-APC"`, `"RWD"`,
#'   `"ARIMA-perSeries"`, `"ETS-Damped"`.
#' @param mx           3D array of (smoothed) mortality rates.
#' @param pop          3D population array.
#' @param years        Integer vector.
#' @param ages         Integer vector of age-group starts.
#' @param groups       data.frame with `sex` and `quintile`.
#' @param h            Forecast horizon.
#' @param calib_coda   Tucker3-CoDa calibration result.
#' @param calib_direct Tucker3-Direct calibration result.
#' @param width        Age-group width.
#' @param calib_all    Optional calibration for non-Tucker3 models.
#' @param return_fitted Return fitted in-sample values?
#' @return List: `mx_forecast`, `mx_fitted`, `full_fc`, `method`.
#' @export
dispatch_forecast <- function(model_name, mx, pop, years, ages, groups,
                              h, calib_coda, calib_direct,
                              width = 5L, calib_all = NULL,
                              return_fitted = TRUE) {
  n_a <- length(ages); n_t <- length(years); n_g <- nrow(groups)

  if (model_name == "Tucker3-CoDa") {
    cp <- calib_coda
    fc <- tucker3_coda_forecast(mx, years, ages, h,
                                cp$R_A, cp$R_B, cp$R_C, cp$ts_method,
                                fast = FALSE, return_full = TRUE, width = width)
    fitted_arr <- NULL
    if (return_fitted && !is.null(fc$X_clr_full)) {
      fitted_arr <- array(NA_real_, dim = c(n_a, n_t, n_g))
      for (kk in seq_len(n_g)) for (jj in seq_len(n_t)) {
        dx_f <- inv_clr_to_dx(fc$X_clr_full[, jj, kk])
        fitted_arr[, jj, kk] <- dx_to_mx(dx_f, ages, width)
      }
    }
    return(list(mx_forecast = fc$mx_forecast, mx_fitted = fitted_arr,
                full_fc = fc, method = "coda"))

  } else if (model_name == "Tucker3-Direct") {
    cp <- calib_direct
    fc <- tucker3_direct_forecast(mx, years, h,
                                  cp$R_A, cp$R_B, cp$R_C, cp$ts_method,
                                  fast = FALSE, return_full = TRUE)
    fitted_arr <- NULL
    if (return_fitted && !is.null(fc$X_log_full))
      fitted_arr <- exp(fc$X_log_full[, 1:n_t, , drop = FALSE])
    return(list(mx_forecast = fc$mx_forecast, mx_fitted = fitted_arr,
                full_fc = fc, method = "direct"))

  } else if (model_name == "Lee-Carter") {
    adj     <- if (!is.null(calib_all$lc$adjust))       calib_all$lc$adjust       else "dt"
    choosep <- if (!is.null(calib_all$lc$chooseperiod)) calib_all$lc$chooseperiod else FALSE
    fc <- lca_forecast(mx, pop, years, ages, groups, h,
                       return_fitted = return_fitted,
                       adjust = adj, chooseperiod = choosep)
    return(list(mx_forecast = fc$mx_forecast, mx_fitted = fc$mx_fitted,
                full_fc = list(mx_forecast = fc$mx_forecast, fit_pct = NA),
                method = "lca"))

  } else if (model_name == "FDM") {
    ord  <- if (!is.null(calib_all$fdm$order))  calib_all$fdm$order  else 6L
    meth <- if (!is.null(calib_all$fdm$method)) calib_all$fdm$method else "arima"
    fc <- fdm_forecast(mx, pop, years, ages, groups, h,
                       return_fitted = return_fitted,
                       order = ord, method = meth)
    return(list(mx_forecast = fc$mx_forecast, mx_fitted = fc$mx_fitted,
                full_fc = list(mx_forecast = fc$mx_forecast, fit_pct = NA),
                method = "fdm"))

  } else if (model_name == "StMoMo-APC") {
    mtype <- if (!is.null(calib_all$stmomo$model_type))
      calib_all$stmomo$model_type else "APC"
    fc <- stmomo_forecast(mx, pop, years, ages, groups, h, model_type = mtype)
    fitted_arr <- NULL
    if (return_fitted && requireNamespace("StMoMo", quietly = TRUE)) {
      fitted_arr <- array(NA_real_, dim = c(n_a, n_t, n_g))
      for (kk in seq_len(n_g)) {
        tryCatch({
          rate_mat <- pmax(mx[, , kk], 1e-12)
          pop_mat  <- pmax(pop[, , kk], 1)
          dxt <- round(rate_mat * pop_mat)
          model_s <- StMoMo::apc(link = "log")
          sfit <- StMoMo::fit(model_s, Dxt = dxt, Ext = pop_mat,
                               ages = seq_len(n_a), years = years,
                               wxt = 1 * (dxt > 0))
          fitted_arr[, , kk] <- fitted(sfit, type = "rates")
        }, error = function(e) NULL)
      }
      if (all(is.na(fitted_arr))) fitted_arr <- NULL
    }
    return(list(mx_forecast = fc$mx_forecast, mx_fitted = fitted_arr,
                full_fc = list(mx_forecast = fc$mx_forecast, fit_pct = NA),
                method = "apc"))

  } else if (model_name == "RWD") {
    dw <- if (!is.null(calib_all$rwd$drift_window))
      calib_all$rwd$drift_window else 0
    fc <- rwd_forecast(mx, years, h, drift_window = dw)
    rwd_fit <- if (return_fitted) mx else NULL
    return(list(mx_forecast = fc$mx_forecast, mx_fitted = rwd_fit,
                full_fc = list(mx_forecast = fc$mx_forecast, fit_pct = NA),
                method = "rwd"))

  } else if (model_name == "ARIMA-perSeries") {
    mp <- if (!is.null(calib_all$arima$max.p)) calib_all$arima$max.p else 2L
    mq <- if (!is.null(calib_all$arima$max.q)) calib_all$arima$max.q else 1L
    md <- if (!is.null(calib_all$arima$max.d)) calib_all$arima$max.d else 1L
    fc <- arima_forecast(mx, years, h, max.p = mp, max.q = mq, max.d = md)
    return(list(mx_forecast = fc$mx_forecast, mx_fitted = fc$mx_fitted,
                full_fc = list(mx_forecast = fc$mx_forecast, fit_pct = NA),
                method = "arima_ps"))

  } else if (model_name == "ETS-Damped") {
    ets_m <- if (!is.null(calib_all$ets$ets_model))
      calib_all$ets$ets_model else "AAN"
    ets_d <- if (!is.null(calib_all$ets$damped))
      calib_all$ets$damped else TRUE
    fc <- ets_forecast(mx, years, h, ets_model = ets_m, damped = ets_d)
    return(list(mx_forecast = fc$mx_forecast, mx_fitted = fc$mx_fitted,
                full_fc = list(mx_forecast = fc$mx_forecast, fit_pct = NA),
                method = "ets_damped"))
  }

  # Fallback
  cp <- calib_direct
  fc <- tucker3_direct_forecast(mx, years, h,
                                cp$R_A, cp$R_B, cp$R_C, cp$ts_method,
                                fast = FALSE, return_full = TRUE)
  list(mx_forecast = fc$mx_forecast, mx_fitted = NULL,
       full_fc = fc, method = "direct")
}

#' Generate bootstrap simulations for any model
#'
#' Routes to the appropriate bootstrap function based on the model's
#' method type.
#'
#' @param model_name Model name string.
#' @param mx         3D smoothed mortality array.
#' @param pop        3D population array.
#' @param years,ages Integer vectors.
#' @param groups     data.frame of groups.
#' @param h          Forecast horizon.
#' @param n_sim      Number of simulations.
#' @param fc_result  Output of [dispatch_forecast()].
#' @param width      Age-group width.
#' @param calib_all  Calibration result for non-Tucker3 models.
#' @param mx_raw     Optional raw (unsmoothed) rates.
#' @return 4D array \[n_a x h x n_g x n_sim\].
#' @export
dispatch_model_sims <- function(model_name, mx, pop, years, ages, groups,
                                h, n_sim, fc_result,
                                width = 5L, calib_all = NULL,
                                mx_raw = NULL) {
  n_a <- length(ages); n_t <- length(years); n_g <- nrow(groups)
  full_fc <- fc_result$full_fc
  method  <- fc_result$method

  if (method %in% c("coda", "direct")) {
    # Tucker3 bootstrap
    A <- full_fc$A; B <- full_fc$B; C <- full_fc$C; G_mat <- full_fc$G_mat
    if (is.null(A) || is.null(B) || is.null(C) || is.null(G_mat))
      stop("Tucker3 bootstrap requires A, B, C, G_mat in fc_result$full_fc. ",
           "Ensure forecast was produced with a Tucker3 model (return_full = TRUE).")
    n_t_b <- full_fc$n_t; R_B <- ncol(B)

    if (!is.null(mx_raw) && method == "coda") {
      X_clr_raw <- array(NA_real_, dim = c(n_a, n_t, n_g))
      for (jj in seq_len(n_t))
        for (kk in seq_len(n_g)) {
          dx_r <- compute_dx(mx_raw[, jj, kk], ages, width)
          X_clr_raw[, jj, kk] <- as.numeric(compositions::clr(dx_r))
        }
      resid <- X_clr_raw - reconstruct_array(A, B, C, G_mat)
    } else if (!is.null(mx_raw) && method == "direct") {
      resid <- log(pmax(mx_raw, 1e-12)) - reconstruct_array(A, B, C, G_mat)
    } else {
      resid <- full_fc$resid_array
    }

    # Winsorize
    for (kk in seq_len(n_g)) for (ii in seq_len(n_a)) {
      rv <- resid[ii, , kk]; rv <- rv[is.finite(rv)]
      if (length(rv) > 2) {
        mu_r <- mean(rv); sd_r <- sd(rv)
        resid[ii, , kk] <- pmax(pmin(resid[ii, , kk],
                                     mu_r + 2 * sd_r), mu_r - 2 * sd_r)
      }
    }

    sims <- array(NA_real_, dim = c(n_a, h, n_g, n_sim))
    for (s in seq_len(n_sim)) {
      B_sim <- matrix(NA_real_, h, R_B)
      for (r in seq_len(R_B))
        B_sim[, r] <- tryCatch(
          simulate(full_fc$arima_fits[[r]], nsim = h,
                   future = TRUE, bootstrap = TRUE),
          error = function(e)
            as.numeric(forecast::forecast(full_fc$arima_fits[[r]], h = h)$mean) +
              stats::rnorm(h, 0, 0.1))
      B_full_s <- rbind(B, B_sim)
      X_sim <- reconstruct_array(A, B_full_s, C, G_mat)
      for (j in seq_len(h)) {
        ry <- sample(1:n_t_b, 1)
        if (method == "coda") {
          for (k in seq_len(n_g)) {
            clr_n <- X_sim[, n_t_b + j, k] + resid[, ry, k]
            dx_s <- inv_clr_to_dx(clr_n)
            sims[, j, k, s] <- dx_to_mx(dx_s, ages, width)
          }
        } else {
          for (k in seq_len(n_g))
            sims[, j, k, s] <- exp(X_sim[, n_t_b + j, k] + resid[, ry, k])
        }
      }
    }
    return(sims)

  } else if (method == "lca") {
    adj     <- if (!is.null(calib_all$lc$adjust))
      calib_all$lc$adjust else "dt"
    choosep <- if (!is.null(calib_all$lc$chooseperiod))
      calib_all$lc$chooseperiod else FALSE
    return(bootstrap_lc_fdm(mx, pop, years, ages, groups, h, n_sim,
                            model_type = "lca", adjust = adj,
                            chooseperiod = choosep, mx_raw = mx_raw))

  } else if (method == "fdm") {
    ord  <- if (!is.null(calib_all$fdm$order))  calib_all$fdm$order  else 6L
    meth <- if (!is.null(calib_all$fdm$method)) calib_all$fdm$method else "arima"
    return(bootstrap_lc_fdm(mx, pop, years, ages, groups, h, n_sim,
                            model_type = "fdm", order = ord, method = meth,
                            mx_raw = mx_raw))

  } else if (method == "apc") {
    sims <- bootstrap_stmomo(mx, pop, years, ages, groups, h, n_sim,
                             mx_raw = mx_raw)
    if (!is.null(sims)) return(sims)
    # Fallback: parametric
    pt <- fc_result$mx_forecast; log_pt <- log(pmax(pt, 1e-12))
    sims <- array(NA_real_, dim = c(n_a, h, n_g, n_sim))
    for (s in seq_len(n_sim))
      sims[, , , s] <- exp(log_pt + stats::rnorm(length(log_pt), 0, 0.1))
    return(sims)

  } else {
    # rwd / arima_ps / ets_damped: per-series residual bootstrap
    sims <- array(NA_real_, dim = c(n_a, h, n_g, n_sim))
    for (i in seq_len(n_a)) for (k in seq_len(n_g)) {
      ls <- log(pmax(mx[i, , k], 1e-12))
      ts_r <- stats::ts(ls, start = min(years), frequency = 1)
      if (method == "rwd") {
        fit_r <- tryCatch(
          forecast::Arima(ts_r, order = c(0, 1, 0), include.drift = TRUE),
          error = function(e) forecast::auto.arima(ts_r))
      } else if (method == "ets_damped") {
        ets_m <- if (!is.null(calib_all$ets$ets_model))
          calib_all$ets$ets_model else "AAN"
        ets_d <- if (!is.null(calib_all$ets$damped))
          calib_all$ets$damped else TRUE
        fit_r <- tryCatch(
          forecast::ets(ts_r, model = ets_m, damped = ets_d),
          error = function(e) forecast::auto.arima(ts_r))
      } else {
        mp <- if (!is.null(calib_all$arima$max.p)) calib_all$arima$max.p else 2L
        mq <- if (!is.null(calib_all$arima$max.q)) calib_all$arima$max.q else 1L
        md <- if (!is.null(calib_all$arima$max.d)) calib_all$arima$max.d else 1L
        fit_r <- forecast::auto.arima(ts_r, max.p = mp, max.q = mq, max.d = md)
      }
      fc_r <- forecast::forecast(fit_r, h = h)
      pt_log <- as.numeric(fc_r$mean)
      if (!is.null(mx_raw)) {
        resid_r <- log(pmax(mx_raw[i, , k], 1e-12)) - as.numeric(fitted(fit_r))
      } else {
        resid_r <- as.numeric(stats::residuals(fit_r))
      }
      resid_r <- resid_r[is.finite(resid_r)]
      if (length(resid_r) == 0) resid_r <- 0
      if (length(resid_r) > 2) {
        mu_r <- mean(resid_r); sd_r <- sd(resid_r)
        resid_r <- pmax(pmin(resid_r, mu_r + 2 * sd_r), mu_r - 2 * sd_r)
      }
      for (s in seq_len(n_sim)) {
        ry <- sample(resid_r, h, replace = TRUE)
        sims[i, , k, s] <- exp(pt_log + ry)
      }
    }
    return(sims)
  }
}
