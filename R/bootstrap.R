#' Bootstrap prediction intervals for Tucker3 forecasts
#'
#' Simulates future B-mode paths via ARIMA bootstrap, adds shuffled
#' historical residuals, and reconstructs mortality rates.
#'
#' @param fc_obj  Tucker3 forecast object (with `return_full = TRUE`).
#' @param n_sim   Number of bootstrap simulations.
#' @param h       Forecast horizon.
#' @param ages    Integer vector of age-group starts.
#' @param n_groups Number of groups (dim 3).
#' @param method  `"coda"` or `"direct"`.
#' @param width   Age-group width.
#' @return List with `mx_lower`, `mx_upper`, `mx_median` (3D arrays).
#' @export
bootstrap_forecast <- function(fc_obj, n_sim, h, ages, n_groups,
                               method = "coda", width = 5) {
  A <- fc_obj$A; B <- fc_obj$B; C <- fc_obj$C; G_mat <- fc_obj$G_mat
  n_a <- nrow(A); n_t <- fc_obj$n_t; R_B <- ncol(B)
  resid <- fc_obj$resid_array

  # Winsorize residuals: clip to +/-2 sigma per age-group series

  for (kk in seq_len(n_groups)) {
    for (ii in seq_len(n_a)) {
      rv <- resid[ii, , kk]
      rv <- rv[is.finite(rv)]
      if (length(rv) > 2) {
        mu_r <- mean(rv); sd_r <- sd(rv)
        resid[ii, , kk] <- pmax(pmin(resid[ii, , kk],
                                     mu_r + 2 * sd_r), mu_r - 2 * sd_r)
      }
    }
  }

  mx_sims <- array(NA_real_, dim = c(n_a, h, n_groups, n_sim))
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

    if (method == "coda") {
      X_sim <- reconstruct_array(A, B_full_sim, C, G_mat)
      for (j in seq_len(h)) {
        ry <- sample(1:n_t, 1)
        for (k in seq_len(n_groups)) {
          clr_n <- X_sim[, n_t + j, k] + resid[, ry, k]
          dx_s  <- inv_clr_to_dx(clr_n)
          mx_sims[, j, k, s] <- dx_to_mx(dx_s, ages, width)
        }
      }
    } else {
      X_sim <- reconstruct_array(A, B_full_sim, C, G_mat)
      for (j in seq_len(h)) {
        ry <- sample(1:n_t, 1)
        for (k in seq_len(n_groups))
          mx_sims[, j, k, s] <- exp(X_sim[, n_t + j, k] + resid[, ry, k])
      }
    }
  }
  list(
    mx_lower  = apply(mx_sims, c(1, 2, 3), stats::quantile, 0.025, na.rm = TRUE),
    mx_upper  = apply(mx_sims, c(1, 2, 3), stats::quantile, 0.975, na.rm = TRUE),
    mx_median = apply(mx_sims, c(1, 2, 3), stats::quantile, 0.50,  na.rm = TRUE)
  )
}

#' Bootstrap PIs for Lee-Carter / FDM models
#'
#' Uses the demography package's point forecast plus residual bootstrap.
#'
#' @param mx_array   3D array of (smoothed) mortality rates.
#' @param pop_array  3D population array.
#' @param years      Integer vector.
#' @param ages       Integer vector.
#' @param groups_dt  data.frame of groups.
#' @param h          Forecast horizon.
#' @param n_sim      Number of simulations.
#' @param model_type `"lca"` or `"fdm"`.
#' @param adjust     LC adjustment type.
#' @param chooseperiod LC choose-period flag.
#' @param order      FDM order.
#' @param method     Forecast method for time series.
#' @param mx_raw     Optional raw (unsmoothed) rates for residual computation.
#' @return 4D array \[n_a x h x n_g x n_sim\].
#' @export
bootstrap_lc_fdm <- function(mx_array, pop_array, years, ages, groups_dt, h,
                              n_sim, model_type = "lca",
                              adjust = "dt", chooseperiod = FALSE,
                              order = 6, method = "arima",
                              mx_raw = NULL) {
  if (!requireNamespace("demography", quietly = TRUE))
    stop("Package 'demography' is required for LC/FDM bootstrap.")

  n_a <- dim(mx_array)[1]; n_g <- dim(mx_array)[3]
  mx_sims <- array(NA_real_, dim = c(n_a, h, n_g, n_sim))

  for (k in seq_len(n_g)) {
    rate_mat <- pmax(mx_array[, , k], 1e-12)
    pop_mat  <- pmax(pop_array[, , k], 1)
    lbl <- tolower(paste0(groups_dt$sex[k], "_q", groups_dt$quintile[k]))
    ddata <- tryCatch(
      demography::demogdata(data = rate_mat, pop = pop_mat,
                            ages = ages, years = years,
                            type = "mortality", label = lbl, name = "total"),
      error = function(e) NULL)
    if (is.null(ddata)) next

    if (model_type == "lca") {
      m_fit <- tryCatch(
        suppressWarnings(demography::lca(ddata, series = "total",
                                         max.age = max(ages),
                                         adjust = adjust,
                                         chooseperiod = chooseperiod)),
        error = function(e) NULL)
    } else {
      order_fdm <- min(order, n_a - 1)
      m_fit <- tryCatch(
        suppressWarnings(demography::fdm(ddata, series = "total",
                                         order = order_fdm,
                                         max.age = max(ages))),
        error = function(e) NULL)
    }
    if (is.null(m_fit)) next

    fc_obj <- tryCatch(forecast::forecast(m_fit, h = h, method = method),
                       error = function(e) NULL)
    if (is.null(fc_obj)) next
    point_rates <- as.matrix(fc_obj$rate[[1]])
    point_log   <- log(pmax(point_rates, 1e-12))

    fitted_log <- m_fit$fitted$y
    if (is.list(fitted_log)) fitted_log <- as.matrix(fitted_log[[1]])
    if (!is.null(mx_raw)) {
      obs_log <- log(pmax(mx_raw[, , k], 1e-12))
    } else {
      obs_log <- log(pmax(rate_mat, 1e-12))
    }
    resid_mat <- obs_log - fitted_log
    n_yr <- ncol(resid_mat)

    # Winsorize residuals
    for (ii in seq_len(nrow(resid_mat))) {
      rv <- resid_mat[ii, ]; rv <- rv[is.finite(rv)]
      if (length(rv) > 2) {
        mu_r <- mean(rv); sd_r <- sd(rv)
        resid_mat[ii, ] <- pmax(pmin(resid_mat[ii, ],
                                     mu_r + 2 * sd_r), mu_r - 2 * sd_r)
      }
    }

    for (s in seq_len(n_sim)) {
      ry <- sample(1:n_yr, h, replace = TRUE)
      mx_sims[, , k, s] <- exp(point_log + resid_mat[, ry])
    }
  }
  mx_sims
}

#' Bootstrap PIs for StMoMo-APC models
#'
#' Uses StMoMo's built-in simulate() or falls back to residual bootstrap.
#'
#' @param mx_array  3D array of mortality rates.
#' @param pop_array 3D population array.
#' @param years     Integer vector.
#' @param ages      Integer vector.
#' @param groups_dt data.frame of groups.
#' @param h         Forecast horizon.
#' @param n_sim     Number of simulations.
#' @param mx_raw    Optional raw rates for residuals.
#' @return 4D array \[n_a x h x n_g x n_sim\] or NULL.
#' @export
bootstrap_stmomo <- function(mx_array, pop_array, years, ages, groups_dt, h,
                              n_sim, mx_raw = NULL) {
  if (!requireNamespace("StMoMo", quietly = TRUE)) return(NULL)
  n_a <- dim(mx_array)[1]; n_g <- dim(mx_array)[3]
  age_idx <- seq_len(n_a)
  mx_sims <- array(NA_real_, dim = c(n_a, h, n_g, n_sim))

  for (k in seq_len(n_g)) {
    tryCatch({
      rate_mat <- pmax(mx_array[, , k], 1e-12)
      pop_mat  <- pmax(pop_array[, , k], 1)
      dxt <- round(rate_mat * pop_mat)
      model <- StMoMo::apc(link = "log")
      sfit <- StMoMo::fit(model, Dxt = dxt, Ext = pop_mat,
                           ages = age_idx, years = years,
                           wxt = 1 * (dxt > 0))
      sim_obj <- tryCatch(simulate(sfit, nsim = n_sim, h = h),
                          error = function(e) NULL)
      if (!is.null(sim_obj) && !is.null(sim_obj$rates)) {
        for (s in seq_len(n_sim))
          mx_sims[, , k, s] <- sim_obj$rates[, , s]
      } else {
        fc <- forecast::forecast(sfit, h = h)
        point_rates <- fc$rates
        fitted_rates <- fitted(sfit, type = "rates")
        obs_for_resid <- if (!is.null(mx_raw)) pmax(mx_raw[, , k], 1e-12) else rate_mat
        resid_log <- log(pmax(obs_for_resid, 1e-12)) -
                     log(pmax(fitted_rates, 1e-12))
        n_yr <- ncol(resid_log)
        for (s in seq_len(n_sim)) {
          ry <- sample(1:n_yr, h, replace = TRUE)
          mx_sims[, , k, s] <- exp(log(pmax(point_rates, 1e-12)) + resid_log[, ry])
        }
      }
    }, error = function(e) NULL)
  }
  mx_sims
}

#' MASE-weighted ensemble bootstrap
#'
#' Combines per-model simulation arrays via weighted averaging.
#'
#' @param model_names Character vector of ensemble models.
#' @param weights     Numeric weight vector (sum to 1).
#' @param sim_arrays  Named list of 4D arrays \[n_a x h x n_g x n_sim\].
#' @param n_a,h,n_g   Dimension sizes.
#' @return List with `mx_lower`, `mx_upper`, `mx_median`, `mx_sims`.
#' @export
ensemble_bootstrap <- function(model_names, weights, sim_arrays, n_a, h, n_g) {
  n_sim <- dim(sim_arrays[[1]])[4]
  mx_ens_sims <- array(0, dim = c(n_a, h, n_g, n_sim))
  w_total <- 0
  for (ii in seq_along(model_names)) {
    nm <- model_names[ii]
    if (!is.null(sim_arrays[[nm]])) {
      sarr <- sim_arrays[[nm]]
      sarr[is.na(sarr)] <- 0
      mx_ens_sims <- mx_ens_sims + weights[ii] * sarr
      w_total <- w_total + weights[ii]
    }
  }
  if (w_total > 0 && w_total < 1) mx_ens_sims <- mx_ens_sims / w_total

  list(
    mx_lower  = apply(mx_ens_sims, c(1, 2, 3), stats::quantile, 0.025, na.rm = TRUE),
    mx_upper  = apply(mx_ens_sims, c(1, 2, 3), stats::quantile, 0.975, na.rm = TRUE),
    mx_median = apply(mx_ens_sims, c(1, 2, 3), stats::quantile, 0.50,  na.rm = TRUE),
    mx_sims   = mx_ens_sims
  )
}

#' GAM Bayesian posterior propagation prediction intervals
#'
#' Hybrid approach: calibration is done once on the posterior-mean smooth
#' (standard `mx_smooth`). Per simulation, a new smooth surface is drawn
#' from the GAM posterior and all ensemble models are re-fit on that draw,
#' yielding a weighted-average forecast. Quantile PIs are computed from
#' the resulting simulation cube.
#'
#' @param gam_fits     List of [mgcv::gam] objects, one per group.
#' @param ages         Numeric vector of age-group starts.
#' @param years        Integer vector of observation years.
#' @param pop_arr      3D population array \[n_a x n_t x n_g\].
#' @param groups       data.table with group information.
#' @param h            Forecast horizon.
#' @param n_sim        Number of posterior draws.
#' @param ens_names    Character vector of ensemble model names.
#' @param ew           Named numeric vector of ensemble weights (sum to 1).
#' @param calib_coda   Calibrated Tucker3-CoDa parameters.
#' @param calib_direct Calibrated Tucker3-Direct parameters.
#' @param calib_all    Calibrated parameters for all models.
#' @param width        Age-group width (default 5).
#' @param dispatch_fn  Forecast dispatch function (default [dispatch_forecast]).
#' @return List with `mx_lower`, `mx_upper`, `mx_median` (3D arrays).
#' @export
bootstrap_gam_bayes <- function(gam_fits, ages, years, pop_arr, groups,
                                h, n_sim, ens_names, ew,
                                calib_coda, calib_direct, calib_all,
                                width = 5L, dispatch_fn = dispatch_forecast) {
  n_a <- length(ages); n_t <- length(years); n_g <- nrow(groups)

  # Pre-compute design matrices and posterior moments per group
  grid_df <- expand.grid(age = ages, year = years)
  Xp_list <- vector("list", n_g)
  Vp_list <- vector("list", n_g)
  mu_list <- vector("list", n_g)

  for (k in seq_len(n_g)) {
    if (!is.null(gam_fits[[k]])) {
      Xp_list[[k]] <- stats::predict(gam_fits[[k]], newdata = grid_df, type = "lpmatrix")
      Vp_list[[k]] <- stats::vcov(gam_fits[[k]])
      mu_list[[k]] <- stats::coef(gam_fits[[k]])
    }
  }

  mx_ens_sims <- array(NA_real_, dim = c(n_a, h, n_g, n_sim))
  cat("  GAM Bayes: running", n_sim, "posterior-propagation draws\n")
  pb <- utils::txtProgressBar(min = 0, max = n_sim, style = 3)

  for (s in seq_len(n_sim)) {
    mx_draw <- array(NA_real_, dim = c(n_a, n_t, n_g))
    for (k in seq_len(n_g)) {
      if (!is.null(gam_fits[[k]])) {
        beta_s <- mgcv::rmvn(1, mu_list[[k]], Vp_list[[k]])
        log_mx_s <- Xp_list[[k]] %*% as.vector(beta_s)
        mx_draw[, , k] <- matrix(exp(as.numeric(log_mx_s)), nrow = n_a, ncol = n_t)
      }
    }

    fc_draw <- array(0, dim = c(n_a, h, n_g))
    w_total <- 0
    for (ii in seq_along(ens_names)) {
      enm <- ens_names[ii]
      fc_i <- tryCatch(
        dispatch_fn(enm, mx_draw, pop_arr, years, ages, groups,
                    h, calib_coda, calib_direct,
                    width = width, calib_all = calib_all,
                    return_fitted = FALSE),
        error = function(e) NULL)
      if (!is.null(fc_i) && !is.null(fc_i$mx_forecast)) {
        fa <- fc_i$mx_forecast
        fa[is.na(fa)] <- 0
        fc_draw <- fc_draw + ew[ii] * fa
        w_total <- w_total + ew[ii]
      }
    }
    if (w_total > 0 && w_total < 1) fc_draw <- fc_draw / w_total

    mx_ens_sims[, , , s] <- fc_draw
    utils::setTxtProgressBar(pb, s)
  }
  close(pb)
  cat("  GAM Bayes: done.\n")

  list(
    mx_lower  = apply(mx_ens_sims, c(1, 2, 3), stats::quantile, 0.025, na.rm = TRUE),
    mx_upper  = apply(mx_ens_sims, c(1, 2, 3), stats::quantile, 0.975, na.rm = TRUE),
    mx_median = apply(mx_ens_sims, c(1, 2, 3), stats::quantile, 0.50,  na.rm = TRUE)
  )
}
