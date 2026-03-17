#' Fit Kannisto logistic model for old-age mortality extrapolation
#'
#' Fits mu(x) = a*exp(b*x) / (1 + a*exp(b*x)) to grouped mortality
#' rates via least-squares on the log scale.
#'
#' @param mx   Numeric vector of grouped mortality rates.
#' @param ages Numeric vector of age-group midpoints.
#' @return Named list with `a`, `b`, and `predict` function.
#' @keywords internal
fit_kannisto <- function(mx, ages) {
  log_mx <- log(pmax(mx, 1e-12))
  start_b <- max(0.05, diff(range(log_mx)) / diff(range(ages)))
  start_a <- exp(mean(log_mx) - start_b * mean(ages))
  tryCatch({
    fit <- stats::optim(
      par = c(log_a = log(start_a), b = start_b),
      fn = function(par) {
        a <- exp(par[1]); b <- par[2]
        mu <- a * exp(b * ages) / (1 + a * exp(b * ages))
        sum((log(mu) - log_mx)^2)
      },
      method = "L-BFGS-B",
      lower = c(-30, 0.01), upper = c(5, 0.5))
    a <- exp(fit$par[1]); b <- fit$par[2]
    list(a = a, b = b,
         predict = function(x) a * exp(b * x) / (1 + a * exp(b * x)))
  }, error = function(e) {
    fit_lm <- stats::lm(log_mx ~ ages)
    list(predict = function(x)
      exp(stats::predict(fit_lm, newdata = data.frame(ages = x))))
  })
}

#' Disaggregate grouped mortality to single-year ages (Hyman + Kannisto)
#'
#' Primary disaggregation method using a Hyman monotone spline for closed
#' age groups and Kannisto logistic extrapolation for open-ended 90+ group.
#' Iterative smooth correction ensures grouped means are preserved.
#'
#' @param mx_grouped  Numeric vector of grouped mortality rates.
#' @param age_starts  Integer vector of group lower bounds (e.g., 30, 35, ..., 90).
#' @param target_ages Integer vector of single-year target ages (default 30:100).
#' @param width       Width of closed age groups (default 5).
#' @return Numeric vector of single-year mortality rates.
#' @export
disaggregate_hyman_kannisto <- function(mx_grouped, age_starts,
                                       target_ages = 30:100, width = 5) {
  n_grp   <- length(mx_grouped)
  n_closed <- n_grp - 1
  log_mx  <- log(pmax(mx_grouped, 1e-12))
  midpoints_closed <- age_starts[1:n_closed] + width / 2

  # Fit Kannisto on ages 60+
  kann_idx <- which(age_starts[1:n_closed] >= 60)
  if (length(kann_idx) < 3) kann_idx <- max(1, n_closed - 3):n_closed
  kann <- fit_kannisto(mx_grouped[kann_idx], midpoints_closed[kann_idx])

  # Virtual Kannisto points for 90+ region

  virt_ages <- c(92.5, 97.5, 102.5)
  virt_mx   <- kann$predict(virt_ages)
  virt_mx   <- pmax(virt_mx, mx_grouped[n_grp])

  all_mids  <- c(midpoints_closed, virt_ages)
  all_logmx <- c(log_mx[1:n_closed], log(virt_mx))

  # Hyman monotone spline if approximately monotone
  if (all(diff(all_logmx) >= -0.05)) {
    for (ii in 2:length(all_logmx))
      if (all_logmx[ii] <= all_logmx[ii - 1])
        all_logmx[ii] <- all_logmx[ii - 1] + 1e-6
    spline_fn <- stats::splinefun(all_mids, all_logmx, method = "hyman")
  } else {
    spline_fn <- stats::splinefun(all_mids, all_logmx, method = "natural")
  }
  mx_single <- exp(spline_fn(target_ages + 0.5))

  # Smooth iterative correction to match group means
  all_group_mids <- c(midpoints_closed, age_starts[n_grp] + 5)
  for (iter in 1:10) {
    factors <- numeric(n_grp)
    for (g in seq_along(age_starts)) {
      lower <- age_starts[g]
      upper <- if (g < n_grp) age_starts[g + 1] - 1 else max(target_ages)
      idx   <- which(target_ages >= lower & target_ages <= upper)
      if (length(idx) > 0 && mx_grouped[g] > 0)
        factors[g] <- mx_grouped[g] / mean(mx_single[idx])
      else
        factors[g] <- 1
    }
    if (max(abs(factors - 1)) < 1e-8) break
    log_f   <- log(factors)
    corr_fn <- stats::splinefun(all_group_mids, log_f, method = "natural")
    corr    <- exp(corr_fn(target_ages + 0.5))
    mx_single <- mx_single * corr
  }
  pmax(mx_single, 1e-12)
}

#' Disaggregate grouped mortality to single-year ages (Natural + Gompertz)
#'
#' Alternative disaggregation using natural cubic spline interpolation
#' with Gompertz virtual extension for old ages.
#'
#' @inheritParams disaggregate_hyman_kannisto
#' @return Numeric vector of single-year mortality rates.
#' @export
disaggregate_natural_gompertz <- function(mx_grouped, age_starts,
                                         target_ages = 30:100, width = 5) {
  n_grp   <- length(mx_grouped)
  n_closed <- n_grp - 1
  log_mx  <- log(pmax(mx_grouped, 1e-12))
  midpoints <- age_starts[1:n_closed] + width / 2

  if (n_closed >= 3) {
    slope <- log_mx[n_closed] - log_mx[n_closed - 1]
    midpoints <- c(midpoints, midpoints[n_closed] + width,
                   midpoints[n_closed] + 2 * width)
    log_mx_ext <- c(log_mx[1:n_closed], log_mx[n_closed] + slope,
                    log_mx[n_closed] + 2 * slope)
  } else {
    log_mx_ext <- log_mx[1:n_closed]
  }
  spline_fn <- stats::splinefun(midpoints, log_mx_ext, method = "natural")
  mx_single <- exp(spline_fn(target_ages + 0.5))

  # Smooth iterative correction
  all_group_mids <- c(age_starts[1:n_closed] + width / 2, age_starts[n_grp] + 5)
  for (iter in 1:10) {
    factors <- numeric(n_grp)
    for (g in seq_along(age_starts)) {
      lower <- age_starts[g]
      upper <- if (g < n_grp) age_starts[g + 1] - 1 else max(target_ages)
      idx   <- which(target_ages >= lower & target_ages <= upper)
      if (length(idx) > 0 && mx_grouped[g] > 0)
        factors[g] <- mx_grouped[g] / mean(mx_single[idx])
      else
        factors[g] <- 1
    }
    if (max(abs(factors - 1)) < 1e-8) break
    log_f   <- log(factors)
    corr_fn <- stats::splinefun(all_group_mids, log_f, method = "natural")
    corr    <- exp(corr_fn(target_ages + 0.5))
    mx_single <- mx_single * corr
  }
  pmax(mx_single, 1e-12)
}

#' Check disaggregation quality
#'
#' Verifies group-mean preservation, monotonicity from age 50+,
#' and smoothness (second differences).
#'
#' @param mx_grouped Grouped rates.
#' @param mx_single  Single-year rates.
#' @param age_starts Group start ages.
#' @param target_ages Single-year ages.
#' @return Named list with `max_rel_error`, `monotone`, `max_jump`.
#' @export
check_disaggregation <- function(mx_grouped, mx_single, age_starts, target_ages) {
  n_grp <- length(mx_grouped)
  errors <- numeric(n_grp)
  for (i in seq_along(age_starts)) {
    lower <- age_starts[i]
    upper <- if (i < n_grp) age_starts[i + 1] - 1 else max(target_ages)
    idx <- which(target_ages >= lower & target_ages <= upper)
    if (length(idx) > 0)
      errors[i] <- abs(mean(mx_single[idx]) - mx_grouped[i]) /
                   max(mx_grouped[i], 1e-12)
  }
  idx50 <- which(target_ages >= 50)
  mono_ok <- all(diff(log(mx_single[idx50])) > -0.1)
  d2 <- diff(log(mx_single), differences = 2)
  list(max_rel_error = max(errors),
       monotone = mono_ok,
       max_jump = max(abs(d2)))
}


#' Gompertz extrapolation: log-linear in age
#'
#' @param anchor_mx   Numeric vector of mortality rates at anchor ages.
#' @param anchor_ages Numeric vector of anchor ages.
#' @param extrap_ages Numeric vector of ages to extrapolate to.
#' @return Numeric vector of extrapolated mortality rates.
#' @export
extrapolate_gompertz <- function(anchor_mx, anchor_ages, extrap_ages) {
  log_mx <- log(pmax(anchor_mx, 1e-12))
  fit <- stats::lm(log_mx ~ anchor_ages)
  exp(stats::predict(fit, newdata = data.frame(anchor_ages = extrap_ages + 0.5)))
}


#' Blend GAM predictions with a parametric tail in log-space
#'
#' The GAM predicts all target ages (including beyond the open-interval
#' midpoint).
#' A Kannisto or Gompertz model fitted on GAM predictions at ages 70+
#' provides stable old-age behaviour.
#' Logistic weights produce a C-inf smooth transition:
#' \eqn{w(x) = 1 / (1 + exp(-(x - open\_mid) / \tau))}.
#'
#' @param mx_gam      Numeric vector of GAM-predicted mx for all target ages.
#' @param target_ages Integer vector of single-year ages.
#' @param open_mid    Centre of the blend (open-interval midpoint, e.g.\
#'   95).
#' @param tau         Logistic blend width (default 2; larger = wider
#'   transition).
#' @param method      \code{"kannisto"} or \code{"gompertz"}.
#' @return Named list with \code{mx_blended}, \code{mx_param},
#'   \code{blend_weights}.
#' @export
blend_parametric_tail <- function(mx_gam, target_ages, open_mid,
                                  tau = 2, method = "kannisto") {
  anchor_idx <- which(target_ages >= 70 & target_ages + 0.5 <= open_mid + 1)
  if (length(anchor_idx) < 5) anchor_idx <- utils::tail(seq_along(target_ages), 15)
  anchor_ages <- target_ages[anchor_idx]
  anchor_mx   <- mx_gam[anchor_idx]

  if (method == "kannisto") {
    kann <- tryCatch(fit_kannisto(anchor_mx, anchor_ages + 0.5),
                     error = function(e) NULL)
    if (!is.null(kann)) {
      mx_param <- kann$predict(target_ages + 0.5)
    } else {
      mx_param <- extrapolate_gompertz(anchor_mx, anchor_ages, target_ages)
    }
  } else {
    mx_param <- extrapolate_gompertz(anchor_mx, anchor_ages, target_ages)
  }
  mx_param <- pmax(mx_param, 1e-12)

  w <- 1 / (1 + exp(-(target_ages + 0.5 - open_mid) / tau))

  log_blended <- (1 - w) * log(pmax(mx_gam, 1e-12)) + w * log(mx_param)
  mx_blended <- exp(log_blended)

  d_log <- diff(log(mx_blended))
  old_idx <- which(target_ages[-1] >= 80)
  if (length(old_idx) > 0 && any(d_log[old_idx] < -0.05))
    message("  Note: log(mx) decrease > 0.05 at ages 80+ in blended curve")

  list(mx_blended = mx_blended, mx_param = mx_param, blend_weights = w)
}


#' Smooth iterative group-mean correction
#'
#' Adjusts single-year rates so that the mean within each 5-year interval
#' matches the original grouped rate.
#' Uses natural spline interpolation of log-correction factors to avoid
#' discontinuous jumps at group boundaries.
#'
#' @param mx_single   Numeric vector of single-year mx.
#' @param mx_grouped  Numeric vector of grouped mx.
#' @param target_ages Integer vector of single-year ages.
#' @param age_starts  Integer vector of age-group starts.
#' @param width       Age-group width (default 5).
#' @param max_age     Maximum target age.
#' @return Corrected single-year mx vector.
#' @export
correct_group_means <- function(mx_single, mx_grouped, target_ages,
                                age_starts, width, max_age) {
  n_grp <- length(age_starts)
  grp_mids <- c(age_starts[-n_grp] + width / 2, age_starts[n_grp] + 5)

  for (iter in 1:15) {
    factors <- numeric(n_grp)
    for (g in seq_along(age_starts)) {
      lower <- age_starts[g]
      upper <- if (g < n_grp) age_starts[g + 1] - 1 else max_age
      idx <- which(target_ages >= lower & target_ages <= upper)
      if (length(idx) > 0 && mx_grouped[g] > 0)
        factors[g] <- mx_grouped[g] / mean(mx_single[idx])
      else
        factors[g] <- 1
    }
    max_err <- max(abs(factors - 1))
    if (max_err < 1e-8) break

    log_f   <- log(factors)
    corr_fn <- stats::splinefun(grp_mids, log_f, method = "natural")
    corr    <- exp(corr_fn(target_ages + 0.5))
    mx_single <- mx_single * corr
  }
  pmax(mx_single, 1e-12)
}


#' Disaggregate prediction interval bounds using ratio interpolation
#'
#' Computes ratio of lower/upper PI to point forecast at grouped level,
#' interpolates to single-year ages via natural splines, and applies to
#' single-year point forecast.
#'
#' @param mx_1yr_point  Numeric vector of single-year point forecast.
#' @param mx_grp_point  Numeric vector of grouped point forecast.
#' @param mx_grp_lo     Numeric vector of grouped lower 95 percent PI.
#' @param mx_grp_hi     Numeric vector of grouped upper 95 percent PI.
#' @param mx_grp_med    Numeric vector of grouped median.
#' @param target_ages   Integer vector of single-year ages.
#' @param age_starts    Integer vector of age-group starts.
#' @param width         Age-group width.
#' @return Named list with \code{mx_1yr}, \code{mx_1yr_lo95},
#'   \code{mx_1yr_hi95}, \code{mx_1yr_median}.
#' @export
disaggregate_pi_ratios <- function(mx_1yr_point, mx_grp_point,
                                   mx_grp_lo, mx_grp_hi, mx_grp_med,
                                   target_ages, age_starts, width) {
  grp_pt  <- pmax(mx_grp_point, 1e-15)
  r_lo    <- mx_grp_lo  / grp_pt
  r_hi    <- mx_grp_hi  / grp_pt
  r_med   <- mx_grp_med / grp_pt
  n_grp   <- length(age_starts)
  grp_mids <- c(age_starts[-n_grp] + width / 2, age_starts[n_grp] + 5)

  fn_lo  <- stats::splinefun(grp_mids, r_lo,  method = "natural")
  fn_hi  <- stats::splinefun(grp_mids, r_hi,  method = "natural")
  fn_med <- stats::splinefun(grp_mids, r_med, method = "natural")

  mx_1yr_lo  <- mx_1yr_point * pmin(fn_lo(target_ages + 0.5),  1)
  mx_1yr_hi  <- mx_1yr_point * pmax(fn_hi(target_ages + 0.5),  1)
  mx_1yr_med <- mx_1yr_point * fn_med(target_ages + 0.5)
  mx_1yr_med <- pmin(pmax(mx_1yr_med, mx_1yr_lo), mx_1yr_hi)

  list(mx_1yr = mx_1yr_point, mx_1yr_lo95 = mx_1yr_lo,
       mx_1yr_hi95 = mx_1yr_hi, mx_1yr_median = mx_1yr_med)
}


#' Temporally-consistent GAM-based disaggregation with logistic blend
#'
#' Fits a 2D GAM on observed grouped data (age-midpoint x year), predicts
#' single-year ages for all target ages, blends with a parametric tail
#' (Kannisto/Gompertz) via logistic weights in log-space, then applies
#' group-mean correction.
#'
#' @param mx_grouped   3D array \code{[n_a x n_t x n_g]} of grouped mx
#'   (observed).
#' @param ages         Integer vector of age-group starts.
#' @param years        Integer vector of observed calendar years.
#' @param mx_forecast  3D array \code{[n_a x h x n_g]} of grouped mx
#'   (forecast).
#' @param forecast_years Integer vector of forecast years.
#' @param config       Config list with elements \code{age_width},
#'   \code{open_interval_midpoint}, \code{extrapolation_method},
#'   \code{covid_years}, \code{output_max_age}, \code{blend_width}.
#' @return Named list with:
#' \describe{
#'   \item{mx_obs_1yr}{Array \code{[n_target x n_t x n_g]}}
#'   \item{mx_fc_1yr}{Array \code{[n_target x h x n_g]}}
#'   \item{target_ages}{Integer vector of single-year ages}
#'   \item{gam_fits}{List of fitted GAM objects}
#'   \item{diag}{List of diagnostic arrays}
#' }
#' @export
disaggregate_forecast_gam <- function(mx_grouped, ages, years,
                                      mx_forecast, forecast_years,
                                      config) {
  n_a <- dim(mx_grouped)[1]
  n_t <- dim(mx_grouped)[2]
  n_g <- dim(mx_grouped)[3]
  h   <- length(forecast_years)
  width <- config$age_width %||% 5L
  max_age <- config$output_max_age %||% 100L
  target_ages <- seq(min(ages), max_age)
  extrap_method <- config$extrapolation_method %||% "kannisto"
  open_mid <- config$open_interval_midpoint %||% 95
  blend_tau <- config$blend_width %||% 2

  n_closed <- n_a - 1
  midpoints <- c(ages[1:n_closed] + width / 2, open_mid)

  k_age  <- min(10, length(midpoints) - 1)
  k_year <- min(10, n_t - 1)
  n_target <- length(target_ages)

  obs_grid <- expand.grid(age = midpoints, year = years)
  pred_obs_grid <- expand.grid(age = target_ages + 0.5, year = years)
  pred_fc_grid  <- expand.grid(age = target_ages + 0.5, year = forecast_years)

  mx_obs_1yr <- array(NA_real_, dim = c(n_target, n_t, n_g))
  mx_fc_1yr  <- array(NA_real_, dim = c(n_target, h, n_g))
  gam_fits   <- vector("list", n_g)
  covid_years <- config$covid_years %||% integer(0)

  diag_gam_raw  <- array(NA_real_, dim = c(n_target, h, n_g))
  diag_param    <- array(NA_real_, dim = c(n_target, h, n_g))
  diag_pre_corr <- array(NA_real_, dim = c(n_target, h, n_g))
  blend_w <- 1 / (1 + exp(-(target_ages + 0.5 - open_mid) / blend_tau))

  for (k in seq_len(n_g)) {
    log_mx <- log(pmax(as.vector(mx_grouped[, , k]), 1e-10))

    w_vec <- rep(1, n_a * n_t)
    for (j in seq_len(n_t)) {
      if (years[j] %in% covid_years)
        w_vec[((j - 1) * n_a + 1):(j * n_a)] <- 0
    }
    is_floor <- as.vector(mx_grouped[, , k]) < 1e-7
    w_vec[is_floor] <- 0

    df <- data.frame(age = obs_grid$age, year = obs_grid$year,
                     log_mx = log_mx, w = w_vec)

    gam_fit <- tryCatch(
      mgcv::gam(log_mx ~ te(age, year, bs = "ps", k = c(k_age, k_year)),
                data = df, weights = w, method = "REML"),
      error = function(e) NULL)
    gam_fits[[k]] <- gam_fit

    if (is.null(gam_fit)) {
      cat("    WARNING: GAM failed for group", k,
          "- using legacy disaggregation\n")
      for (j in seq_len(n_t))
        mx_obs_1yr[, j, k] <- disaggregate_hyman_kannisto(
          mx_grouped[, j, k], ages, target_ages, width)
      for (j in seq_len(h))
        mx_fc_1yr[, j, k] <- disaggregate_hyman_kannisto(
          mx_forecast[, j, k], ages, target_ages, width)
      next
    }

    mx_obs_gam <- matrix(
      exp(stats::predict(gam_fit, newdata = pred_obs_grid)),
      nrow = n_target, ncol = n_t)
    mx_fc_gam  <- matrix(
      exp(stats::predict(gam_fit, newdata = pred_fc_grid)),
      nrow = n_target, ncol = h)

    for (j in seq_len(n_t)) {
      bl <- blend_parametric_tail(mx_obs_gam[, j], target_ages,
                                  open_mid, blend_tau, extrap_method)
      mx_obs_1yr[, j, k] <- bl$mx_blended
    }
    for (j in seq_len(h)) {
      bl <- blend_parametric_tail(mx_fc_gam[, j], target_ages,
                                  open_mid, blend_tau, extrap_method)
      mx_fc_1yr[, j, k] <- bl$mx_blended
      diag_gam_raw[, j, k]  <- mx_fc_gam[, j]
      diag_param[, j, k]    <- bl$mx_param
      diag_pre_corr[, j, k] <- bl$mx_blended
    }

    for (j in seq_len(n_t))
      mx_obs_1yr[, j, k] <- correct_group_means(
        mx_obs_1yr[, j, k], mx_grouped[, j, k],
        target_ages, ages, width, max_age)
    for (j in seq_len(h))
      mx_fc_1yr[, j, k] <- correct_group_means(
        mx_fc_1yr[, j, k], mx_forecast[, j, k],
        target_ages, ages, width, max_age)
  }

  n_floor_total <- sum(mx_grouped < 1e-7)
  if (n_floor_total > 0)
    cat("    GAM disaggregation: ", n_floor_total, " floor cells (weight=0)\n")
  cat("    GAM disaggregation: ", n_target, " single-year ages (",
      min(target_ages), "-", max(target_ages), "), ",
      n_g, " groups, extrap=", extrap_method,
      ", blend tau=", blend_tau, "\n")

  list(mx_obs_1yr  = mx_obs_1yr,
       mx_fc_1yr   = mx_fc_1yr,
       target_ages = target_ages,
       gam_fits    = gam_fits,
       diag = list(gam_raw  = diag_gam_raw,
                   param_raw = diag_param,
                   blend_weights = blend_w,
                   mx_pre_correction = diag_pre_corr))
}
