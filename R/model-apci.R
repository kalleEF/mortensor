#' Age-period-cohort-improvement (APCI) model
#'
#' The stochastic implementation of the APCI model of Richards, Currie,
#' Kleinow and Ritchie (2019), fitted per stratum:
#' `log m[x, y] = alpha[x] + beta[x] (y - ybar) + kappa[y] + gamma[y - x]`
#' (eq. 8), with `ybar` the mean of the years of the data. The model needs
#' at least five single years of age, without an open interval, and at
#' least ten years (package-specific, for the choice of the ARIMA models).
#'
#' Following the paper:
#'
#' * the deaths are Poisson with mean `m[x, y]` times the central exposure
#'   (eq. 4), and the parameters are estimated by maximum likelihood
#'   (Section 3);
#' * cohorts with four or fewer observations get `gamma = 0` and are not
#'   estimated (Appendix B);
#' * five identifiability constraints (Appendix C): `sum(kappa) = 0` (eq.
#'   17), `sum((y - ymin) kappa) = 0` (eq. 21), and `sum(w gamma) = 0`,
#'   `sum(w (c - cmin + 1) gamma) = 0` and `sum(w (c - cmin + 1)^2 gamma) = 0`
#'   (eq. 19, 20 and 22), with `w` the number of cells of cohort `c`. With
#'   the corner cohorts fixed at 0, only the constraints on `kappa` are
#'   needed; the constraints on `gamma` restrict the fit (the paper's
#'   over-constrained model, Section 6). The estimates maximise the
#'   likelihood subject to all five constraints, which hold exactly. For
#'   the paper's 55 ages and 45 years, the unsmoothed model then has its
#'   effective dimension of 241 (Table 2);
#' * with `smooth = TRUE` (the paper's APCI(S), Table 1), `alpha` and `beta`
#'   are P-splines over age (Eilers and Marx 1996) with penalties on the
#'   second-order differences of their coefficients (Appendix D; Appendix G
#'   calls the penalties quadratic), and the two smoothing parameters
#'   minimise the BIC, deviance plus `log(n)` times the effective dimension,
#'   with `n` the number of cells (eq. 24). `kappa` and `gamma`
#'   are not smoothed. With `smooth = FALSE`, every age has its own `alpha`
#'   and `beta`;
#' * `kappa` and `gamma` are projected by ARIMA models without a mean,
#'   fitted with [stats::arima()] to the estimates as if they were observed
#'   (Section 5, Appendix E). By default the orders are chosen as in Appendix
#'   E: for `kappa`, the ARMA orders `p, q` in 0 to 3 with the smallest
#'   AICc for `d = 1` and for `d = 2`, then of these two the model with the
#'   smaller root mean squared error (Solo 1984); for `gamma`, the orders `p,
#'   q` in 0 to 3 with the smallest AICc for `d = 1`. `gamma` is projected
#'   for the cohorts born after the last estimated one, which include the
#'   youngest cohorts with four or fewer observations (Appendix B). The two
#'   processes are independent;
#' * draws are simulated paths of both ARIMA processes, with their
#'   coefficients at the estimates. The paper adds parameter risk by a
#'   bootstrap of the mean of the process for `kappa` only, and the APCI
#'   processes have no mean (Appendix E). `alpha`, `beta` and the estimated
#'   `kappa` and `gamma` are fixed, so the deterministic part `beta[x] (y -
#'   ybar)` continues linearly. Predictive draws
#'   (`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
#'   given the rates.
#'
#' Package-specific (reference unknown): the B-splines are cubic with
#' equally spaced knots about five years of age apart, on the range of the
#' ages widened by 1% at each end (as [p_splines()]; the paper does not give
#' the knots); the smoothing parameters are searched between 1e-4 and 1e10
#' with the search of [p_splines()]; the estimates come from
#' penalised iteratively reweighted least squares with the constraints as
#' Lagrange conditions (the paper: the algorithm of Currie 2013, with the
#' constraints as rows of an augmented model matrix); the AICc counts the
#' ARMA coefficients and the innovation variance as parameters, on the
#' differenced series, as `forecast::Arima()`; the root mean squared error
#' that chooses `d` is taken over the one-step residuals from the third
#' year, which both orders of differencing have; ARIMA paths are simulated
#' from the state-space form of the fitted model, from the filtered state at
#' the last estimate. The paper's value-at-risk framework and annuity
#' factors are not implemented.
#'
#' @param smooth Smooth `alpha` and `beta` over age with P-splines, the
#'   smoothing parameters chosen by BIC (Richards et al. 2019, Table 1).
#' @param kappa_order,gamma_order ARIMA orders `c(p, d, q)` of the models of
#'   `kappa` and `gamma` (without a mean). `NULL` (default) chooses them as
#'   in Richards et al. (2019, Appendix E).
#' @return A model specification.
#' @references
#' Richards SJ, Currie ID, Kleinow T, Ritchie GP (2019). A stochastic
#' implementation of the APCI model for mortality projections. *British
#' Actuarial Journal* 24: e13. \doi{10.1017/S1357321718000260}
#'
#' Eilers PHC, Marx BD (1996). Flexible smoothing with B-splines and
#' penalties. *Statistical Science* 11(2): 89-121.
#' \doi{10.1214/ss/1038425655}
#'
#' Currie ID (2013). Smoothing constrained generalized linear models with an
#' application to the Lee-Carter model. *Statistical Modelling* 13(1):
#' 69-93. \doi{10.1177/1471082X12471373}
#'
#' Solo V (1984). The order of differencing in ARIMA models. *Journal of the
#' American Statistical Association* 79(388): 916-921.
#' \doi{10.1080/01621459.1984.10477111}
#' @seealso [stmomo()] for the age-period-cohort and Renshaw-Haberman
#'   models.
#' @examples
#' \donttest{
#' if (requireNamespace("StMoMo", quietly = TRUE)) {
#'   ew <- mt_subset(as_mt_data(StMoMo::EWMaleData), ages = c(55, 90), periods = c(1971, 2011))
#'   fit <- mt_fit(apci(), ew)
#'   p <- mt_parameters(fit)
#'   p[parameter %in% c("order_kappa", "order_gamma")]
#'   fc <- mt_forecast(fit, h = 10, n_draws = 200, seed = 1)
#'   life_expectancy(fc, at = 65)
#' }
#' }
#' @export
apci <- function(smooth = TRUE, kappa_order = NULL, gamma_order = NULL) {
  check_flag(smooth, "smooth")
  check_order <- function(o, name) {
    if (is.null(o)) return(NULL)
    if (!is.numeric(o) || length(o) != 3 || anyNA(o) || any(o < 0) || any(o != round(o))) {
      mt_abort(sprintf("`%s` must be an ARIMA order c(p, d, q) of non-negative whole numbers.", name),
               class = "argument", call = NULL)
    }
    as.integer(o)
  }
  kappa_order <- check_order(kappa_order, "kappa_order")
  gamma_order <- check_order(gamma_order, "gamma_order")
  new_mt_spec("apci", list(smooth = smooth, kappa_order = kappa_order, gamma_order = gamma_order),
              min_periods = 10, label = paste0("APCI", if (!smooth) " (unsmoothed)" else ""))
}

#' @export
fit_model.mt_spec_apci <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  if (any(ages$upper - ages$lower != 1)) {
    mt_abort(paste("apci() needs single years of age, without an open interval; regroup the ages, for example",
                   "with mt_subset(ages = ) or mt_regroup()."), class = "model_data", call = NULL)
  }
  if (length(ages) < 5) {
    mt_abort("apci() needs at least five single years of age (cohorts with five or more observations).",
             class = "model_data", call = NULL)
  }
  est <- apci_estimate(D, E, ages$lower, periods, smooth = p$smooth)
  est$kappa_ts <- apci_arima(est$kappa, p$kappa_order, d = 1:2, what = "kappa")
  est$gamma_ts <- apci_arima(est$gamma[est$estimated], p$gamma_order, d = 1, what = "gamma")
  est
}

#' @export
simulate_model.mt_spec_apci <- function(spec, fit, h, n_draws, ...) {
  nx <- length(fit$ages)
  yrs <- max(fit$years) + seq_len(h)
  c_last <- max(fit$cohorts[fit$estimated])
  coh <- outer(-fit$ages, yrs, `+`)
  known <- coh <= c_last
  # cohorts after the last estimated one, including the youngest cohorts with
  # four or fewer observations, take the projected gamma
  steps <- coh[!known] - c_last
  kp <- arima_paths(fit$kappa_ts, h, n_draws)
  gp <- arima_paths(fit$gamma_ts, max(coh) - c_last, n_draws)
  base <- fit$alpha + outer(fit$beta, yrs - fit$ybar)
  base[known] <- base[known] + fit$gamma[match(coh[known], fit$cohorts)]
  out <- array(NA_real_, c(nx, h, n_draws))
  for (j in seq_len(n_draws)) {
    eta <- base + rep(kp[, j], each = nx)
    eta[!known] <- eta[!known] + gp[steps, j]
    out[, , j] <- exp(eta)
  }
  out
}

# Model matrix (sparse, cells [age, year] with age fastest) and constraint
# matrix of the APCI model (Richards et al. 2019, eq. 8 and Appendix C).
# Coefficients: alpha (B-spline coefficients with `smooth`, else one per
# age), beta (the same), kappa (one per year), gamma (cohorts with five or
# more observations, Appendix B). Penalties: second-order differences of the
# alpha and beta coefficients.
apci_design <- function(ages, years, smooth) {
  nx <- length(ages)
  ny <- length(years)
  n <- nx * ny
  ix <- rep(seq_len(nx), ny)
  iy <- rep(seq_len(ny), each = nx)
  ybar <- mean(years)
  tt <- years[iy] - ybar
  cohorts <- seq(min(years) - max(ages), max(years) - min(ages))
  coh <- years[iy] - ages[ix]
  ic <- coh - cohorts[1] + 1L
  w <- tabulate(ic, length(cohorts))
  estimated <- w > 4
  basis <- if (smooth) ps_basis(ages, ps_ndx(ages), 3, 0.01)$B else diag(nx)
  ka <- ncol(basis)
  nc <- sum(estimated)
  cpos <- cumsum(estimated)
  sp <- function(rows, cols, vals, ncol) Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = c(n, ncol))
  nz <- which(basis[ix, , drop = FALSE] != 0, arr.ind = TRUE)
  bvals <- basis[ix, , drop = FALSE][nz]
  in_c <- estimated[ic]
  X <- cbind(sp(nz[, 1], nz[, 2], bvals, ka), sp(nz[, 1], nz[, 2], bvals * tt[nz[, 1]], ka),
             sp(seq_len(n), iy, rep(1, n), ny),
             sp(which(in_c), cpos[ic[in_c]], rep(1, sum(in_c)), nc))
  p <- ncol(X)
  idx <- list(alpha = seq_len(ka), beta = ka + seq_len(ka), kappa = 2 * ka + seq_len(ny),
              gamma = 2 * ka + ny + seq_len(nc))
  # constraints: eq. 17, 21 (kappa), 19, 20, 22 (gamma; w = cells of each
  # cohort, c - cmin + 1 counted from the oldest cohort of the data)
  cc <- seq_along(cohorts)[estimated]
  wc <- w[estimated]
  H <- matrix(0, 5, p)
  H[1, idx$kappa] <- 1
  H[2, idx$kappa] <- years - min(years)
  H[3, idx$gamma] <- wc
  H[4, idx$gamma] <- wc * cc
  H[5, idx$gamma] <- wc * cc^2
  pen <- if (smooth) diff_penalty(ka, 2) else matrix(0, ka, ka)
  list(X = X, H = H, idx = idx, basis = basis, pen = pen, ages = ages, years = years, ybar = ybar,
       cohorts = cohorts, w = w, estimated = estimated, smooth = smooth)
}

# Penalty matrix of the coefficients for smoothing parameters lambda = c(alpha, beta)
apci_penalty <- function(des, lambda) {
  p <- ncol(des$X)
  P <- matrix(0, p, p)
  if (des$smooth) {
    P[des$idx$alpha, des$idx$alpha] <- lambda[1] * des$pen
    P[des$idx$beta, des$idx$beta] <- lambda[2] * des$pen
  }
  P
}

# Penalised Poisson IRLS with the linear constraints H theta = 0 as Lagrange
# conditions: each step solves [X'WX + P, H'; H, 0] [theta; omega] = [X'Wz;
# 0], with step halving when the penalised deviance increases. `theta`:
# starting coefficients (NULL: from mu = deaths + 0.1, as glm()). Returns the
# coefficients, fitted deaths, deviance and effective dimension tr(G X'WX),
# G the top left block of the inverse of the system matrix (the trace of the
# hat matrix of the constrained fit).
apci_irls <- function(des, d, off, P, theta = NULL, maxit = 100, tol = 1e-11) {
  X <- des$X
  H <- des$H
  p <- ncol(X)
  nh <- nrow(H)
  pdev <- function(mu, th) poisson_dev(d, mu) + sum(th * (P %*% th))
  system <- function(mu) {
    A <- as.matrix(Matrix::crossprod(X, X * mu))
    list(A = A, M = rbind(cbind(A + P, t(H)), cbind(H, matrix(0, nh, nh))))
  }
  if (is.null(theta)) {
    mu <- d + 0.1
    eta <- log(mu) - off
  } else {
    eta <- as.vector(X %*% theta)
    mu <- exp(off + eta)
  }
  obj <- if (is.null(theta)) Inf else pdev(mu, theta)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    z <- eta + (d - mu) / mu
    s <- system(mu)
    rhs <- c(as.vector(Matrix::crossprod(X, mu * z)), numeric(nh))
    new <- tryCatch(solve(s$M, rhs)[seq_len(p)], error = function(e) NULL)
    if (is.null(new) || anyNA(new)) mt_abort("The APCI fit failed (singular system).", class = "model_fit")
    step <- if (is.null(theta)) new else new - theta
    base <- if (is.null(theta)) numeric(p) else theta
    for (half in 0:30) {
      th <- base + step
      eta_new <- as.vector(X %*% th)
      mu_new <- exp(off + eta_new)
      obj_new <- if (all(is.finite(mu_new))) pdev(mu_new, th) else Inf
      if (obj_new <= obj + 1e-10 * abs(obj) || !is.finite(obj)) break
      step <- step / 2
    }
    if (!is.finite(obj_new)) mt_abort("The APCI fit diverged.", class = "model_fit")
    change <- abs(obj - obj_new)
    theta <- th
    eta <- eta_new
    mu <- mu_new
    done <- is.finite(obj) && change <= tol * (abs(obj_new) + 0.1)
    obj <- obj_new
    if (done) {
      converged <- TRUE
      break
    }
  }
  if (!converged) mt_abort("The APCI fit did not converge.", class = "model_fit")
  s <- system(mu)
  G <- solve(s$M)[seq_len(p), seq_len(p)]
  list(theta = theta, mu = mu, dev = poisson_dev(d, mu), ed = sum(G * s$A), iterations = it)
}

# APCI estimates of one stratum: coefficients by constrained (penalised)
# maximum likelihood; with `smooth`, smoothing parameters minimising the BIC
# (Richards et al. 2019, eq. 24) unless `lambda` is given
apci_estimate <- function(D, E, ages, years, smooth = TRUE, lambda = NULL, lambda_range = c(1e-4, 1e10)) {
  des <- apci_design(ages, years, smooth)
  d <- as.vector(D)
  off <- log(as.vector(E))
  n <- length(d)
  state <- new.env(parent = emptyenv())
  state$theta <- NULL
  fit_at <- function(lam) {
    f <- apci_irls(des, d, off, apci_penalty(des, lam), state$theta)
    state$theta <- f$theta
    f$bic <- f$dev + log(n) * f$ed
    f$lambda <- lam
    f
  }
  if (!smooth) {
    f <- fit_at(c(0, 0))
    f$lambda <- NULL
  } else if (!is.null(lambda)) {
    f <- fit_at(lambda)
  } else {
    sel <- ps_select_lambda(function(ll) tryCatch(fit_at(10^ll)$bic, mt_error_model_fit = function(e) Inf),
                            log10(lambda_range))
    f <- fit_at(10^sel$ll)
  }
  th <- f$theta
  idx <- des$idx
  gamma <- numeric(length(des$cohorts))
  gamma[des$estimated] <- th[idx$gamma]
  list(alpha = as.vector(des$basis %*% th[idx$alpha]), beta = as.vector(des$basis %*% th[idx$beta]),
       kappa = th[idx$kappa], gamma = gamma, cohorts = des$cohorts, cohort_cells = des$w,
       estimated = des$estimated, ages = ages, years = years, ybar = des$ybar, smooth = smooth,
       lambda = f$lambda, dev = f$dev, ed = f$ed, bic = f$dev + log(n) * f$ed,
       fitted = matrix(f$mu, nrow(D)) / E)
}

# ARIMA model without a mean of an estimated index (Richards et al. 2019,
# Appendix E): the given `order`, or for each d in `d` the orders p, q in
# 0:3 with the smallest AICc, then the d with the smallest root mean squared
# error of the one-step residuals (Solo 1984), taken from the time point
# max(d) + 1 that all candidates have
apci_arima <- function(x, order = NULL, d = 1:2, what = "kappa") {
  fit_one <- function(ord) {
    f <- tryCatch(suppressWarnings(stats::arima(x, order = ord, include.mean = FALSE)), error = function(e) NULL)
    if (is.null(f) || !is.finite(f$aic) || !is.finite(f$sigma2)) return(NULL)
    k <- ord[1] + ord[3] + 1
    f$aicc <- if (f$nobs - k - 1 > 0) f$aic + 2 * k * (k + 1) / (f$nobs - k - 1) else Inf
    f$order <- ord
    f
  }
  if (!is.null(order)) {
    f <- fit_one(order)
    if (is.null(f)) {
      mt_abort(sprintf("The ARIMA(%s) model of %s could not be fitted.", paste(order, collapse = ","), what),
               class = "model_fit", call = NULL)
    }
    return(f)
  }
  best <- lapply(d, function(dd) {
    cand <- Filter(Negate(is.null), lapply(seq_len(16) - 1L, function(i) fit_one(c(i %/% 4, dd, i %% 4))))
    aicc <- vapply(cand, `[[`, 0, "aicc")
    if (!length(cand) || !any(is.finite(aicc))) return(NULL)
    cand[[which.min(aicc)]]
  })
  best <- Filter(Negate(is.null), best)
  if (!length(best)) {
    mt_abort(sprintf("No ARIMA model of %s could be fitted (series of %d values).", what, length(x)),
             class = "model_fit", call = NULL)
  }
  if (length(best) == 1L) return(best[[1]])
  from <- max(d) + 1L
  rmse <- vapply(best, function(f) sqrt(mean(stats::residuals(f)[from:length(x)]^2)), 0)
  best[[which.min(rmse)]]
}

# Simulated future paths [h, n_draws] of a fitted stats::arima() model,
# conditional on the series: the state-space form of the fit (state
# transition T, observation vector Z, innovation loading R with V = R R')
# from the filtered state at the end of the series and its covariance
# (`model$a`, sigma2 `model$P`, as predict.Arima() uses them), with normal
# innovations of variance sigma2. Means and variances agree with
# predict.Arima().
arima_paths <- function(fit, h, n_draws) {
  mod <- fit$model
  s2 <- fit$sigma2
  rd <- length(mod$a)
  P <- s2 * (mod$P + t(mod$P)) / 2
  ev <- eigen(P, symmetric = TRUE)
  L <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), rd)
  R <- mod$V[, 1]
  A <- mod$a + L %*% matrix(stats::rnorm(rd * n_draws), rd)
  out <- matrix(NA_real_, h, n_draws)
  for (j in seq_len(h)) {
    A <- mod$T %*% A + outer(R, stats::rnorm(n_draws, 0, sqrt(s2)))
    out[j, ] <- crossprod(mod$Z, A)
  }
  out
}

# Parameters of APCI fits, for mt_parameters(): alpha, beta (ages), kappa
# (years), gamma (years of birth; 0 for the cohorts with four or fewer
# observations), the mean year, the smoothing parameters, effective
# dimension and BIC, and the ARIMA models of kappa and gamma
apci_parameters <- function(fit) {
  lab <- fit$strata$label
  ages <- fit$age$labels
  row <- param_rows
  arima_rows <- function(s, f, name) {
    cf <- stats::coef(f)
    rbind(row(s, paste0("order_", name), c("p", "d", "q"), f$order),
          if (length(cf)) row(s, paste0("arima_", name), names(cf), cf),
          row(s, paste0("sd_", name), "", sqrt(f$sigma2)))
  }
  rbindlist(lapply(seq_along(fit$fits), function(s) {
    f <- fit$fits[[s]]
    if (is.null(f)) return(NULL)
    rbind(row(lab[s], "alpha", ages, f$alpha), row(lab[s], "beta", ages, f$beta),
          row(lab[s], "kappa", f$years, f$kappa), row(lab[s], "gamma", f$cohorts, f$gamma),
          row(lab[s], "ybar", "", f$ybar),
          if (!is.null(f$lambda)) row(lab[s], c("lambda_alpha", "lambda_beta"), "", f$lambda),
          row(lab[s], c("ed", "bic"), "", c(f$ed, f$bic)),
          arima_rows(lab[s], f$kappa_ts, "kappa"), arima_rows(lab[s], f$gamma_ts, "gamma"))
  }))
}
