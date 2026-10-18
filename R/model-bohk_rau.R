#' Bohk-Ewald and Rau model of age-specific rates of mortality improvement
#'
#' The probabilistic model of Bohk-Ewald and Rau (2017) forecasts rates of
#' mortality improvement (ROMIs) by age, whose trends over time differ
#' between ages, so that the age pattern of mortality decline changes over
#' the forecast. For each age `x` and year `t` the ROMI is
#'
#' `rho[x, t] = 1 - m[x, t + 1] / m[x, t]`
#'
#' (equation 1 of the paper; positive when mortality falls), computed from
#' death rates smoothed over age and time with two-dimensional P-splines
#' (Currie, Durban and Eilers 2004). One of two Bayesian core models
#' describes the ROMIs of every age over the years of the base period:
#'
#' * `core = "exponential"` (the core model of the paper's applications):
#'   `rho[x, t] ~ N(exp(b1[x] + b2[x] log t), sigma^2)`, `b1[x] ~
#'   N(theta1[x], sigma1^2)`, `b2[x] ~ N(theta2[x], sigma2^2)`, with
#'   `theta1[x]` and `theta2[x]` the least-squares estimates of `log rho[x,
#'   t] = theta1[x] + theta2[x] log t` for each age and `sigma`, `sigma1`,
#'   `sigma2` uniform on (0, 1). Improvements can level off; mortality
#'   increases cannot be forecast.
#' * `core = "linear"`: `rho[x, t] ~ N(b1[x] + b2[x] t, sigma^2)` with
#'   `(b1[x], b2[x])` bivariate normal with means `mu1`, `mu2` and covariance
#'   from standard deviations `omega1`, `omega2` and correlation `r`;
#'   uniform priors `sigma`, `omega1`, `omega2` on (0, 1), `mu1`, `mu2` on
#'   (-0.1, 0.1) and `r` on (-1, 1). Trends of the improvements continue
#'   linearly.
#'
#' The posterior and the future ROMIs (each one the curve of a posterior
#' sample plus normal noise with that sample's `sigma`) are sampled by
#' Markov chain Monte Carlo. Forecast ROMIs outside `limits` are set to the
#' nearest limit (the paper uses 0.005 and 0.035: mortality declines by at
#' least 0.5% and at most 3.5% a year; it excludes implausible forecasts "by
#' setting thresholds", which this package reads as setting the values to
#' the thresholds). Death rates follow from the last
#' smoothed rates of the base period by `m[x, t] = m[x, t - 1] (1 - rho[x,
#' t - 1])` (equation 20). With `draws = "quantiles"` (the paper's method)
#' draw `j` of `n` applies the quantile `(j - 0.5) / n` of the forecast ROMIs
#' at every age and year, so that the draws give the quantiles of the death
#' rates as in the paper. These draws are comonotone across ages and years;
#' their order is permuted at random, so that the draws of different strata
#' are independent. With `draws = "samples"` (package-specific) every draw
#' is the path of one posterior sample with its own noise. Such paths give
#' much narrower intervals of summary measures such as life expectancy,
#' because the yearly noise of the improvements averages out over ages and
#' years.
#'
#' The defaults follow the paper's application (its Table 3): the
#' exponential core model, five chains of 5,200 iterations with the first
#' 200 discarded and every fifth kept (5,000 samples), and the limits 0.005
#' and 0.035. The paper recommends a base period at least as long as the
#' forecast horizon, and uses single years of age from 0.
#'
#' **What differs from the paper.** The paper samples with JAGS; here a
#' sampler written for this package is used (Metropolis steps for the curve
#' parameters of each age, with proposals scaled by the curvature at the
#' least-squares estimates and adapted during the burn-in; conjugate draws of
#' the variances truncated to the prior range; moves that scale `sigma1` or
#' `sigma2` together with the deviations of the curve parameters from
#' `theta`; Gibbs and slice sampling steps for the linear core model). The
#' paper does not state the following, which are choices of this package
#' (package-specific, reference unknown):
#'
#' * time `t` counts the years of the base period from 1 (the change from
#'   the first to the second year) and continues into the forecast;
#' * ROMIs below the lower limit (0.001 with `limits = NULL` or a lower limit
#'   of 0 or below) are set to that limit in the least-squares fit of
#'   `theta` (the paper sets negative values "to a minimum that is greater
#'   than zero");
#' * the P-spline surface has cubic B-splines with knots every 5 years of age
#'   and time and second-order penalties, with both smoothing parameters
#'   chosen by BIC (the paper uses the package MortalitySmooth); when the
#'   grid starts with age 0 as its own interval, the rates of age 0 are
#'   smoothed over time only (as [mt_smooth()] leaves out the first year of
#'   life);
#' * forecast ROMIs are kept below 1 so that rates stay positive.
#'
#' The model's option to let a population converge to the mortality trend of
#' reference populations, by weights that shift from the population to the
#' references over the forecast years, is not implemented: the paper leaves
#' the weighting scheme open. Each stratum is forecast from its own data.
#'
#' Draws describe the expected (smoothed) rates; predictive draws
#' (`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
#' given the rates. Fitted rates are the smoothed rates (the observed rates
#' with `smooth = FALSE`). The Gelman-Rubin statistics of the scalar
#' parameters and the largest one of the curve parameters are reported by
#' [mt_parameters()] (`rhat`); values above 1.1 give a warning. With the
#' default settings a fit takes about 6 seconds per stratum for 71 ages and 36
#' years.
#'
#' @param core `"exponential"` or `"linear"` core model.
#' @param limits Lower and upper limit of the forecast ROMIs, or `NULL` for
#'   none.
#' @param smooth Smooth the death rates with two-dimensional P-splines before
#'   computing the ROMIs (`FALSE`: the observed rates, with 0.5 deaths in
#'   cells without deaths, for example when the rates were smoothed in a
#'   preprocessing step).
#' @param draws `"quantiles"` (equation 20 of the paper) or `"samples"`.
#' @param mcmc Settings of the sampler, a named vector: `chains`,
#'   `iterations` per chain (burn-in included), `burn_in` and `thin`.
#'   Elements left out keep their defaults.
#' @return A model specification.
#' @references
#' Bohk-Ewald C, Rau R (2017). Probabilistic mortality forecasting with
#' varying age-specific survival improvements. *Genus* 73(1): 1.
#' \doi{10.1186/s41118-016-0017-8}
#'
#' Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
#' mortality rates. *Statistical Modelling* 4(4): 279-298.
#' @examples
#' x <- mt_subset(mt_example(deaths = "single", years = 1995:2019), sex == "female")
#' # shorter chains than the defaults of the paper, for a quick example
#' spec <- bohk_rau(mcmc = c(chains = 2, iterations = 3200, burn_in = 200, thin = 3))
#' set.seed(1)
#' fit <- mt_fit(spec, x)
#' fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
#' life_expectancy(fc)
#' @export
bohk_rau <- function(core = c("exponential", "linear"), limits = c(0.005, 0.035), smooth = TRUE,
                     draws = c("quantiles", "samples"),
                     mcmc = c(chains = 5, iterations = 5200, burn_in = 200, thin = 5)) {
  core <- match.arg(core)
  draws <- match.arg(draws)
  check_flag(smooth, "smooth")
  if (!is.null(limits) && (!is.numeric(limits) || length(limits) != 2L || anyNA(limits) || limits[1] >= limits[2] ||
                           limits[2] >= 1)) {
    mt_abort("`limits` must be NULL or two increasing numbers below 1.", class = "argument")
  }
  mcmc <- bohk_rau_mcmc(mcmc)
  new_mt_spec("bohk_rau", list(core = core, limits = limits, smooth = smooth, draws = draws, mcmc = mcmc),
              min_periods = 6,
              label = paste0("Bohk-Ewald-Rau", if (core == "linear") " (linear)" else "",
                             if (draws == "samples") ", sample paths" else ""))
}

# Settings of the sampler: defaults completed by the elements given
bohk_rau_mcmc <- function(mcmc) {
  out <- c(chains = 5, iterations = 5200, burn_in = 200, thin = 5)
  if (is.list(mcmc)) mcmc <- unlist(mcmc)
  nm <- names(mcmc)
  if (!is.numeric(mcmc) || (length(mcmc) && (is.null(nm) || !all(nm %in% names(out)) || anyDuplicated(nm)))) {
    mt_abort("`mcmc` must be a named numeric vector with elements `chains`, `iterations`, `burn_in` and `thin`.",
             class = "argument", call = NULL)
  }
  out[nm] <- mcmc
  if (anyNA(out) || any(abs(out - round(out)) > 1e-8) || out[["chains"]] < 1 || out[["thin"]] < 1 ||
      out[["burn_in"]] < 0 || out[["iterations"]] - out[["burn_in"]] < out[["thin"]]) {
    mt_abort("`mcmc` needs whole numbers with at least one chain and more iterations than burn-in plus `thin`.",
             class = "argument", call = NULL)
  }
  round(out)
}

#' @export
fit_model.mt_spec_bohk_rau <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  if (nrow(D) < 3L) mt_abort("bohk_rau() needs at least three age groups.", class = "model_data")
  m <- if (p$smooth) br_smooth_rates(D, E, ages, periods) else pmax(D, 0.5) / E
  nt <- ncol(m)
  R <- 1 - m[, -1, drop = FALSE] / m[, -nt, drop = FALSE]
  tt <- seq_len(nt - 1L)
  romi_floor <- if (is.null(p$limits) || p$limits[1] <= 0) 0.001 else p$limits[1]
  smp <- if (p$core == "exponential") br_mcmc_exp(R, tt, p$mcmc, romi_floor) else br_mcmc_lin(R, tt, p$mcmc)
  if (any(smp$rhat > 1.1, na.rm = TRUE)) {
    mt_warn(sprintf(paste("The chains of the Bohk-Ewald-Rau model may not have converged (largest Gelman-Rubin",
                          "statistic %.2f); use longer chains (`mcmc`)."), max(smp$rhat, na.rm = TRUE)),
            class = "convergence", call = NULL)
  }
  dimnames(m) <- NULL
  list(core = p$core, samples = smp, t_last = nt - 1L, jump_off = log(m[, nt]), fitted = m, romi = R,
       limits = p$limits, draws = p$draws, age_labels = ages$labels)
}

# Smoothed death rates: a P-spline surface over age and time; when the grid
# starts with age 0 as its own interval, age 0 is smoothed over time only and
# the surface starts at the next age (as mt_smooth() leaves out the first
# year of life)
br_smooth_rates <- function(D, E, ages, periods) {
  if (!(ages$lower[1] == 0 && ages$upper[1] == 1 && nrow(D) > 4)) return(br_smooth(D, E, ages, periods))
  rest <- list(lower = ages$lower[-1], upper = ages$upper[-1])
  m <- matrix(NA_real_, nrow(D), ncol(D))
  m[-1, ] <- br_smooth(D[-1, , drop = FALSE], E[-1, , drop = FALSE], rest, periods)
  xp <- periods + 0.5
  bp <- pclm_basis(xp, min(xp), max(xp), kr = 5)
  f <- tryCatch(pclm_fit(D[1, ], Matrix::Diagonal(ncol(D)), E[1, ], bp$B, list(bp$P)), error = function(e) e)
  if (inherits(f, "error") || !isTRUE(f$converged) || any(!is.finite(f$gamma)) || any(f$gamma <= 0)) {
    mt_abort("The P-spline smoothing of the death rates failed; try smooth = FALSE.", class = "model_fit",
             call = NULL)
  }
  m[1, ] <- f$gamma
  m
}

# Death rates smoothed over age and time: Poisson P-spline surface with the
# exposures as offset (Currie, Durban and Eilers 2004), cubic B-splines with
# knots every `kr` years, second-order penalties, both smoothing parameters
# chosen by BIC (two rounds of searches along each parameter on a log scale,
# as pclm_fit()). The products with the two-dimensional basis are computed
# from the bases of age and time (row tensors), without forming the full
# basis, so that a surface takes about a second.
br_smooth <- function(D, E, ages, periods, kr = 5, lambda_range = c(1e-2, 1e8)) {
  xa <- cell_ages(ages$lower, ages$upper)
  xp <- periods + 0.5
  ba <- pclm_basis(xa, min(xa), max(xa), kr = kr)
  bp <- pclm_basis(xp, min(xp), max(xp), kr = kr)
  Ba <- ba$B
  Bp <- bp$B
  ca <- ncol(Ba)
  cp <- ncol(Bp)
  # row tensors: row i holds the products of all pairs of the basis values of row i
  row_tensor <- function(X) {
    X[, rep(seq_len(ncol(X)), ncol(X)), drop = FALSE] * X[, rep(seq_len(ncol(X)), each = ncol(X)), drop = FALSE]
  }
  Ga <- row_tensor(Ba)
  Gp <- row_tensor(Bp)
  # B' W B for B = Bp %x% Ba and weights W [age, period]
  btwb <- function(W) {
    A <- array(crossprod(Ga, W %*% Gp), c(ca, ca, cp, cp))
    matrix(aperm(A, c(1, 3, 2, 4)), ca * cp)
  }
  bt <- function(Z) as.vector(crossprod(Ba, Z) %*% Bp)
  linpred <- function(th) Ba %*% matrix(th, ca) %*% t(Bp)
  P1 <- kronecker(diag(cp), ba$P)
  P2 <- kronecker(bp$P, diag(ca))
  n <- length(D)
  th0 <- solve(btwb(matrix(1, nrow(D), ncol(D))) + diag(1e-6, ca * cp), bt(log(pmax(D, 0.5) / E)))
  irls <- function(lambda, th) {
    P <- lambda[1] * P1 + lambda[2] * P2
    pdev <- function(th) {
      eta <- linpred(th)
      if (any(!is.finite(eta)) || max(eta) > 50) return(Inf)
      poisson_dev(D, E * exp(eta)) + sum(th * (P %*% th))
    }
    obj <- pdev(th)
    conv <- FALSE
    for (it in 1:100) {
      eta <- linpred(th)
      mu <- E * exp(eta)
      W <- btwb(mu)
      A <- W + P
      th_new <- tryCatch(solve(A, bt(mu * eta + D - mu)), error = function(e) NULL)
      if (is.null(th_new) || anyNA(th_new)) break
      step <- th_new - th
      obj_new <- pdev(th + step)
      halvings <- 0L
      while (!(obj_new <= obj + 1e-10 * abs(obj)) && halvings < 25L) {
        step <- step / 2
        obj_new <- pdev(th + step)
        halvings <- halvings + 1L
      }
      if (!is.finite(obj_new)) break
      th <- th + step
      done <- max(abs(step)) < 1e-6 * (1 + max(abs(th))) || abs(obj - obj_new) < 1e-8 * (1 + abs(obj_new))
      obj <- obj_new
      if (done) {
        conv <- TRUE
        break
      }
    }
    mu <- E * exp(linpred(th))
    W <- btwb(mu)
    ed <- tryCatch(sum(diag(solve(W + P, W))), error = function(e) NA_real_)
    dev <- poisson_dev(D, mu)
    list(th = th, rate = mu / E, conv = conv, bic = dev + log(n) * ed)
  }
  state <- new.env()
  state$th <- th0
  obj <- function(ll) {
    f <- irls(10^ll, state$th)
    if (f$conv) state$th <- f$th
    if (is.finite(f$bic)) f$bic else 1e300
  }
  lr <- log10(lambda_range)
  ll <- rep(mean(lr), 2)
  for (round in 1:2) {
    for (k in 1:2) {
      ll[k] <- stats::optimize(function(z) obj(replace(ll, k, z)), lr, tol = 0.1)$minimum
    }
  }
  f <- irls(10^ll, state$th)
  if (!f$conv || any(!is.finite(f$rate)) || any(f$rate <= 0)) {
    mt_abort("The P-spline smoothing of the death rates failed; try smooth = FALSE.", class = "model_fit",
             call = NULL)
  }
  out <- matrix(f$rate, nrow(D))
  attr(out, "lambda") <- 10^ll
  out
}

# Gelman-Rubin potential scale reduction factor of samples [iteration, chain]
br_psrf <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (ncol(x) < 2L || n < 2L) return(NA_real_)
  w <- mean(apply(x, 2, stats::var))
  b <- n * stats::var(colMeans(x))
  if (!is.finite(w) || w <= 0) return(NA_real_)
  sqrt(((n - 1) / n * w + b / n) / w)
}

# Largest PSRF of the rows of samples [parameter, iteration x chain]
br_psrf_max <- function(x, chains) {
  n <- ncol(x) / chains
  max(vapply(seq_len(nrow(x)), function(i) br_psrf(matrix(x[i, ], n, chains)), numeric(1)), na.rm = TRUE)
}

# Variance sigma^2 given the sum of squares `ss` of `n` normal values, with
# sigma uniform on (0, 1): inverse gamma with shape (n - 1) / 2 and rate
# ss / 2, truncated below 1 (inverse distribution function of the precision)
br_variance <- function(n, ss) {
  a <- (n - 1) / 2
  b <- max(ss, 1e-300) / 2
  q1 <- stats::pgamma(1, a, rate = b, lower.tail = FALSE)
  if (!is.finite(q1) || q1 <= 0) return(1)
  1 / stats::qgamma(stats::runif(1) * q1, a, rate = b, lower.tail = FALSE)
}

# Normal truncated to (lo, hi), by the inverse distribution function
br_trunc_norm <- function(m, s, lo, hi) {
  pl <- stats::pnorm(lo, m, s)
  ph <- stats::pnorm(hi, m, s)
  if (!(ph > pl)) return(if (m < lo) lo else hi)
  min(max(stats::qnorm(stats::runif(1, pl, ph), m, s), lo), hi)
}

# Slice sampling of a univariate density on a bounded interval (Neal 2003),
# shrinking the interval from its bounds
br_slice <- function(f, x0, lo, hi) {
  y <- f(x0) - stats::rexp(1)
  for (i in 1:200) {
    x1 <- stats::runif(1, lo, hi)
    if (f(x1) > y) return(x1)
    if (x1 < x0) lo <- x1 else hi <- x1
  }
  x0
}

# Indexes of the kept iterations of a chain
br_keep <- function(mcmc) {
  seq(mcmc[["burn_in"]] + mcmc[["thin"]], mcmc[["iterations"]], by = mcmc[["thin"]])
}

# Exponential core model: rho[x, t] ~ N(exp(b1[x] + b2[x] log t), sigma^2),
# b[x] ~ N(theta[x], diag(sigma1^2, sigma2^2)); sigma, sigma1, sigma2 ~ U(0, 1)
br_mcmc_exp <- function(R, tt, mcmc, romi_floor) {
  na <- nrow(R)
  nt <- ncol(R)
  lt <- log(tt)
  X <- cbind(1, lt)
  theta <- log(pmax(R, romi_floor)) %*% X %*% solve(crossprod(X))
  th1 <- theta[, 1]
  th2 <- theta[, 2]
  # Gauss-Newton information of the curve at theta, for the proposals
  f0 <- exp(th1 + outer(th2, lt))
  j11 <- rowSums(f0^2)
  j12 <- rowSums(f0^2 * rep(lt, each = na))
  j22 <- rowSums(f0^2 * rep(lt^2, each = na))
  ssq <- function(b1, b2) rowSums((R - exp(b1 + outer(b2, lt)))^2)
  keep <- br_keep(mcmc)
  nk <- length(keep)
  nc <- mcmc[["chains"]]
  out <- list(b1 = matrix(0, na, nk * nc), b2 = matrix(0, na, nk * nc), sigma = numeric(nk * nc),
              sigma1 = numeric(nk * nc), sigma2 = numeric(nk * nc))
  s_ols <- sqrt(mean(ssq(th1, th2)))
  for (ch in seq_len(nc)) {
    s <- min(max(s_ols * exp(stats::rnorm(1, 0, 0.3)), 1e-6), 0.99)
    s1 <- stats::runif(1, 0.02, 0.5)
    s2 <- stats::runif(1, 0.02, 0.5)
    b1 <- th1 + stats::rnorm(na, 0, 0.01)
    b2 <- th2 + stats::rnorm(na, 0, 0.01)
    ss <- ssq(b1, b2)
    sc <- rep(2.38^2 / 2, na)
    acc <- numeric(na)
    k <- 0L
    for (it in seq_len(mcmc[["iterations"]])) {
      # curve parameters of each age: random walk Metropolis, proposal
      # covariance sc * H^-1 with H the curvature of the log posterior
      h11 <- j11 / s^2 + 1 / s1^2
      h12 <- j12 / s^2
      h22 <- j22 / s^2 + 1 / s2^2
      det_h <- h11 * h22 - h12^2
      l11 <- sqrt(h22 / det_h)
      l21 <- -h12 / det_h / l11
      l22 <- sqrt(pmax(h11 / det_h - l21^2, 0))
      z1 <- stats::rnorm(na)
      z2 <- stats::rnorm(na)
      p1 <- b1 + sqrt(sc) * l11 * z1
      p2 <- b2 + sqrt(sc) * (l21 * z1 + l22 * z2)
      ss_p <- ssq(p1, p2)
      lp_c <- -ss / (2 * s^2) - (b1 - th1)^2 / (2 * s1^2) - (b2 - th2)^2 / (2 * s2^2)
      lp_p <- -ss_p / (2 * s^2) - (p1 - th1)^2 / (2 * s1^2) - (p2 - th2)^2 / (2 * s2^2)
      ok <- is.finite(lp_p) & log(stats::runif(na)) < lp_p - lp_c
      b1[ok] <- p1[ok]
      b2[ok] <- p2[ok]
      ss[ok] <- ss_p[ok]
      acc <- acc + ok
      s <- sqrt(br_variance(na * nt, sum(ss)))
      s1 <- sqrt(br_variance(na, sum((b1 - th1)^2)))
      s2 <- sqrt(br_variance(na, sum((b2 - th2)^2)))
      # scale moves: sigma1 (sigma2) and the deviations of b1 (b2) from theta
      # multiplied by the same factor, for ages whose curve the data hardly
      # determine (the acceptance ratio is the likelihood ratio times the
      # factor)
      for (v in 1:2) {
        cf <- exp(stats::rnorm(1, 0, 0.3))
        if ((if (v == 1L) s1 else s2) * cf >= 1) next
        q1 <- if (v == 1L) th1 + (b1 - th1) * cf else b1
        q2 <- if (v == 2L) th2 + (b2 - th2) * cf else b2
        ss_q <- ssq(q1, q2)
        lr <- -(sum(ss_q) - sum(ss)) / (2 * s^2) + log(cf)
        if (is.finite(lr) && log(stats::runif(1)) < lr) {
          b1 <- q1
          b2 <- q2
          ss <- ss_q
          if (v == 1L) s1 <- s1 * cf else s2 <- s2 * cf
        }
      }
      # adaptation of the proposal scale during the burn-in
      if (it <= mcmc[["burn_in"]] && it %% 50 == 0) {
        rate <- acc / 50
        sc <- sc * ifelse(rate > 0.45, 2, ifelse(rate < 0.15, 0.5, 1))
        acc[] <- 0
      }
      if (k < nk && it == keep[k + 1L]) {
        k <- k + 1L
        j <- (ch - 1L) * nk + k
        out$b1[, j] <- b1
        out$b2[, j] <- b2
        out$sigma[j] <- s
        out$sigma1[j] <- s1
        out$sigma2[j] <- s2
      }
    }
  }
  out$theta <- theta
  out$rhat <- c(sigma = br_psrf(matrix(out$sigma, nk)), sigma1 = br_psrf(matrix(out$sigma1, nk)),
                sigma2 = br_psrf(matrix(out$sigma2, nk)),
                curve_max = max(br_psrf_max(out$b1, nc), br_psrf_max(out$b2, nc)))
  out
}

# Linear core model: rho[x, t] ~ N(b1[x] + b2[x] t, sigma^2), b[x] ~ N(mu,
# Omega) with Omega from omega1, omega2 and r; sigma, omega ~ U(0, 1),
# mu ~ U(-0.1, 0.1), r ~ U(-1, 1)
br_mcmc_lin <- function(R, tt, mcmc) {
  na <- nrow(R)
  nt <- ncol(R)
  X <- cbind(1, tt)
  XtX <- crossprod(X)
  XtR <- t(R %*% X)
  b_ols <- solve(XtX, XtR)
  keep <- br_keep(mcmc)
  nk <- length(keep)
  nc <- mcmc[["chains"]]
  out <- list(b1 = matrix(0, na, nk * nc), b2 = matrix(0, na, nk * nc), sigma = numeric(nk * nc),
              mu1 = numeric(nk * nc), mu2 = numeric(nk * nc), omega1 = numeric(nk * nc),
              omega2 = numeric(nk * nc), r = numeric(nk * nc))
  ll_omega <- function(w1, w2, r, S) {
    -na / 2 * log(w1^2 * w2^2 * (1 - r^2)) -
      0.5 * (S[1, 1] / w1^2 - 2 * r * S[1, 2] / (w1 * w2) + S[2, 2] / w2^2) / (1 - r^2)
  }
  res0 <- R - b_ols[1, ] - outer(b_ols[2, ], tt)
  s_ols <- sqrt(mean(res0^2))
  sd_b <- pmin(pmax(apply(b_ols, 1, stats::sd), 1e-6), 0.9)
  for (ch in seq_len(nc)) {
    s <- min(max(s_ols * exp(stats::rnorm(1, 0, 0.3)), 1e-6), 0.99)
    mu <- pmin(pmax(rowMeans(b_ols) + stats::rnorm(2, 0, sd_b / 4), -0.099), 0.099)
    w1 <- min(sd_b[1] * exp(stats::rnorm(1, 0, 0.3)), 0.99)
    w2 <- min(sd_b[2] * exp(stats::rnorm(1, 0, 0.3)), 0.99)
    r <- stats::runif(1, -0.5, 0.5)
    k <- 0L
    for (it in seq_len(mcmc[["iterations"]])) {
      Om <- matrix(c(w1^2, r * w1 * w2, r * w1 * w2, w2^2), 2)
      Oi <- solve(Om)
      U <- chol(XtX / s^2 + Oi)
      M <- backsolve(U, forwardsolve(t(U), XtR / s^2 + as.vector(Oi %*% mu)))
      B <- M + backsolve(U, matrix(stats::rnorm(2 * na), 2))
      ss <- sum((R - B[1, ] - outer(B[2, ], tt))^2)
      s <- sqrt(br_variance(na * nt, ss))
      bbar <- rowMeans(B)
      mu[1] <- br_trunc_norm(bbar[1] + Om[1, 2] / Om[2, 2] * (mu[2] - bbar[2]),
                             sqrt((Om[1, 1] - Om[1, 2]^2 / Om[2, 2]) / na), -0.1, 0.1)
      mu[2] <- br_trunc_norm(bbar[2] + Om[1, 2] / Om[1, 1] * (mu[1] - bbar[1]),
                             sqrt((Om[2, 2] - Om[1, 2]^2 / Om[1, 1]) / na), -0.1, 0.1)
      S <- tcrossprod(B - mu)
      w1 <- br_slice(function(v) ll_omega(v, w2, r, S), w1, 0, 1)
      w2 <- br_slice(function(v) ll_omega(w1, v, r, S), w2, 0, 1)
      r <- br_slice(function(v) ll_omega(w1, w2, v, S), r, -1, 1)
      if (k < nk && it == keep[k + 1L]) {
        k <- k + 1L
        j <- (ch - 1L) * nk + k
        out$b1[, j] <- B[1, ]
        out$b2[, j] <- B[2, ]
        out$sigma[j] <- s
        out$mu1[j] <- mu[1]
        out$mu2[j] <- mu[2]
        out$omega1[j] <- w1
        out$omega2[j] <- w2
        out$r[j] <- r
      }
    }
  }
  sc <- c("sigma", "mu1", "mu2", "omega1", "omega2", "r")
  out$rhat <- c(vapply(sc, function(v) br_psrf(matrix(out[[v]], nk)), numeric(1)),
                curve_max = max(br_psrf_max(out$b1, nc), br_psrf_max(out$b2, nc)))
  out
}

# Mean ROMI curves [age, sample] at time t
br_curve <- function(fit, t, idx = NULL) {
  smp <- fit$samples
  b1 <- if (is.null(idx)) smp$b1 else smp$b1[, idx, drop = FALSE]
  b2 <- if (is.null(idx)) smp$b2 else smp$b2[, idx, drop = FALSE]
  if (fit$core == "exponential") exp(b1 + b2 * log(t)) else b1 + b2 * t
}

# Forecast ROMIs: limits, and below 1 so that rates stay positive
br_limit <- function(rho, limits) {
  if (!is.null(limits)) rho <- pmin(pmax(rho, limits[1]), limits[2])
  pmin(rho, 1 - 1e-6)
}

#' @export
simulate_model.mt_spec_bohk_rau <- function(spec, fit, h, n_draws, ...) {
  smp <- fit$samples
  ns <- length(smp$sigma)
  na <- length(fit$jump_off)
  logm <- matrix(fit$jump_off, na, n_draws)
  out <- array(NA_real_, c(na, h, n_draws))
  if (identical(fit$draws, "samples")) {
    idx <- if (n_draws <= ns) sample.int(ns, n_draws) else sample.int(ns, n_draws, replace = TRUE)
    for (j in seq_len(h)) {
      mu <- br_curve(fit, fit$t_last + j, idx)
      rho <- br_limit(mu + rep(smp$sigma[idx], each = na) * stats::rnorm(na * n_draws), fit$limits)
      logm <- logm + log1p(-rho)
      out[, j, ] <- exp(logm)
    }
    return(out)
  }
  # quantiles of the forecast ROMIs at every age and year (equation 20)
  pos <- 1 + (ns - 1) * (seq_len(n_draws) - 0.5) / n_draws
  lo <- floor(pos)
  hi <- pmin(lo + 1L, ns)
  w <- rep(pos - lo, each = na)
  for (j in seq_len(h)) {
    mu <- br_curve(fit, fit$t_last + j)
    rho <- br_limit(mu + rep(smp$sigma, each = na) * stats::rnorm(na * ns), fit$limits)
    srt <- t(apply(matrix(rho, na), 1, sort))
    q <- srt[, lo, drop = FALSE] + w * (srt[, hi, drop = FALSE] - srt[, lo, drop = FALSE])
    logm <- logm + log1p(-q)
    out[, j, ] <- exp(logm)
  }
  out[, , sample.int(n_draws), drop = FALSE]
}
