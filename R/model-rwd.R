#' Random walk with drift on log death rates
#'
#' Benchmark model: the log death rate of every age follows a random walk
#' with drift, `y[t] = y[t-1] + mu + e[t]`. The drift is estimated from the
#' last `drift_window` years. Draws include drift uncertainty and
#' innovations:
#'
#' * `"bootstrap"` (default): innovations are age vectors resampled from the
#'   observed year-to-year changes, which keeps the correlation between ages;
#'   the drift of each path is the mean of a bootstrap sample of those
#'   changes;
#' * `"gaussian"`: independent normal innovations per age with the
#'   estimated variance, and normal drift uncertainty.
#'
#' Zero death counts are replaced by 0.5 before taking logs.
#'
#' With `observation_noise = TRUE` the log rate is a random walk with drift
#' observed with noise (a local level model with drift). Year-to-year
#' changes of observed rates then contain the noise of two years, and a
#' plain random walk accumulates it at every step, which makes its
#' intervals for observed rates too wide when deaths are few or rates
#' fluctuate. The variances are estimated per age from the variance and
#' first autocovariance of the changes (the noise variance is the Poisson
#' part, one over the deaths, plus an extra part common to all ages), the
#' level is Kalman filtered, and the forecast draws are paths of the level;
#' predictive draws add the observation noise. The variant is more accurate
#' than the plain random walk where deaths are few, but over 5 to 10 years
#' its intervals can be too narrow, because changes of trend exceed the
#' uncertainty of a constant drift; the plain random walk is wider. With a
#' `drift_window`, the model is fitted to the years of the window and fitted
#' rates of earlier years are missing.
#'
#' @param drift_window Number of most recent years used for drift and
#'   innovations (`NULL`: all years).
#' @param innovations `"bootstrap"` or `"gaussian"`.
#' @param observation_noise Model the log rates as a random walk plus noise.
#' @return A model specification.
#' @examples
#' x <- mt_example(deaths = "single")
#' fc <- mt_forecast(mt_fit(rwd(drift_window = 15), x), h = 5, n_draws = 100, seed = 1)
#' @export
rwd <- function(drift_window = NULL, innovations = c("bootstrap", "gaussian"), observation_noise = FALSE) {
  innovations <- match.arg(innovations)
  if (!is.null(drift_window)) check_number(drift_window, "drift_window", lower = 3, integerish = TRUE)
  check_flag(observation_noise, "observation_noise")
  new_mt_spec("rwd", list(drift_window = drift_window, innovations = innovations,
                          observation_noise = observation_noise),
              min_periods = if (observation_noise) 6 else 3,
              label = if (observation_noise) "Random walk with drift plus noise" else "Random walk with drift",
              smooth = observation_noise)
}

log_rates <- function(D, E) log(pmax(D, 0.5) / E)

#' @export
fit_model.mt_spec_rwd <- function(spec, D, E, ages, periods, ...) {
  y <- log_rates(D, E)
  nt <- ncol(y)
  w <- min(spec$params$drift_window %||% nt, nt)
  y <- y[, (nt - w + 1):nt, drop = FALSE]
  dy <- y[, -1, drop = FALSE] - y[, -w, drop = FALSE]
  mu <- rowMeans(dy)
  out <- list(last = y[, w], mu = mu, resid = dy - mu, sigma = sqrt(rowSums((dy - mu)^2) / max(w - 2, 1)),
              innovations = spec$params$innovations)
  if (isTRUE(spec$params$observation_noise)) {
    Dw <- D[, (nt - w + 1):nt, drop = FALSE]
    out <- c(out, local_level_fit(y, dy, mu, rowMeans(1 / pmax(Dw, 0.5)), weights = rowSums(Dw)))
    # fitted rates of the years before the window are missing
    out$fitted <- cbind(matrix(NA_real_, nrow(y), nt - w), out$fitted)
  }
  out
}

#' @export
simulate_model.mt_spec_rwd <- function(spec, fit, h, n_draws, ...) {
  if (!is.null(fit$level)) return(local_level_paths(fit, h, n_draws))
  rw_paths(fit$last, fit$mu, fit$resid, h, n_draws, fit$innovations, sigma = fit$sigma)
}

# Random walk with drift observed with noise (local level model with drift),
# per age: y[t] = L[t] + e[t], L[t] = L[t-1] + mu + eta[t]. Moments of the
# differences: Var(dy) = s2_eta + 2 s2_e, Cov(dy[t], dy[t-1]) = -s2_e. The
# noise variance is the Poisson part (1 / deaths) plus a common extra part
# (median over ages of -Cov - Poisson part); the level is Kalman filtered.
local_level_fit <- function(y, dy, mu, pois, weights = NULL) {
  na <- nrow(y)
  w <- ncol(y)
  r <- dy - mu
  g0 <- rowSums(r^2) / max(ncol(r) - 1, 1)
  g1 <- if (ncol(r) > 2) rowSums(r[, -1, drop = FALSE] * r[, -ncol(r), drop = FALSE]) / (ncol(r) - 1) else -pois
  extra <- max(0, stats::median(-g1 - pois))
  s2e <- pmin(pois + extra, 0.45 * g0)
  s2eta <- pmax(g0 - 2 * s2e, 0.1 * g0)
  # partial pooling over ages: where deaths are few, the moment estimate is
  # dominated by sampling noise (variance about 2 g0^2 / n) and is shrunk
  # towards the median of the ages, weighted by their deaths; the prior is
  # weak (standard deviation three times the median), so well observed ages
  # keep their own estimate
  if (!is.null(weights) && na > 2) {
    pooled <- weighted_median(s2eta, weights)
    v_age <- 2 * g0^2 / max(ncol(r) - 1, 1)
    v_pool <- (3 * pooled)^2
    s2eta <- (s2eta / v_age + pooled / v_pool) / (1 / v_age + 1 / v_pool)
  }
  L <- matrix(NA_real_, na, w)
  P <- matrix(NA_real_, na, w)
  Lp <- L
  Pp <- P
  L[, 1] <- y[, 1]
  P[, 1] <- s2e
  for (t in 2:w) {
    Lp[, t] <- L[, t - 1] + mu
    Pp[, t] <- P[, t - 1] + s2eta
    K <- Pp[, t] / (Pp[, t] + s2e)
    L[, t] <- Lp[, t] + K * (y[, t] - Lp[, t])
    P[, t] <- (1 - K) * Pp[, t]
  }
  S <- L
  for (t in (w - 1):1) {
    J <- P[, t] / Pp[, t + 1]
    S[, t] <- L[, t] + J * (S[, t + 1] - Lp[, t + 1])
  }
  list(level = L[, w], level_var = P[, w], s2eta = s2eta, s2e = s2e,
       mu_var = (s2eta * (w - 1) + 2 * s2e) / (w - 1)^2, scale = sqrt(s2eta / pmax(g0, 1e-12)), fitted = exp(S))
}

local_level_paths <- function(fit, h, n_draws) {
  na <- length(fit$level)
  nr <- ncol(fit$resid)
  out <- array(NA_real_, c(na, h, n_draws))
  for (j in seq_len(n_draws)) {
    start <- fit$level + stats::rnorm(na, 0, sqrt(fit$level_var))
    drift <- fit$mu + stats::rnorm(na, 0, sqrt(fit$mu_var))
    eta <- if (fit$innovations == "bootstrap") {
      fit$resid[, sample.int(nr, h, replace = TRUE), drop = FALSE] * fit$scale
    } else {
      matrix(stats::rnorm(na * h, 0, sqrt(fit$s2eta)), na, h)
    }
    path <- start + t(apply(eta + drift, 1, cumsum))
    if (h == 1) path <- matrix(path, na, 1)
    out[, , j] <- path
  }
  exp(out)
}

weighted_median <- function(x, w) {
  o <- order(x)
  cw <- cumsum(w[o]) / sum(w)
  x[o][which(cw >= 0.5)[1]]
}

# Multivariate random walk with drift: returns exp(paths) as [dim, h, draw]
rw_paths <- function(last, mu, resid, h, n_draws, innovations = "bootstrap", sigma = NULL, exp_out = TRUE) {
  na <- length(last)
  nr <- ncol(resid)
  out <- array(NA_real_, c(na, h, n_draws))
  for (j in seq_len(n_draws)) {
    if (innovations == "bootstrap") {
      drift <- mu + rowMeans(resid[, sample.int(nr, nr, replace = TRUE), drop = FALSE])
      eps <- resid[, sample.int(nr, h, replace = TRUE), drop = FALSE]
    } else {
      drift <- mu + stats::rnorm(na, 0, sigma / sqrt(nr))
      eps <- matrix(stats::rnorm(na * h, 0, sigma), na, h)
    }
    path <- last + t(apply(eps + drift, 1, cumsum))
    if (h == 1) path <- matrix(path, na, 1)
    out[, , j] <- path
  }
  if (exp_out) exp(out) else out
}
