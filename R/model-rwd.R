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
#' Zero death counts are replaced by 0.5 before taking logs. The draws
#' describe observed rates, so predictive draws add nothing.
#'
#' A random walk with drift observed with noise (the local level model with
#' drift, Harvey 1989) has ARIMA(0,1,1) with drift as its reduced form for
#' the observed log rates: `arima_series(order = c(0, 1, 1), drift = TRUE)`
#' fits it by maximum likelihood.
#'
#' @param drift_window Number of most recent years used for drift and
#'   innovations (`NULL`: all years).
#' @param innovations `"bootstrap"` or `"gaussian"`.
#' @return A model specification.
#' @references
#' Harvey AC (1989). *Forecasting, Structural Time Series Models and the
#' Kalman Filter*. Cambridge University Press.
#' @seealso [arima_series()].
#' @examples
#' x <- mt_example(deaths = "single")
#' fc <- mt_forecast(mt_fit(rwd(drift_window = 15), x), h = 5, n_draws = 100, seed = 1)
#' @export
rwd <- function(drift_window = NULL, innovations = c("bootstrap", "gaussian")) {
  innovations <- match.arg(innovations)
  if (!is.null(drift_window)) check_number(drift_window, "drift_window", lower = 3, integerish = TRUE)
  new_mt_spec("rwd", list(drift_window = drift_window, innovations = innovations),
              min_periods = 3, label = "Random walk with drift", smooth = FALSE)
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
  list(last = y[, w], mu = mu, resid = dy - mu, sigma = sqrt(rowSums((dy - mu)^2) / max(w - 2, 1)),
       innovations = spec$params$innovations)
}

#' @export
simulate_model.mt_spec_rwd <- function(spec, fit, h, n_draws, ...) {
  rw_paths(fit$last, fit$mu, fit$resid, h, n_draws, fit$innovations, sigma = fit$sigma)
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
