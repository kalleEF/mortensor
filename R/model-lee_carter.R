#' Lee-Carter model
#'
#' `log m[x, t] = a[x] + b[x] k[t]` with the constraints `sum(b) = 1` and
#' `sum(k) = 0`, fitted per stratum:
#'
#' * `method = "poisson"`: Poisson log-bilinear model (Brouhns, Denuit and
#'   Vermunt 2002), maximum likelihood by Newton updates;
#' * `method = "svd"`: original method (Lee and Carter 1992): singular value
#'   decomposition of centred log rates, with `k[t]` re-estimated so that the
#'   fitted deaths match the observed total deaths of each year.
#'
#' The period index follows a random walk with drift. Draws include the
#' innovations and the uncertainty of the drift. With `bootstrap = TRUE` the
#' parameters are re-estimated on `n_boot` Poisson resamples of the fitted
#' deaths (semi-parametric bootstrap of Brouhns et al. 2005), which adds the
#' uncertainty of `a` and `b`.
#'
#' @param method `"poisson"` or `"svd"`.
#' @param jump_off `"fitted"` (forecasts start from the fitted rates of the
#'   last year) or `"observed"` (from the observed rates of the last year).
#' @param bootstrap Semi-parametric bootstrap of the parameters.
#' @param n_boot Number of bootstrap fits.
#' @return A model specification.
#' @references
#' Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
#' *Journal of the American Statistical Association* 87(419): 659-671.
#'
#' Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear regression
#' approach to the construction of projected lifetables. *Insurance:
#' Mathematics and Economics* 31(3): 373-393.
#' \doi{10.1016/S0167-6687(02)00185-3}
#'
#' Brouhns N, Denuit M, Van Keilegom I (2005). Bootstrapping the Poisson
#' log-bilinear model for mortality forecasting. *Scandinavian Actuarial
#' Journal* 2005(3): 212-224.
#' @examples
#' x <- mt_example(deaths = "single")
#' fit <- mt_fit(lee_carter(), x)
#' fc <- mt_forecast(fit, h = 10, n_draws = 200, seed = 1)
#' @export
lee_carter <- function(method = c("poisson", "svd"), jump_off = c("fitted", "observed"), bootstrap = FALSE,
                       n_boot = 100) {
  method <- match.arg(method)
  jump_off <- match.arg(jump_off)
  new_mt_spec("lee_carter", list(method = method, jump_off = jump_off, bootstrap = bootstrap, n_boot = n_boot),
              min_periods = 5, label = "Lee-Carter")
}

#' @export
fit_model.mt_spec_lee_carter <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  est <- if (p$method == "poisson") lc_poisson(D, E) else lc_svd(D, E)
  est$last_log_rate <- log_rates(D[, ncol(D)], E[, ncol(E)])
  est$fitted <- exp(est$a + outer(est$b, est$k))
  est$jump_off <- p$jump_off
  if (isTRUE(p$bootstrap)) {
    fitted <- E * est$fitted
    est$boot <- lapply(seq_len(p$n_boot), function(i) {
      Db <- matrix(stats::rpois(length(fitted), fitted), nrow(fitted))
      tryCatch(if (p$method == "poisson") lc_poisson(Db, E, init = est) else lc_svd(Db, E), error = function(e) NULL)
    })
    est$boot <- Filter(Negate(is.null), est$boot)
  }
  est
}

#' @export
simulate_model.mt_spec_lee_carter <- function(spec, fit, h, n_draws, ...) {
  sets <- if (length(fit$boot)) fit$boot else list(fit)
  which_set <- rep_len(seq_along(sets), n_draws)
  out <- array(NA_real_, c(length(fit$a), h, n_draws))
  for (i in seq_along(sets)) {
    js <- which(which_set == i)
    if (!length(js)) next
    s <- sets[[i]]
    kp <- rw_paths(s$k[length(s$k)], s$drift, matrix(s$dk - s$drift, nrow = 1), h, length(js), "gaussian",
                   sigma = s$sigma_k, exp_out = FALSE)
    for (m in seq_along(js)) {
      kk <- kp[1, , m]
      base <- if (fit$jump_off == "observed") fit$last_log_rate - s$b * s$k[length(s$k)] else s$a
      out[, , js[m]] <- exp(base + outer(s$b, kk))
    }
  }
  out
}

lc_finish <- function(a, b, k) {
  s <- sum(b)
  kbar <- mean(k)
  a <- a + b * kbar
  k <- (k - kbar) * s
  b <- b / s
  dk <- diff(k)
  drift <- mean(dk)
  list(a = a, b = b, k = k, dk = dk, drift = drift, sigma_k = sqrt(sum((dk - drift)^2) / max(length(dk) - 1, 1)))
}

lc_svd <- function(D, E) {
  y <- log_rates(D, E)
  a <- rowMeans(y)
  s <- svd(y - a, nu = 1, nv = 1)
  b <- s$u[, 1]
  k <- s$d[1] * s$v[, 1]
  if (sum(b) < 0) {
    b <- -b
    k <- -k
  }
  # second stage: refit k so that fitted total deaths equal observed deaths
  for (t in seq_len(ncol(D))) {
    f <- function(kt) sum(E[, t] * exp(a + b * kt)) - sum(D[, t])
    rng <- k[t] + c(-1, 1) * (10 + 10 * abs(k[t]))
    k[t] <- tryCatch(stats::uniroot(f, rng, extendInt = "yes", tol = 1e-12)$root, error = function(e) k[t])
  }
  out <- lc_finish(a, b, k)
  out$method <- "svd"
  out
}

lc_poisson <- function(D, E, init = NULL, maxit = 500, tol = 1e-8) {
  if (is.null(init)) init <- lc_svd(D, E)
  a <- init$a
  b <- init$b
  k <- init$k
  ll_old <- -Inf
  for (it in seq_len(maxit)) {
    Dh <- E * exp(a + outer(b, k))
    a <- a + rowSums(D - Dh) / rowSums(Dh)
    Dh <- E * exp(a + outer(b, k))
    k <- k + colSums((D - Dh) * b) / colSums(Dh * b^2)
    a <- a + b * mean(k)
    k <- k - mean(k)
    Dh <- E * exp(a + outer(b, k))
    kk <- matrix(k, nrow(D), ncol(D), byrow = TRUE)
    b <- b + rowSums((D - Dh) * kk) / rowSums(Dh * kk^2)
    Dh <- E * exp(a + outer(b, k))
    ll <- sum(D * log(Dh) - Dh)
    if (!is.finite(ll)) mt_abort("Lee-Carter fit diverged.", class = "model_fit")
    if (abs(ll - ll_old) < tol * abs(ll)) break
    ll_old <- ll
  }
  out <- lc_finish(a, b, k)
  out$method <- "poisson"
  out$loglik <- ll
  out$iterations <- it
  out
}

#' Parameters of a fitted Lee-Carter model
#'
#' @param fit An `mt_fit` object from [lee_carter()].
#' @return A `data.table` with `stratum`, `parameter` (`a`, `b` or `k`),
#'   `index` (age label or year) and `value`.
#' @examples
#' x <- mt_example(deaths = "single")
#' head(lc_parameters(mt_fit(lee_carter(), x)))
#' @export
lc_parameters <- function(fit) {
  check_class(fit, "mt_fit", "fit")
  if (!inherits(fit$spec, "mt_spec_lee_carter")) mt_abort("Not a Lee-Carter fit.", class = "argument")
  rbindlist(lapply(seq_along(fit$fits), function(s) {
    f <- fit$fits[[s]]
    if (is.null(f)) return(NULL)
    lab <- fit$strata$label[s]
    rbind(data.table(stratum = lab, parameter = "a", index = fit$age$labels, value = f$a),
          data.table(stratum = lab, parameter = "b", index = fit$age$labels, value = f$b),
          data.table(stratum = lab, parameter = "k", index = as.character(fit$periods), value = f$k))
  }))
}
