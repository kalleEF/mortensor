#' Lee-Carter model
#'
#' `log m[x, t] = a[x] + b[x] k[t]` with the constraints `sum(b) = 1` and
#' `sum(k) = 0`, fitted per stratum:
#'
#' * `method = "poisson"`: Poisson log-bilinear model (Brouhns, Denuit and
#'   Vermunt 2002), maximum likelihood by Newton updates;
#' * `method = "svd"`: original method (Lee and Carter 1992): singular value
#'   decomposition of centred log rates, with `k[t]` re-estimated in a
#'   second stage (`adjust`, as in `demography::lca()`): so that the fitted
#'   deaths match the observed total deaths of each year (`"dt"`, Lee and
#'   Carter 1992), by Poisson maximum likelihood of the deaths by age given
#'   `a` and `b` (`"dxt"`, Booth, Maindonald and Smith 2002), so that the
#'   fitted life expectancy at the first age matches the observed one
#'   (`"e0"`, Lee and Miller 2001), or not at all (`"none"`).
#'
#' `choose_period = TRUE` fits the model to the last years only, starting in
#' the year that minimises the ratio of the mean deviance of a fit with a
#' linear `k[t]` to the mean deviance of the Lee-Carter fit (Booth,
#' Maindonald and Smith 2002): the longest recent period over which the
#' decline was close to linear (`demography::lca(chooseperiod = TRUE,
#' breakmethod = "bms")` gives the same period). Fitted rates of earlier
#' years are missing.
#'
#' Compared with `demography::lca()`, `k[t]` is centred again after the
#' second stage (a constant shift of `k` and `a`, with identical fitted
#' rates), and life tables use `a_x = n/2` at age 0 (as [life_table()]), so
#' `adjust = "e0"` matches a slightly different life expectancy than
#' demography's life table where age 0 is in the data.
#'
#' The period index follows a random walk with drift. Draws include the
#' innovations and the uncertainty of the drift. With `bootstrap = TRUE` the
#' parameters are re-estimated on `n_boot` Poisson resamples of the fitted
#' deaths of the fitting period (semi-parametric bootstrap of Brouhns et al.
#' 2005), which adds the uncertainty of `a` and `b`.
#'
#' Age groups without deaths in the fitting period get a constant rate (0.5
#' deaths over the period) and no trend; with `adjust = "e0"` they enter the
#' life tables at that rate.
#'
#' @param method `"poisson"` or `"svd"`.
#' @param adjust Re-estimation of `k[t]` for `method = "svd"`: `"dt"`,
#'   `"dxt"`, `"e0"` or `"none"`.
#' @param choose_period Choose the fitting period (Booth, Maindonald and
#'   Smith 2002).
#' @param min_years `choose_period` considers fitting periods of more than
#'   `min_years` years (`minperiod` of `demography::lca()`).
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
#' Booth H, Maindonald J, Smith L (2002). Applying Lee-Carter under
#' conditions of variable mortality decline. *Population Studies* 56(3):
#' 325-336. \doi{10.1080/00324720215935}
#'
#' Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
#' method for forecasting mortality. *Demography* 38(4): 537-549.
#' \doi{10.1353/dem.2001.0036}
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
                       n_boot = 100, adjust = c("dt", "dxt", "e0", "none"), choose_period = FALSE,
                       min_years = 20) {
  method <- match.arg(method)
  jump_off <- match.arg(jump_off)
  adjust <- match.arg(adjust)
  check_flag(choose_period, "choose_period")
  check_number(min_years, "min_years", lower = 5, integerish = TRUE)
  new_mt_spec("lee_carter", list(method = method, jump_off = jump_off, bootstrap = bootstrap, n_boot = n_boot,
                                 adjust = adjust, choose_period = choose_period, min_years = min_years),
              min_periods = if (choose_period) min_years else 5,
              label = paste0("Lee-Carter", if (method == "svd") sprintf(" (SVD, %s)", adjust) else "",
                             if (choose_period) ", chosen period" else "",
                             if (jump_off == "observed") ", observed jump-off" else "",
                             if (bootstrap) ", bootstrap" else ""))
}

#' @export
fit_model.mt_spec_lee_carter <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  zero <- rowSums(D) == 0
  if (sum(!zero) < 2) mt_abort("Fewer than two age groups with deaths.", class = "model_fit")
  fallback <- FALSE
  # Poisson estimates, or SVD estimates where the Poisson likelihood has no
  # finite maximum (very sparse counts). `fixed`: constant rates of the age
  # groups left out of D (NA for the rows of D), for life tables over all ages
  est_one <- function(D, E, fixed = NULL) {
    if (p$method == "svd") return(lc_svd(D, E, adjust = p$adjust %||% "dt", ages = ages, fixed = fixed))
    tryCatch(lc_poisson(D, E), mt_error_model_fit = function(e) {
      fallback <<- TRUE
      tryCatch(lc_poisson(D, E, ridge = 1), mt_error_model_fit = function(e) lc_svd(D, E))
    })
  }
  fit_lc <- function(D, E) {
    if (!any(zero)) return(est_one(D, E))
    # age groups without any death: constant rate of 0.5 deaths over the
    # training years, no trend (the Poisson estimate would be minus infinity)
    fixed <- rep(NA_real_, nrow(D))
    fixed[zero] <- 0.5 / rowSums(E[zero, , drop = FALSE])
    e <- est_one(D[!zero, , drop = FALSE], E[!zero, , drop = FALSE], fixed)
    a <- b <- numeric(nrow(D))
    a[!zero] <- e$a
    b[!zero] <- e$b
    a[zero] <- log(fixed[zero])
    e$a <- a
    e$b <- b
    e
  }
  first <- 1L
  if (isTRUE(p$choose_period)) {
    first <- bms_start(D, E, est_one, p$min_years)
    D <- D[, first:ncol(D), drop = FALSE]
    E <- E[, first:ncol(E), drop = FALSE]
    zero <- rowSums(D) == 0
  }
  est <- fit_lc(D, E)
  if (any(zero)) {
    mt_warn(sprintf(paste("No deaths in age group(s) %s in the training years: constant rate (0.5 deaths in all",
                          "years), no trend."), paste(ages$labels[zero], collapse = ", ")), class = "sparse")
  }
  if (fallback) {
    mt_warn(paste("The Poisson likelihood has no finite maximum (sparse counts); the age pattern of change b",
                  "is shrunk towards equal shares (ridge penalty)."), class = "sparse")
  }
  est$last_log_rate <- log_rates(D[, ncol(D)], E[, ncol(E)])
  est$fitted <- exp(est$a + outer(est$b, est$k))
  est$jump_off <- p$jump_off
  if (isTRUE(p$bootstrap)) {
    # resamples of the (chosen) fitting period
    fitted <- E * est$fitted
    est$boot <- lapply(seq_len(p$n_boot), function(i) {
      Db <- matrix(stats::rpois(length(fitted), fitted), nrow(fitted))
      tryCatch(fit_lc(Db, E), error = function(e) NULL)
    })
    est$boot <- Filter(Negate(is.null), est$boot)
  }
  if (first > 1L) {
    # fitted rates only for the chosen period
    est$fitted <- cbind(matrix(NA_real_, nrow(D), first - 1L), est$fitted)
    est$first_year <- periods[first]
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

# `fixed`: for adjust = "e0" with rows of `ages` left out of D, their constant
# rates (NA for the rows of D), so that life tables cover the whole age range
lc_svd <- function(D, E, adjust = "dt", ages = NULL, fixed = NULL) {
  y <- log_rates(D, E)
  a <- rowMeans(y)
  s <- svd(y - a, nu = 1, nv = 1)
  b <- s$u[, 1]
  k <- s$d[1] * s$v[, 1]
  if (sum(b) < 0) {
    b <- -b
    k <- -k
  }
  # second stage: re-estimate k (Lee and Carter: total deaths; Booth et al.:
  # deaths by age; Lee and Miller: life expectancy)
  if (adjust == "e0" && is.null(ages)) adjust <- "dt"
  all_ages <- function(m) {
    if (is.null(fixed)) return(array(m, c(nrow(m), ncol(m), 1)))
    out <- matrix(fixed, length(fixed), ncol(m))
    out[is.na(fixed), ] <- m
    array(out, c(nrow(out), ncol(out), 1))
  }
  if (adjust == "e0") e_obs <- life_expectancy(all_ages(pmax(D, 0.5) / E), ages)[, 1]
  for (t in seq_len(ncol(D))) {
    if (adjust == "none") break
    if (adjust == "dxt") {
      kt <- k[t]
      for (it in 1:50) {
        mu <- E[, t] * exp(a + b * kt)
        step <- sum((D[, t] - mu) * b) / sum(mu * b^2)
        kt <- kt + step
        if (abs(step) < 1e-10) break
      }
      if (is.finite(kt)) k[t] <- kt
      next
    }
    f <- if (adjust == "dt") {
      function(kt) sum(E[, t] * exp(a + b * kt)) - sum(D[, t])
    } else {
      function(kt) life_expectancy(all_ages(matrix(exp(a + b * kt))), ages)[1, 1] - e_obs[t]
    }
    rng <- k[t] + c(-1, 1) * (10 + 10 * abs(k[t]))
    # without a root in reach, total deaths keep the SVD index (with `b` of
    # both signs the nearest root of total deaths can lie far off and makes
    # the index erratic); life expectancy searches for the nearest root
    k[t] <- tryCatch(stats::uniroot(f, rng, extendInt = "yes", tol = 1e-10)$root, error = function(e) {
      if (adjust == "e0") nearest_root(f, k[t], 0.1 * (1 + abs(k[t]))) else k[t]
    })
  }
  out <- lc_finish(a, b, k)
  out$method <- "svd"
  out
}

# Root of f nearest to x0, for f that is not monotone (life expectancy as a
# function of k when b has both signs, with sparse counts): steps of growing
# size in both directions until the sign changes; x0 if there is none
nearest_root <- function(f, x0, step) {
  f0 <- f(x0)
  if (!is.finite(f0) || f0 == 0) return(x0)
  prev <- c(x0, x0)
  for (d in step * 2^(-4:12)) {
    for (i in 1:2) {
      x <- x0 + c(-1, 1)[i] * d
      fx <- f(x)
      if (is.finite(fx) && sign(fx) != sign(f0)) {
        return(tryCatch(stats::uniroot(f, sort(c(prev[i], x)), tol = 1e-10)$root, error = function(e) x0))
      }
      prev[i] <- x
    }
  }
  x0
}

lc_poisson <- function(D, E, init = NULL, maxit = 500, tol = 1e-8, ridge = 0) {
  if (is.null(init)) init <- lc_svd(D, E)
  a <- init$a
  b <- if (ridge > 0) rep(1 / nrow(D), nrow(D)) else init$b
  k <- init$k
  b0 <- 1 / nrow(D)
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
    # a ridge penalty shrinks b towards equal shares where the counts are
    # too sparse for a finite maximum likelihood estimate
    b <- b + (rowSums((D - Dh) * kk) - ridge * (b - b0)) / (rowSums(Dh * kk^2) + ridge)
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

# First year of the fitting period with the smallest ratio of the mean
# deviance of a linear-k fit to the mean deviance of the Lee-Carter fit
# (Booth, Maindonald and Smith 2002)
bms_start <- function(D, E, est_one, min_years) {
  nt <- ncol(D)
  na <- nrow(D)
  if (nt <= min_years + 1) return(1L)
  starts <- seq_len(nt - min_years)
  ratio <- vapply(starts, function(s) {
    cols <- s:nt
    Ds <- D[, cols, drop = FALSE]
    Es <- E[, cols, drop = FALSE]
    est <- tryCatch(est_one(Ds, Es), error = function(e) NULL)
    if (is.null(est)) return(Inf)
    n <- length(cols)
    # linear k through its mean with the mean annual change as slope
    klin <- mean(est$k) + mean(diff(est$k)) * (seq_len(n) - (n + 1) / 2)
    dev_base <- poisson_dev(Ds, Es * exp(est$a + outer(est$b, est$k)))
    dev_lin <- poisson_dev(Ds, Es * exp(est$a + outer(est$b, klin)))
    (dev_lin / (na * (n - 2))) / (dev_base / ((na - 1) * (n - 2)))
  }, 0)
  starts[which.min(ratio)]
}
