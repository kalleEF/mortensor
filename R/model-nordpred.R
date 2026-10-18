#' Nordpred age-period-cohort model with a power link
#'
#' The Nordpred method of Møller et al. (2003) for projections of cancer
#' incidence and mortality: an age-period-cohort model for 5-year age groups
#' and 5-year periods, fitted by Poisson maximum likelihood with a power link,
#' whose linear trend (the drift) is projected and damped in later periods.
#' The deaths `D` of age group `a` in period `p` with exposure `E` are
#' Poisson with mean
#'
#' `E (A[a] + d p + P[p] + C[c])^5`
#'
#' (`link = "power5"`), with age effects `A`, drift `d`, non-linear period
#' effects `P` and cohort effects `C` (cohort `c = p - a` on the scale of
#' 5-year steps), or `E exp(A[a] + d p + P[p] + C[c])` (`link = "log"`, the
#' multiplicative model). The power 5 levels off the exponential growth of
#' the multiplicative model (Møller et al. 2003). The model is fitted as a
#' generalised linear model with a Poisson variance function and the power
#' link `mu^(1/5)` on the rates, written for this package.
#'
#' The specification follows Møller et al. (2003) and the R package nordpred
#' (version 1.1, by Møller and Weedon-Fekjær, GPL-2), whose documentation and
#' code define the details below. This implementation is written
#' independently of nordpred's code. Its point projections agree with those
#' of nordpred on the same data to a relative difference below `1e-7` (see
#' the unit tests).
#'
#' * **Identification.** The period effects of the first and the last period
#'   and the cohort effects of the oldest and the youngest cohort are zero,
#'   so that the drift carries the linear trend (the identification of
#'   nordpred).
#' * **Prediction base.** `n_periods` lists the candidate numbers of recent
#'   5-year periods to fit. The model is fitted to the largest; while its
#'   deviance rejects the fit at the 1% level (chi-squared goodness of fit),
#'   the oldest period is dropped, down to the smallest candidate.
#' * **Trend.** With `recent = NULL` the projected drift is the slope of the
#'   last 10 years (the change of the period part of the linear predictor
#'   between the last two periods) if a quadratic period term is significant
#'   at the 5% level in a Poisson model with log link (age, drift, quadratic
#'   term and cohort effects), and the average drift `d` otherwise.
#' * **Cut trend.** In future period `j` the drift is multiplied by `1 -
#'   cut_trend[j]` (cumulated over the periods): the default reduces it by
#'   0%, 25%, 50%, 75% and 75% in the five 5-year periods after the last
#'   observation.
#' * **Cohorts.** A future cell keeps the effect of its cohort if that cohort
#'   is older than the cohort of age group `use_age` in the last observed
#'   period, and takes the effect of that cohort otherwise.
#' * **Young ages.** Age groups below `use_age` are projected at the mean of
#'   the observed rates of the last two 5-year periods; age groups below
#'   `est_age` do not enter the regression.
#'
#' **Annual data (package-specific).** mortensor models work on single
#' calendar years. The deaths and exposures of the training years are summed
#' to 5-year periods that end with the last training year; up to four
#' earlier years that do not fill a 5-year period are left out. The model is
#' fitted and projected on these periods. Annual rates are obtained by linear
#' interpolation of the linear predictor (the fifth root of the rate, or its
#' logarithm with `link = "log"`) between the midpoints of the periods: the
#' fitted value of the last observed period at its midpoint (two years
#' before the last training year) and the projections of the future periods
#' at theirs. Where the linear predictor changes linearly over the periods,
#' its mean over the five years of a future period equals the Nordpred
#' projection of that period. Fitted rates are the fitted period rates in
#' each year of the periods fitted, and missing for earlier years and for
#' age groups below `est_age`. The age groups must be 5 years wide; the last
#' may be open or wider (it counts as one 5-year step for the cohorts, as
#' nordpred's 85+). Use [mt_regroup()] or `step_regroup()` for single years
#' of age.
#'
#' **Uncertainty (package-specific).** Nordpred gives point projections and
#' documents no method for prediction intervals. Here each draw takes the
#' coefficients from their asymptotic normal distribution (estimates and the
#' inverse Fisher information of the Poisson model) and projects them as
#' above. Draws describe the uncertainty of the coefficients only, with no
#' variation of future period effects and no overdispersion, so intervals are
#' narrow when the data are large; check coverage in a backtest and widen
#' them with [recalibrated()] if needed. With the power link, linear
#' predictors below 0.001 of a draw are set to 0.001 (a safeguard). Age
#' groups below `use_age` have no uncertainty. Predictive draws
#' (`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
#' given the rates.
#'
#' Beyond five future periods (25 years) the last value of `cut_trend` is
#' used (nordpred projects at most five periods). The candidate number of
#' periods is chosen with the link of the model (nordpred uses the power
#' link for this choice also when `link = "log"`).
#'
#' @param link `"power5"` (Nordpred's power link) or `"log"` (multiplicative
#'   Poisson age-period-cohort model).
#' @param n_periods Candidate numbers of 5-year periods in the prediction
#'   base, at least 3 each (nordpred's default `4:6`). A single number fixes
#'   the base.
#' @param recent Project the slope of the last 10 years (`TRUE`), the
#'   average drift (`FALSE`), or choose by the test described above
#'   (`NULL`).
#' @param cut_trend Proportions of the drift cut in the future 5-year
#'   periods.
#' @param est_age Lower bound of the youngest age group in the regression
#'   (nordpred's `startestage`). Default: the first age group of the grid.
#' @param use_age Lower bound of the youngest age group projected with the
#'   regression (nordpred's `startuseage`). Default: `est_age`.
#' @return A model specification.
#' @references
#' Møller B, Fekjær H, Hakulinen T, Sigvaldason H, Storm HH, Talbäck M,
#' Haldorsen T (2003). Prediction of cancer incidence in the Nordic
#' countries: empirical comparison of different approaches. *Statistics in
#' Medicine* 22(17): 2751-2766. \doi{10.1002/sim.1481}
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2019)
#' x5 <- mt_regroup(x, age = seq(30, 85, 5))
#' fit <- mt_fit(nordpred_apc(), x5)
#' fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
#' @export
nordpred_apc <- function(link = c("power5", "log"), n_periods = 4:6, recent = NULL,
                         cut_trend = c(0, 0.25, 0.5, 0.75, 0.75), est_age = NULL, use_age = NULL) {
  link <- match.arg(link)
  if (!is.numeric(n_periods) || !length(n_periods) || anyNA(n_periods) || any(n_periods < 3) ||
      any(abs(n_periods - round(n_periods)) > 1e-8)) {
    mt_abort("`n_periods` must be whole numbers of at least 3.", class = "argument")
  }
  if (!is.null(recent)) check_flag(recent, "recent")
  if (!is.numeric(cut_trend) || !length(cut_trend) || anyNA(cut_trend) || any(cut_trend < 0 | cut_trend > 1)) {
    mt_abort("`cut_trend` must be proportions between 0 and 1.", class = "argument")
  }
  if (!is.null(est_age)) check_number(est_age, "est_age")
  if (!is.null(use_age)) check_number(use_age, "use_age")
  if (!is.null(est_age) && !is.null(use_age) && use_age < est_age) {
    mt_abort("`use_age` must not be below `est_age`.", class = "argument")
  }
  new_mt_spec("nordpred", list(link = link, n_periods = sort(unique(round(n_periods))), recent = recent,
                               cut_trend = cut_trend, est_age = est_age, use_age = use_age),
              min_periods = 15, label = paste0("Nordpred", if (link == "log") " (log link)" else ""))
}

#' @export
fit_model.mt_spec_nordpred <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  lo <- ages$lower
  na <- length(lo)
  w <- ages$upper - lo
  if (na < 3 || any(abs(utils::head(w, -1) - 5) > 1e-9)) {
    mt_abort(paste("nordpred_apc() needs 5-year age groups (the last group may be open or wider); regroup with",
                   "mt_regroup(x, age = seq(<first age>, <last age>, 5)) or step_regroup()."),
             class = "model_data", call = NULL)
  }
  nt <- length(periods)
  nb <- nt %/% 5L
  if (nb < 3L) mt_abort("nordpred_apc() needs at least 15 years (three 5-year periods).", class = "model_data")
  keep <- (nt - 5L * nb + 1L):nt
  blk <- rep(seq_len(nb), each = 5L)
  Db <- t(rowsum(t(D[, keep, drop = FALSE]), blk, reorder = FALSE))
  Eb <- t(rowsum(t(E[, keep, drop = FALSE]), blk, reorder = FALSE))
  ie <- nordpred_age_index(p$est_age %||% lo[1], lo, "est_age")
  iu <- nordpred_age_index(p$use_age %||% lo[ie], lo, "use_age")
  if (iu < ie) mt_abort("`use_age` must not be below `est_age`.", class = "argument", call = NULL)
  rows <- ie:na
  if (length(rows) < 3L) {
    mt_abort("nordpred_apc() needs at least three age groups from `est_age`.", class = "model_data")
  }
  # prediction base: drop the oldest period while the fit is rejected
  cand <- sort(unique(pmin(p$n_periods, nb)))
  while (length(cand) > 1L) {
    cols <- (nb - max(cand) + 1L):nb
    f <- nordpred_glm(Db[rows, cols, drop = FALSE], Eb[rows, cols, drop = FALSE], p$link)
    cand <- if (f$gof_p < 0.01) utils::head(cand, -1L) else max(cand)
  }
  np <- cand
  cols <- (nb - np + 1L):nb
  f <- nordpred_glm(Db[rows, cols, drop = FALSE], Eb[rows, cols, drop = FALSE], p$link)
  p_recent <- nordpred_recent_test(Db[rows, cols, drop = FALSE], Eb[rows, cols, drop = FALSE])
  recent <- p$recent %||% (p_recent < 0.05)
  # age groups below use_age: mean of the observed rates of the last two periods
  last2 <- (nb - 1L):nb
  young <- rowMeans(Db[, last2, drop = FALSE] / Eb[, last2, drop = FALSE])
  young <- pmax(young, 0.5 / rowSums(Eb[, last2, drop = FALSE]))
  # fitted rates: fitted period rates in each year of the periods fitted
  fitted <- matrix(NA_real_, na, nt)
  yrs <- keep[blk %in% cols]
  fitted[rows, yrs] <- f$fitted[, match(blk[blk %in% cols], cols)]
  first <- min(periods[keep]) + 5L * (cols - 1L)
  list(coef = f$coef, vcov = f$vcov, link = p$link, n_age = length(rows), n_periods = as.integer(np), rows = rows,
       use = iu - ie + 1L, young = young, recent = recent, p_recent = p_recent, gof_p = f$gof_p,
       deviance = f$deviance, df = f$df, cut_trend = p$cut_trend, fitted = fitted,
       period_labels = paste0(first, "-", first + 4L), first_years = first, last_year = max(periods),
       age_labels = ages$labels, age_lower = lo)
}

nordpred_age_index <- function(age, lower, name) {
  i <- which(abs(lower - age) < 1e-9)
  if (length(i) != 1L) {
    mt_abort(sprintf("`%s` (%s) must be the lower bound of an age group of the grid.", name, format(age)),
             class = "argument", call = NULL)
  }
  i
}

# Design matrix of the age-period-cohort model for `na` age groups and `np`
# periods, cells ordered with age fastest: age dummies, linear period index
# (the drift), period dummies of periods 2 to np - 1 and cohort dummies of
# cohorts 2 to nk - 1 (cohort k = p - a + na, 1 the oldest). The effects of
# the first and last period and cohort are zero. `periods = FALSE` leaves out
# the period dummies, `quadratic = TRUE` adds the squared period index.
nordpred_design <- function(na, np, periods = TRUE, quadratic = FALSE) {
  a <- rep(seq_len(na), np)
  pp <- rep(seq_len(np), each = na)
  k <- pp - a + na
  nk <- na + np - 1L
  dummies <- function(v, levels, prefix) {
    if (!length(levels)) return(NULL)
    m <- outer(v, levels, "==") * 1
    colnames(m) <- paste0(prefix, levels)
    m
  }
  X <- cbind(dummies(a, seq_len(na), "age"), drift = pp,
             if (periods) dummies(pp, seq_len(np)[-c(1L, np)], "period"),
             if (quadratic) cbind(drift2 = pp^2),
             dummies(k, seq_len(nk)[-c(1L, nk)], "cohort"))
  X
}

# Poisson family for rates y = D / E with prior weights E: the deviance and
# the score are those of Poisson deaths D with mean E * mu. Link mu^(1/5)
# (power5) or log.
nordpred_family <- function(link) {
  lk <- if (link == "power5") stats::power(1 / 5) else stats::make.link("log")
  structure(list(
    family = "poisson", link = lk$name, linkfun = lk$linkfun, linkinv = lk$linkinv, mu.eta = lk$mu.eta,
    valideta = lk$valideta, variance = function(mu) mu, validmu = function(mu) all(is.finite(mu)) && all(mu > 0),
    dev.resids = function(y, mu, wt) 2 * wt * (ifelse(y > 0, y * log(y / mu), 0) - (y - mu)),
    aic = function(y, n, mu, wt, dev) {
      d <- wt * y
      m <- wt * mu
      -2 * sum(ifelse(d > 0, d * log(m), 0) - m - lgamma(d + 1))
    },
    initialize = expression({
      if (any(y < 0)) stop("negative rates")
      n <- rep.int(1, nobs)
      mustart <- y + 0.5 / pmax(weights, 1e-8)
    })), class = "family")
}

# Fit of the model to deaths and exposures [age, period]
nordpred_glm <- function(D, E, link, X = nordpred_design(nrow(D), ncol(D))) {
  y <- as.vector(D) / as.vector(E)
  wt <- as.vector(E)
  f <- stats::glm.fit(X, y, weights = wt, family = nordpred_family(link),
                      control = stats::glm.control(epsilon = 1e-10, maxit = 100))
  if (!isTRUE(f$converged) || anyNA(f$coefficients)) {
    mt_abort("The Nordpred model did not converge.", class = "model_fit", call = NULL)
  }
  V <- tryCatch(chol2inv(chol(crossprod(X * sqrt(f$weights)))), error = function(e) NULL)
  if (is.null(V)) mt_abort("The information matrix of the Nordpred model is singular.", class = "model_fit")
  dimnames(V) <- list(colnames(X), colnames(X))
  list(coef = f$coefficients, vcov = V, deviance = f$deviance, df = f$df.residual,
       gof_p = stats::pchisq(f$deviance, f$df.residual, lower.tail = FALSE),
       fitted = matrix(f$fitted.values, nrow(D)))
}

# p-value of the test for departure from a linear trend: Poisson models with
# log link, age, drift and cohort effects, without and with a quadratic
# period term (likelihood ratio test)
nordpred_recent_test <- function(D, E) {
  na <- nrow(D)
  np <- ncol(D)
  f1 <- nordpred_glm(D, E, "log", nordpred_design(na, np, periods = FALSE))
  f2 <- nordpred_glm(D, E, "log", nordpred_design(na, np, periods = FALSE, quadratic = TRUE))
  stats::pchisq(max(f1$deviance - f2$deviance, 0), 1, lower.tail = FALSE)
}

# Coefficient weights of the linear predictor of the age groups projected by
# the regression (from use_age), for the last observed period (j = 0) and
# future periods j = 1..J: rows ordered with age fastest
nordpred_lincomb <- function(fit, J) {
  nm <- names(fit$coef)
  na <- fit$n_age
  np <- fit$n_periods
  ages <- fit$use:na
  col <- function(name) match(name, nm)
  cut <- fit$cut_trend[pmin(seq_len(J), length(fit$cut_trend))]
  steps <- c(0, cumsum(1 - cut))
  # cohort of age group use_age in the last observed period; younger cohorts
  # take its effect
  k_use <- np - fit$use + na
  L <- matrix(0, length(ages) * (J + 1L), length(nm))
  r <- 0L
  for (j in 0:J) {
    for (a in ages) {
      r <- r + 1L
      L[r, col(paste0("age", a))] <- 1
      L[r, col("drift")] <- np + steps[j + 1L]
      if (j > 0 && fit$recent) L[r, col(paste0("period", np - 1L))] <- -steps[j + 1L]
      k <- np + j - a + na
      if (j > 0 && k >= k_use) k <- k_use
      ck <- col(paste0("cohort", k))
      if (!is.na(ck)) L[r, ck] <- 1
    }
  }
  L
}

# Weights [h, J + 1] of the linear interpolation between the midpoints of the
# periods (j = 0: last observed period, 2 years before the last year; period
# j at 5 j - 2 years after it) for years 1..h after the last year
nordpred_interp <- function(h, J) {
  mid <- 5 * (0:J) - 2
  W <- matrix(0, h, J + 1L)
  for (y in seq_len(h)) {
    j <- findInterval(y, mid, rightmost.closed = TRUE)
    u <- (y - mid[j]) / 5
    W[y, j] <- 1 - u
    if (u > 0) W[y, j + 1L] <- u
  }
  W
}

nordpred_inverse <- function(eta, link) if (link == "power5") pmax(eta, 1e-3)^5 else exp(eta)

# Point projections of the 5-year periods 1..J after the last observed one,
# [age, period]: the Nordpred projections on the period scale
nordpred_projection <- function(fit, J = 5) {
  L <- nordpred_lincomb(fit, J)
  nu <- fit$n_age - fit$use + 1L
  eta <- matrix(L %*% fit$coef, nu)[, -1L, drop = FALSE]
  out <- matrix(fit$young, length(fit$young), J)
  out[fit$rows[fit$use:fit$n_age], ] <- nordpred_inverse(eta, fit$link)
  dimnames(out) <- list(fit$age_labels, paste0("+", seq_len(J)))
  out
}

#' @export
simulate_model.mt_spec_nordpred <- function(spec, fit, h, n_draws, ...) {
  J <- ceiling((h + 2) / 5)
  L <- nordpred_lincomb(fit, J)
  nu <- fit$n_age - fit$use + 1L
  B <- draw_mvnorm(n_draws, fit$coef, fit$vcov)
  eta <- array(L %*% B, c(nu, J + 1L, n_draws))
  W <- nordpred_interp(h, J)
  eta_y <- array(apply(eta, 3, function(m) m %*% t(W)), c(nu, h, n_draws))
  out <- array(fit$young, c(length(fit$young), h, n_draws))
  out[fit$rows[fit$use:fit$n_age], , ] <- nordpred_inverse(eta_y, fit$link)
  out
}

# n draws [length(mu), n] from the normal distribution with mean mu and
# covariance V (eigen decomposition, negative eigenvalues of rounding set to 0)
draw_mvnorm <- function(n, mu, V) {
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  R <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
  mu + R %*% matrix(stats::rnorm(length(mu) * n), length(mu))
}
