# Negative binomial Lee-Carter model with a dispersion by age and the
# forecast origin correction (Li, Hardy and Tan 2009). The paper's data (US
# and Canada, Human Mortality Database, ages 0-99, 1950-2004) are not in the
# demography or StMoMo packages, so its estimates cannot be reproduced; the
# tests check the likelihood, recovery on simulated data, limiting cases and
# exact identities.

# Deaths from the negative binomial Lee-Carter model with dispersion alpha[x]
# by row (variance mu + alpha[x] mu^2); sum(b) = 1, sum(k) = 0
sim_lc_disp <- function(seed, nx = 30, ny = 50, exposure = 2e5, alpha = rep(c(0.005, 0.01, 0.02), each = 10)) {
  set.seed(seed)
  a <- -5 + 0.1 * (seq_len(nx) - 1)
  b <- 1 + 0.5 * sin(seq_len(nx) / 5)
  b <- b / sum(b)
  k <- seq(20, -20, length.out = ny)
  k <- k - mean(k)
  E <- matrix(exposure, nx, ny)
  mu <- E * exp(a + outer(b, k))
  D <- matrix(stats::rnbinom(nx * ny, size = 1 / rep_len(alpha, nx * ny), mu = mu), nx)
  list(D = D, E = E, mu = mu, a = a, b = b, k = k, alpha = alpha)
}

# Log-likelihood of deaths D with means mu and one dispersion per row, from
# dnbinom() and dpois() (dispersion 0)
ll_ref <- function(D, mu, alpha) {
  alpha <- rep_len(alpha, length(D))
  pois <- alpha == 0
  sum(stats::dnbinom(D[!pois], size = 1 / alpha[!pois], mu = mu[!pois], log = TRUE)) +
    sum(stats::dpois(D[pois], mu[pois], log = TRUE))
}

# Central differences of f at v
num_grad <- function(f, v, eps = 1e-6) {
  vapply(seq_along(v), function(i) {
    e <- replace(numeric(length(v)), i, eps)
    (f(v + e) - f(v - e)) / (2 * eps)
  }, 0)
}

test_that("the log-likelihood with one dispersion per row is that of eq. (13)", {
  set.seed(1)
  mu <- matrix(stats::runif(200, 0.5, 300), 4)
  alpha <- c(0.3, 0, 1e-6, 0.02)
  D <- matrix(stats::rnbinom(200, size = 20, mu = mu), 4)
  expect_equal(mortensor:::nb_loglik(D, mu, alpha), ll_ref(D, mu, alpha), tolerance = 1e-10)
  # eq. (13) term by term, with the constant c = -sum(log(D!))
  eq13 <- function(D, lam, al) {
    sum(vapply(seq_along(D), function(i) {
      d <- D[i]
      s <- if (d > 0) sum(log((1 + al[i] * (seq_len(d) - 1)) / al[i])) else 0
      s + d * log(al[i] * lam[i]) - (d + 1 / al[i]) * log1p(al[i] * lam[i]) - lgamma(d + 1)
    }, 0))
  }
  pos <- c(0.3, 0.05, 0.002, 0.02)
  expect_equal(mortensor:::nb_loglik(D, mu, pos), eq13(D, mu, rep_len(pos, length(D))), tolerance = 1e-10)
  # one value for all rows: the global dispersion
  expect_equal(mortensor:::nb_loglik(D, mu, rep(0.05, 4)), mortensor:::nb_loglik(D, mu, 0.05), tolerance = 1e-12)
})

test_that("the fit with a dispersion per group of ages solves the likelihood equations of eq. (13)", {
  x <- sim_lc_disp(1, nx = 12, ny = 30, alpha = rep(c(0.004, 0.02, 0.01), each = 4))
  g <- rep(1:3, each = 4)
  f <- mortensor:::lc_nb(x$D, x$E, groups = g)
  expect_length(f$dispersion, 12)
  expect_identical(f$age_nb_size, 1 / f$dispersion)
  expect_null(f$nb_size)
  expect_equal(f$loglik, ll_ref(x$D, x$E * exp(f$a + outer(f$b, f$k)), f$dispersion), tolerance = 1e-10)
  ll <- function(a = f$a, b = f$b, k = f$k, al = f$dispersion) ll_ref(x$D, x$E * exp(a + outer(b, k)), al)
  scale <- sum(abs(x$D - x$E * exp(f$a + outer(f$b, f$k))))
  expect_lt(max(abs(num_grad(function(v) ll(a = v), f$a))) / scale, 1e-6)
  expect_lt(max(abs(num_grad(function(v) ll(b = v), f$b))) / scale, 1e-6)
  expect_lt(max(abs(num_grad(function(v) ll(k = v), f$k))) / scale, 1e-6)
  # each dispersion maximises the likelihood of its group (relative change
  # of the dispersion)
  for (j in 1:3) {
    al <- f$dispersion[g == j][1]
    dl <- (ll(al = ifelse(g == j, al * (1 + 1e-5), f$dispersion)) -
             ll(al = ifelse(g == j, al * (1 - 1e-5), f$dispersion))) / 2e-5
    expect_lt(abs(dl), 1e-3)
  }
  expect_equal(sum(f$b), 1, tolerance = 1e-12)
  expect_equal(sum(f$k), 0, tolerance = 1e-9)
})

test_that("the forecast origin correction fits the observed rates of the last year (Appendix 2)", {
  x <- sim_lc_disp(2, nx = 12, ny = 30, alpha = rep(c(0.004, 0.02, 0.01), each = 4))
  g <- rep(1:3, each = 4)
  f <- mortensor:::lc_nb(x$D, x$E, groups = g, origin = TRUE)
  free <- mortensor:::lc_nb(x$D, x$E, groups = g)
  # a = log of the observed rates of the last year, k[T] = 0, sum(b) = 1;
  # 0.5 deaths where there are none
  expect_equal(f$a, log(x$D[, 30] / x$E[, 30]), tolerance = 1e-14)
  x0 <- x
  x0$D[1, 30] <- 0
  x0$E[1, ] <- 200
  expect_equal(mortensor:::lc_nb(x0$D, x0$E, groups = g, origin = TRUE)$a[1], log(0.5 / 200), tolerance = 1e-14)
  expect_identical(f$k[30], 0)
  expect_equal(sum(f$b), 1, tolerance = 1e-12)
  # b, k[1..T-1] and the dispersions maximise the likelihood under these
  # constraints
  ll <- function(b = f$b, k = f$k, al = f$dispersion) ll_ref(x$D, x$E * exp(f$a + outer(b, k)), al)
  scale <- sum(abs(x$D - x$E * exp(f$a + outer(f$b, f$k))))
  expect_lt(max(abs(num_grad(function(v) ll(b = v), f$b))) / scale, 1e-6)
  expect_lt(max(abs(num_grad(function(v) ll(k = c(v, 0)), f$k[-30]))) / scale, 1e-6)
  for (j in 1:3) {
    al <- f$dispersion[g == j][1]
    dl <- (ll(al = ifelse(g == j, al * (1 + 1e-5), f$dispersion)) -
             ll(al = ifelse(g == j, al * (1 - 1e-5), f$dispersion))) / 2e-5
    expect_lt(abs(dl), 1e-3)
  }
  # the constraints lower the log-likelihood (footnote 3)
  expect_lt(f$loglik, free$loglik)
  # Poisson model with the correction (lee_carter(jump_off = "constrained")):
  # the Poisson likelihood equations under the same constraints
  long <- data.frame(age = rep(60:71, 30), year = rep(1981:2010, each = 12), deaths = as.vector(x$D),
                     exposure = as.vector(x$E))
  fp <- mt_fit(lee_carter(jump_off = "constrained"), mt_data(deaths = long, exposure = long))$fits[[1]]
  expect_null(fp$dispersion)
  expect_identical(fp$method, "poisson")
  expect_equal(unname(fp$a), f$a, tolerance = 1e-14)
  llp <- function(b = fp$b, k = fp$k) ll_ref(x$D, x$E * exp(fp$a + outer(b, k)), 0)
  expect_lt(max(abs(num_grad(function(v) llp(b = v), fp$b))) / scale, 1e-6)
  expect_lt(max(abs(num_grad(function(v) llp(k = c(v, 0)), fp$k[-30]))) / scale, 1e-6)
  expect_equal(unname(fp$fitted[, 30]), x$D[, 30] / x$E[, 30], tolerance = 1e-12)
})

test_that("the dispersions by group of ages and by age are recovered from simulated data", {
  # a, b and k within sampling error: standard errors from the expected
  # information (weights mu / (1 + alpha[x] mu)) under sum(b) = 1, sum(k) = 0
  x <- sim_lc_disp(3)
  g <- rep(1:3, each = 10)
  f <- mortensor:::lc_nb(x$D, x$E, groups = g)
  nx <- length(f$a)
  ny <- length(f$k)
  mu <- x$E * exp(f$a + outer(f$b, f$k))
  w <- as.vector(mu / (1 + f$dispersion * mu))
  ix <- rep(seq_len(nx), ny)
  iy <- rep(seq_len(ny), each = nx)
  J <- cbind(outer(ix, seq_len(nx), "==") * 1, outer(ix, seq_len(nx), "==") * f$k[iy],
             outer(iy, seq_len(ny), "==") * f$b[ix])
  info <- crossprod(J * sqrt(w))
  H <- rbind(c(rep(0, nx), rep(1, nx), rep(0, ny)), c(rep(0, 2 * nx), rep(1, ny)))
  cov <- solve(rbind(cbind(info, t(H)), cbind(H, matrix(0, 2, 2))))[seq_len(ncol(info)), seq_len(ncol(info))]
  z <- (c(f$a, f$b, f$k) - c(x$a, x$b, x$k)) / sqrt(diag(cov))
  expect_lt(max(abs(z)), 4.5)
  expect_gt(stats::sd(z), 0.7)
  expect_lt(stats::sd(z), 1.4)
  # the dispersions over 20 independent data sets. Maximum likelihood
  # estimates of a variance are biased downwards by the share of the
  # information used for the other parameters: a[x] and b[x] are estimated
  # from the 50 years of each age, and k[t] mostly from the ages with the
  # largest weight mu / (1 + alpha mu) (here the first group). The mean
  # estimate lies within 18% below the truth
  truth <- c(0.005, 0.01, 0.02)
  reps <- 20
  sets <- lapply(seq_len(reps), function(i) sim_lc_disp(100 + i))
  est <- t(vapply(sets, function(s) {
    tapply(mortensor:::lc_nb(s$D, s$E, groups = g)$dispersion, g, mean)
  }, numeric(3)))
  se <- apply(est, 2, stats::sd) / sqrt(reps)
  expect_true(all(colMeans(est) > 0.82 * truth))
  expect_true(all(colMeans(est) < truth + 3 * se))
  # given the true expected deaths, the estimates have no such bias
  est0 <- t(vapply(sets, function(s) {
    vapply(1:3, function(j) mortensor:::nb_dispersion(s$D[g == j, ], s$mu[g == j, ]), 0)
  }, numeric(3)))
  se0 <- apply(est0, 2, stats::sd) / sqrt(reps)
  expect_true(all(abs(colMeans(est0) - truth) < 4 * se0))
  # one dispersion per age: the mean ratio to the truth over ages and data
  # sets, and the ordering of the groups
  ratio <- t(vapply(sets[1:10], function(s) {
    mortensor:::lc_nb(s$D, s$E, groups = seq_len(30))$dispersion / s$alpha
  }, numeric(30)))
  expect_gt(mean(ratio), 0.82)
  expect_lt(mean(ratio), 1.02)
  by_group <- tapply(colMeans(ratio) * x$alpha, g, mean)
  expect_true(all(diff(by_group) > 0))
})

test_that("one dispersion group gives the fit of family = \"nb\" (Delwarde, Denuit and Partrat 2007)", {
  x <- sim_lc_disp(4, nx = 12, ny = 30, alpha = rep(c(0.004, 0.02, 0.01), each = 4))
  f1 <- mortensor:::lc_nb(x$D, x$E, groups = rep(1L, 12))
  f0 <- mortensor:::lc_nb(x$D, x$E)
  expect_equal(f1$a, f0$a, tolerance = 1e-12)
  expect_equal(f1$b, f0$b, tolerance = 1e-12)
  expect_equal(f1$k, f0$k, tolerance = 1e-12)
  expect_equal(f1$dispersion, rep(f0$dispersion, 12), tolerance = 1e-12)
  expect_equal(f1$loglik, f0$loglik, tolerance = 1e-12)
})

test_that("as the dispersions vanish the fit is the Poisson Lee-Carter fit", {
  x <- sim_lc_disp(5, nx = 12, ny = 30, alpha = rep(c(0.004, 0.02, 0.01), each = 4))
  g <- rep(1:3, each = 4)
  ref <- mortensor:::lc_poisson(x$D, x$E, tol = 1e-15, maxit = 5000)
  gap <- vapply(10^c(-3, -5, -7, -9), function(eps) {
    f <- mortensor:::lc_nb(x$D, x$E, disp = eps * c(1, 3, 2)[g], groups = g)
    max(abs(f$b - ref$b), abs(f$a - ref$a), abs(f$k - ref$k) / max(abs(ref$k)))
  }, 0)
  expect_true(all(diff(gap) < 0))
  expect_lt(gap[4], 1e-6)
  # without overdispersion every dispersion is estimated as 0 and the fit is
  # the Poisson maximum likelihood fit
  D0 <- round(x$mu)
  f0 <- mortensor:::lc_nb(D0, x$E, groups = seq_len(12))
  p0 <- mortensor:::lc_poisson(D0, x$E, tol = 1e-15, maxit = 5000)
  expect_identical(f0$dispersion, rep(0, 12))
  expect_identical(f0$age_nb_size, rep(Inf, 12))
  expect_equal(f0$a, p0$a, tolerance = 1e-7)
  expect_equal(f0$b, p0$b, tolerance = 1e-7)
  expect_equal(f0$k, p0$k, tolerance = 1e-7)
  # with the forecast origin correction, the limit is the Poisson fit with
  # the same constraints
  pc <- mortensor:::lc_nb(x$D, x$E, disp = 0, origin = TRUE)
  fc <- mortensor:::lc_nb(x$D, x$E, disp = rep(1e-10, 12), groups = g, origin = TRUE)
  expect_equal(fc$b, pc$b, tolerance = 1e-7)
  expect_equal(fc$k, pc$k, tolerance = 1e-7)
})

test_that("the settings of the dispersion and the jump-off are checked", {
  expect_identical(lee_carter(family = "nb", dispersion = "age")$params$dispersion, "age")
  expect_identical(lee_carter(family = "nb")$params$dispersion, "global")
  expect_error(lee_carter(dispersion = "age"), class = "mt_error_argument")
  expect_error(lee_carter(method = "svd", jump_off = "constrained"), class = "mt_error_argument")
  expect_error(lee_carter(jump_off = "constrained", choose_period = TRUE), class = "mt_error_argument")
})

test_that("lee_carter(dispersion = \"age\") fits through mt_fit() with one dispersion per age", {
  x <- mt_example(deaths = "single", years = 1995:2015)
  spec <- lee_carter(family = "nb", dispersion = "age", jump_off = "constrained")
  expect_match(spec$label, "dispersion by age")
  expect_match(spec$label, "constrained jump-off")
  fit <- mt_fit(spec, x)
  expect_true(all(fit$diagnostics$status == "ok"))
  p <- mt_parameters(fit)
  for (s in unique(p$stratum)) {
    disp <- p[stratum == s & parameter == "nb_dispersion"]
    expect_identical(disp$index, fit$age$labels)
    expect_true(all(disp$value >= 0))
    expect_equal(sum(p[stratum == s & parameter == "b"]$value), 1, tolerance = 1e-10)
    k <- p[stratum == s & parameter == "k"]$value
    expect_identical(k[length(k)], 0)
  }
  # the fitted rates of the last year are the observed ones (0.5 deaths
  # where there are none)
  f <- fit$fits[[1]]
  Dl <- fit$train_D[, 21, 1]
  expect_equal(f$fitted[, 21], pmax(Dl, 0.5) / fit$train_E[, 21, 1], tolerance = 1e-12)
  a <- summary(fit)$assumptions
  expect_match(a, "Li, Hardy and Tan 2009")
  expect_match(a, "alpha\\(x\\)")
  expect_match(a, "constraints of the fit")
  # defaults unchanged
  expect_identical(lee_carter()$params$jump_off, "fitted")
  expect_identical(lee_carter()$label, "Lee-Carter")
})

test_that("predictive draws are negative binomial with the dispersion of each age", {
  x <- sim_lc_disp(6, nx = 10, ny = 25, alpha = rep(c(0.002, 0.03), each = 5))
  long <- data.frame(age = rep(60:69, 25), year = rep(1991:2015, each = 10), deaths = as.vector(x$D),
                     exposure = as.vector(x$E))
  d <- mt_data(deaths = long, exposure = long)
  fit <- mt_fit(lee_carter(family = "nb", dispersion = "age"), d)
  alpha <- fit$fits[[1]]$dispersion
  expect_length(alpha, 10)
  expect_gt(mean(alpha[6:10]), 5 * mean(alpha[1:5]))
  n <- 4000
  fc <- mt_forecast(fit, h = 1, n_draws = n, seed = 1, predictive = TRUE)
  expect_equal(dim(fc$nb_size), c(10, 1, n))
  expect_equal(fc$nb_size[, 1, 1], 1 / alpha)
  # the rates of the draws are those without `predictive` (same seed): the
  # deaths given them have variance mu + alpha[x] mu^2
  rates <- mt_forecast(fit, h = 1, n_draws = n, seed = 1)$rates
  E <- x$E[, 25]
  mu <- rates[, 1, 1, ] * E
  z <- (fc$rates[, 1, 1, ] * E - mu) / sqrt(mu + alpha * mu^2)
  expect_lt(max(abs(apply(z, 1, stats::var) - 1)), 0.08)
  expect_lt(max(abs(rowMeans(z))), 0.08)
  # the same through a wrapper that adds the observation noise itself
  fa <- mt_fit(anchored(lee_carter(family = "nb", dispersion = "age")), d)
  pa <- mt_forecast(fa, h = 1, n_draws = n, seed = 2, predictive = TRUE)
  ra <- mt_forecast(fa, h = 1, n_draws = n, seed = 2)$rates
  mu <- ra[, 1, 1, ] * E
  z <- (pa$rates[, 1, 1, ] * E - mu) / sqrt(mu + alpha * mu^2)
  expect_lt(max(abs(apply(z, 1, stats::var) - 1)), 0.08)
})

test_that("the parametric bootstrap re-estimates a, b and k from negative binomial resamples (Section 3)", {
  x <- sim_lc_disp(7, nx = 10, ny = 25, alpha = rep(c(0.002, 0.03), each = 5))
  long <- data.frame(age = rep(60:69, 25), year = rep(1991:2015, each = 10), deaths = as.vector(x$D),
                     exposure = as.vector(x$E))
  d <- mt_data(deaths = long, exposure = long)
  nb <- 150
  fit <- with_seed_test(1, mt_fit(lee_carter(family = "nb", dispersion = "age", jump_off = "constrained",
                                             bootstrap = TRUE, n_boot = nb), d))
  f <- fit$fits[[1]]
  expect_length(f$boot, nb)
  # the refits keep the dispersions (step 2a), have k[T] = 0 and take a from
  # the resampled deaths of the last year: its variance across refits is
  # that of the log of negative binomial deaths, about 1 / mu + alpha[x]
  expect_true(all(vapply(f$boot, function(s) identical(s$dispersion, f$dispersion), logical(1))))
  expect_true(all(vapply(f$boot, function(s) s$k[length(s$k)], 0) == 0))
  A <- vapply(f$boot, `[[`, numeric(10), "a")
  mu <- exp(f$a) * x$E[, 25]
  ratio <- apply(A, 1, stats::var) / (1 / mu + f$dispersion)
  expect_gt(mean(ratio), 0.85)
  expect_lt(mean(ratio), 1.15)
  # draws of the first year start from the a of their refit
  fc <- mt_forecast(fit, h = 1, n_draws = nb, seed = 1)
  expect_equal(dim(fc$nb_size), c(10, 1, nb))
  expect_equal(as.vector(fc$nb_size[, 1, ]), rep(1 / f$dispersion, nb))
  kdraw <- log(fc$rates[, 1, 1, ]) - A
  B <- vapply(f$boot, `[[`, numeric(10), "b")
  expect_equal(unname(kdraw / B), matrix(colMeans(kdraw / B), 10, nb, byrow = TRUE), tolerance = 1e-10)
})

test_that("sizes by age follow the strata and draws of reordered forecasts", {
  set.seed(1)
  x <- array(stats::runif(3 * 4 * 5), c(3, 4, 5))
  perm <- t(replicate(4, sample.int(5)))
  out <- mortensor:::permute_rows(x, perm)
  for (s in 1:4) expect_identical(out[, s, ], x[, s, perm[s, ]])
  hs <- list(S = matrix(1, 6, 4), bottom = c(2, 3, 5, 6))
  nr <- mortensor:::node_rows(x, hs)
  expect_identical(nr[, hs$bottom, ], x)
  expect_true(all(is.na(nr[, c(1, 4), ])))
  expect_identical(mortensor:::nb_size_rows(x, c(4, 2)), x[, c(4, 2), , drop = FALSE])
  expect_identical(mortensor:::nb_size_stratum(x, 3), x[, 3, ])
  # the size of each cell [age, period, stratum, draw]
  d <- c(3, 2, 4, 5)
  cells <- array(mortensor:::nb_cell_size(x, d), d)
  for (t in 1:2) expect_identical(cells[, t, , ], x)
  m <- matrix(stats::runif(20), 4)
  cells <- array(mortensor:::nb_cell_size(m, d), d)
  for (a in 1:3) for (t in 1:2) expect_identical(cells[a, t, , ], m)
  # draws of an ensemble take the sizes of their component: by age for the
  # component with a dispersion by age, Inf (Poisson) for the other
  x <- mt_example(deaths = "single", years = 1995:2015)
  x <- mt_subset(x, sex == "male", ages = c(60, 70))
  comps <- list(mt_fit(lee_carter(family = "nb", dispersion = "age"), x), mt_fit(lee_carter(), x))
  alpha <- comps[[1]]$fits[[1]]$dispersion
  nz <- mortensor:::draw_noise(list(components = comps), list(component = rep(1:2, 3)), 6)
  sz <- nz$size(1)
  expect_identical(dim(sz), c(10L, 6L))
  expect_equal(sz[, c(1, 3, 5)], matrix(1 / alpha, 10, 3))
  expect_true(all(sz[, c(2, 4, 6)] == Inf))
})
