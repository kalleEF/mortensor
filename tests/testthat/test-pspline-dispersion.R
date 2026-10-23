# p_splines(overdispersion = "age"): a dispersion that varies smoothly over
# age (Djeundje and Currie 2011). The paper's data (CMI assured lives and
# pensioners) are not public; the tests check exact identities, the limiting
# cases, mgcv as a reference implementation of both penalised
# quasi-likelihood fits, and the recovery of known dispersions.

# Overdispersed deaths with mean E m and variance phi[x] E m (negative
# binomial with size mu / (phi - 1)), ages 40-89, years 1971-2010
disp_sim <- function(seed, ages = 40:89, years = 1971:2010, exposure = 20000) {
  xa <- ages + 0.5
  phi <- 1.5 + 2 * exp(-((xa - 65) / 12)^2)
  t <- years - 1990
  eta <- outer(-10 + 0.09 * xa + 2e-4 * (xa - 65)^2, rep(1, length(t))) +
    outer(rep(1, length(xa)), -0.015 * t - 2e-4 * t^2)
  E <- matrix(exposure, length(xa), length(t))
  mu <- E * exp(eta)
  set.seed(seed)
  D <- matrix(stats::rnbinom(length(mu), size = mu / (phi - 1), mu = mu), nrow(mu))
  list(D = D, E = E, xa = xa, years = years, phi = phi, age = mortensor:::mt_axis(ages, ages + 1, kind = "age"))
}

ps_age_fit <- function(s, lambda = NULL, overdispersion = "age", ndx = NULL) {
  ndx <- ndx %||% c(mortensor:::ps_ndx(s$xa), mortensor:::ps_ndx(s$years))
  mortensor:::ps2_model(s$D, s$E, s$xa, s$years, ndx, c(3, 3), c(2, 2), 0.01, FALSE, lambda, "bic", c(1e-4, 1e6),
                        overdispersion)
}

test_that("dispersions 1 give the Poisson fit, a constant dispersion the fit of overdispersion = TRUE", {
  s <- disp_sim(1, ages = 50:79, years = 1986:2010)
  lam <- c(1, 1e3)
  f0 <- ps_age_fit(s, lam, FALSE)
  f1 <- ps_age_fit(s, lam, TRUE)
  dat <- f0$data
  start <- mortensor:::ps_start(dat$Y, dat$E, dat$wt, f0$ma$B, f0$my$B, s$xa, s$years)
  P <- lam[1] * f0$pen$Pa + lam[2] * f0$pen$Py
  ones <- matrix(1, nrow(dat$Y), ncol(dat$Y))
  # the scoring algorithm with dispersions 1 is the Poisson one
  g0 <- mortensor:::ps2_irls(dat$Y, dat$O, dat$wt, f0$ma, f0$my, P, start, phi = ones)
  expect_equal(g0$alpha, f0$alpha, tolerance = 1e-9)
  expect_equal(mortensor:::ps2_ed(g0, P), f0$ed, tolerance = 1e-9)
  # a constant dispersion psi2 in the weights (Djeundje and Currie 2011, eq.
  # 3.13) is psi2 in the penalty (Mort2Dsmooth(overdispersion = TRUE)): the
  # same coefficients, effective dimension and criterion (eq. 3.16)
  psi2 <- f1$psi2
  g1 <- mortensor:::ps2_irls(dat$Y, dat$O, dat$wt, f0$ma, f0$my, P, start, phi = psi2 * ones)
  expect_equal(g1$alpha, f1$alpha, tolerance = 1e-9)
  ed1 <- mortensor:::ps2_ed(g1, P)
  expect_equal(ed1, f1$ed, tolerance = 1e-9)
  obs <- dat$wt > 0
  ic <- sum(mortensor:::poisson_dev_cells(dat$Y[obs], g1$mu[obs]) / psi2) + log(sum(obs)) * ed1
  expect_equal(ic, f1$ic, tolerance = 1e-9)
  # and the same forecast and draws: psi2 (B'WB + psi2 P)^-1 = (B'(W / psi2)B + P)^-1
  fa <- f1
  fa$psi2 <- 1
  fa$phi <- rep(psi2, nrow(dat$Y))
  fx1 <- mortensor:::ps2_extend_fit(f1, 10, c(2, 2))
  fxa <- mortensor:::ps2_extend_fit(fa, 10, c(2, 2))
  expect_equal(fxa$eta, fx1$eta, tolerance = 1e-8)
  d1 <- with_seed_test(3, mortensor:::ps_draws(fx1, f1$ma$B, psi2, 20))
  da <- with_seed_test(3, mortensor:::ps_draws(fxa, f1$ma$B, 1, 20))
  expect_equal(da, d1, tolerance = 1e-8)
})

test_that("the fit solves the estimating equations of both steps (Djeundje and Currie 2011, sections 3-4)", {
  s <- disp_sim(2, ages = 50:79, years = 1986:2010)
  f <- ps_age_fit(s, c(1, 1e3))
  dat <- f$data
  obs <- dat$wt > 0
  P <- f$lambda[1] * f$pen$Pa + f$lambda[2] * f$pen$Py
  X <- kronecker(f$my$B, f$ma$B)
  # rates: X'((y - mu) / phi) = P alpha (penalised eq. 2.9, log link)
  score <- crossprod(X, as.vector((dat$Y - f$mu) / f$phi)) - P %*% f$alpha
  expect_lt(max(abs(score)), 1e-6 * max(abs(P %*% f$alpha)))
  # leverages: the diagonal of the hat matrix X (X'WX + P)^-1 X'W (eq. 3.17),
  # summing to the effective dimension
  W <- as.vector(f$mu / f$phi)
  lev <- mortensor:::ps2_leverage(f, P, f$ma$B, f$my$B, matrix(W, nrow(dat$Y)))
  H <- X %*% solve(crossprod(X, W * X) + P, t(X * W))
  expect_equal(as.vector(lev), diag(H), tolerance = 1e-10)
  expect_equal(sum(lev), f$ed, tolerance = 1e-10)
  # bias-corrected mean Pearson statistics by age (eq. 4.20, 4.21), from the
  # converged fit
  n <- rowSums(obs)
  raw <- rowSums((dat$Y - f$mu)^2 / f$mu) / (n - rowSums(lev))
  expect_equal(f$dispersion_raw, raw, tolerance = 1e-4)
  # dispersions: penalised quasi-likelihood with E[d] = phi, var(d) = tau
  # phi^2 / n, link -1 / phi (eq. 4.22): B'(n (d - phi)) / tau = lambda2 Dp beta
  Dp <- mortensor:::ps_penalties(ncol(f$ma$B), 1, c(2, 2))$Pa
  g <- mortensor:::ps_disp_irls(f$dispersion_raw, n, f$dispersion_tau, f$ma$B, f$lambda_dispersion * Dp,
                                rep(-0.5, ncol(f$ma$B)))
  expect_equal(g$phi, f$phi, tolerance = 1e-6)
  sc <- crossprod(f$ma$B, n * (f$dispersion_raw - g$phi)) / f$dispersion_tau - f$lambda_dispersion * Dp %*% g$beta
  expect_lt(max(abs(sc)), 1e-6 * max(n))
  # tau: the mean Pearson statistic of the smooth (package-specific)
  tau <- sum(n * (f$dispersion_raw - g$phi)^2 / g$phi^2) / (length(n) - g$ed)
  expect_equal(f$dispersion_tau, tau, tolerance = 1e-2)
  # the dispersion of an age is the same in every year and divides the
  # weights of the forecast fit (the forecast years have weight 0)
  fx <- mortensor:::ps2_extend_fit(f, 5, c(2, 2))
  Xx <- kronecker(fx$my$B, f$ma$B)
  Wx <- cbind(dat$wt, matrix(0, nrow(dat$Y), 5)) * fx$mu / f$phi
  expect_equal(fx$BWB, crossprod(Xx, as.vector(Wx) * Xx), tolerance = 1e-8)
})

test_that("both fits agree with mgcv (penalised quasi-Poisson and gamma fits)", {
  skip_if_not_installed("mgcv")
  s <- disp_sim(3, ages = 55:74, years = 1991:2010)
  f <- ps_age_fit(s, c(1, 1e3))
  dat <- f$data
  X <- kronecker(f$my$B, f$ma$B)
  # rates: quasi-Poisson with prior weights 1 / phi and the smoothing
  # parameters fixed; mgcv minimises sum(w dev) + sum(sp b'Sb)
  dd <- data.frame(y = as.vector(dat$Y), off = as.vector(dat$O), w = rep(1 / f$phi, ncol(dat$Y)))
  dd$X <- X
  g <- mgcv::gam(y ~ X - 1 + offset(off), family = stats::quasipoisson(), weights = w, data = dd,
                 paraPen = list(X = list(f$pen$Pa, f$pen$Py, sp = f$lambda)))
  expect_equal(unname(stats::coef(g)), f$alpha, tolerance = 1e-8)
  expect_equal(sum(g$edf), f$ed, tolerance = 1e-8)
  lev <- mortensor:::ps2_leverage(f, f$lambda[1] * f$pen$Pa + f$lambda[2] * f$pen$Py, f$ma$B, f$my$B,
                                  dat$wt * f$mu / f$phi)
  expect_equal(g$hat, as.vector(lev), tolerance = 1e-8)
  # dispersions: gamma with the inverse link (1 / phi, the negative of the
  # paper's -1 / phi), prior weights n and sp = tau lambda2
  n <- rowSums(dat$wt > 0)
  Dp <- mortensor:::ps_penalties(ncol(f$ma$B), 1, c(2, 2))$Pa
  sm <- mortensor:::ps_disp_irls(f$dispersion_raw, n, f$dispersion_tau, f$ma$B, f$lambda_dispersion * Dp,
                                 rep(-0.5, ncol(f$ma$B)))
  d2 <- data.frame(d = f$dispersion_raw, n = n)
  d2$B <- f$ma$B
  g2 <- mgcv::gam(d ~ B - 1, family = stats::Gamma(link = "inverse"), weights = n, data = d2,
                  paraPen = list(B = list(Dp, sp = f$dispersion_tau * f$lambda_dispersion)))
  expect_equal(unname(stats::fitted(g2)), sm$phi, tolerance = 1e-8)
  expect_equal(sum(g2$edf), sm$ed, tolerance = 1e-8)
})

test_that("one age gives the bias-corrected Pearson estimate and recovers the dispersion of duplicates", {
  # Djeundje and Currie (2011, section 5.1.2): males aged 75, 1950-2005,
  # 1000 policies a year held by 200, 150, 100 and 50 lives with 1, 2, 3
  # and 4 policies, so that var(C) = 2500 m and phi = 2.5; log m quadratic
  # (read off Figure 1); 13 cubic B-splines over years
  years <- 1950:2005
  t <- years - 1950
  m <- exp(-2.71 - 0.0014122 * t - 0.00027515 * t^2)
  set.seed(4)
  res <- t(vapply(1:40, function(i) {
    C <- vapply(m, function(mt) sum((1:4) * stats::rpois(4, c(200, 150, 100, 50) * mt)), 0)
    s <- list(D = matrix(C, 1), E = matrix(1000, 1, length(years)), xa = 75.5, years = years)
    f0 <- ps_age_fit(s, ndx = c(1, 10), overdispersion = FALSE)
    f <- ps_age_fit(s, ndx = c(1, 10))
    # one class: the dispersion is the bias-corrected mean Pearson
    # statistic of the fit (eq. 4.20 with v = ed)
    pearson <- sum((s$D - f$mu)^2 / f$mu) / (length(years) - f$ed)
    c(phi = f$phi, pearson = pearson, ed0 = f0$ed, ed = f$ed)
  }, numeric(4)))
  expect_equal(ncol(mortensor:::ps_basis(years, 10)$B), 13)
  expect_equal(res[, "phi"], res[, "pearson"], tolerance = 1e-4)
  expect_lt(abs(mean(res[, "phi"]) - 2.5), 0.25)
  # the dispersion lowers the effective dimension (section 5.1.2)
  expect_lt(mean(res[, "ed"]), mean(res[, "ed0"]))
})

test_that("a dispersion that varies smoothly over age is recovered", {
  s <- disp_sim(5)
  f <- ps_age_fit(s)
  rel <- f$phi / s$phi - 1
  expect_lt(mean(abs(rel)), 0.1)
  expect_lt(max(abs(rel)), 0.25)
  # the bias-corrected statistics are about unbiased: over three data sets,
  # their mean relative to the true dispersions is about 1
  ratio <- vapply(6:8, function(seed) {
    g <- ps_age_fit(disp_sim(seed), lambda = f$lambda)
    mean(g$dispersion_raw / s$phi)
  }, 0)
  expect_lt(abs(mean(ratio) - 1), 0.05)
})

test_that("p_splines(overdispersion = \"age\") fits, forecasts and reports the dispersion by age", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  spec <- p_splines(lambda = c(1, 1e3), overdispersion = "age")
  out <- expect_mt_model(spec, x)
  expect_equal(spec$label, "P-splines, dispersion by age")
  f <- out$fit$fits[[1]]
  expect_length(f$phi, 71)
  expect_identical(f$psi2, 1)
  p <- mt_parameters(out$fit)
  d <- p[parameter == "dispersion" & stratum == out$fit$strata$label[1]]
  expect_equal(d$index, out$fit$age$labels)
  expect_equal(d$value, f$phi)
  expect_true(all(c("dispersion_raw", "lambda_dispersion", "dispersion_tau") %in% p$parameter))
  expect_match(mortensor:::model_assumptions(spec), "Djeundje and Currie 2011")
  # draws follow the normal approximation with the covariance (B'WB + P)^-1,
  # W the fitted deaths divided by their dispersion
  fx <- mortensor:::ps2_extend_fit(f, 8, f$pord)
  eta <- with_seed_test(1, log(simulate_model(spec, f, 8, 4000)))
  Ainv <- chol2inv(fx$R)
  for (cell in list(c(10, 1), c(40, 8))) {
    v <- as.vector(kronecker(fx$my$B[length(f$xy) + cell[2], ], f$ma$B[cell[1], ]))
    expect_equal(stats::sd(eta[cell[1], cell[2], ]), sqrt(sum(v * (Ainv %*% v))), tolerance = 0.06)
  }
  # predictive draws add Poisson variation, as with overdispersion = TRUE
  fcp <- mt_forecast(out$fit, h = 2, n_draws = 5, seed = 1, predictive = TRUE)
  expect_true(all(is.finite(fcp$rates)))
})

test_that("the dispersion by age works with age groups and missing cells", {
  x <- mt_regroup(mt_example(deaths = "single", years = 1995:2012), age = c(30, seq(35, 95, 5)))
  d <- x$obs$deaths
  x$obs$deaths <- d[!(d$stratum == 1 & d$age_lo == 50 & d$period_lo %in% c(2000, 2001))]
  fit <- mt_fit(p_splines(lambda = c(1, 1e3), overdispersion = "age"), x)
  f <- fit$fits[[1]]
  expect_equal(sum(f$data$wt == 0), 2)
  expect_length(f$phi, 14)
  expect_true(all(f$phi > 0))
  fc <- mt_forecast(fit, h = 3, n_draws = 10, seed = 1)
  expect_true(all(is.finite(fc$rates)))
})

test_that("overdispersion is checked", {
  expect_error(p_splines(overdispersion = "year"), class = "mt_error_argument")
  expect_error(p_splines(overdispersion = NA), class = "mt_error_argument")
  expect_error(p_splines(overdispersion = c(TRUE, FALSE)), class = "mt_error_argument")
  expect_s3_class(p_splines(overdispersion = "age"), "mt_spec_p_splines")
})
