# Negative binomial Lee-Carter model (Delwarde, Denuit and Partrat 2007)

# Deaths from the negative binomial Lee-Carter model with dispersion phi
# (variance mu + phi mu^2); sum(b) = 1, sum(k) = 0
sim_lc_nb <- function(seed, nx = 20, ny = 30, exposure = 1e6, phi = 0.005) {
  set.seed(seed)
  a <- -6 + 0.08 * (seq_len(nx) - 1)
  b <- 1 + 0.5 * sin(seq_len(nx) / 4)
  b <- b / sum(b)
  k <- seq(12, -12, length.out = ny)
  k <- k - mean(k)
  E <- matrix(exposure, nx, ny)
  mu <- E * exp(a + outer(b, k))
  D <- matrix(if (phi > 0) stats::rnbinom(nx * ny, size = 1 / phi, mu = mu) else stats::rpois(nx * ny, mu), nx)
  list(D = D, E = E, mu = mu, a = a, b = b, k = k, phi = phi)
}

test_that("the negative binomial log-likelihood is that of eq. (5)", {
  set.seed(1)
  mu <- stats::runif(300, 0.2, 2000)
  D <- stats::rnbinom(300, size = 30, mu = mu)
  for (phi in c(0.5, 0.02, 1e-4, 1e-7)) {
    expect_equal(mortensor:::nb_loglik(D, mu, phi), sum(stats::dnbinom(D, size = 1 / phi, mu = mu, log = TRUE)),
                 tolerance = 1e-10)
  }
  expect_equal(mortensor:::nb_loglik(D, mu, 0), sum(stats::dpois(D, mu, log = TRUE)), tolerance = 1e-12)
})

test_that("the dispersion is the maximum likelihood estimate given the means", {
  set.seed(2)
  mu <- stats::runif(400, 5, 3000)
  D <- stats::rnbinom(400, size = 1 / 0.01, mu = mu)
  ll <- function(lp) sum(stats::dnbinom(D, size = exp(-lp), mu = mu, log = TRUE))
  ref <- exp(stats::optimize(ll, c(-12, 2), maximum = TRUE, tol = 1e-12)$maximum)
  expect_equal(mortensor:::nb_dispersion(D, mu), ref, tolerance = 1e-6)
  # no overdispersion: the estimate is 0, the Poisson model
  expect_identical(mortensor:::nb_dispersion(round(mu), mu), 0)
})

test_that("the fit solves the likelihood equations of eq. (5)", {
  x <- sim_lc_nb(1)
  f <- mortensor:::lc_nb(x$D, x$E)
  # log-likelihood from dnbinom(), differentiated numerically
  ll <- function(a, b, k, phi) sum(stats::dnbinom(x$D, size = 1 / phi, mu = x$E * exp(a + outer(b, k)), log = TRUE))
  grad <- function(v, set) {
    vapply(seq_along(v), function(i) {
      e <- replace(numeric(length(v)), i, 1e-6)
      args <- list(f$a, f$b, f$k, f$dispersion)
      lo <- hi <- args
      hi[[set]] <- v + e
      lo[[set]] <- v - e
      (do.call(ll, hi) - do.call(ll, lo)) / 2e-6
    }, 0)
  }
  # gradients relative to the sums of the absolute terms of the score
  mu <- x$E * exp(f$a + outer(f$b, f$k))
  scale <- sum(abs(x$D - mu))
  expect_lt(max(abs(grad(f$a, 1))) / scale, 1e-6)
  expect_lt(max(abs(grad(f$b, 2))) / scale, 1e-6)
  expect_lt(max(abs(grad(f$k, 3))) / scale, 1e-6)
  dphi <- (ll(f$a, f$b, f$k, f$dispersion * (1 + 1e-5)) - ll(f$a, f$b, f$k, f$dispersion * (1 - 1e-5))) / 2e-5
  expect_lt(abs(dphi), 1e-3)
  # identification constraints
  expect_equal(sum(f$b), 1, tolerance = 1e-12)
  expect_equal(sum(f$k), 0, tolerance = 1e-9)
})

test_that("negative binomial Lee-Carter recovers the parameters of simulated data", {
  # a, b and k within sampling error: standard errors from the expected
  # information of the negative binomial likelihood (weights mu / (1 + phi
  # mu) of the log link) under the constraints sum(b) = 1 and sum(k) = 0
  x <- sim_lc_nb(1)
  f <- mortensor:::lc_nb(x$D, x$E)
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
  # the dispersion, over 20 independent data sets: maximum likelihood
  # estimates of a variance are biased downwards by about the share of
  # estimated parameters among the cells (68 of 600 here), so the mean
  # estimate lies within 15% below the truth
  reps <- 20
  truth <- sim_lc_nb(0)
  fits <- lapply(seq_len(reps), function(i) {
    x <- sim_lc_nb(100 + i)
    mortensor:::lc_nb(x$D, x$E)
  })
  phi <- vapply(fits, `[[`, 0, "dispersion")
  expect_gt(mean(phi), 0.85 * truth$phi)
  expect_lt(mean(phi), truth$phi + 3 * stats::sd(phi) / sqrt(reps))
  # given the true expected deaths, the dispersion estimate has no such bias
  phi_true_mu <- vapply(seq_len(reps), function(i) {
    x <- sim_lc_nb(100 + i)
    mortensor:::nb_dispersion(x$D, x$mu)
  }, 0)
  expect_lt(abs(mean(phi_true_mu) - truth$phi), 4 * stats::sd(phi_true_mu) / sqrt(reps))
})

test_that("as the dispersion vanishes the fit is the Poisson Lee-Carter fit", {
  x <- sim_lc_nb(3)
  ref <- mortensor:::lc_poisson(x$D, x$E, tol = 1e-15, maxit = 5000)
  gap <- vapply(10^c(-3, -5, -7, -9), function(phi) {
    f <- mortensor:::lc_nb(x$D, x$E, disp = phi)
    max(abs(f$b - ref$b), abs(f$a - ref$a), abs(f$k - ref$k) / max(abs(ref$k)))
  }, 0)
  expect_true(all(diff(gap) < 0))
  expect_lt(gap[4], 1e-6)
  # data without overdispersion: the estimated dispersion is 0 and the fit
  # is the Poisson maximum likelihood fit
  D0 <- round(x$mu)
  f0 <- mortensor:::lc_nb(D0, x$E)
  p0 <- mortensor:::lc_poisson(D0, x$E, tol = 1e-15, maxit = 5000)
  expect_identical(f0$dispersion, 0)
  expect_identical(f0$nb_size, Inf)
  expect_equal(f0$a, p0$a, tolerance = 1e-7)
  expect_equal(f0$b, p0$b, tolerance = 1e-7)
  expect_equal(f0$k, p0$k, tolerance = 1e-7)
})

test_that("lee_carter(family = \"nb\") fits through mt_fit() with the constraints and the dispersion", {
  x <- mt_example(deaths = "single", years = 1995:2015)
  spec <- lee_carter(family = "nb")
  expect_identical(spec$params$family, "nb")
  expect_match(spec$label, "negative binomial")
  fit <- mt_fit(spec, x)
  expect_true(all(fit$diagnostics$status == "ok"))
  p <- mt_parameters(fit)
  for (s in unique(p$stratum)) {
    expect_equal(sum(p[stratum == s & parameter == "b"]$value), 1, tolerance = 1e-10)
    expect_equal(sum(p[stratum == s & parameter == "k"]$value), 0, tolerance = 1e-8)
    disp <- p[stratum == s & parameter == "nb_dispersion"]$value
    expect_length(disp, 1)
    expect_gte(disp, 0)
  }
  expect_match(summary(fit)$assumptions, "negative binomial")
  # the Poisson default is unchanged and has no dispersion
  expect_identical(lee_carter()$params$family, "poisson")
  expect_false("nb_dispersion" %in% mt_parameters(mt_fit(lee_carter(), x))$parameter)
  expect_error(lee_carter(method = "svd", family = "nb"), class = "mt_error_argument")
  expect_error(lee_carter(choose_period = TRUE, family = "nb"), class = "mt_error_argument")
})

test_that("predictive draws are negative binomial with the estimated dispersion", {
  x <- sim_lc_nb(4, nx = 10, ny = 25, exposure = 2e5, phi = 0.01)
  ages <- 60:69
  long <- data.frame(age = rep(ages, 25), year = rep(1991:2015, each = 10), deaths = as.vector(x$D),
                     exposure = as.vector(x$E))
  d <- mt_data(deaths = long, exposure = long)
  fit <- mt_fit(lee_carter(family = "nb"), d)
  expect_identical(dim(fit$train_D), c(10L, 25L, 1L))
  phi <- fit$fits[[1]]$dispersion
  expect_gt(phi, 0)
  fc <- mt_forecast(fit, h = 1, n_draws = 4000, seed = 1, predictive = TRUE)
  expect_equal(as.vector(fc$nb_size), rep(1 / phi, 4000))
  rates <- mt_forecast(fit, h = 1, n_draws = 4000, seed = 1)$rates
  # variance of the predictive deaths: variance of the expected deaths plus
  # the negative binomial variance mu + phi mu^2
  E <- x$E[, 25]
  Dp <- fc$rates[, 1, 1, ] * E
  mu <- rates[, 1, 1, ] * E
  expected <- apply(mu, 1, stats::var) + rowMeans(mu + phi * mu^2)
  expect_lt(max(abs(apply(Dp, 1, stats::var) / expected - 1)), 0.12)
  # the negative binomial variance is several times the Poisson variance
  expect_gt(min(rowMeans(mu + phi * mu^2) / rowMeans(mu)), 2)
})

test_that("the bootstrap resamples negative binomial deaths and keeps the dispersion of each refit", {
  x <- mt_example(deaths = "single", years = 1995:2015)
  x <- mt_subset(x, sex == "male", ages = c(60, 80))
  fit <- with_seed_test(1, mt_fit(lee_carter(family = "nb", bootstrap = TRUE, n_boot = 6), x))
  f <- fit$fits[[1]]
  expect_length(f$boot, 6)
  sizes <- vapply(f$boot, `[[`, 0, "nb_size")
  fc <- mt_forecast(fit, h = 2, n_draws = 12, seed = 1)
  expect_equal(as.vector(fc$nb_size), rep_len(sizes, 12))
})

test_that("the paper's estimates for France are reproduced from the HMD data", {
  skip_if_not_installed("demography")
  # Delwarde, Denuit and Partrat (2007) use French data of INED and INSEE,
  # ages 0-99, 1950-2002; demography::fr.mort holds the Human Mortality
  # Database series for France (deaths = rates x exposures). Section 3.3:
  # dispersion 0.004737 (males) and 0.002343 (females); Section 3.5: drift
  # of k (theta) for the Poisson and negative binomial models; the column
  # labelled sigma^2 there holds the standard deviation of the innovations
  # (maximum likelihood, divisor T - 1); Section 3.4: AIC.
  fr <- mortensor:::pkg_data("fr.mort", "demography")
  ages <- 0:99
  paper <- list(male = c(phi = 0.004737, theta_p = -1.456, sd_p = 1.928, theta_nb = -1.510, sd_nb = 2.017,
                         aic_p = 93647.96, aic_nb = 65259.80),
                female = c(phi = 0.002343, theta_p = -1.968, sd_p = 2.851, theta_nb = -2.049, sd_nb = 2.575,
                           aic_p = 71312.47, aic_nb = 60496.70))
  for (sex in names(paper)) {
    yrs <- 1950:2002
    D <- fr$rate[[sex]][as.character(ages), as.character(yrs)] * fr$pop[[sex]][as.character(ages), as.character(yrs)]
    E <- fr$pop[[sex]][as.character(ages), as.character(yrs)]
    nb <- mortensor:::lc_nb(D, E)
    po <- mortensor:::lc_poisson(D, E, tol = 1e-12, maxit = 5000)
    ml_sd <- function(f) f$sigma_k * sqrt((length(yrs) - 2) / (length(yrs) - 1))
    np <- 2 * length(ages) + length(yrs) - 2
    aic_p <- -2 * mortensor:::nb_loglik(D, E * exp(po$a + outer(po$b, po$k)), 0) + 2 * np
    aic_nb <- -2 * nb$loglik + 2 * (np + 1)
    ref <- paper[[sex]]
    expect_equal(nb$dispersion, ref[["phi"]], tolerance = 0.01)
    expect_equal(po$drift, ref[["theta_p"]], tolerance = 0.005)
    expect_equal(nb$drift, ref[["theta_nb"]], tolerance = 0.005)
    expect_equal(ml_sd(po), ref[["sd_p"]], tolerance = 0.01)
    expect_equal(ml_sd(nb), ref[["sd_nb"]], tolerance = 0.01)
    expect_equal(aic_p, ref[["aic_p"]], tolerance = 0.005)
    expect_equal(aic_nb, ref[["aic_nb"]], tolerance = 0.005)
  }
  # Section 3.6, training period 1950-1980: drift and standard deviation
  yrs <- 1950:1980
  D <- fr$rate$male[as.character(ages), as.character(yrs)] * fr$pop$male[as.character(ages), as.character(yrs)]
  E <- fr$pop$male[as.character(ages), as.character(yrs)]
  nb <- mortensor:::lc_nb(D, E)
  expect_equal(nb$drift, -1.01685612, tolerance = 0.005)
  expect_equal(nb$sigma_k * sqrt((length(yrs) - 2) / (length(yrs) - 1)), 1.726567609, tolerance = 0.01)
})
