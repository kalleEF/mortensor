test_that("stick-breaking shares add up to one", {
  set.seed(1)
  p <- array(runif(2 * 3 * 3 * 4), c(2, 3, 3, 4))
  s <- stick_shares(p)
  expect_equal(dim(s), c(2, 3, 4, 4))
  expect_equal(apply(s, c(1, 2, 4), sum), array(1, c(2, 3, 4)))
  expect_equal(s[, , 1, ], p[, , 1, ])
  expect_equal(s[, , 2, ], (1 - p[, , 1, ]) * p[, , 2, ])
})

test_that("the forecast of coupled period effects has the variance of the multivariate random walk", {
  # one posterior draw with known precisions; the period effects of two
  # fields at horizon h have covariance h (s Omega)^-1
  set.seed(2)
  Y <- array(rbinom(4 * 6 * 2, 1000, 0.1), c(4, 6, 2))
  st <- mc_setup(Y, array(1000, c(4, 6, 2)), 1, c(1L, 1L, 1L))
  Om <- matrix(c(400, -150, -150, 300), 2)
  fit <- list(setup = st, x = matrix(0, 1, st$P), Om = array(Om, c(1, 2, 2)), Oc = array(diag(1e8, 2), c(1, 2, 2)))
  h <- 4
  sims <- replicate(4000, {
    p <- mc_draws(fit, h, 1L)
    stats::qlogis(p[1, st$J + h, , 1])
  })
  expect_equal(stats::cov(t(sims)), h * solve(st$Kp$scale * Om), tolerance = 0.06)
})

test_that("apc_causes() forecasts causes that add up to all causes", {
  skip_on_cran()
  x <- mt_example_causes(years = 2005:2019)
  x <- mt_regroup(x, age = seq(30, 100, 10))
  fit <- mt_fit(apc_causes(mcmc = mcmc_control("test"), seed = 1), x)
  expect_true(all(fit$diagnostics$status == "ok"))
  expect_lt(max(fit$diagnostics$max_rhat), 1.1)
  fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
  r <- fc$rates
  lab <- fit$strata$label
  for (sx in c("female", "male")) {
    causes <- match(paste0(sx, ".", c("cvd", "cancer", "other")), lab)
    expect_equal(apply(r[, , causes, ], c(1, 2, 4), sum), r[, , match(paste0(sx, ".all"), lab), ])
  }
  f <- fitted(fit)
  expect_false(anyNA(f))
  causes <- match(paste0("female.", c("cvd", "cancer", "other")), lab)
  expect_equal(apply(f[, , causes], c(1, 2), sum), f[, , match("female.all", lab)], tolerance = 0.02)
  p <- mt_parameters(fit)
  expect_true(all(c("period_correlation", "cohort_correlation", "rhat", "period") %in% p$parameter))
})

test_that("apc_causes() checks its arguments and data", {
  expect_error(apc_causes(var = c("a", "b")), class = "mt_error_argument")
  expect_error(apc_causes(residual = "yes"), class = "mt_error_argument")
  expect_error(apc_causes(period = "rw3"), class = "mt_error_argument")
  x <- mt_regroup(mt_example_causes(years = 2010:2019), age = seq(30, 100, 10))
  expect_error(mt_fit(apc_causes(var = "sex", mcmc = mcmc_control("test")), x), class = "mt_error_model_data")
  expect_error(mt_fit(apc_causes(var = "none", mcmc = mcmc_control("test")), x), class = "mt_error_argument")
})
