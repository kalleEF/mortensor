test_that("apc_bayes validates its arguments", {
  skip_if_not_installed("bamp")
  expect_error(apc_bayes(age = "rw2+het"), class = "mt_error_argument")
  expect_equal(mcmc_control("fast")$options$number_of_iterations, 5000)
  expect_equal(mcmc_control("test", iterations = 10)$options$number_of_iterations, 10)
})

test_that("apc_bayes converts data at the death resolution", {
  x <- mt_example(years = 2001:2008)
  cells <- mortensor:::apc_cells(x, 2001:2008)
  expect_equal(dim(cells$D), c(13, 8, 2))
  expect_equal(cells$ppa, 5)
  expect_equal(sum(cells$E), sum(mt_obs(x, "exposure")$value))
  expect_equal(sum(cells$D), sum(mt_obs(x, "deaths")$value))
  u <- mt_ungroup(x)
  md <- mortensor:::model_inputs(u, lee_carter())
  sh <- mortensor:::apc_shape(md, cells$age, 3)
  # shape factors average to 1 within each group (weighted by last-year exposures)
  e <- md$E[, 8, 1]
  expect_equal(as.numeric(tapply(e * sh$r[, 1], sh$group, sum) / tapply(e, sh$group, sum)), rep(1, 13))
  bad <- suppressWarnings(mt_regroup(x, age = c(30, 40, 45, 90)))
  expect_error(mortensor:::apc_cells(bad, 2001:2008), class = "mt_error_model_data")
})

test_that("apc_bayes fits, converges and forecasts through the common interface", {
  skip_if_not_installed("bamp")
  skip_on_cran()
  x <- mt_example(years = 2004:2016)
  x$obs <- lapply(x$obs, function(o) o[o$stratum == 1])
  x$strata <- x$strata[1]
  fit <- mt_fit(apc_bayes(age = "rw1", period = "rw1", cohort = "rw1", mcmc = mcmc_control("test")), x,
                prep = mt_prep(step_ungroup()))
  expect_equal(fit$diagnostics$status, "ok")
  fc <- mt_forecast(fit, h = 3, n_draws = 50, seed = 1)
  expect_equal(dim(fc$rates), c(71, 3, 1, 50))
  expect_true(all(is.finite(fc$rates) & fc$rates > 0))
  truth <- attr(x, "truth")[, "2016", 1]
  med <- apply(fc$rates[, 1, 1, ], 1, median)
  expect_lt(median(abs(log(med / truth))), 0.15)
  e <- apc_effects(fit)
  expect_setequal(unique(e$effect), c("age", "period", "cohort"))
})
