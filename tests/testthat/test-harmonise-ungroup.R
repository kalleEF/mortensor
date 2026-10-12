noise_free_fixture <- function(omega = 101) {
  ages <- 30:(omega - 1)
  m <- 2e-5 * exp(0.1 * ages) + 3e-4
  e <- 1e5 * exp(-0.03 * (ages - 30)) * (1 + 0.2 * sin(ages / 3))
  d <- m * e
  lo <- c(seq(30, 85, 5), 90)
  hi <- c(seq(35, 90, 5), omega)
  dg <- data.frame(year = 2000, age_lo = lo, age_hi = hi, deaths = tapply(d, findInterval(ages, lo), sum))
  ex <- data.frame(year = 2000, age_lo = ages, age_hi = ages + 1, exposure = e)
  list(m = m, x = suppressWarnings(mt_data(deaths = dg, exposure = ex, age = c("age_lo", "age_hi"))))
}

test_that("rate-based PCLM recovers a smooth schedule without a level shift", {
  f <- noise_free_fixture()
  # fixed, moderate smoothing: no level shift and group totals fitted without rescaling
  u <- mt_ungroup(f$x, lambda = 1, preserve_totals = FALSE)
  est <- mt_array(u, "deaths")[, 1, 1] / mt_array(u, "exposure")[, 1, 1]
  rel <- est / f$m - 1
  expect_lt(max(abs(rel)), 0.01)
  expect_lt(abs(mean(rel)), 0.002)
  dg <- attr(u, "ungroup")
  expect_equal(dg$method, "pclm_rate_1d")
  expect_lt(dg$max_rel_diff, 0.01)
  # smoothing selected by BIC, totals preserved
  u2 <- mt_ungroup(f$x)
  est2 <- mt_array(u2, "deaths")[, 1, 1] / mt_array(u2, "exposure")[, 1, 1]
  # BIC on noise-free data prefers strong smoothing (the deviance is close to 0)
  expect_lt(max(abs(est2 / f$m - 1)), 0.05)
  expect_equal(sum(mt_obs(u2, "deaths")$value), sum(mt_obs(f$x, "deaths")$value), tolerance = 1e-8)
})

test_that("counts-only PCLM without offset matches group totals and warns about the open interval", {
  f <- noise_free_fixture()
  x <- f$x
  x$obs$exposure <- NULL
  expect_warning(u <- mt_ungroup(x), class = "mt_warning_ungroup_counts")
  expect_equal(sum(mt_obs(u, "deaths")$value), sum(mt_obs(x, "deaths")$value), tolerance = 1e-8)
  expect_error(suppressWarnings(mt_ungroup(x, dims = 2)), class = "mt_error_ungroup")
})

test_that("open intervals are expanded to omega and exposures can be ungrouped first", {
  x <- mt_example(years = 2001:2002)
  xr <- suppressWarnings(mt_regroup(x, age = c(30, seq(35, 90, 5))))
  xr <- mt_set_grid(xr, age = single_ages(30, 100))
  expect_warning(u <- mt_ungroup(xr, exposure_method = "monotone"), class = "mt_warning_ungroup_open_counts")
  e_true <- mt_array(x, "exposure")
  e_est <- mt_array(u, "exposure")
  expect_equal(sum(e_est), sum(e_true), tolerance = 1e-8)
  expect_lt(median(abs(log(e_est[1:60, , ] / e_true[1:60, , ]))), 0.05)
  expect_false(anyNA(mt_array(u, "deaths")))
})

test_that("ungrouped rates are close to the simulated truth", {
  x <- mt_example(years = 2001:2003)
  u <- mt_ungroup(x)
  est <- mt_array(u, "deaths") / mt_array(u, "exposure")
  truth <- attr(x, "truth")
  le <- abs(log(est / truth))
  expect_lt(median(le[1:60, , ]), 0.05)
  prov <- unique(mt_obs(u, "deaths")$provenance)
  expect_true("pclm" %in% prov)
})

test_that("2D ungrouping runs and keeps totals", {
  x <- mt_example(years = 2001:2004)
  u <- mt_ungroup(x, dims = 2, lambda = c(100, 10))
  expect_equal(sum(mt_obs(u, "deaths")$value), sum(mt_obs(x, "deaths")$value), tolerance = 1e-8)
  expect_equal(attr(u, "ungroup")$method, rep("pclm_rate_2d", 2))
})

test_that("England & Wales benchmark (HMD) meets the accuracy thresholds", {
  skip_if_not_installed("StMoMo")
  skip_if(!identical(Sys.getenv("MORTENSOR_SLOW_TESTS"), "true"), "slow benchmark")
  res <- benchmark_ungroup_ew(years = 2000:2011)
  expect_lt(res$male_30_89, 0.02)
  expect_lt(res$male_90_100, 0.035)
})
