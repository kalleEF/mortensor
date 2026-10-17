test_that("mt_smooth reduces noise and keeps exposures", {
  x <- mt_example(years = 2001:2004, deaths = "single", size = 2e4)
  s <- mt_smooth(x)
  truth <- attr(x, "truth")
  raw <- mt_array(x, "deaths") / mt_array(x, "exposure")
  sm <- mt_array(s, "deaths") / mt_array(s, "exposure")
  expect_lt(median(abs(log(sm / truth))), median(abs(log(pmax(raw, 1e-6) / truth))))
  expect_equal(mt_array(s, "exposure"), mt_array(x, "exposure"))
  expect_equal(unique(mt_obs(s, "deaths")$provenance), "smoothed")
  expect_error(mt_smooth(mt_example(years = 2001:2002)), class = "mt_error_resolution")
})

test_that("step_shocks interpolates inside the data and excludes at the edges", {
  x <- mt_example(years = 2001:2008, deaths = "single")
  y <- mt_prep_apply(mt_prep(step_shocks(2005)), x)
  m <- mt_array(y, "deaths") / mt_array(y, "exposure")
  m0 <- log(pmax(mt_array(x, "deaths"), 0.5) / mt_array(x, "exposure"))
  expect_equal(log(m[, "2005", ]), (m0[, "2004", ] + m0[, "2006", ]) / 2, tolerance = 1e-10)
  expect_equal(sort(unique(mt_obs(y)$provenance)), c("interpolated", "observed"))
  z <- mt_prep_apply(mt_prep(step_shocks(2007:2008, action = "exclude")), x)
  expect_equal(max(mt_grid(z)$period$lower), 2006)
  expect_error(mt_prep_apply(mt_prep(step_shocks(2005, action = "exclude")), x), class = "mt_error_prep")
  # interpolation has no neighbour after the last year: the shock year is removed
  expect_equal(max(mt_grid(mt_prep_apply(mt_prep(step_shocks(2008)), x))$period$lower), 2007)
  # the input object is not modified by reference
  expect_equal(unique(mt_obs(x)$provenance), "observed")
})

test_that("preprocessing pipelines validate and print their steps", {
  expect_error(mt_prep(1), class = "mt_error_prep")
  p <- mt_prep(step_regroup(age = age_groups(seq(30, 90, 10))), step_grid())
  expect_output(print(p), "step_regroup")
  x <- mt_prep_apply(p, mt_example(years = 2001:2002))
  expect_equal(mt_grid(x)$age$labels[1], "30-39")
  expect_identical(mt_prep_apply(NULL, x), x)
})

test_that("steps fitted inside a backtest never see data after the origin", {
  x <- mt_example(years = 1995:2010)
  y <- x
  y$obs$deaths <- data.table::copy(y$obs$deaths)
  y$obs$deaths[period_lo > 2005, value := 0]
  prep <- mt_prep(step_ungroup(dims = 2, lambda = c(100, 10)), step_smooth())
  f1 <- mt_fit(rwd(), mt_window(x, to = 2005), prep = prep)
  f2 <- mt_fit(rwd(), mt_window(y, to = 2005), prep = prep)
  expect_equal(f1$observed, f2$observed)
})
