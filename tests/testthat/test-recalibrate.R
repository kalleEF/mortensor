test_that("the calibration estimator recovers a known spread factor", {
  set.seed(1)
  n <- 200
  G <- 400
  L <- matrix(stats::rnorm(G * n, 0, 0.05), G)
  pieces <- list(list(label = "a", L = L, m = rep(0, G), h = rep(1:4, G / 4), E = rep(1, G),
                      C = Matrix::Diagonal(G), Ecell = rep(1, G), y = stats::rnorm(G, 0, 0.1),
                      noise = rep(FALSE, n), size = rep(Inf, n)))
  P <- with_seed_test(2, mortensor:::stack_pieces(pieces, n))
  r <- mortensor:::calibration_fit(P, "constant", FALSE)
  expect_equal(r$table$scale_h1, 2, tolerance = 0.15)
  expect_lt(r$table$crps_after, r$table$crps_before)
  expect_gt(r$table$cov95_after, 0.9)
  ri <- mortensor:::calibration_fit(P, "constant", FALSE, target = "interval", level = 0.9)
  expect_equal(ri$table$scale_h1, 2, tolerance = 0.2)
  rh <- mortensor:::calibration_fit(P, "horizon", TRUE)
  expect_equal(abs(rh$table$slope) < 0.3, TRUE)
  expect_lt(abs(rh$table$drift), 0.02)
})

test_that("the transformation keeps the median and the ranks of the draws", {
  sh <- mortensor:::calibration_shift(c(log(2), 0), 1:3, "horizon", FALSE)
  expect_equal(sh$s, rep(2, 3))
  shc <- mortensor:::calibration_shift(c(0, 0.5), c(1, 4, 16), "horizon", FALSE, h_max = 4)
  expect_equal(shc$s, c(1, 2, 2))
  sh0 <- mortensor:::calibration_shift(0, 1:3, "constant", FALSE)
  expect_equal(sh0$s, 1)
  expect_equal(sh0$b, 0)
})

test_that("recalibrated() wraps any model and forecasts on the grid", {
  x <- mt_example(deaths = "single", years = 1995:2016)
  sp <- recalibrated(lee_carter(), years = 6, n_origins = 3, n_draws = 30)
  expect_s3_class(sp, "mt_spec_recalibrated")
  expect_equal(sp$min_periods, lee_carter()$min_periods + 6)
  r <- expect_mt_model(sp, data = x, n_draws = 20)
  cal <- r$fit$calibration
  expect_true(all(c("group", "scale_h1", "slope", "drift", "crps_before", "crps_after", "cov95_after") %in%
                    names(cal)))
  expect_lte(cal$crps_after, cal$crps_before + 1e-12)
  expect_equal(r$fit$calibration_origins, c(2010, 2012, 2014))
  expect_identical(r$fc$model, "Lee-Carter, recalibrated")
  pr <- mt_forecast(r$fit, h = 3, n_draws = 20, seed = 1, predictive = TRUE)
  expect_true(all(is.finite(pr$rates) & pr$rates >= 0))
  # summary and parameters of the base model
  expect_output(print(summary(r$fit)), "calibration \\(all\\): spread factor")
  pa <- mt_parameters(r$fit)
  expect_true(all(c("a", "b", "k", "scale_h1") %in% pa$parameter))
  # separate parameters by sex, pooled ones as well
  fs <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 2, by = "sex", n_draws = 30), x)
  expect_setequal(fs$calibration$group, c("all", "female", "male"))
  expect_identical(fs$spec$label, "Lee-Carter, recalibrated by sex")
  expect_identical(recalibrated(rwd(), target = "interval", level = 0.9)$label,
                   "Random walk with drift, recalibrated (90% interval)")
  expect_error(recalibrated(lee_carter(), years = 1), class = "mt_error_argument")
  expect_error(mt_fit(recalibrated(lee_carter(), years = 6, by = "region", n_draws = 30), x),
               class = "mt_error_argument")
  expect_warning(mt_fit(recalibrated(lee_carter(), years = 20, n_draws = 20), x), class = "mt_warning_calibration")
  short <- mt_subset(x, periods = c(1995, 1999))
  expect_error(mt_fit(recalibrated(lee_carter(), years = 3), short), class = "mt_error_model_data")
})

test_that("recalibrated ensembles, reconciled models and backtests work", {
  x <- mt_example(deaths = "single", years = 1995:2014)
  ens <- recalibrated(ensemble(lc = lee_carter(), rw = rwd()), years = 5, n_origins = 2, n_draws = 30)
  r <- expect_mt_model(ens, data = x, n_draws = 20)
  expect_equal(nrow(r$fit$calibration), 1)
  bt <- mt_backtest(recalibrated(rwd(), years = 5, n_origins = 2, n_draws = 30), x, origins = 2011, h = 3,
                    n_draws = 30)
  expect_equal(unique(bt$model), "Random walk with drift, recalibrated")
  xr <- mt_example_regions(years = 2000:2014, regions = 2)
  h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
  fit <- mt_fit(recalibrated(reconciled(lee_carter(), h, method = "wls"), years = 5, n_origins = 2, n_draws = 30),
                mt_aggregate_strata(xr, h))
  fc <- mt_forecast(fit, h = 2, n_draws = 20, seed = 1)
  expect_lt(coherence_gap(fc, h), 1e-8)
  pr <- mt_forecast(fit, h = 2, n_draws = 20, seed = 1, predictive = TRUE)
  expect_lt(coherence_gap(pr, h), 1e-8)
})

test_that("recalibrated() estimates factors by age band and on life expectancy", {
  x <- mt_example(deaths = "single", years = 1985:2015)
  fb <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 3, n_draws = 30, ages = c(60, 80)), x)
  expect_identical(fb$calibration$band, c("<60", "60-80", "80+"))
  fa <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 3, n_draws = 30), x)
  expect_equal(sum(fb$calibration$cells), fa$calibration$cells)
  fc <- mt_forecast(fb, h = 12, n_draws = 40, seed = 1)
  expect_true(all(is.finite(fc$rates)))
  # each band gets its own factor in the forecast
  base <- mt_forecast(mt_fit(lee_carter(), x), h = 12, n_draws = 40, seed = 1)
  ratio <- function(a) {
    w <- function(r) stats::IQR(log(r))
    w(fc$rates[a, 1, 1, ]) / w(base$rates[a, 1, 1, ])
  }
  expect_equal(ratio(10), fb$calibration$scale_h1[1], tolerance = 0.05)
  expect_equal(ratio(60), fb$calibration$scale_h1[3], tolerance = 0.05)
  expect_error(recalibrated(rwd(), ages = c(80, 60)), class = "mt_error_argument")
  expect_error(recalibrated(rwd(), ages = 60, quantity = "asr"), class = "mt_error_argument")
  fe <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 3, n_draws = 30, quantity = "life_expectancy"), x)
  # one life expectancy per stratum and test year of each origin (6 + 4 + 2 years)
  expect_equal(fe$calibration$scored, 2 * 12)
  expect_match(fe$spec$label, "on life expectancy")
})

test_that("recalibrated() reports intervals of the factors and damps the drift", {
  x <- mt_example(deaths = "single", years = 1985:2015)
  fit <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 3, n_draws = 30, n_boot = 5), x)
  cal <- fit$calibration
  expect_true(all(c("scale_h1_lo", "scale_h1_hi", "slope_lo", "slope_hi") %in% names(cal)))
  expect_lte(cal$scale_h1_lo, cal$scale_h1_hi)
  expect_output(print(summary(fit)), "spread factor .* \\(")
  sh <- mortensor:::calibration_shift(c(0, 0.02), 1:20, "constant", TRUE, h_max = 8)
  expect_equal(sh$b[8:20], rep(0.16, 13))
  expect_warning(mt_fit(recalibrated(rwd(), years = 3, n_origins = 1, n_draws = 20), x),
                 class = "mt_warning_calibration")
})
