test_that("mt_tune selects parameters with an inner backtest", {
  x <- mt_example(deaths = "single", years = 1990:2010)
  tn <- mt_tune(function(drift_window) rwd(drift_window = drift_window), data.frame(drift_window = c(5, 15)), x,
                inner_origins = 2, inner_h = 3, n_draws = 50)
  expect_equal(nrow(tn$scores), 2)
  expect_s3_class(tn$best, "mt_spec_rwd")
  expect_equal(tn$best$params$drift_window, tn$scores$drift_window[which.min(tn$scores$crps)])
  expect_equal(tn$origins, c(2006, 2007))
  expect_error(mt_tune(1, data.frame(a = 1), x), class = "mt_error_argument")
})

test_that("mt_tuned tunes inside each fit using only the training data", {
  x <- mt_example(deaths = "single", years = 1990:2010)
  sp <- mt_tuned(function(drift_window) rwd(drift_window = drift_window), data.frame(drift_window = c(5, 15)),
                 inner_origins = 2, inner_h = 3, n_draws = 30)
  fit <- mt_fit(sp, mt_window(x, to = 2005))
  expect_equal(nrow(fit$tuning), 2)
  # changing data after 2005 does not change the tuning at origin 2005
  y <- x
  y$obs$deaths <- data.table::copy(y$obs$deaths)
  y$obs$deaths[period_lo > 2005, value := value * 2]
  fit2 <- mt_fit(sp, mt_window(y, to = 2005))
  expect_equal(fit$tuning, fit2$tuning)
  bt <- mt_backtest(sp, x, origins = 2005, h = 3, n_draws = 30)
  expect_s3_class(bt, "mt_backtest")
})
