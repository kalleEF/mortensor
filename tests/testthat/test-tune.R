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

test_that("tuned tunes inside each fit using only the training data", {
  x <- mt_example(deaths = "single", years = 1990:2010)
  sp <- tuned(function(drift_window) rwd(drift_window = drift_window), data.frame(drift_window = c(5, 15)),
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

test_that("mt_tune() supports a hold-out, point metrics and scores by horizon", {
  x <- mt_example(deaths = "single", years = 1990:2012)
  tn <- mt_tune(function(drift_window) rwd(drift_window = drift_window), data.frame(drift_window = c(5, 15)), x,
                inner_origins = 1, inner_h = 4, metric = "rmse", n_draws = 30)
  expect_equal(tn$origins, 2008)
  expect_named(tn$scores, c("drift_window", "rmse", "se", "selected"))
  expect_equal(sum(tn$scores$selected), 1)
  expect_equal(sort(unique(tn$by_horizon$horizon)), 1:4)
  expect_equal(tn$chosen$drift_window, tn$scores$drift_window[tn$scores$selected])
  expect_true("rmse" %in% names(mt_score(mt_backtest(rwd(), x, origins = 2008, h = 2, n_draws = 20))))
  expect_error(mt_tune(function(w) rwd(), data.frame(w = 1), x, metric = "mase", target = "life_expectancy"),
               class = "mt_error_argument")
})

test_that("mt_tune() tunes preprocessing together with the model", {
  x <- mt_example(deaths = "single", years = 1990:2012)
  f <- function(lambda) list(spec = lee_carter(), prep = mt_prep(step_smooth(lambda = lambda)))
  tn <- mt_tune(f, data.frame(lambda = c(1, 1e4)), x, inner_origins = 1, inner_h = 3, n_draws = 30,
                target = "life_expectancy")
  expect_true(all(is.finite(tn$scores$crps)))
  expect_s3_class(tn$best_prep, "mt_prep")
  expect_identical(tn$best_prep$steps[[1]]$args$lambda, tn$chosen$lambda)
  expect_error(mt_tune(f, data.frame(lambda = 1), x, prep = mt_prep(step_smooth())), class = "mt_error_argument")
  fit <- mt_fit(tuned(f, data.frame(lambda = c(1, 1e4)), inner_origins = 1, inner_h = 3, n_draws = 30), x)
  expect_s3_class(fit$prep, "mt_prep")
  expect_match(fit$spec$label, "lambda =")
})

test_that("the one-standard-error rule takes the simplest setting near the best", {
  sc <- data.table::data.table(row = 1:3, value = c(0.105, 0.101, 0.100), se = c(0.004, 0.004, 0.004))
  expect_equal(mortensor:::tune_select(sc, "min"), 3)
  expect_equal(mortensor:::tune_select(sc, "one_se"), 2)
  sc$se <- 0.01
  expect_equal(mortensor:::tune_select(sc, "one_se"), 1)
})

test_that("tuning by group uses only the strata of each group", {
  x <- mt_example(deaths = "single", years = 1990:2012)
  f <- function(drift_window) rwd(drift_window = drift_window)
  g <- data.frame(drift_window = c(5, 15))
  tn <- mt_tune(f, g, x, inner_origins = 2, inner_h = 3, n_draws = 30, by = "sex")
  expect_named(tn$best, c("female", "male"))
  alone <- mt_tune(f, g, mt_subset(x, sex == "female"), inner_origins = 2, inner_h = 3, n_draws = 30)
  expect_equal(tn$scores[group == "female", crps], alone$scores$crps)
  sp <- tuned(f, g, inner_origins = 2, inner_h = 3, n_draws = 30, by = "sex")
  fit <- mt_fit(sp, x)
  expect_s3_class(fit, "mt_fit_split")
  expect_identical(dimnames(fitted(fit))[[3]], mt_strata(x)$label)
  fc <- mt_forecast(fit, h = 3, n_draws = 20, seed = 1, predictive = TRUE)
  expect_identical(dim(fc$rates)[3], 2L)
  expect_true(all(is.finite(fc$rates)))
  expect_identical(mt_forecast(fit, h = 3, n_draws = 20, seed = 1)$rates,
                   mt_forecast(fit, h = 3, n_draws = 20, seed = 1)$rates)
  expect_output(print(fit), "female: ")
  expect_s3_class(mt_backtest(sp, x, origins = 2008, h = 2, n_draws = 20), "mt_backtest")
})

test_that("wrappers around tuned() forecast with the chosen model", {
  x <- mt_example(deaths = "single", years = 1985:2012)
  tu <- tuned(function(drift_window) rwd(drift_window = drift_window), data.frame(drift_window = c(5, 10)),
              inner_origins = 2, inner_h = 3, n_draws = 30)
  fit <- mt_fit(calibrated(anchored(tu), years = 6, n_origins = 2, n_draws = 30), x)
  fc <- mt_forecast(fit, h = 3, n_draws = 20, seed = 1)
  expect_s3_class(fc, "mt_draws")
  expect_true(all(is.finite(fc$rates)))
  expect_s3_class(mt_parameters(fit), "data.table")
})
