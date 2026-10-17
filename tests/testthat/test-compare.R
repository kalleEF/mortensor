test_that("mt_compare() backtests named models and mt_score_table() puts targets side by side", {
  x <- mt_example(deaths = "single", years = 1990:2010)
  bt <- mt_compare(list(rw = rwd(), lc = lee_carter(), lcs = list(spec = lee_carter(), prep = mt_prep(step_smooth()))),
                   x, origins = 2006, h = 3, n_draws = 30, quantities = "life_expectancy")
  expect_s3_class(bt, "mt_backtest")
  expect_identical(levels(bt$model), c("rw", "lc", "lcs"))
  tab <- mt_score_table(bt, metrics = c("crps", "cov95"))
  expect_named(tab, c("model", "crps", "cov95", "crps_le", "cov95_le"))
  expect_equal(tab$crps[2], signif(mt_score(bt)[model == "lc", crps], 3))
  expect_error(mt_score_table(bt, metrics = "nonsense"), class = "mt_error_argument")
  expect_error(mt_compare(list(a = rwd(), a = rwd()), x, origins = 2006, h = 2), class = "mt_error_argument")
  expect_error(mt_compare(list(1), x), class = "mt_error_argument")
})

test_that("mt_compare_forecasts() tabulates one quantity of several forecasts", {
  x <- mt_example(deaths = "single", years = 1990:2010)
  fc <- list(rw = mt_forecast(mt_fit(rwd(), x), h = 5, n_draws = 50, seed = 1),
             lc = mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 50, seed = 1))
  tab <- mt_compare_forecasts(fc, periods = 2015)
  expect_named(tab, c("label", "period", "rw", "lc"))
  e <- life_expectancy(fc$lc)[period == 2015 & label == "female"]
  expect_match(tab[label == "female", lc], sprintf("^%.1f \\(", e$q0.5))
  obs <- mt_compare_forecasts(list(observed = x), periods = 2010, digits = 2)
  expect_equal(nrow(obs), 2)
  d <- mt_compare_forecasts(fc, quantity = "deaths", periods = 2015, wide = FALSE)
  expect_true(all(d$lo <= d$med & d$med <= d$hi))
  r <- mt_compare_forecasts(fc, quantity = "rate", ages = "70", periods = 2015, digits = 4)
  expect_true("age" %in% names(r))
  expect_error(mt_compare_forecasts(list(obs = x), quantity = "deaths"), class = "mt_error_argument")
})

test_that("comparison and design plots are ggplot objects", {
  skip_if_not_installed("ggplot2")
  x <- mt_example(deaths = "single", years = 1990:2010)
  fc <- list(rw = mt_forecast(mt_fit(rwd(), x), h = 5, n_draws = 30, seed = 1),
             lc = mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 30, seed = 1))
  expect_s3_class(mt_plot_compare(fc, data = x), "ggplot")
  expect_s3_class(mt_plot_compare(fc, quantity = "rate", ages = "70", periods = c(2005, 2015)), "ggplot")
  expect_s3_class(mt_plot_backtest_design(x, origins = 2003:2006, h = 4), "ggplot")
  tn <- mt_tune(function(w) rwd(drift_window = w), data.frame(w = c(5, 10)), x, inner_origins = 2, inner_h = 2,
                n_draws = 20)
  expect_s3_class(mt_plot_tuning(tn), "ggplot")
})
