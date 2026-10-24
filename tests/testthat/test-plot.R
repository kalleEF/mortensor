test_that("plot functions return ggplot objects", {
  skip_if_not_installed("ggplot2")
  x <- mt_example(years = 2001:2012)
  u <- mt_ungroup(x)
  fc <- mt_forecast(mt_fit(lee_carter(), u), h = 3, n_draws = 30, seed = 1)
  bt <- mt_backtest(rwd(), x, origins = 2008, h = 3, n_draws = 30, prep = mt_prep(step_ungroup()))
  plots <- list(mt_plot_forecast(fc), mt_plot_age(fc), mt_plot_life_expectancy(fc, at = 65),
                mt_plot_rates(u), mt_plot_ungroup(x, u, 2012), mt_plot_coverage(bt))
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
})

test_that("Tucker, share and score plots build", {
  skip_if_not_installed("ggplot2")
  x <- mt_example(deaths = "single", years = 2005:2019)
  p <- mt_plot_tucker(mt_fit(tucker3(ranks = c(2, 2, 2)), x))
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
  expect_error(mt_plot_tucker(mt_fit(rwd(), x)), class = "mt_error_argument")
  bt <- rbind(mt_backtest(rwd(), x, origins = 2014, h = 3, n_draws = 20),
              mt_backtest(lee_carter(), x, origins = 2014, h = 3, n_draws = 20))
  for (by in c("horizon", "age")) expect_silent(ggplot2::ggplot_build(mt_plot_scores(bt, "wis", by)))
  expect_error(mt_plot_scores(bt, "nope"), class = "mt_error_argument")
  xc <- mt_example_causes()
  hc <- mt_hierarchy(cause = list(total = "all", additive = "rates"))
  fc <- mt_forecast(mt_fit(reconciled(lee_carter(), hc, method = "bottom_up"), xc), h = 3, n_draws = 20, seed = 1)
  p <- mt_plot_shares(fc)
  expect_silent(ggplot2::ggplot_build(p))
  dd <- apply(mt_deaths(fc, noise = FALSE), c(2, 3, 4), sum)
  st <- fc$strata
  parts <- which(st$cause != "all" & st$sex == "female")
  i <- parts[st$cause[parts] == "cvd"]
  manual <- apply(dd[, i, ] / apply(dd[, parts, ], c(1, 3), sum), 1, stats::median)
  expect_equal(p$data[kind == "forecast" & group == "female" & cause == "cvd"]$med, unname(manual))
  expect_error(mt_plot_shares(fc, var = "region"), class = "mt_error_argument")
})

test_that("derived tables of a backtest print as data tables", {
  x <- mt_example(deaths = "single", years = 2005:2019)
  bt <- mt_backtest(rwd(), x, origins = 2014, h = 2, n_draws = 20)
  tab <- bt[, .N, by = horizon]
  expect_output(print(tab), "horizon")
  expect_s3_class(bt, "mt_backtest")
})
