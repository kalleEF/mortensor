test_that("plot functions return ggplot objects", {
  skip_if_not_installed("ggplot2")
  x <- mt_example(years = 2001:2012)
  u <- mt_ungroup(x)
  fc <- mt_forecast(mt_fit(lee_carter(), u), h = 3, n_draws = 30, seed = 1)
  bt <- mt_backtest(rwd(), u, origins = 2008, h = 3, n_draws = 30)
  plots <- list(mt_plot_forecast(fc), mt_plot_age(fc), mt_plot_life_expectancy(fc, at = 65),
                mt_plot_rates(u), mt_plot_ungroup(x, u, 2012), mt_plot_coverage(bt))
  for (p in plots) {
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
})
