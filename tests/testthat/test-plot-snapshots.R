# Visual regression tests (vdiffr). Inputs are small and seeded; Tucker
# factor signs are fixed in hooi(), so the plots are reproducible across
# platforms. Update snapshots with testthat::snapshot_review() after an
# intended change.

snapshot_inputs <- function() {
  x <- mt_example(deaths = "single", years = 2006:2019)
  xg <- mt_example(years = 2016:2019)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 100, seed = 1)
  bt <- rbind(mt_backtest(rwd(), x, origins = c(2014, 2016), h = 3, n_draws = 50, seed = 1),
              mt_backtest(lee_carter(), x, origins = c(2014, 2016), h = 3, n_draws = 50, seed = 1))
  list(x = x, xg = xg, u = mt_ungroup(xg), fc = fc, bt = bt)
}

test_that("plots match their snapshots", {
  skip_if_not_installed("vdiffr")
  skip_if_not_installed("ggplot2")
  skip_on_cran()
  s <- snapshot_inputs()
  vdiffr::expect_doppelganger("forecast fan", mt_plot_forecast(s$fc, ages = c("50", "80")))
  vdiffr::expect_doppelganger("age profile", mt_plot_age(s$fc))
  vdiffr::expect_doppelganger("life expectancy", mt_plot_life_expectancy(s$fc, at = 65))
  vdiffr::expect_doppelganger("lexis surface", mt_plot_rates(s$x))
  vdiffr::expect_doppelganger("ungrouping", mt_plot_ungroup(s$xg, s$u, period = 2019))
  vdiffr::expect_doppelganger("coverage", mt_plot_coverage(s$bt))
  vdiffr::expect_doppelganger("scores by horizon", mt_plot_scores(s$bt))
  vdiffr::expect_doppelganger("tucker factors", mt_plot_tucker(mt_fit(tucker3(ranks = c(2, 2, 2)), s$x)))
  xc <- mt_example_causes(years = 2008:2019)
  fcc <- mt_forecast(mt_fit(cause_shares(ranks = c(2, 2, 2, 2)), xc), h = 3, n_draws = 50, seed = 1)
  vdiffr::expect_doppelganger("cause shares", mt_plot_shares(fcc))
})
