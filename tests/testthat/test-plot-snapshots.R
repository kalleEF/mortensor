# Visual regression tests (vdiffr). Inputs are small and seeded; Tucker
# factor signs are fixed in hooi(), so the plots are reproducible across
# platforms. The snapshots were made with ggplot2 4.0; older versions draw
# differently and skip. Update snapshots with testthat::snapshot_review()
# after an intended change or a new ggplot2 release.

# Deaths by cause with a fixed split of the deaths of mt_example(): the
# expected counts of two causes rounded down, the rest to "other".
# mt_example_causes() splits by rmultinom(), whose binomial draws differed
# on macOS (arm64) for some cells and moved the shares in the SVG
causes_fixed <- function(years) {
  x <- mt_example(deaths = "single", years = years)
  D <- mt_obs(x, "deaths")
  E <- mt_obs(x, "exposure")
  t0 <- (D$period_lo - min(years)) / 10
  a <- (pmin(D$age_lo, 100) - 30) / 10
  p <- cbind(exp(-1.5 + 0.35 * a - 0.25 * t0), exp(0.3 - 0.15 * a + 0.05 * t0), 1)
  p <- p / rowSums(p)
  cvd <- floor(D$value * p[, 1])
  cancer <- floor(D$value * p[, 2])
  split <- list(cvd = cvd, cancer = cancer, other = D$value - cvd - cancer)
  base <- D[, list(sex, period_lo, age_lo, age_hi)]
  long <- rbind(cbind(base, cause = "all", deaths = D$value),
                rbindlist(lapply(names(split), function(k) cbind(base, cause = k, deaths = split[[k]]))))
  ex <- E[, list(sex, period_lo, age_lo, age_hi, exposure = value)]
  ex <- rbindlist(lapply(c("all", names(split)), function(k) cbind(ex, cause = k)))
  mt_data(deaths = long, exposure = ex, by = c("sex", "cause"), age = c("age_lo", "age_hi"), period = "period_lo")
}

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
  skip_if_not_installed("ggplot2", minimum_version = "4.0.0")
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
  xc <- causes_fixed(2008:2019)
  hc <- mt_hierarchy(cause = list(total = "all", additive = "rates"))
  fcc <- mt_forecast(mt_fit(reconciled(rwd(), hc, method = "bottom_up"), xc), h = 3, n_draws = 50, seed = 1)
  vdiffr::expect_doppelganger("cause shares", mt_plot_shares(fcc))
})
