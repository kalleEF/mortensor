# Simulated data with a known shock: deaths at ages `from` and above
# multiplied by `factor` in `years`; attribute "added" holds the injected
# deaths by sex and year
shocked_example <- function(years = 2000:2019, shock = 2016:2017, factor = 1.1, from = 70, seed = 1,
                            deaths = "single") {
  x <- mt_example(years = years, deaths = deaths, seed = seed)
  d <- x$obs$deaths
  hit <- d$period_lo %in% shock & d$age_lo >= from
  add <- round(factor * d$value[hit]) - d$value[hit]
  x$obs$deaths[hit, value := value + add]
  inj <- stats::aggregate(list(added = add), list(sex = x$strata$sex[d$stratum[hit]], period = d$period_lo[hit]),
                          sum)
  attr(x, "added") <- data.table::as.data.table(inj)
  x
}

test_that("mt_excess sums draws and computes excess and P-scores path by path", {
  x <- shocked_example()
  ex <- mt_excess(lee_carter(), x, baseline = c(2000, 2014), ages = c(30, 70), n_draws = 201, seed = 1)
  expect_equal(nrow(ex), 2 * 2 * 5)
  expect_true(all(c("sex", "age", "age_lo", "age_hi", "period", "observed", "expected_mean", "expected_q0.025",
                    "expected_q0.5", "expected_q0.975", "excess_mean", "excess_q0.5", "p_score_q0.025",
                    "p_score_q0.975", "p_positive", "cells", "skipped") %in% names(ex)))
  expect_equal(unique(ex$age), c("30-69", "70+"))
  expect_equal(ex$excess_mean, ex$observed - ex$expected_mean)
  # with 201 draws the quantiles are order statistics, so the P-score
  # quantiles are those of the expected deaths
  expect_equal(ex$p_score_q0.025, 100 * (ex$observed / ex$expected_q0.975 - 1))
  expect_equal(ex$excess_q0.5, ex$observed - ex$expected_q0.5)
  # observed deaths are the sums of the data
  d <- mt_obs(x, "deaths")
  expect_equal(ex[sex == "male" & age == "70+" & period == 2016]$observed,
               sum(d[sex == "male" & age_lo >= 70 & period_lo == 2016]$value))
  # totals: the same draws summed over ages and strata
  tot <- mt_excess(lee_carter(), x, baseline = c(2000, 2014), by = character(), n_draws = 201, seed = 1)
  expect_equal(tot$observed, ex[, sum(observed), by = "period"]$V1)
  expect_equal(tot$expected_mean, ex[, sum(expected_mean), by = "period"]$V1)
  expect_equal(tot$age, rep("30+", 5))
  expect_equal(tot$cells, ex[, sum(cells), by = "period"]$V1)
  # the shock is found, and the injected excess is inside the interval
  added <- attr(x, "added")[, list(added = sum(added)), by = "period"]
  sh <- merge(tot, added, by = "period")
  expect_true(all(sh$excess_q0.025 <= sh$added & sh$added <= sh$excess_q0.975))
  expect_true(all(sh$p_positive > 0.95))
  expect_true(all(tot[period %in% c(2015, 2018, 2019)]$excess_q0.025 < 0))
  # cumulative sums over the years, draw by draw
  cu <- mt_excess(lee_carter(), x, baseline = c(2000, 2014), by = character(), cumulative = TRUE, n_draws = 201,
                  seed = 1)
  expect_equal(cu$observed, cumsum(tot$observed))
  expect_equal(cu$expected_mean, cumsum(tot$expected_mean))
  expect_equal(cu$from, rep(2015, 5))
  expect_true(all(diff(cu$excess_q0.975 - cu$excess_q0.025) > 0))
})

test_that("mt_excess takes a fit or a specification and is reproducible", {
  x <- shocked_example(years = 2005:2019)
  a <- mt_excess(rwd(), x, baseline = c(2005, 2015), periods = 2017, n_draws = 50, seed = 3)
  fit <- mt_fit(rwd(), mt_subset(x, periods = c(2005, 2015)))
  b <- mt_excess(fit, x, periods = 2017, n_draws = 50, seed = 3)
  expect_equal(a, b, ignore_attr = TRUE)
  expect_equal(attr(a, "excess")$baseline, c(2005, 2015))
  expect_equal(a$period, c(2017, 2017))
  expect_error(mt_excess(fit, x, baseline = c(2005, 2014)), class = "mt_error_argument")
  expect_error(mt_excess(rwd(), x), class = "mt_error_argument")
  expect_error(mt_excess(fit, x, periods = 2015), class = "mt_error_argument")
  expect_error(mt_excess(fit, x, periods = 2025), class = "mt_error_argument")
  expect_error(mt_excess(fit, x, by = "region"), class = "mt_error_argument")
  expect_error(mt_excess(fit, x, ages = c(30, 72)), NA)
  expect_error(mt_excess(fit, x, quantity = "life_expectancy", ages = c(30, 70)), class = "mt_error_argument")
  expect_error(mt_excess(fit, x, probs = c(0, 0.5)), class = "mt_error_argument")
})

test_that("mt_excess compares at the resolution of the observed deaths", {
  single <- mt_example(years = 2000:2019, deaths = "single")
  grouped <- mt_example(years = 2000:2019, deaths = "grouped")
  fit <- mt_fit(lee_carter(), mt_subset(single, periods = c(2000, 2015)))
  es <- mt_excess(fit, single, ages = c(30, 70), n_draws = 100, seed = 1)
  eg <- mt_excess(fit, grouped, ages = c(30, 70), n_draws = 100, seed = 1)
  expect_equal(eg$observed, es$observed)
  expect_equal(eg$expected_mean, es$expected_mean)
  ob <- mt_excess(fit, grouped, periods = 2016, ages = "observed", n_draws = 100, seed = 1)
  expect_equal(unique(ob$age), c(paste0(seq(30, 85, 5), "-", seq(34, 89, 5)), "90+"))
  expect_true(all(ob$cells == 1))
  # observed groups that cross the limits of `ages`
  expect_error(mt_excess(fit, grouped, ages = c(30, 72), n_draws = 10), class = "mt_error_argument")
})

test_that("mt_excess leaves out estimated cells", {
  x <- mt_example(years = 2005:2019, deaths = "single")
  fit <- mt_fit(rwd(), mt_subset(x, periods = c(2005, 2015)))
  x$obs$deaths[stratum == 1 & age_lo == 50 & period_lo == 2017, provenance := "estimated"]
  ex <- mt_excess(fit, x, n_draws = 20, seed = 1)
  expect_equal(ex[sex == "female" & period == 2017]$skipped, 1L)
  expect_equal(ex[sex == "female" & period == 2017]$cells, ex[sex == "female" & period == 2016]$cells - 1L)
  expect_equal(sum(ex$skipped), 1L)
  inc <- mt_excess(fit, x, n_draws = 20, seed = 1, include_estimated = TRUE)
  expect_equal(sum(inc$skipped), 0L)
  expect_gt(inc[sex == "female" & period == 2017]$observed, ex[sex == "female" & period == 2017]$observed)
  # life expectancy needs every age: the year is skipped with a warning
  expect_warning(le <- mt_excess(fit, x, quantity = "life_expectancy", n_draws = 20, seed = 1),
                 class = "mt_warning_excess")
  expect_equal(nrow(le), 2 * 4 - 1)
})

test_that("mt_excess leaves out the totals of a hierarchy when it sums strata", {
  x <- mt_example_regions(years = 2005:2019, regions = 2)
  xa <- mt_aggregate_strata(x, mt_hierarchy(sex = "both"))
  ex <- mt_excess(rwd(), xa, baseline = c(2005, 2015), by = character(), n_draws = 20, seed = 1)
  d <- mt_obs(x, "deaths")
  expect_equal(ex$observed, d[period_lo >= 2016, sum(value), by = "period_lo"]$V1)
  bs <- mt_excess(rwd(), xa, baseline = c(2005, 2015), by = "sex", periods = 2016, n_draws = 20, seed = 1)
  expect_equal(bs[sex == "both"]$observed, sum(bs[sex != "both"]$observed))
})

test_that("mt_excess of life expectancy", {
  x <- shocked_example(years = 2005:2019, shock = 2017)
  ex <- mt_excess(lee_carter(), x, baseline = c(2005, 2015), quantity = "life_expectancy", at = 65, n_draws = 100,
                  seed = 1)
  expect_equal(nrow(ex), 2 * 4)
  expect_equal(ex$at, rep(65, 8))
  expect_false("p_score_mean" %in% names(ex))
  m <- mt_array(x, "deaths") / mt_array(x, "exposure")
  e_obs <- life_expectancy(m, mt_grid(x)$age, at = 65)
  expect_equal(ex[sex == "male" & period == 2017]$observed, unname(e_obs["2017", "male"]))
  expect_equal(ex$excess_mean, ex$observed - ex$expected_mean)
  # the shock lowers life expectancy
  expect_true(all(ex[period == 2017]$excess_q0.975 < 0))
  # both sexes together: deaths and exposures summed
  both <- mt_excess(lee_carter(), x, baseline = c(2005, 2015), quantity = "life_expectancy", by = character(),
                    n_draws = 100, seed = 1)
  D <- apply(mt_array(x, "deaths"), c(1, 2), sum)
  E <- apply(mt_array(x, "exposure"), c(1, 2), sum)
  expect_equal(both[period == 2016]$observed, unname(life_expectancy(D / E, mt_grid(x)$age)["2016"]))
})

test_that("mt_plot_excess returns ggplot objects", {
  skip_if_not_installed("ggplot2")
  x <- shocked_example(years = 2005:2019)
  ex <- mt_excess(rwd(), x, baseline = c(2005, 2014), ages = c(30, 70), n_draws = 50, seed = 1)
  for (type in c("expected", "excess", "p_score")) {
    p <- mt_plot_excess(ex, type = type)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
  one <- mt_excess(rwd(), x, baseline = c(2005, 2014), periods = 2016, ages = seq(30, 90, 10), n_draws = 50, seed = 1)
  expect_s3_class(mt_plot_excess(one), "ggplot")
  expect_s3_class(mt_plot_excess(one, type = "excess", colour = "sex"), "ggplot")
  le <- mt_excess(rwd(), x, baseline = c(2005, 2014), quantity = "life_expectancy", by = character(), n_draws = 50,
                  seed = 1)
  expect_s3_class(mt_plot_excess(le), "ggplot")
  expect_error(mt_plot_excess(le, type = "p_score"), class = "mt_error_argument")
  expect_error(mt_plot_excess(ex, level = 0.8), class = "mt_error_argument")
  expect_error(mt_plot_excess(data.frame(a = 1)), class = "mt_error_argument")
})

test_that("excess intervals have nominal coverage on simulated shocks", {
  skip_if(!identical(Sys.getenv("MORTENSOR_SLOW_TESTS"), "true"), "slow benchmark")
  res <- data.table::rbindlist(lapply(1:20, function(r) {
    x <- shocked_example(years = 1990:2019, seed = r)
    ex <- mt_excess(lee_carter(), x, baseline = c(1990, 2014), n_draws = 300, seed = r)
    ex <- merge(ex, attr(x, "added"), by = c("sex", "period"), all.x = TRUE)
    ex[is.na(added), added := 0]
    ex[, list(added, covered = excess_q0.025 <= added & added <= excess_q0.975, med = excess_q0.5)]
  }))
  expect_gt(mean(res[added == 0]$covered), 0.85)
  expect_gt(mean(res[added > 0]$covered), 0.85)
  expect_lt(abs(mean(res[added > 0, (med - added) / added])), 0.05)
})
