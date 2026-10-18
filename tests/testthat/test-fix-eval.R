# Regression tests for fixes in backtests, comparisons, derived quantities
# and plots

test_that("predictive deaths are floored after summing them to observed groups", {
  # grouped deaths, ungrouped inside the fit: few deaths at young ages
  x <- mt_example(deaths = "grouped", years = 1995:2019, size = 150)
  bt <- suppressWarnings(mt_backtest(lee_carter(), x, origins = 2014, h = 5, n_draws = 100,
                                     prep = mt_prep(step_ungroup())))
  young <- bt[target == "rate" & age_lo < 45]
  expect_setequal(unique(young$age_hi - young$age_lo), 5)
  # a floor of 0.5 deaths per single year of age put every observed group
  # below all draws (mean PIT and coverage 0)
  expect_gt(mean(young$pit), 0.15)
  expect_gt(mean(young$cov95), 0.7)
})

test_that("derived quantities use the same (unfloored) deaths for draws and observations", {
  E <- c(2000, 1500, 800, 300)
  o <- data.table::data.table(period_lo = 2020, age_lo = c(60, 70, 80, 90), age_hi = c(70, 80, 90, Inf),
                              value = c(0, 3, 0, 0))
  # every draw equals the observed deaths (zero cells included)
  D <- matrix(o$value, 4, 50)
  q <- mortensor:::quantity_scores(D, o, E, c("life_expectancy", "asr"), mortensor:::wis_probs(), c(60, Inf))
  expect_equal(q$target, c("life_expectancy", "asr"))
  expect_true(all(is.finite(q$observed)))
  expect_equal(q$scores$ae, c(0, 0))
  expect_equal(q$scores$crps, c(0, 0))
  ax <- mt_axis(o$age_lo, o$age_hi, kind = "age")
  expect_equal(q$observed[1], life_expectancy(o$value / E, ax))
  # no deaths at any age: a finite (floored) log ASR, the same for draws and data
  q0 <- mortensor:::quantity_scores(D * 0, data.table::copy(o)[, value := 0], E, "asr", mortensor:::wis_probs(),
                                    c(60, Inf))
  expect_true(is.finite(q0$observed))
  expect_equal(q0$draws[[1]], rep(q0$observed, 50))
})

test_that("derived quantities need the full age range in a year", {
  E <- rep(1000, 4)
  o <- data.table::data.table(period_lo = rep(2020:2021, each = 4), age_lo = rep(c(60, 70, 80, 90), 2),
                              age_hi = rep(c(70, 80, 90, Inf), 2), value = 10)
  D <- matrix(10, 8, 20)
  probs <- mortensor:::wis_probs()
  full <- mortensor:::quantity_scores(D, o, rep(E, 2), "life_expectancy", probs, c(60, Inf))
  expect_equal(full$period, 2020:2021)
  # the youngest cell of 2020 and the open cell of 2021 are missing (suppressed or estimated)
  keep <- -c(1, 8)
  part <- mortensor:::quantity_scores(D[keep, ], o[keep], rep(E, 2)[keep], "life_expectancy", probs, c(60, Inf))
  expect_null(part)
  # gaps between strata only for the same first age
  qd <- list(`1` = list(period = 2020, target = "life_expectancy", at = 60, observed = 20, draws = list(rnorm(20))),
             `2` = list(period = 2020, target = "life_expectancy", at = 70, observed = 15, draws = list(rnorm(20))))
  st <- data.table::data.table(stratum = 1:2, sex = c("female", "male"), label = c("female", "male"))
  expect_length(mortensor:::contrast_quantity_scores(qd, st, c(var = "sex", ref = "female"), probs), 0)
  qd$`2`$at <- 60
  expect_length(mortensor:::contrast_quantity_scores(qd, st, c(var = "sex", ref = "female"), probs), 1)
})

test_that("backtests skip years whose youngest cell is estimated, also in contrasts", {
  x <- mt_example_regions(years = 2005:2019, regions = 1)
  x$obs$deaths <- data.table::copy(x$obs$deaths)
  first <- min(x$obs$deaths$age_lo)
  x$obs$deaths[stratum == 1 & age_lo == first & period_lo == 2016, provenance := "estimated"]
  bt <- mt_backtest(rwd(), x, origins = 2014, h = 4, n_draws = 50, quantities = "life_expectancy",
                    contrast = c(var = "quintile", ref = "Q1"))
  le <- bt[target == "life_expectancy"]
  expect_equal(nrow(le[stratum == 1 & period == 2016]), 0)
  expect_equal(nrow(le[stratum == 1]), 3)
  gap <- bt[target == "life_expectancy_gap" & sex == "female"]
  expect_false(2016 %in% gap$period)
  expect_true(all(c(2015, 2017, 2018) %in% gap$period))
})

test_that("mt_inequality() leaves the totals of a hierarchy out of the default order", {
  x <- mt_example_regions(years = 2017:2019, regions = 2)
  xa <- suppressWarnings(mt_aggregate_strata(x, mt_hierarchy(sex = "both", region = "all", quintile = "all")))
  i_def <- mt_inequality(xa, "quintile", quantity = "asr")
  i_ord <- mt_inequality(xa, "quintile", order = paste0("Q", 1:5), quantity = "asr")
  expect_equal(i_def, i_ord)
})

test_that("backtests with step_regroup() regroup the test data to the fit", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  pp <- mt_prep(step_regroup(age = age_groups(seq(30, 100, 5))))
  expect_no_warning(bt <- mt_backtest(rwd(), x, origins = 2008, h = 3, prep = pp, n_draws = 30,
                                      quantities = "life_expectancy"))
  r <- bt[target == "rate"]
  expect_setequal(unique(r$age_hi - r$age_lo), c(5, Inf))
  expect_equal(sort(unique(r$horizon)), 1:3)
  expect_equal(nrow(bt[target == "life_expectancy"]), 3 * 2)
  # a group containing an estimated single year counts as estimated
  test <- mt_window(x, from = 2009)
  test$obs$deaths <- data.table::copy(test$obs$deaths)
  test$obs$deaths[stratum == 1 & age_lo == 32 & period_lo == 2010, provenance := "estimated"]
  ag <- mt_fit(rwd(), mt_window(x, to = 2008), prep = pp)$age
  tg <- mortensor:::test_on_age_axis(test, ag)
  d <- tg$obs$deaths
  expect_identical(d[stratum == 1 & age_lo == 30 & period_lo == 2010]$provenance, "estimated")
  expect_equal(nrow(d[provenance == "estimated"]), 1)
})

test_that("mt_compare() uses the same default origins for every model", {
  x <- mt_example(deaths = "single", years = 2006:2015)
  rwn <- suppressWarnings(rwd(observation_noise = TRUE))
  expect_warning(bt <- mt_compare(list(rw = rwd(), rwn = rwn), x, h = 3, n_draws = 20),
                 class = "mt_warning_backtest_origin")
  o <- bt[, list(origins = paste(sort(unique(origin)), collapse = ",")), by = model]
  expect_equal(o$origins[1], o$origins[2])
  expect_equal(o$origins[1], "2011,2012")
})

test_that("mt_score_table() keeps the preprocessing pipelines apart", {
  x <- mt_example(deaths = "single", years = 1990:2012)
  b1 <- mt_backtest(rwd(), x, origins = c(2005, 2007), h = 3, n_draws = 30, quantities = "life_expectancy")
  b2 <- mt_backtest(rwd(), x, origins = c(2005, 2007), h = 3, n_draws = 30, quantities = "life_expectancy",
                    prep = mt_prep(step_window(10)))
  bt <- rbind(b1, b2)
  tab <- mt_score_table(bt, digits = 6)
  expect_equal(nrow(tab), 2)
  expect_named(tab, c("model", "prep", "crps", "cov95", "crps_le", "cov95_le"))
  s_le <- mt_score(bt, target = "life_expectancy")
  expect_equal(tab$crps_le, signif(s_le$crps[match(tab$prep, s_le$prep)], 6))
  s <- mt_score(bt)
  expect_equal(tab$crps, signif(s$crps[match(tab$prep, s$prep)], 6))
})

test_that("mt_backtest() works without a seed", {
  x <- mt_example(deaths = "single", years = 1995:2010)
  expect_no_warning(bt <- mt_backtest(rwd(), x, origins = c(2005, 2007), h = 2, n_draws = 20, seed = NULL))
  expect_setequal(unique(bt$origin), c(2005, 2007))
})

test_that("mt_backtest() refuses data with ungrouped deaths", {
  x <- mt_example(years = 2001:2012)
  u <- mt_ungroup(x)
  expect_error(mt_backtest(rwd(), u, origins = 2008, h = 3, n_draws = 20), class = "mt_error_ungrouped_data")
  expect_error(mt_compare(list(rw = rwd()), u, origins = 2008, h = 3, n_draws = 20),
               class = "mt_error_ungrouped_data")
  bt <- mt_backtest(rwd(), x, origins = 2008, h = 3, n_draws = 20, prep = mt_prep(step_ungroup()))
  expect_setequal(unique(bt$age_hi - bt$age_lo), c(5, Inf))
})

test_that("mt_compare_forecasts() keeps a column for a forecast without values", {
  x <- mt_example(deaths = "single", years = 2000:2015)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 30, seed = 1)
  tab <- mt_compare_forecasts(list(obs = x, lc = fc), periods = 2020)
  expect_named(tab, c("label", "period", "obs", "lc"))
  expect_true(all(is.na(tab$obs)))
  expect_false(anyNA(tab$lc))
})

test_that("comparison, contrast and age plots handle one year, other quantiles and several levels", {
  skip_if_not_installed("ggplot2")
  x <- mt_example(deaths = "single", years = 2000:2015)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 30, seed = 1)
  fc2 <- mt_forecast(mt_fit(rwd(), x), h = 5, n_draws = 30, seed = 1)
  p <- mt_plot_compare(list(a = fc, b = fc2), quantity = "rate", periods = 2020, strata = "female")
  expect_false(anyNA(p$data$row))
  expect_silent(ggplot2::ggplot_build(p))
  # quantiles other than 2.5% and 97.5%
  cf <- mt_contrast(fc, var = "sex", ref = "male", probs = c(0.1, 0.5, 0.9))
  b <- ggplot2::ggplot_build(mt_plot_contrast(cf))
  expect_equal(sort(b$data[[2]]$ymin), sort(cf$q0.1))
  expect_equal(sort(b$data[[2]]$ymax), sort(cf$q0.9))
  ob <- mt_contrast(x, var = "sex", ref = "male", probs = c(0.1, 0.5, 0.9))
  expect_silent(ggplot2::ggplot_build(mt_plot_contrast(observed = ob)))
  # several interval levels in the age profile
  expect_no_warning(pa <- mt_plot_age(fc, levels = c(0.8, 0.95)))
  ba <- ggplot2::ggplot_build(pa)
  s <- mt_summarise(fc, c(0.025, 0.1, 0.5, 0.9, 0.975))[period == 2020]
  s <- s[order(label, age_lo)]
  med <- ba$data[[3]]
  expect_equal(sort(med$y), sort(log10(s$q0.5)))
  expect_equal(sort(ba$data[[1]]$ymin), sort(log10(s$q0.025)))
  expect_equal(sort(ba$data[[2]]$ymin), sort(log10(s$q0.1)))
  expect_length(ba$plot$labels$subtitle, 1)
})

test_that("a zero rate in the open interval gives a finite life expectancy", {
  m <- 2e-4 + 3e-5 * exp(0.1 * (0:100))
  m[101] <- 0
  lt <- life_table(m)
  expect_true(all(is.finite(lt$ex)))
  m2 <- m
  m2[101] <- m[100]
  expect_equal(lt$ex, life_table(m2)$ex)
  expect_equal(lt$mx[101], 0)
  M <- cbind(m, m2)
  expect_equal(as.vector(life_expectancy(M, single_ages(0, 100))), rep(life_table(m2)$ex[1], 2))
})

test_that("a reconciled model is backtested on bottom-level data like on aggregated data", {
  x <- mt_subset(mt_example_regions(years = 2000:2019, regions = 2), quintile %in% c("Q1", "Q5"), sex == "female")
  h <- mt_hierarchy(region = "national")
  spec <- reconciled(lee_carter(), h, method = "wls")
  b1 <- suppressWarnings(mt_backtest(spec, x, origins = 2014, h = 3, n_draws = 30, seed = 1))
  b2 <- suppressWarnings(mt_backtest(spec, mt_aggregate_strata(x, h), origins = 2014, h = 3, n_draws = 30, seed = 1))
  m <- merge(b1[, .(label = mt_strata(x)$label[stratum], age_lo, period, crps)],
             b2[, .(label = mt_strata(mt_aggregate_strata(x, h))$label[stratum], age_lo, period, crps)],
             by = c("label", "age_lo", "period"))
  expect_equal(nrow(m), nrow(b1))
  expect_equal(m$crps.x, m$crps.y)
})

test_that("mt_tune() refuses ungrouped data with the reason", {
  x <- mt_example(deaths = "grouped", years = 1995:2012)
  u <- suppressWarnings(mt_ungroup(x))
  expect_error(mt_tune(function(w) rwd(drift_window = w), data.frame(w = c(5, 10)), u),
               class = "mt_error_ungrouped_data")
})
