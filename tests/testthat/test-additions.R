# Smoothing over time, raw residuals, jump-off, fixed reconciliation,
# quantity scores, Lee-Carter variants, tuning, disaggregation

test_that("smoothing over time and the age-period GAM work and are checked", {
  skip_if_not_installed("mgcv")
  x <- mt_example(deaths = "single", years = 2000:2012)
  st <- mt_smooth(x, along = "period", method = "whittaker")
  expect_equal(attr(st, "smooth")$along[1], "period")
  sg <- mt_smooth(x, dims = 2, method = "gam")
  expect_true(all(attr(sg, "smooth")$dispersion > 0))
  ck <- mt_smooth_check(x, st)
  expect_true(all(c("mad_log", "dispersion", "trend_ratio", "correlation") %in% names(ck)))
  expect_true(all(ck$trend_ratio < 1))
  expect_error(mt_smooth(x, dims = 2, method = "whittaker"), class = "mt_error_argument")
  expect_error(mt_smooth(x, along = "period"), class = "mt_error_argument")
  skip_if_not_installed("ggplot2")
  expect_s3_class(mt_plot_smooth(x, st), "ggplot")
})

test_that("models fitted to smoothed rates take observation noise from the raw rates", {
  x <- mt_example(deaths = "single", years = 2000:2015)
  # variation beyond Poisson noise in every cell
  set.seed(3)
  x$obs$deaths[, value := round(value * exp(stats::rnorm(.N, 0, 0.1)))]
  raw <- mt_fit(lee_carter(), x)
  sm <- mt_fit(lee_carter(), x, prep = mt_prep(step_smooth(along = "period", method = "whittaker")))
  expect_false(is.null(sm$raw_D))
  expect_equal(sm$observed, raw$observed)
  expect_gt(stats::median(raw$noise_sd), 0.01)
  expect_gt(stats::median(sm$noise_sd), 0.5 * stats::median(raw$noise_sd))
  on_smoothed <- mortensor:::extra_noise_sd(sm$spec, sm$fits, sm$train_D, sm$train_E)
  expect_gt(stats::median(sm$noise_sd), stats::median(on_smoothed))
  # residuals refer to the raw deaths
  expect_equal(sum(stats::residuals(sm)$deaths), sum(sm$raw_D))
})

test_that("anchored() moves the jump-off and keeps wrappers coherent", {
  x <- mt_example(deaths = "single", years = 1995:2015)
  a1 <- mt_forecast(mt_fit(anchored(lee_carter(), years = 1, decay = Inf, shrink = FALSE), x), h = 3,
                    n_draws = 40, seed = 1)
  lo <- mt_forecast(mt_fit(lee_carter(jump_off = "observed"), x), h = 3, n_draws = 40, seed = 1)
  expect_equal(apply(a1$rates, 1:3, stats::median), apply(lo$rates, 1:3, stats::median), tolerance = 1e-8)
  a0 <- mt_forecast(mt_fit(anchored(lee_carter(), decay = 0), x), h = 3, n_draws = 40, seed = 1)
  b0 <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 40, seed = 1)
  expect_equal(a0$rates, b0$rates, ignore_attr = TRUE)
  fe <- mt_fit(anchored(ensemble(lc = lee_carter(), rw = rwd())), x)
  expect_equal(fe$anchor$components$rw, matrix(0, length(fe$age), nrow(fe$strata)))
  expect_s3_class(mt_forecast(fe, h = 2, n_draws = 20, seed = 1, predictive = TRUE), "mt_draws")
  jo <- mt_jump_off(a1)
  expect_true(all(c("observed", "fitted", "forecast", "step_observed", "step_fitted") %in% names(jo)))
  expect_equal(nrow(mt_jump_off(a1, quantity = "rate")), length(a1$age) * nrow(a1$strata))
  xr <- mt_example_regions(years = 2000:2014, regions = 2)
  h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
  fr <- mt_forecast(mt_fit(anchored(reconciled(lee_carter(), h, method = "wls")), mt_aggregate_strata(xr, h)),
                    h = 2, n_draws = 20, seed = 1)
  expect_lt(coherence_gap(fr, h), 1e-8)
})

test_that("fixed nodes are kept by reconciliation", {
  x <- mt_example_regions(years = 2005:2019, regions = 2)
  h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
  xa <- mt_aggregate_strata(x, h)
  fc <- mt_forecast(mt_fit(lee_carter(), xa), h = 3, n_draws = 30, seed = 1)
  rw <- mt_reconcile(fc, h, method = "wls", fixed = "total", dependence = "none")
  tot <- which(fc$strata$label == "both.all.all")
  expect_lt(max(abs(rw$rates[, , tot, ] / fc$rates[, , tot, ] - 1)), 1e-6)
  expect_lt(rw$reconciliation$max_relative_change_fixed, 1e-6)
  co <- mt_coherence(rw, h)
  expect_lt(max(abs(c(co$q0.025, co$q0.975) - 1)), 1e-8)
  expect_gt(max(abs(mt_coherence(fc, h)$q0.5 - 1)), 1e-4)
  expect_error(mt_reconcile(fc, h, method = "bottom_up", fixed = "total"), class = "mt_error_argument")
  expect_error(mt_reconcile(fc, h, fixed = "nonsense"), class = "mt_error_argument")
  fit <- mt_fit(reconciled(lee_carter(), h, method = "mint", fixed = "total"), xa)
  pr <- mt_forecast(fit, h = 2, n_draws = 20, seed = 1, predictive = TRUE)
  expect_lt(coherence_gap(pr, h), 1e-8)
  # a total far below the sum of its parts: bottoms are held at zero, the total is kept
  low <- fc
  low$rates[, , tot, ] <- 0.3 * low$rates[, , tot, ]
  for (m in c("wls", "mint", "ols")) {
    rl <- mt_reconcile(low, h, method = m, fixed = "total", dependence = "none")
    expect_lt(rl$reconciliation$max_relative_change_fixed, 1e-6)
    expect_gt(rl$reconciliation$share_set_to_zero, 0)
    # rates are floored at 1e-12, which matters only for aggregates of zeros
    expect_lt(coherence_gap(rl, h), 1e-4)
    expect_true(all(rl$rates >= 0))
  }
  skip_if_not_installed("ggplot2")
  expect_s3_class(mt_plot_coherence(fc, h), "ggplot")
})

test_that("backtests score derived quantities", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  bt <- mt_backtest(lee_carter(), x, origins = 2007, h = 3, n_draws = 50, quantities = c("life_expectancy", "asr"))
  expect_setequal(unique(bt$target), c("rate", "life_expectancy", "asr"))
  s_rate <- mt_score(bt)
  s_e0 <- mt_score(bt, target = "life_expectancy")
  expect_equal(s_e0$n, 3 * 2)
  expect_true(is.na(s_e0$mase))
  expect_false(isTRUE(all.equal(s_rate$crps, s_e0$crps)))
  expect_error(mt_backtest(lee_carter(), x, origins = 2007, h = 3, quantities = "e0"), class = "mt_error_argument")
  plain <- mt_backtest(rwd(), x, origins = 2007, h = 3, n_draws = 20)
  expect_error(mt_score(plain, target = "asr"), class = "mt_error_argument")
  skip_if_not_installed("ggplot2")
  expect_s3_class(mt_plot_scores(rbind(bt, plain, fill = TRUE), by = "model"), "ggplot")
})

test_that("Lee-Carter adjustments and the choice of the fitting period follow their definitions", {
  x <- mt_example(deaths = "single", years = 1980:2015)
  f <- mt_fit(lee_carter(method = "svd", adjust = "e0"), x)
  e_fit <- life_expectancy(stats::fitted(f), f$age)
  e_obs <- life_expectancy(f$observed, f$age)
  expect_equal(e_fit, e_obs, tolerance = 1e-6, ignore_attr = TRUE)
  fd <- mt_fit(lee_carter(method = "svd", adjust = "dxt"), x)
  s1 <- fd$fits[[1]]
  D <- fd$train_D[, , 1]
  mu <- fd$train_E[, , 1] * s1$fitted
  expect_lt(max(abs(colSums((D - mu) * s1$b)) / colSums(mu * s1$b^2)), 1e-6)
  fc <- mt_fit(lee_carter(choose_period = TRUE, min_years = 10), x)
  expect_true(fc$fits[[1]]$first_year >= 1980 && fc$fits[[1]]$first_year <= 2005)
  expect_true(anyNA(stats::fitted(fc)))
  expect_s3_class(mt_forecast(fc, h = 3, n_draws = 10, seed = 1), "mt_draws")
})

test_that("ensembles keep the best components and tuning works in two stages", {
  x <- mt_example(deaths = "single", years = 1990:2012)
  e <- mt_fit(ensemble(lc = lee_carter(), rw = rwd(), tu = tucker3(), top_k = 2, inner_origins = 2, inner_h = 3), x)
  expect_equal(sum(e$weights > 0), 2)
  g <- data.frame(drift_window = c(3, 5, 8, 10, 15, 20))
  tn <- mt_tune(function(drift_window) rwd(drift_window = drift_window), g, x, inner_origins = 2, inner_h = 3,
                strategy = "two_stage")
  expect_lt(nrow(tn$scores), nrow(g))
  expect_s3_class(tn$best, "mt_spec")
})

test_that("disaggregation options, shapes and checks", {
  x <- mt_example(years = 2000:2019, deaths = "single")
  xg <- mt_regroup(x, age = c(seq(30, 85, 5), 90))
  fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 3, n_draws = 20, seed = 1)
  e <- mt_array(x, "exposure")[, 20, ]
  sh <- mt_shape_rates(x, periods = 2020:2022)
  expect_equal(dim(sh), c(length(mt_grid(x)$age), 3, 2))
  fs <- mt_disaggregate(fc, mt_grid(x)$age, exposure = e, method = "shape", shape_rates = sh)
  ck <- mt_disaggregate_check(fc, fs, e)
  expect_lt(max(ck$max_group_error), 1e-8)
  expect_error(mt_disaggregate(fc, mt_grid(x)$age, exposure = e, method = "shape", shape_rates = sh[1:3, , ]),
               class = "mt_error_argument")
  skip_if_not_installed("mgcv")
  u <- mt_ungroup(mt_example(years = 2012:2019), method = "gam", tail = "gompertz", tail_width = 3)
  expect_s3_class(u, "mt_data")
  expect_s3_class(mt_ungroup(mt_example(years = 2019), method = "gam"), "mt_data")
  shg <- mt_shape_rates(x, periods = 2020:2022, method = "gam")
  expect_true(all(is.finite(shg) & shg > 0))
  skip_if_not_installed("ggplot2")
  expect_s3_class(mt_plot_disaggregation(fc, fs), "ggplot")
  expect_s3_class(mt_plot_disaggregation(fc, fs, type = "smoothness"), "ggplot")
})
