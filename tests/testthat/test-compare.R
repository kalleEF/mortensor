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

test_that("mt_plot_smooth() compares smoothers and mt_plot_pit() shows calibration", {
  skip_if_not_installed("ggplot2")
  x <- mt_example(deaths = "single", years = 1995:2012)
  sm <- list(age = mt_smooth(x), time = mt_smooth(x, along = "period", method = "whittaker"))
  expect_s3_class(mt_plot_smooth(x, sm, strata = "female"), "ggplot")
  expect_s3_class(mt_plot_smooth(x, sm, along = "age", periods = c(1995, 2012)), "ggplot")
  expect_s3_class(mt_plot_smooth(x, sm$time, along = "age"), "ggplot")
  expect_error(mt_plot_smooth(x, unname(sm)), class = "mt_error_argument")
  bt <- mt_backtest(rwd(), x, origins = 2008, h = 3, n_draws = 30)
  expect_s3_class(mt_plot_pit(bt), "ggplot")
})

test_that("smoothing over time is the weighted Whittaker smoother with GCV per age", {
  x <- mt_example(deaths = "single", years = 1995:2015)
  s <- mt_smooth(x, along = "period", method = "whittaker")
  dg <- attr(s, "smooth")
  # one smoothing parameter per age and stratum
  expect_identical(nrow(dg), 2L * length(mt_grid(x)$age))
  expect_true(all(c("age", "lambda_period", "edf", "at_bound") %in% names(dg)))
  D <- mt_array(x, "deaths")[, , 1]
  E <- mt_array(x, "exposure")[, , 1]
  S <- mt_array(s, "deaths")[, , 1]
  # each age: (W + lambda P) z = W y with the deaths as weights and lambda by GCV
  P <- crossprod(diff(diag(21), differences = 2))
  gcv <- function(l, y, w) {
    A <- diag(w) + l * P
    z <- solve(A, w * y)
    21 * sum(w * (y - z)^2) / (21 - sum(diag(solve(A, diag(w)))))^2
  }
  for (a in c(1, 30, 60)) {
    y <- log(pmax(D[a, ], 0.5) / E[a, ])
    w <- pmax(D[a, ], 0.5)
    l <- dg[stratum == 1]$lambda_period[a]
    expect_equal(as.vector(log(S[a, ] / E[a, ])), as.vector(solve(diag(w) + l * P, w * y)), tolerance = 1e-8)
    grid <- 10^seq(-2, 8, 0.25)
    expect_lte(gcv(l, y, w), min(vapply(grid, gcv, 0, y = y, w = w)) * (1 + 1e-4))
  }
  # lambda_range is honoured in both directions
  sr <- suppressWarnings(mt_smooth(x, along = "period", lambda_range = c(10, 100)), classes = "mt_warning_lambda_bound")
  expect_true(all(attr(sr, "smooth")$lambda_period >= 10 - 1e-6 & attr(sr, "smooth")$lambda_period <= 100 + 1e-6))
  sa <- suppressWarnings(mt_smooth(mt_window(x, from = 2013), method = "whittaker", lambda_range = c(10, 100)),
                         classes = "mt_warning_lambda_bound")
  expect_true(all(attr(sa, "smooth")$lambda_age >= 10 - 1e-6 & attr(sa, "smooth")$lambda_age <= 100 + 1e-6))
})

test_that("mt_plot_contrast() draws gaps of one or several forecasts", {
  skip_if_not_installed("ggplot2")
  x <- mt_example(deaths = "single", years = 1995:2012)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 50, seed = 1)
  cf <- mt_contrast(fc, var = "sex", ref = "male", at = 30)
  expect_s3_class(mt_plot_contrast(cf, observed = mt_contrast(x, var = "sex", ref = "male", at = 30)), "ggplot")
  expect_s3_class(mt_plot_contrast(list(a = cf, b = cf)), "ggplot")
  expect_error(mt_plot_contrast(list(cf, cf)), class = "mt_error_argument")
  expect_error(mt_plot_contrast(data.frame(a = 1)), class = "mt_error_argument")
  # observed gaps alone, optionally for some levels
  obs <- mt_contrast(x, var = "sex", ref = "male")
  p <- mt_plot_contrast(observed = obs, levels = "female")
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(as.character(p$data$level)), "female")
  expect_error(mt_plot_contrast(), class = "mt_error_argument")
})

test_that("backtests score differences between strata path by path", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  bt <- mt_backtest(lee_carter(), x, origins = 2007, h = 3, n_draws = 50, quantities = "life_expectancy",
                    contrast = c(var = "sex", ref = "male"))
  gap <- bt[target == "life_expectancy_gap"]
  le <- bt[target == "life_expectancy"]
  expect_identical(unique(gap$sex), "female")
  expect_equal(nrow(gap), 3L)
  expect_equal(gap$observed, le[sex == "female"]$observed - le[sex == "male"]$observed)
  s <- mt_score(bt, target = "life_expectancy_gap")
  expect_identical(s$unit, "years (difference)")
  expect_true("crps_le_gap" %in% names(mt_score_table(bt)))
  expect_error(mt_backtest(lee_carter(), x, origins = 2007, h = 3, contrast = c(var = "sex", ref = "male")),
               class = "mt_error_argument")
  expect_error(mt_backtest(lee_carter(), x, origins = 2007, h = 3, quantities = "life_expectancy",
                           contrast = c(var = "sex", ref = "both")), class = "mt_error_argument")
})

test_that("mt_score_table() separates the levels of a hierarchy", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  h <- mt_hierarchy(sex = "both")
  xa <- mt_aggregate_strata(x, h)
  bt <- mt_backtest(reconciled(lee_carter(), h, method = "bottom_up"), xa, origins = 2009, h = 3, n_draws = 50)
  tab <- mt_score_table(bt, hierarchy = h)
  expect_identical(as.character(tab$level), c("bottom", "total"))
  expect_false("level" %in% names(bt))
})
