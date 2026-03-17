test_that("mort_calibration constructor and print work", {
  # Minimal mort_calibration object
  md <- mort_data(
    mx = array(runif(120, 0.001, 0.1), dim = c(3, 4, 10)),
    pop = array(round(runif(120, 1000, 50000)), dim = c(3, 4, 10)),
    ages = c(30, 35, 40), years = 2001:2004,
    groups = paste0("G", 1:10)
  )

  comp <- data.table::data.table(
    model = c("Tucker3-CoDa", "RWD"),
    mase_combined = c(0.8, 1.2)
  )

  cal <- mort_calibration(
    data = md,
    models = list(),
    comparison = comp,
    best_model = "Tucker3-CoDa",
    ensemble_weights = c(0.6, 0.4)
  )

  expect_s3_class(cal, "mort_calibration")
  expect_output(print(cal), "mort_calibration")
})

test_that("mort_shared_config works", {
  sc <- mort_shared_config(output_dir = tempdir(), forecast_h = 5)
  expect_s3_class(sc, "mort_shared_config")
  expect_equal(sc$forecast_h, 5)
})

test_that("life_expectancy.mort_data returns data.table", {
  md <- mort_data(
    mx = array(exp(-10.5 + 0.085 * rep(seq(30, 90, 5), 4 * 2)),
               dim = c(13, 4, 2)),
    pop = array(10000, dim = c(13, 4, 2)),
    ages = seq(30, 90, 5), years = 2001:2004,
    groups = c("A", "B")
  )

  le <- life_expectancy(md)
  expect_s3_class(le, "data.table")
  expect_true("ex" %in% names(le))
  expect_true(all(le$ex > 0))
  expect_equal(nrow(le), 4 * 2)  # years * groups
})

test_that("life_expectancy.mort_forecast returns data.table with PI", {
  md <- mort_data(
    mx = array(exp(-10.5 + 0.085 * rep(seq(30, 90, 5), 4 * 2)),
               dim = c(13, 4, 2)),
    pop = array(10000, dim = c(13, 4, 2)),
    ages = seq(30, 90, 5), years = 2001:2004,
    groups = c("A", "B")
  )

  mx_fc <- array(exp(-10.5 + 0.085 * rep(seq(30, 90, 5), 3 * 2)),
                 dim = c(13, 3, 2))
  mx_lo <- mx_fc * 0.8
  mx_hi <- mx_fc * 1.2

  fc <- mort_forecast(
    mx_point = mx_fc, train = md,
    forecast_years = 2005:2007, model = "test",
    mx_lower = mx_lo, mx_upper = mx_hi
  )

  le <- life_expectancy(fc)
  expect_s3_class(le, "data.table")
  expect_true(all(c("ex", "ex_lower", "ex_upper") %in% names(le)))
  # Lower mortality -> higher life expectancy
  expect_true(all(le$ex_upper >= le$ex_lower))
})

test_that("tucker3_info_criterion returns valid output", {
  skip_if_not_installed("rTensor")
  skip_if_not_installed("compositions")

  set.seed(123)
  mx <- array(exp(-10 + 0.08 * rep(seq(30, 50, 5), 10 * 3)),
              dim = c(5, 10, 3))
  # Add noise

  mx <- mx * exp(array(rnorm(150, 0, 0.05), dim = c(5, 10, 3)))

  res <- tucker3_info_criterion(mx, R_A = 2, R_B = 2, R_C = 2,
                                 method = "direct", criterion = "bic")
  expect_true(is.finite(res$criterion))
  expect_true(res$fit_pct > 0)
  expect_equal(length(res$ranks), 3)

  # AIC should differ from BIC
  res_aic <- tucker3_info_criterion(mx, R_A = 2, R_B = 2, R_C = 2,
                                     method = "direct", criterion = "aic")
  expect_true(res$criterion != res_aic$criterion)
})

test_that("compute_wis works with symmetric quantiles", {
  n <- 10
  obs <- rnorm(n)
  # Create quantile forecasts: 2.5%, 25%, 50%, 75%, 97.5%
  probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  qmat <- matrix(NA, n, 5)
  for (i in 1:n) {
    qmat[i, ] <- qnorm(probs, mean = obs[i], sd = 1)
  }

  wis <- compute_wis(obs, qmat, probs)
  expect_true(is.finite(wis))
  expect_true(wis > 0)
})

test_that("compute_interval_score_array works on 3D arrays", {
  obs <- array(exp(rnorm(24)), dim = c(3, 4, 2))
  lower <- obs * 0.8
  upper <- obs * 1.2

  # Obs inside interval
  is_val <- compute_interval_score_array(obs, lower, upper)
  expect_true(is.finite(is_val))
  expect_true(is_val > 0)
})

# ---- Tests for unified API (Phases 1-4) ----

test_that("plot.mort_calibration works", {
  md <- mort_data(
    mx = array(runif(120, 0.001, 0.1), dim = c(3, 4, 10)),
    ages = c(30, 35, 40), years = 2001:2004,
    groups = paste0("G", 1:10)
  )
  comp <- data.table::data.table(
    model = c("Tucker3-CoDa", "RWD", "ETS"),
    mase_combined = c(0.8, 1.2, 1.0)
  )
  cal <- mort_calibration(
    data = md, models = c("Tucker3-CoDa", "RWD", "ETS"),
    comparison = comp, best_model = "Tucker3-CoDa"
  )

  # comparison plot should produce a ggplot
  p <- plot(cal, type = "comparison")
  expect_s3_class(p, "gg")

  # ranks plot gracefully returns NULL when no rank data
  expect_message(plot(cal, type = "ranks"), "No rank search")
})

test_that("plot.mort_forecast extended types work", {
  md <- mort_data(
    mx = array(exp(-10 + 0.08 * rep(seq(30, 50, 5), 4 * 2)),
               dim = c(5, 4, 2)),
    ages = seq(30, 50, 5), years = 2001:2004,
    groups = c("A", "B")
  )
  mx_fc <- array(exp(-10 + 0.08 * rep(seq(30, 50, 5), 3 * 2)),
                 dim = c(5, 3, 2))

  fc <- mort_forecast(
    mx_point = mx_fc, train = md,
    forecast_years = 2005:2007, model = "test",
    mx_raw = md$mx,
    extras = list(mx_smooth = md$mx * 0.99)
  )

  # smoothing plot should work
  p <- plot(fc, type = "smoothing")
  expect_s3_class(p, "gg")

  # components gracefully returns NULL without B matrix
  expect_message(plot(fc, type = "components"), "No Tucker3 B-mode")

  # model_comparison gracefully returns NULL without fit_info
  expect_message(plot(fc, type = "model_comparison"), "No model comparison")

  # disagg gracefully returns NULL without mx_1yr_dt
  expect_message(plot(fc, type = "disagg"), "No single-year")
})

test_that("export.mort_forecast writes CSV", {
  md <- mort_data(
    mx = array(0.01, dim = c(3, 4, 2)),
    ages = c(30, 35, 40), years = 2001:2004,
    groups = c("A", "B")
  )
  mx_fc <- array(0.01, dim = c(3, 2, 2))
  fc <- mort_forecast(
    mx_point = mx_fc, train = md,
    forecast_years = 2005:2006, model = "test"
  )

  tmp <- tempdir()
  export(fc, dir = tmp, prefix = "test_export")
  expect_true(file.exists(file.path(tmp, "test_export_forecast_grouped.csv")))
  unlink(file.path(tmp, "test_export_forecast_grouped.csv"))
})

test_that("export.mort_calibration writes CSV", {
  md <- mort_data(
    mx = array(0.01, dim = c(3, 4, 2)),
    ages = c(30, 35, 40), years = 2001:2004,
    groups = c("A", "B")
  )
  comp <- data.table::data.table(
    model = c("A", "B"), mase_combined = c(0.5, 0.8)
  )
  cal <- mort_calibration(
    data = md, models = c("A", "B"),
    comparison = comp, best_model = "A"
  )

  tmp <- tempdir()
  export(cal, dir = tmp, prefix = "test_calib")
  expect_true(file.exists(file.path(tmp, "test_calib_model_comparison.csv")))
  unlink(file.path(tmp, "test_calib_model_comparison.csv"))
})

test_that("reconcile produces mort_forecast objects", {
  md_a <- mort_data(
    mx = array(0.01, dim = c(3, 4, 2)),
    ages = c(30, 35, 40), years = 2001:2004,
    groups = c("A", "B"), cause = "cancer"
  )
  md_b <- mort_data(
    mx = array(0.005, dim = c(3, 4, 2)),
    ages = c(30, 35, 40), years = 2001:2004,
    groups = c("A", "B"), cause = "cvd"
  )
  md_all <- mort_data(
    mx = array(0.02, dim = c(3, 4, 2)),
    ages = c(30, 35, 40), years = 2001:2004,
    groups = c("A", "B"), cause = "all"
  )

  fc_a <- mort_forecast(
    mx_point = array(0.012, dim = c(3, 2, 2)), train = md_a,
    forecast_years = 2005:2006, model = "test",
    mx_lower = array(0.010, dim = c(3, 2, 2)),
    mx_upper = array(0.014, dim = c(3, 2, 2))
  )
  fc_b <- mort_forecast(
    mx_point = array(0.006, dim = c(3, 2, 2)), train = md_b,
    forecast_years = 2005:2006, model = "test",
    mx_lower = array(0.004, dim = c(3, 2, 2)),
    mx_upper = array(0.008, dim = c(3, 2, 2))
  )
  fc_all <- mort_forecast(
    mx_point = array(0.02, dim = c(3, 2, 2)), train = md_all,
    forecast_years = 2005:2006, model = "test"
  )

  recon <- reconcile(
    list(cancer = fc_a, cvd = fc_b),
    all_cause = fc_all
  )

  expect_type(recon, "list")
  expect_length(recon, 2)
  expect_s3_class(recon$cancer, "mort_forecast")
  expect_s3_class(recon$cvd, "mort_forecast")
  expect_true(grepl("reconciled", recon$cancer$model))

  # Reconciled point forecasts should differ from originals
  # (unless they already summed perfectly)
  total_recon <- recon$cancer$mx_point + recon$cvd$mx_point
  total_orig <- fc_a$mx_point + fc_b$mx_point
  expect_false(isTRUE(all.equal(total_recon, total_orig)))
})
