test_that("compute_dx and dx_to_mx roundtrip", {
  ages <- seq(30, 90, by = 5)
  # Gompertz-like mortality schedule
  mx_orig <- exp(-10.5 + 0.085 * ages)

  dx <- compute_dx(mx_orig, ages, width = 5)
  expect_true(all(dx > 0))
  expect_equal(sum(dx), 1, tolerance = 1e-10)

  mx_back <- dx_to_mx(dx, ages, width = 5)
  # Roundtrip should be close for closed intervals (not open)
  n <- length(ages)
  expect_equal(mx_back[1:(n - 2)], mx_orig[1:(n - 2)], tolerance = 0.05)
})

test_that("inv_clr_to_dx returns valid composition", {
  clr_vec <- rnorm(10)
  dx <- inv_clr_to_dx(clr_vec)
  expect_true(all(dx > 0))
  expect_equal(sum(dx), 1, tolerance = 1e-10)
})

test_that("life_table produces valid output", {
  ages <- seq(30, 90, by = 5)
  mx <- exp(-10.5 + 0.085 * ages)

  lt <- life_table(mx, ages, width = 5)
  expect_s3_class(lt, "data.frame")
  expect_equal(nrow(lt), length(ages))
  expect_true(all(c("age", "mx", "qx", "px", "lx", "dx", "Lx", "Tx", "ex") %in% names(lt)))

  # qx in [0, 1]
  expect_true(all(lt$qx >= 0 & lt$qx <= 1))
  # lx monotonically non-increasing
  expect_true(all(diff(lt$lx) <= 0))
  # ex positive
  expect_true(all(lt$ex >= 0))
  # Last qx is 1 (open interval)
  expect_equal(lt$qx[nrow(lt)], 1)
  # lx starts at radix
  expect_equal(lt$lx[1], 100000)
})

test_that("life_table handles different open_age_methods", {
  ages <- seq(30, 90, by = 5)
  mx <- exp(-10.5 + 0.085 * ages)

  lt_rec <- life_table(mx, ages, open_age_method = "reciprocal")
  lt_kan <- life_table(mx, ages, open_age_method = "kannisto")
  lt_gom <- life_table(mx, ages, open_age_method = "gompertz")

  # All should produce valid life tables
  for (lt in list(lt_rec, lt_kan, lt_gom)) {
    expect_true(all(lt$ex >= 0))
    expect_true(all(lt$lx >= 0))
  }
})

test_that("mort_data constructor validates dimensions", {
  mx <- array(runif(60), dim = c(3, 4, 5))
  pop <- array(runif(60), dim = c(3, 4, 5))

  md <- mort_data(mx, pop, ages = 1:3, years = 2001:2004,
                  groups = paste0("G", 1:5))
  expect_s3_class(md, "mort_data")

  # Mismatched ages
  expect_error(mort_data(mx, pop, ages = 1:4, years = 2001:2004,
                         groups = paste0("G", 1:5)))
})

test_that("mort_config creates valid config with defaults", {
  cfg <- mort_config()
  expect_s3_class(cfg, "mort_config")
  expect_null(cfg$cause_cols)
  expect_equal(cfg$forecast_h, 10)
  expect_equal(cfg$disagg_method, "natural_gompertz")
  expect_equal(cfg$pi_method, "bootstrap")
})

test_that("mort_forecast constructor works", {
  mx_point <- array(runif(30), dim = c(3, 5, 2))
  md <- mort_data(array(runif(24), dim = c(3, 4, 2)),
                  array(runif(24), dim = c(3, 4, 2)),
                  ages = 1:3, years = 2001:2004, groups = c("A", "B"))

  fc <- mort_forecast(mx_point, train = md,
                      forecast_years = 2005:2009, model = "test")
  expect_s3_class(fc, "mort_forecast")
  expect_equal(fc$model, "test")
  expect_equal(fc$forecast_years, 2005:2009)
})

test_that("compute_mase returns expected values", {
  train <- c(1, 2, 3, 4, 5)
  # Perfect forecast
  actual_test <- c(6, 7)
  predicted_test <- c(6, 7)
  expect_equal(compute_mase(actual_test, predicted_test, train), 0)

  # Naive forecast (repeat last)
  predicted_naive <- c(5, 5)
  mase_val <- compute_mase(actual_test, predicted_naive, train)
  expect_true(is.finite(mase_val))
  expect_true(mase_val > 0)
})

test_that("compute_interval_score works correctly", {
  # Observation inside interval: score = width
  expect_equal(compute_interval_score(5, 3, 7, alpha = 0.05), 4)

  # Observation below interval
  score_below <- compute_interval_score(1, 3, 7, alpha = 0.05)
  expect_true(score_below > 4)  # width + penalty

  # Observation above interval
  score_above <- compute_interval_score(9, 3, 7, alpha = 0.05)
  expect_true(score_above > 4)
})

test_that("compute_rmse_log handles edge cases", {
  expect_true(is.na(compute_rmse_log(c(-1, 0), c(1, 1))))
  expect_equal(compute_rmse_log(c(1, 1), c(1, 1)), 0)
})

test_that("compute_mae works", {
  expect_equal(compute_mae(c(1, 2, 3), c(1, 2, 3)), 0)
  expect_equal(compute_mae(c(1, 2, 3), c(2, 3, 4)), 1)
})

test_that("build_array constructs correct dimensions", {
  # build_array(dt, value_col, ages, years, groups) - uses pre-specified vectors
  ages  <- c(30, 35)
  years <- 2001:2002
  groups <- data.frame(group = c("A", "B"))
  dt <- data.table::data.table(
    age_group = rep(paste0(ages, "-", ages + 4), 4),
    year      = rep(rep(years, each = 2), 2),
    group     = rep(c("A", "B"), 4),
    value     = runif(8)
  )
  arr <- build_array(dt, value_col = "value",
                     ages = ages, years = years, groups = groups)
  expect_equal(dim(arr), c(2, 2, 2))
})
