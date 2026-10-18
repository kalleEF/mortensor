# Single years of age 55-89, as in the examples of StMoMo for the CBD family
old_ages <- function(years = 1990:2019) {
  x <- mt_example(deaths = "single", years = years)
  a <- single_ages(55, 89, open = FALSE)
  mt_regroup(mt_set_grid(x, age = a), age = a)
}

test_that("stmomo() M7 and M8 agree with StMoMo and forecast plausible rates", {
  skip_if_not_installed("StMoMo")
  x <- old_ages()
  D <- mt_array(x, "deaths")
  E <- mt_array(x, "exposure")
  ages <- 55:89
  years <- 1990:2019
  for (m in c("m7", "m8")) {
    fit <- mt_fit(stmomo(m), x)
    expect_true(all(fit$diagnostics$status == "ok"))
    f <- stats::fitted(fit)
    expect_identical(dim(f), c(35L, 30L, 2L))
    expect_false(anyNA(f))
    for (s in 1:2) {
      # StMoMo called directly: logit link, initial exposures E + D / 2, cohorts clipped
      Dx <- round(D[, , s])
      Ex <- E[, , s] + D[, , s] / 2
      dimnames(Dx) <- dimnames(Ex) <- list(ages, years)
      mod <- if (m == "m7") StMoMo::m7() else StMoMo::m8(xc = 89)
      ref <- StMoMo::fit(mod, Dxt = Dx, Ext = Ex, ages = ages, years = years,
                         wxt = StMoMo::genWeightMat(ages, years, clip = 3), verbose = FALSE)
      q <- stats::fitted(ref, type = "rates")
      ok <- !is.na(q)
      expect_gt(mean(ok), 0.9)
      expect_equal(f[, , s][ok], -log(1 - q[ok]), tolerance = 1e-8)
      # fitted rates are close to the observed ones
      expect_lt(stats::median(abs(log(f[, , s] / (pmax(D[, , s], 0.5) / E[, , s])))), 0.05)
    }
    fc <- mt_forecast(fit, h = 10, n_draws = 50, seed = 1)
    expect_identical(dim(fc$rates), c(35L, 10L, 2L, 50L))
    expect_true(all(is.finite(fc$rates) & fc$rates > 0 & fc$rates < 1))
    expect_equal(fc$periods, 2020:2029)
    # medians stay near the last fitted rates and keep rising with age
    med <- apply(fc$rates, 1:3, stats::median)
    expect_true(all(abs(log(med[, 1, ] / f[, 30, ])) < 0.2))
    expect_true(all(abs(log(med[, 10, ] / f[, 30, ])) < 0.6))
    expect_true(all(diff(log(med[, 10, 1])) > -0.1))
    expect_true(all(apply(fc$rates[, 10, , ], c(1, 2), stats::sd) > 0))
  }
})

test_that("stmomo() M8 uses `xc` and M7 forecasts its cohort index with ARIMA(2,0,0)", {
  skip_if_not_installed("StMoMo")
  x <- old_ages(1995:2019)
  s100 <- stmomo("m8", xc = 100)
  expect_match(s100$label, "xc = 100")
  expect_identical(stmomo("m8")$params$xc, NULL)
  f100 <- mt_fit(s100, x)
  f89 <- mt_fit(stmomo("m8"), x)
  st <- f100$fits[[1]]$fit$model
  expect_equal(st$cohortAgeFun(60, 55:89), 40)
  expect_equal(f89$fits[[1]]$fit$model$cohortAgeFun(60, 55:89), 29)
  expect_false(isTRUE(all.equal(stats::fitted(f100), stats::fitted(f89))))
  expect_identical(stmomo("m7")$params$gc_order, c(2, 0, 0))
  expect_identical(stmomo("m6")$params$gc_order, c(1, 1, 0))
  expect_identical(stmomo("m7", gc_order = c(1, 1, 0))$params$gc_order, c(1, 1, 0))
  expect_error(stmomo("m7", xc = 90), class = "mt_error_argument")
  expect_error(stmomo("m8", xc = "a"), class = "mt_error_argument")
  expect_error(stmomo("m7", gc_order = c(1, 1)), class = "mt_error_argument")
})

test_that("stmomo() M7 and M8 need single years of age and run in backtests", {
  skip_if_not_installed("StMoMo")
  x <- old_ages(1990:2012)
  xg <- mt_regroup(x, age = seq(55, 85, 5))
  for (m in c("m7", "m8")) {
    fit <- tryCatch(mt_fit(stmomo(m), xg), error = function(e) e)
    expect_true(inherits(fit, "mt_error_model_fit") || all(fit$diagnostics$status == "failed"))
  }
  bt <- mt_compare(list(cbd = stmomo("cbd"), m7 = stmomo("m7"), m8 = stmomo("m8")), x, origins = 2005:2007, h = 3,
                   n_draws = 30)
  s <- mt_score(bt)
  expect_identical(as.character(s$model), c("cbd", "m7", "m8"))
  expect_true(all(s$crps > 0 & s$crps < 0.2))
})
