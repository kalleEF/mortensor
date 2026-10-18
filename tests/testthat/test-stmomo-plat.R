plat_ages <- function(years = 1990:2019) {
  x <- mt_example(deaths = "single", years = years)
  a <- single_ages(55, 89, open = FALSE)
  mt_regroup(mt_set_grid(x, age = a), age = a)
}

test_that("stmomo(\"plat\") is the Plat model of the StMoMo vignette", {
  skip_if_not_installed("StMoMo")
  # the model as defined in the code of the vignette (Villegas, Kaishev and
  # Millossovich 2018, Section 4)
  f2 <- function(x, ages) mean(ages) - x
  const_plat <- function(ax, bx, kt, b0x, gc, wxt, ages) {
    n_years <- dim(wxt)[2]
    x <- ages
    t <- 1:n_years
    c <- (1 - tail(ages, 1)):(n_years - ages[1])
    xbar <- mean(x)
    phi <- coef(lm(gc ~ 1 + c + I(c^2), na.action = na.omit))
    gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2
    kt[2, ] <- kt[2, ] + 2 * phi[3] * t
    kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t^2 - 2 * xbar * t)
    ax <- ax + phi[1] - phi[2] * x + phi[3] * x^2
    ci <- rowMeans(kt, na.rm = TRUE)
    ax <- ax + ci[1] + ci[2] * (xbar - x)
    kt[1, ] <- kt[1, ] - ci[1]
    kt[2, ] <- kt[2, ] - ci[2]
    list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
  }
  vignette_plat <- StMoMo::StMoMo(link = "logit", staticAgeFun = TRUE, periodAgeFun = c("1", f2),
                                  cohortAgeFun = "1", constFun = const_plat)
  x <- plat_ages()
  fit <- mt_fit(stmomo("plat"), x)
  expect_true(all(fit$diagnostics$status == "ok"))
  ages <- 55:89
  years <- 1990:2019
  D <- mt_array(x, "deaths")[, , 1]
  E <- mt_array(x, "exposure")[, , 1]
  Dx <- round(D)
  Ex <- E + D / 2
  dimnames(Dx) <- dimnames(Ex) <- list(ages, years)
  ref <- StMoMo::fit(vignette_plat, Dxt = Dx, Ext = Ex, ages = ages, years = years,
                     wxt = StMoMo::genWeightMat(ages, years, clip = 3), verbose = FALSE)
  own <- fit$fits[[1]]$fit
  expect_equal(own$kt, ref$kt, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(own$gc, ref$gc, tolerance = 1e-10, ignore_attr = TRUE)
  q <- stats::fitted(ref, type = "rates")
  ok <- !is.na(q)
  expect_equal(stats::fitted(fit)[, , 1][ok], -log(1 - q[ok]), tolerance = 1e-10)
  # constraints: centred period indexes, cohort effect without level, linear and quadratic trend
  expect_lt(max(abs(rowMeans(own$kt))), 1e-8)
  cc <- seq_along(own$gc)
  ok_g <- !is.na(own$gc)
  expect_lt(max(abs(stats::coef(stats::lm(own$gc[ok_g] ~ cc[ok_g] + I(cc[ok_g]^2))))), 1e-8)
  # forecasts: simulate() with ARIMA(2,0,0) for the cohort index, as in the vignette
  fc <- mt_forecast(fit, h = 10, n_draws = 50, seed = 4)
  set.seed(4, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  utils::capture.output(sim <- suppressWarnings(stats::simulate(ref, nsim = 50, h = 10, gc.order = c(2, 0, 0))))
  expect_equal(fc$rates[, , 1, ], -log(1 - sim$rates), tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("stmomo(\"plat\") uses ARIMA(2,0,0) by default and needs single years of age", {
  skip_if_not_installed("StMoMo")
  expect_identical(stmomo("plat")$params$gc_order, c(2, 0, 0))
  expect_identical(stmomo("plat", gc_order = c(1, 1, 0))$params$gc_order, c(1, 1, 0))
  expect_error(stmomo("plat", xc = 89), class = "mt_error_argument")
  x <- plat_ages(1995:2012)
  xg <- mt_regroup(x, age = seq(55, 85, 5))
  fit <- tryCatch(mt_fit(stmomo("plat"), xg), error = function(e) e)
  expect_true(inherits(fit, "mt_error_model_fit") || all(fit$diagnostics$status == "failed"))
  expect_match(summary(mt_fit(stmomo("plat"), x))$assumptions, "Plat")
  fc <- mt_forecast(mt_fit(stmomo("plat"), x), h = 5, n_draws = 20, seed = 1)
  expect_true(all(is.finite(fc$rates) & fc$rates > 0 & fc$rates < 1))
})
