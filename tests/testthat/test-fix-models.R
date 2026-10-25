# Regression tests for fixes in the base models

test_that("Lee-Carter with a chosen period keeps its bootstrap fits", {
  x <- mt_example(deaths = "single", years = 2000:2019)
  set.seed(1)
  fit <- mt_fit(lee_carter(choose_period = TRUE, min_years = 10, bootstrap = TRUE, n_boot = 4), x)
  expect_true(all(fit$diagnostics$status == "ok"))
  for (f in fit$fits) {
    expect_length(f$boot, 4)
    expect_length(f$boot[[1]]$k, length(f$k))
  }
  fc <- mt_forecast(fit, h = 2, n_draws = 8, seed = 1)
  expect_false(anyNA(fc$rates))
})

test_that("SVD Lee-Carter matches life expectancy over all ages when age groups have no deaths", {
  x <- mt_subset(mt_example(deaths = "single", years = 2000:2019), sex == "female")
  D <- mt_array(x, "deaths")[, , 1]
  E <- mt_array(x, "exposure")[, , 1]
  ages <- mt_grid(x)$age
  D[c(2, 5), ] <- 0
  expect_warning(est <- mortensor:::fit_model(lee_carter(method = "svd", adjust = "e0"), D, E, ages, 2000:2019),
                 class = "mt_warning_sparse")
  obs <- pmax(D, 0.5) / E
  obs[c(2, 5), ] <- est$fitted[c(2, 5), ]
  expect_equal(as.vector(life_expectancy(est$fitted, ages)), as.vector(life_expectancy(obs, ages)), tolerance = 1e-6)
  # sparse counts: life expectancy is not monotone in k; the root nearest to the SVD estimate is used
  xs <- mt_example(deaths = "single", years = 2000:2019, size = 200, seed = 1)
  D <- mt_array(xs, "deaths")[, , 1]
  E <- mt_array(xs, "exposure")[, , 1]
  e <- mortensor:::lc_svd(D, E, adjust = "e0", ages = ages)
  expect_equal(as.vector(life_expectancy(exp(e$a + outer(e$b, e$k)), ages)),
               as.vector(life_expectancy(pmax(D, 0.5) / E, ages)), tolerance = 1e-6)
})

test_that("StMoMo forecasts one year ahead", {
  skip_if_not_installed("StMoMo")
  x <- mt_example(deaths = "single", years = 2000:2019)
  fit <- suppressWarnings(mt_fit(stmomo("lc"), x))
  fc <- mt_forecast(fit, h = 1, n_draws = 5, seed = 1)
  expect_equal(dim(fc$rates), c(length(fit$age), 1, 2, 5))
  expect_true(all(is.finite(fc$rates) & fc$rates > 0))
})

test_that("joint models run with a single age group", {
  x <- mt_regroup(mt_example(deaths = "single", years = 2000:2019), age = 30)
  for (spec in list(tucker3(ranks = c(1, 2, 2)), li_lee())) {
    fit <- mt_fit(spec, x)
    expect_true(all(fit$diagnostics$status == "ok"))
    fc <- mt_forecast(fit, h = 2, n_draws = 5, seed = 1, predictive = TRUE)
    expect_equal(dim(fc$rates), c(1, 2, 2, 5))
  }
  # the default ranks are reduced to the one age group
  expect_warning(mt_fit(tucker3(), x), "c\\(1, 2, 2\\) are used", class = "mt_warning_argument")
})

test_that("fits with random numbers follow set.seed() and do not depend on cores", {
  skip_on_os("windows")
  x <- mt_example(deaths = "single", years = 2000:2019)
  one <- function(cores) {
    set.seed(42)
    fit <- mt_fit(lee_carter(bootstrap = TRUE, n_boot = 2), x, cores = cores)
    list(k = lapply(fit$fits, function(f) f$boot[[2]]$k), after = stats::runif(1))
  }
  r1 <- one(1)
  expect_identical(one(1), r1)
  expect_identical(one(2), r1)
  expect_false(identical(r1$k[[1]], r1$k[[2]]))
})

test_that("coherent_fdm() needs two strata per group", {
  skip_if_not_installed("demography")
  x <- mt_subset(mt_example(deaths = "single", years = 2000:2019), sex == "female")
  expect_error(mt_fit(coherent_fdm(order1 = 2, order2 = 1), x), class = "mt_error_model_fit", regexp = "two strata")
})

test_that("models of observed rates fitted to smoothed rates add observation noise", {
  x <- mt_regroup(mt_example(deaths = "grouped", years = 1995:2019), age = "deaths")
  raw <- mt_fit(rwd(), x)
  expect_false(raw$smoothed_input)
  # observed rates: predictive draws add nothing
  a <- mt_forecast(raw, h = 2, n_draws = 20, seed = 1)
  b <- mt_forecast(raw, h = 2, n_draws = 20, seed = 1, predictive = TRUE)
  expect_equal(b$rates, a$rates)
  fit <- mt_fit(rwd(), x, prep = mt_prep(step_smooth(along = "period", method = "whittaker")))
  expect_true(fit$smoothed_input)
  d <- mt_forecast(fit, h = 3, n_draws = 200, seed = 1)
  p <- mt_forecast(fit, h = 3, n_draws = 200, seed = 1, predictive = TRUE)
  # predictive draws are deaths given the exposures: wider than the draws of the smoothed rates
  sd_log <- function(a) mean(apply(log(pmax(a, 1e-12)), 1:3, stats::sd))
  expect_gt(sd_log(p$rates), sd_log(d$rates))
})
