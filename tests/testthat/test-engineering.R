test_that("mt_parameters() of ensembles gives each component and the weights", {
  x <- mt_example(deaths = "single", years = 1995:2015)
  fit <- mt_fit(ensemble(lc = lee_carter(), rw = rwd()), x)
  p <- mt_parameters(fit)
  expect_setequal(unique(p$component), c("lc", "rw"))
  w <- p[parameter == "weight"]
  expect_identical(w$index, c("lc", "rw"))
  expect_equal(sum(w$value), 1)
  expect_true(all(c("a", "b", "k") %in% p[component == "lc", parameter]))
  expect_true(all(c("drift", "sd") %in% p[component == "rw", parameter]))
  # a wrapper around an ensemble keeps the rows of the ensemble
  rc <- mt_fit(recalibrated(ensemble(lc = lee_carter(), rw = rwd())), x)
  pr <- mt_parameters(rc)
  expect_true(all(c("weight", "scale_h1") %in% pr$parameter))
})

test_that("mt_parameters() of fits by groups gives each group", {
  x <- mt_example(deaths = "single", years = 1995:2015)
  fit <- mt_fit(stratified("sex", female = lee_carter(), male = rwd()), x)
  p <- mt_parameters(fit)
  expect_setequal(unique(p$group), c("female", "male"))
  expect_true(all(p[group == "female", stratum] == "female"))
  expect_true("k" %in% p[group == "female", parameter])
  expect_true("drift" %in% p[group == "male", parameter])
  tf <- mt_fit(tuned(function(w) rwd(drift_window = w), data.frame(w = c(5, 10)), by = "sex", inner_origins = 2,
                     inner_h = 3), x)
  pt <- mt_parameters(tf)
  expect_setequal(unique(pt$group), c("female", "male"))
})

test_that("backtests and recalibration do not depend on the number of cores", {
  skip_on_os("windows")
  x <- mt_example(deaths = "single", years = 1995:2015)
  spec <- lee_carter(bootstrap = TRUE, n_boot = 3)
  bt <- function(cores) {
    set.seed(3)
    mt_backtest(spec, x, origins = c(2008, 2010), h = 3, n_draws = 20, cores = cores)
  }
  b1 <- bt(1)
  expect_identical(bt(1), b1)
  expect_equal(bt(2), b1)
  rc <- function(cores) {
    set.seed(3)
    f <- mt_fit(recalibrated(spec, years = 6, n_origins = 2), x, cores = cores)
    f$calibration
  }
  expect_equal(rc(2), rc(1))
})

test_that("par_lapply() with random number streams passes the elements, not their positions", {
  set.seed(1)
  expect_identical(par_lapply(c(10, 20, 30), function(v) v, streams = TRUE), list(10, 20, 30))
  set.seed(1)
  a <- par_lapply(1:3, function(i) stats::runif(1), streams = TRUE)
  set.seed(1)
  b <- par_lapply(1:3, function(i) stats::runif(1), cores = 2, streams = TRUE)
  expect_identical(a, b)
})

test_that("the prepared ARFIMA models of coherent_fdm() give the draws of forecast's simulate()", {
  skip_if_not_installed("demography")
  skip_if_not_installed("fracdiff")
  x <- mt_example(deaths = "single", years = 1990:2015)
  fit <- suppressWarnings(mt_fit(coherent_fdm(order1 = 2, order2 = 2), x))
  r <- fit$fits$ratios[[1]]
  fc <- quiet_try(forecast::forecast(r, h = 4, method = "arfima", estim = "mle", drange = fit$fits$drange))
  inner <- function(f) lapply(f$coeff, function(cf) if (inherits(cf$model, "forecast")) cf$model$model else cf$model)
  expect_true(any(vapply(inner(fc), inherits, logical(1), what = "fracdiff")))
  prepared <- fracdiff_sim_models(fc)
  expect_false(any(vapply(inner(prepared), inherits, logical(1), what = "fracdiff")))
  a <- withr::with_seed(5, stats::simulate(fc, nsim = 6))
  b <- withr::with_seed(5, stats::simulate(prepared, nsim = 6))
  expect_identical(b, a)
})
