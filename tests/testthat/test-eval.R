test_that("WIS matches the definition of Bracher et al. (2021) on hand-computed cases", {
  # one central 80% interval [1, 3] and median 2: K = 1
  probs <- c(0.1, 0.5, 0.9)
  Q <- matrix(c(1, 2, 3), 1)
  # y inside: (0.5 * |y - m| + 0.1 * (u - l)) / 1.5
  expect_equal(mortensor:::wis_from_quantiles(Q, 2.5, probs), (0.5 * 0.5 + 0.1 * 2) / 1.5)
  # y above: interval score = (u - l) + 2 / 0.2 * (y - u)
  expect_equal(mortensor:::wis_from_quantiles(Q, 4, probs), (0.5 * 2 + 0.1 * (2 + 10 * 1)) / 1.5)
})

test_that("sample CRPS agrees with the closed form for a normal distribution", {
  set.seed(1)
  X <- matrix(rnorm(20000), 1)
  y <- 0.7
  closed <- y * (2 * pnorm(y) - 1) + 2 * dnorm(y) - 1 / sqrt(pi)
  sc <- mortensor:::score_samples(X, y)
  expect_equal(sc$crps, closed, tolerance = 0.03)
  expect_equal(sc$cov95, 1)
  expect_gt(sc$pit, 0.7)
  expect_lt(sc$pit, 0.8)
})

test_that("backtests refit on the training window only (no leakage)", {
  x <- mt_example(years = 1990:2012)
  prep <- mt_prep(step_ungroup())
  b1 <- mt_backtest(lee_carter(), x, origins = 2005, h = 3, prep = prep, n_draws = 50)
  # change deaths after the origin: forecasts must not change
  y <- x
  y$obs$deaths <- data.table::copy(y$obs$deaths)
  y$obs$deaths[period_lo > 2005, value := value * 3]
  b2 <- mt_backtest(lee_carter(), y, origins = 2005, h = 3, prep = prep, n_draws = 50)
  expect_equal(b1$median, b2$median)
  expect_false(isTRUE(all.equal(b1$observed, b2$observed)))
  # scored at the observed resolution (5-year groups and 90+)
  expect_setequal(unique(b1$age_hi - b1$age_lo), c(5, Inf))
  expect_equal(sort(unique(b1$horizon)), 1:3)
})

test_that("mt_score summarises by model and horizon", {
  x <- mt_example(deaths = "single", years = 1990:2010)
  bt <- rbind(mt_backtest(rwd(), x, origins = c(2004, 2006), h = 3, n_draws = 50),
              mt_backtest(lee_carter(), x, origins = c(2004, 2006), h = 3, n_draws = 50))
  s <- mt_score(bt, by = c("model", "horizon"))
  expect_equal(nrow(s), 6)
  expect_true(all(s$cov95 >= 0 & s$cov95 <= 1))
  expect_true(all(is.finite(s$mase)))
  expect_output(print(structure(bt, class = c("mt_backtest", class(bt)))), "mt_backtest")
})

test_that("interval calibration on the simulation benchmark", {
  skip_if(!identical(Sys.getenv("MORTENSOR_SLOW_TESTS"), "true"), "slow benchmark")
  # the generator is the Lee-Carter model with Poisson deaths and one period
  # index per stratum (10 strata): Tucker3 needs period and stratum ranks 10
  # to describe it (with the default ranks c(2, 2, 2) the 95% coverage is
  # about 0.74)
  for (sp in list(lee_carter(), tucker3(ranks = c(1, 10, 10)))) {
    s <- attr(mortensor:::benchmark_coverage(sp, n_rep = 6), "summary")
    expect_gt(s$cov95_mean, 0.90)
    expect_lt(s$cov95_mean, 0.98)
    expect_gt(s$cov95_min, 0.85)
  }
})
