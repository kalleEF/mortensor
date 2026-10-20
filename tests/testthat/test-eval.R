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

test_that("score_samples() gives the spread, 90% coverage and half-widths of the draws", {
  set.seed(4)
  X <- matrix(stats::rnorm(3 * 4000, sd = 2), 3)
  y <- c(0, 3.5, -5)
  sc <- mortensor:::score_samples(X, y)
  expect_equal(sc$sd, apply(X, 1, stats::sd))
  Q <- mortensor:::row_quantiles(X, c(0.025, 0.05, 0.1, 0.9, 0.95, 0.975))
  expect_equal(sc$hw80, (Q[, 4] - Q[, 3]) / 2, ignore_attr = TRUE)
  expect_equal(sc$hw90, (Q[, 5] - Q[, 2]) / 2, ignore_attr = TRUE)
  expect_equal(sc$hw95, (Q[, 6] - Q[, 1]) / 2, ignore_attr = TRUE)
  expect_equal(sc$cov90, as.numeric(y >= Q[, 2] & y <= Q[, 5]))
  # normal quantiles: half-width of the 90% interval 1.645 sd
  expect_equal(sc$hw90, rep(2 * stats::qnorm(0.95), 3), tolerance = 0.05)
})

# A toy backtest table: two strata, three origins, life expectancy rows
toy_backtest <- function() {
  bt <- data.table::CJ(origin = 2001:2003, stratum = 1:2, period = 2002:2005)[period > origin]
  bt[, `:=`(model = "m", prep = "none", horizon = period - origin, age_lo = NA_real_, age_hi = NA_real_,
            target = "life_expectancy", scale = NA_real_)]
  bt[, median := 70 + 0.2 * (period - 2000) + 0.05 * (origin - 2000) + stratum]
  bt[, observed := median + c(0.3, -0.2, 0.5, -0.1)[period - 2001]]
  bt[, `:=`(mean = median, ae = abs(observed - median), sd = 0.4 + 0.1 * horizon, crps = 0, wis = 0,
            cov50 = 1, cov80 = 1, cov95 = 1, pit = 0.5, cov90 = 1, hw80 = 1, hw90 = 1.5, hw95 = 2)]
  bt[]
}

test_that("SAPE of mt_score() is the standardised absolute error of Raftery et al. (2013)", {
  bt <- toy_backtest()
  s <- mt_score(bt, by = c("model", "horizon"), target = "life_expectancy", measures = c("sape", "halfwidth"))
  hand <- bt[, list(sape = mean(sqrt(pi / 2) * abs(observed - median) / sd)), by = horizon]
  expect_equal(s$sape, hand$sape[match(s$horizon, hand$horizon)])
  expect_equal(s$hw90, rep(1.5, 4))
  # calibrated normal forecasts: mean SAPE about 1
  set.seed(9)
  X <- matrix(stats::rnorm(2000 * 1000), 2000)
  sc <- mortensor:::score_samples(X, stats::rnorm(2000))
  bn <- cbind(data.table(model = "normal", target = "life_expectancy", scale = NA_real_), sc)
  expect_equal(mt_score(bn, target = "life_expectancy", measures = "sape")$sape, 1, tolerance = 0.05)
  # intervals half as wide as they should be: about 2
  sc2 <- mortensor:::score_samples(X / 2, stats::rnorm(2000))
  bn2 <- cbind(data.table(model = "narrow", target = "life_expectancy", scale = NA_real_), sc2)
  expect_equal(mt_score(bn2, target = "life_expectancy", measures = "sape")$sape, 2, tolerance = 0.05)
  # defined for life expectancy only; backtests without `sd` are refused
  expect_error(mt_score(toy_backtest()[, target := "rate"], measures = "sape"), class = "mt_error_argument")
  expect_error(mt_score(toy_backtest()[, sd := NULL], target = "life_expectancy", measures = "sape"),
               class = "mt_error_argument")
  expect_error(mt_score(toy_backtest(), target = "life_expectancy", measures = "mape"), class = "mt_error_argument")
})

test_that("the revision measure is the change of the median between consecutive origins", {
  bt <- toy_backtest()
  s <- mt_score(bt, by = c("model", "horizon"), target = "life_expectancy", measures = c("revision", "cov90"))
  # by hand: the median of origin o minus that of origin o - 1 for the same
  # stratum and year (here 0.05 per origin), by the horizon of origin o
  hand <- merge(bt[, list(origin, stratum, period, horizon, median)],
                bt[, list(origin = origin + 1, stratum, period, prev = median)], by = c("origin", "stratum", "period"))
  hand <- hand[, list(revision = mean(abs(median - prev))), by = horizon]
  expect_equal(s[match(hand$horizon, horizon), revision], hand$revision)
  expect_equal(hand$revision, rep(0.05, 3), tolerance = 1e-10)
  # the longest horizon is forecast only from the first origin: no revision
  expect_true(is.na(s[horizon == 4, revision]))
  expect_equal(s$cov90, rep(1, 4))
  # cells missing at the previous origin are left out; models are separate
  b2 <- rbind(bt, copy(bt)[, `:=`(model = "m2", median = median * 2)])[!(origin == 2002 & stratum == 2)]
  s2 <- mt_score(b2, target = "life_expectancy", measures = "revision")
  m1 <- b2[model == "m"]
  r1 <- c(m1[origin == 2002, median] - m1[origin == 2001 & stratum == 1 & period > 2002, median],
          m1[origin == 2003 & stratum == 1, median] - m1[origin == 2002 & period > 2003, median])
  expect_equal(s2[model == "m", revision], mean(abs(r1)))
  # the default summary is unchanged
  expect_named(mt_score(bt, target = "life_expectancy"),
               c("model", "crps", "wis", "mae", "rmse", "mase", "cov50", "cov80", "cov95", "pit_mean", "n", "unit"))
})
