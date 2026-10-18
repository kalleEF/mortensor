# A backtest table with known scores: models, origins, horizons and ages,
# with an origin effect shared by the cells of an origin
fake_backtest <- function(models = c("a", "b"), origins = 2000:2019, horizons = 1:3, n_age = 3, seed = 1) {
  set.seed(seed)
  g <- expand.grid(age_lo = 59 + seq_len(n_age), horizon = horizons, origin = origins, model = models,
                   stringsAsFactors = FALSE)
  g$age_hi <- g$age_lo + 1
  g$period <- g$origin + g$horizon
  g$stratum <- 1L
  g$sex <- "male"
  g$target <- "rate"
  shared <- stats::rexp(length(origins))[match(g$origin, origins)]
  g$crps <- 0.05 * (stats::rexp(nrow(g)) + shared) * ifelse(g$model == models[1], 1, 1.15) * (1 + 0.1 * g$horizon)
  g$wis <- 1.1 * g$crps + 0.01 * stats::runif(nrow(g))
  g$ae <- 1.3 * g$crps + 0.01 * stats::runif(nrow(g))
  data.table::as.data.table(g)
}

# Losses per origin of one model, metric and horizon (all horizons: the mean
# of the per-horizon means)
origin_losses <- function(bt, mod, metric = "crps", h = NULL) {
  b <- as.data.frame(bt)
  b <- b[b$model == mod, ]
  v <- if (metric == "se") b$ae^2 else b[[metric]]
  if (!is.null(h)) return(tapply(v[b$horizon == h], b$origin[b$horizon == h], mean))
  per_h <- tapply(v, list(b$origin, b$horizon), mean)
  rowMeans(per_h)
}

# The EWC test computed by hand (Lazarus et al. 2018): statistic and
# two-sided p-value for the loss differences d
ewc_by_hand <- function(d) {
  n <- length(d)
  nu <- max(1, floor(0.4 * n^(2 / 3)))
  lam <- numeric(nu)
  for (j in seq_len(nu)) {
    for (t in seq_len(n)) lam[j] <- lam[j] + sqrt(2 / n) * (d[t] - mean(d)) * cos(pi * j * (t - 0.5) / n)
  }
  stat <- sqrt(n) * mean(d) / sqrt(mean(lam^2))
  list(nu = nu, statistic = stat, p_value = 2 * stats::pt(-abs(stat), df = nu))
}

test_that("the default EWC test matches a hand computation", {
  bt <- fake_backtest()
  r <- mt_dm_test(bt, metrics = c("crps", "se"))
  expect_true(all(is.na(r$lag)))
  expect_true(all(r$df == 2))  # 20 origins: floor(0.4 * 20^(2/3)) = 2
  for (m in c("crps", "se")) {
    for (h in c(1:3, NA)) {
      e1 <- origin_losses(bt, "b", m, if (is.na(h)) NULL else h)
      e2 <- origin_losses(bt, "a", m, if (is.na(h)) NULL else h)
      ref <- ewc_by_hand(e1 - e2)
      row <- r[metric == m & horizon == if (is.na(h)) "all" else h]
      expect_equal(row$statistic, ref$statistic, tolerance = 1e-10)
      expect_equal(row$p_value, ref$p_value, tolerance = 1e-10)
      expect_equal(row$df, ref$nu)
      expect_equal(row$mean_diff, mean(e1 - e2), tolerance = 1e-12)
    }
  }
  # one-sided p-values and the degrees of freedom by the number of origins
  d <- c(0.3, -0.1, 0.4, 0.2, 0.5, -0.2, 0.1)
  s <- mortensor:::dm_stat(d, 1L, "less", "ewc")
  expect_equal(s$p_value, stats::pt(ewc_by_hand(d)$statistic, df = 1), tolerance = 1e-12)
  expect_identical(vapply(c(2, 7, 11, 12, 20, 21, 31, 32), mortensor:::dm_ewc_df, 0), c(1, 1, 1, 2, 2, 3, 3, 4))
})

test_that("the EWC test keeps its size with few overlapping origins", {
  # loss differences of forecasts 3 years ahead from 16 consecutive origins:
  # a moving average of order 2 with mean zero (no difference in accuracy)
  set.seed(20)
  n <- 16
  h <- 3
  p <- vapply(seq_len(500), function(i) {
    e <- stats::rnorm(n + h - 1)
    d <- e[h:(n + h - 1)] + e[2:(n + 1)] + e[1:n]
    mortensor:::dm_stat(d, h, varestimator = "ewc")$p_value
  }, 0)
  size <- mean(p < 0.05)
  expect_gt(size, 0.02)
  expect_lt(size, 0.10)
})

test_that("mt_dm_test(varestimator = \"acf\") equals forecast::dm.test() on the losses per origin", {
  skip_if_not_installed("forecast")
  bt <- fake_backtest()
  for (alt in c("two.sided", "less", "greater")) {
    for (ve in c("acf", "bartlett")) {
      for (m in c("crps", "se")) {
        r <- mt_dm_test(bt, metrics = m, alternative = alt, varestimator = ve)
        expect_identical(levels(r$horizon), c("1", "2", "3", "all"))
        expect_true(all(r$model == "b" & r$ref == "a"))
        for (h in c(1:3, NA)) {
          e1 <- origin_losses(bt, "b", m, if (is.na(h)) NULL else h)
          e2 <- origin_losses(bt, "a", m, if (is.na(h)) NULL else h)
          ref <- forecast::dm.test(e1, e2, h = if (is.na(h)) 3 else h, power = 1, alternative = alt, varestimator = ve)
          row <- r[r$horizon == if (is.na(h)) "all" else h, ]
          expect_equal(row$statistic, unname(ref$statistic), tolerance = 1e-10)
          expect_equal(row$p_value, unname(ref$p.value), tolerance = 1e-10)
          expect_equal(row$mean_diff, mean(e1 - e2), tolerance = 1e-12)
          expect_equal(row$score, mean(e1), tolerance = 1e-12)
          expect_identical(row$n_origins, 20L)
          expect_identical(row$lag, as.integer(if (is.na(h)) 2 else h - 1))
          expect_equal(row$df, 19)
        }
      }
    }
  }
})

test_that("the acf version uses lags by the spacing of origins and gives no test with too few origins", {
  skip_if_not_installed("forecast")
  # origins 2 years apart: forecasts 3 years ahead overlap at one lag
  bt <- fake_backtest(origins = seq(1970, 2010, by = 2), horizons = c(1, 3))
  r <- mt_dm_test(bt, varestimator = "acf")
  expect_identical(r$lag, c(0L, 1L, 1L))
  e1 <- origin_losses(bt, "b", h = 3)
  e2 <- origin_losses(bt, "a", h = 3)
  expect_equal(r$statistic[2], unname(forecast::dm.test(e1, e2, h = 2, power = 1)$statistic), tolerance = 1e-10)
  # 4 origins and horizons up to 5: no test where the lags cannot be estimated
  bt <- fake_backtest(origins = 2000:2003, horizons = 1:5)
  expect_warning(r <- mt_dm_test(bt, varestimator = "acf"), class = "mt_warning_dm")
  expect_identical(r$lag, c(0L, 1L, 2L, NA, NA, NA))
  expect_true(all(is.na(r$p_value[4:6])))
  expect_false(anyNA(r$p_value[1:3]))
  # the EWC test does not use the lags: a test at every horizon
  expect_no_warning(r <- mt_dm_test(bt))
  expect_false(anyNA(r$p_value))
  expect_true(all(r$df == 1))
  # identical losses: no test
  bt <- fake_backtest(origins = 2000:2009, horizons = 1)
  bt[model == "b", crps := bt[model == "a", crps]]
  for (ve in c("ewc", "acf")) {
    expect_warning(r <- mt_dm_test(bt, varestimator = ve), class = "mt_warning_dm")
    expect_true(all(is.na(r$statistic)))
  }
})

test_that("the acf version gives no test instead of lag 0 when the variance is not positive", {
  skip_if_not_installed("forecast")
  # alternating loss differences: the autocovariance at lag 1 is negative
  # enough to make the variance estimate with lags up to 1 negative
  bt <- fake_backtest(origins = 2000:2011, horizons = 2, n_age = 1)
  alt <- rep(c(0.02, -0.01), 6)
  bt[model == "b", crps := bt[model == "a", crps] + alt]
  e1 <- origin_losses(bt, "b", h = 2)
  e2 <- origin_losses(bt, "a", h = 2)
  expect_warning(forecast::dm.test(e1, e2, h = 2, power = 1), "negative")
  expect_warning(r <- mt_dm_test(bt, varestimator = "acf"), class = "mt_warning_dm")
  expect_true(is.na(r[horizon == "2"]$statistic))
  expect_true(is.na(r[horizon == "2"]$p_value))
  expect_true(is.na(r[horizon == "2"]$lag))
  # the EWC test is defined
  r <- mt_dm_test(bt)
  expect_false(anyNA(r$p_value))
})

test_that("mt_dm_test() pairs cells, splits by strata and identifies models", {
  bt <- fake_backtest(models = c("a", "b", "c"), horizons = 1:2)
  # cells scored for one model only are left out
  drop <- which(bt$model == "b" & bt$age_lo == 60 & bt$origin < 2005)
  r <- mt_dm_test(bt[-drop], ref = "c")
  expect_identical(unique(r$model), c("a", "b"))
  e1 <- origin_losses(bt[-drop], "b", h = 1)
  e2 <- origin_losses(bt[!(age_lo == 60 & origin < 2005)], "c", h = 1)
  expect_equal(r[model == "b" & horizon == "1", statistic], ewc_by_hand(e1 - e2)$statistic, tolerance = 1e-10)
  # by stratum variable
  b2 <- data.table::copy(bt)[, `:=`(sex = "female", stratum = 2L)]
  r2 <- mt_dm_test(rbind(bt, b2), by = "sex", metrics = c("crps", "wis"))
  expect_identical(nrow(r2), 2L * 2L * 2L * 3L)
  expect_equal(r2[model == "b" & sex == "male" & metric == "crps"]$statistic,
               mt_dm_test(bt)[model == "b"]$statistic, tolerance = 1e-12)
  # the same model with two pipelines
  b3 <- fake_backtest(models = c("lc", "lc"))
  b3[, prep := rep(c("none", "smooth"), each = .N / 2)]
  r3 <- mt_dm_test(b3)
  expect_identical(unique(r3$model), "lc (smooth)")
  expect_identical(unique(r3$ref), "lc (none)")
  # errors
  expect_error(mt_dm_test(bt, metrics = "cov95"), class = "mt_error_argument")
  expect_error(mt_dm_test(bt, ref = "z"), class = "mt_error_argument")
  expect_error(mt_dm_test(bt, by = "origin"), class = "mt_error_argument")
  expect_error(mt_dm_test(bt, targets = "asr"), class = "mt_error_argument")
  expect_error(mt_dm_test(bt[model == "a"]), class = "mt_error_argument")
  expect_error(mt_dm_test(rbind(bt, bt)), class = "mt_error_argument")
  expect_error(mt_dm_test(list()), class = "mt_error_argument")
})

test_that("mt_dm_test() works on a backtest of two models", {
  skip_if_not_installed("forecast")
  x <- mt_example(deaths = "single", years = 1990:2019)
  bt <- mt_compare(list(lc = lee_carter(), rw = rwd()), x, origins = 2005:2010, h = 5, n_draws = 50,
                   quantities = "life_expectancy")
  r <- mt_dm_test(bt, metrics = c("crps", "ae"))
  expect_s3_class(r, "data.table")
  expect_named(r, c("model", "ref", "target", "metric", "horizon", "n_origins", "lag", "df", "score", "score_ref",
                    "mean_diff", "statistic", "p_value"))
  expect_identical(nrow(r), 2L * 2L * 6L)
  expect_true(all(r$model == "rw" & r$ref == "lc" & r$n_origins == 6L & r$df == 1))
  expect_true(all(r$p_value >= 0 & r$p_value <= 1))
  # the acf version as forecast::dm.test() on the losses per origin, except
  # where that falls back to lag 0 (variance estimate not positive): no test
  ra <- suppressWarnings(mt_dm_test(bt, metrics = c("crps", "ae"), varestimator = "acf"))
  b <- as.data.frame(bt)
  b <- b[b$target == "rate", ]
  for (h in 1:5) {
    e1 <- origin_losses(b, "rw", h = h)
    e2 <- origin_losses(b, "lc", h = h)
    fell_back <- FALSE
    ref <- withCallingHandlers(forecast::dm.test(e1, e2, h = h, power = 1), warning = function(w) {
      fell_back <<- TRUE
      invokeRestart("muffleWarning")
    })
    row <- ra[target == "rate" & metric == "crps" & horizon == h]
    if (fell_back) {
      expect_true(is.na(row$p_value))
    } else {
      expect_equal(row$statistic, unname(ref$statistic), tolerance = 1e-10)
      expect_equal(row$p_value, unname(ref$p.value), tolerance = 1e-10)
    }
  }
  # mean scores agree with mt_score()
  s <- mt_score(bt, by = c("model", "horizon"))
  expect_equal(r[target == "rate" & metric == "crps" & horizon == "1", score_ref],
               s[model == "lc" & horizon == 1, crps], tolerance = 1e-12)
})
