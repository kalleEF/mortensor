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

test_that("mt_dm_test() equals forecast::dm.test() on the losses per origin", {
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
        }
      }
    }
  }
})

test_that("mt_dm_test() uses lags by the spacing of origins and gives no test with too few origins", {
  skip_if_not_installed("forecast")
  # origins 2 years apart: forecasts 3 years ahead overlap at one lag
  bt <- fake_backtest(origins = seq(1970, 2010, by = 2), horizons = c(1, 3))
  r <- mt_dm_test(bt)
  expect_identical(r$lag, c(0L, 1L, 1L))
  e1 <- origin_losses(bt, "b", h = 3)
  e2 <- origin_losses(bt, "a", h = 3)
  expect_equal(r$statistic[2], unname(forecast::dm.test(e1, e2, h = 2, power = 1)$statistic), tolerance = 1e-10)
  # 4 origins and horizons up to 5: no test where the lags cannot be estimated
  bt <- fake_backtest(origins = 2000:2003, horizons = 1:5)
  expect_warning(r <- mt_dm_test(bt), class = "mt_warning_dm")
  expect_identical(r$lag, c(0L, 1L, 2L, NA, NA, NA))
  expect_true(all(is.na(r$p_value[4:6])))
  expect_false(anyNA(r$p_value[1:3]))
  # identical losses: no test
  bt <- fake_backtest(origins = 2000:2009, horizons = 1)
  bt[model == "b", crps := bt[model == "a", crps]]
  expect_warning(r <- mt_dm_test(bt), class = "mt_warning_dm")
  expect_true(all(is.na(r$statistic)))
})

test_that("mt_dm_test() pairs cells, splits by strata and identifies models", {
  skip_if_not_installed("forecast")
  bt <- fake_backtest(models = c("a", "b", "c"), horizons = 1:2)
  # cells scored for one model only are left out
  drop <- which(bt$model == "b" & bt$age_lo == 60 & bt$origin < 2005)
  r <- mt_dm_test(bt[-drop], ref = "c")
  expect_identical(unique(r$model), c("a", "b"))
  e1 <- origin_losses(bt[-drop], "b", h = 1)
  e2 <- origin_losses(bt[!(age_lo == 60 & origin < 2005)], "c", h = 1)
  expect_equal(r[model == "b" & horizon == "1", statistic], unname(forecast::dm.test(e1, e2, power = 1)$statistic),
               tolerance = 1e-10)
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
  r <- suppressWarnings(mt_dm_test(bt, metrics = c("crps", "ae")))
  expect_s3_class(r, "data.table")
  expect_named(r, c("model", "ref", "target", "metric", "horizon", "n_origins", "lag", "score", "score_ref",
                    "mean_diff", "statistic", "p_value"))
  expect_identical(nrow(r), 2L * 2L * 6L)
  expect_true(all(r$model == "rw" & r$ref == "lc" & r$n_origins == 6L))
  expect_true(all(r$p_value >= 0 & r$p_value <= 1))
  # the same as forecast::dm.test() on the losses per origin (both fall back
  # to h = 1 when the variance estimate is not positive)
  b <- as.data.frame(bt)
  b <- b[b$target == "rate", ]
  for (h in 1:5) {
    e1 <- origin_losses(b, "rw", h = h)
    e2 <- origin_losses(b, "lc", h = h)
    ref <- suppressWarnings(forecast::dm.test(e1, e2, h = h, power = 1))
    row <- r[target == "rate" & metric == "crps" & horizon == h]
    expect_equal(row$statistic, unname(ref$statistic), tolerance = 1e-10)
    expect_equal(row$p_value, unname(ref$p.value), tolerance = 1e-10)
  }
  # mean scores agree with mt_score()
  s <- mt_score(bt, by = c("model", "horizon"))
  expect_equal(r[target == "rate" & metric == "crps" & horizon == "1", score_ref],
               s[model == "lc" & horizon == 1, crps], tolerance = 1e-12)
})
