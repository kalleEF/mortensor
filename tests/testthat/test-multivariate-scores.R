test_that("the energy score equals the CRPS in one dimension", {
  set.seed(1)
  X <- matrix(rnorm(200), 1)
  y <- 0.3
  expect_equal(mortensor:::energy_score(X, y), score_samples(X, y)$crps)
})

test_that("energy and variogram scores agree with scoringRules", {
  skip_if_not_installed("scoringRules")
  set.seed(2)
  X <- matrix(rnorm(4 * 150), 4)
  X[2, ] <- X[1, ] + rnorm(150, 0, 0.2)
  y <- c(0.1, -0.4, 1.2, 0.5)
  expect_equal(mortensor:::energy_score(X, y), scoringRules::es_sample(y, X))
  expect_equal(mortensor:::variogram_score(X, y), scoringRules::vs_sample(y, X, p = 0.5))
})

test_that("the variogram score rewards the observed dependence between strata", {
  # two strata that move together: draws with the same dependence score
  # better on average than independent draws with the same margins
  set.seed(3)
  sc <- replicate(200, {
    z <- rnorm(1)
    y <- c(z, z + rnorm(1, 0, 0.1))
    a <- rnorm(300)
    dep <- rbind(a, a + rnorm(300, 0, 0.1))
    ind <- rbind(rnorm(300), rnorm(300, 0, sqrt(1.01)))
    c(mortensor:::variogram_score(dep, y), mortensor:::variogram_score(ind, y))
  })
  expect_lt(mean(sc[1, ]), mean(sc[2, ]))
})

test_that("mt_backtest(multivariate = TRUE) scores all strata jointly", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  bt <- mt_backtest(lee_carter(), x, origins = 2007, h = 3, n_draws = 50, quantities = c("life_expectancy", "asr"),
                    multivariate = TRUE)
  jt <- bt[target == "life_expectancy_joint"]
  expect_equal(nrow(jt), 3L)
  expect_true(all(is.na(jt$stratum)))
  expect_equal(jt$n_dim, rep(nrow(x$strata), 3))
  expect_true(all(jt$energy > 0) && all(jt$variogram >= 0))
  expect_equal(nrow(bt[target == "asr_joint"]), 3L)
  # the default backtest has no joint rows, and the other rows are unchanged
  bt0 <- mt_backtest(lee_carter(), x, origins = 2007, h = 3, n_draws = 50, quantities = c("life_expectancy", "asr"))
  expect_false(any(grepl("_joint", bt0$target)))
  expect_equal(bt[!grepl("_joint", target)]$crps, bt0$crps)
  s <- mt_score(bt, target = "asr_joint")
  expect_named(s, c("model", "energy", "variogram", "n", "n_dim", "unit"))
  tab <- mt_score_table(bt, metrics = c("crps", "energy"))
  expect_true(all(c("crps", "crps_le", "energy_le_joint", "energy_asr_joint") %in% names(tab)))
  expect_false("crps_le_joint" %in% names(tab))
  expect_false("energy_le_joint" %in% names(mt_score_table(bt)))
  expect_error(mt_score(bt, target = "asr_joint", measures = "cov90"), class = "mt_error_argument")
  expect_error(mt_backtest(lee_carter(), x, origins = 2007, h = 3, multivariate = TRUE), class = "mt_error_argument")
})
