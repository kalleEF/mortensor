h3 <- function() mt_hierarchy(sex = "both", region = "all", quintile = "all")

test_that("hierarchies are defined and printed", {
  h <- mt_hierarchy(sex = "both", region = list(total = "DE", groups = list(East = c("BE", "SN"))))
  expect_s3_class(h, "mt_hierarchy")
  expect_output(print(h), "groups East")
  expect_error(mt_hierarchy("both"), class = "mt_error_hierarchy")
  expect_error(mt_hierarchy(sex = 1), class = "mt_error_hierarchy")
})

test_that("aggregate strata are sums of the bottom strata", {
  x <- mt_example_regions(years = 2010:2013, regions = 3)
  xa <- mt_aggregate_strata(x, h3())
  st <- mt_strata(xa)
  expect_equal(nrow(st), 3 * 4 * 6)
  expect_equal(st$stratum, seq_len(nrow(st)))
  D <- mt_array(xa, "deaths")
  expect_equal(sum(D[, , "both.all.all"]), sum(mt_obs(x, "deaths")$value))
  expect_equal(D[, , "female.all.Q1"], apply(D[, , paste0("female.R", 1:3, ".Q1")], 1:2, sum))
  v <- mt_validate(xa, hierarchy = h3())
  expect_false(any(v$severity == "error"))
  # an inconsistent observed total is flagged
  bad <- xa
  bad$obs$deaths <- data.table::copy(bad$obs$deaths)
  tot <- st$stratum[st$label == "both.all.all"]
  bad$obs$deaths[stratum == tot, value := value * 2]
  vb <- mt_validate(bad, hierarchy = h3())
  expect_true(any(vb$severity == "error" & grepl("median 2.000", vb$detail)))
  # intermediate groups
  hg <- mt_hierarchy(region = list(total = "all", groups = list(north = c("R1", "R2"))))
  xg <- mt_aggregate_strata(x, hg)
  expect_true(all(c("female.north.Q3", "male.all.Q5") %in% mt_strata(xg)$label))
  expect_identical(mt_aggregate_strata(xa, h3())$strata, xa$strata)
})

test_that("reconciliation makes all nodes add up for every method", {
  x <- mt_example_regions(years = 2005:2016, regions = 2)
  xa <- mt_aggregate_strata(x, h3())
  fc <- mt_forecast(mt_fit(lee_carter(), xa), h = 2, n_draws = 40, seed = 1)
  expect_gt(coherence_gap(fc, h3()), 0.01)
  for (m in c("wls", "mint", "structural", "ols", "bottom_up")) {
    r <- mt_reconcile(fc, h3(), method = m, seed = 1)
    expect_lt(coherence_gap(r, h3()), 1e-8)
    expect_true(all(r$rates > 0))
    expect_equal(r$reconciliation$method, m)
  }
  r0 <- mt_reconcile(fc, h3(), dependence = "none")
  expect_lt(coherence_gap(r0, h3()), 1e-8)
  expect_error(mt_reconcile(fc, mt_hierarchy(cause = "all")), class = "mt_error_argument")
})

test_that("the empirical copula keeps each stratum's marginal distribution", {
  x <- mt_example_regions(years = 2005:2016, regions = 2)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 2, n_draws = 50, seed = 1)
  re <- with_seed_test(1, mortensor:::copula_reorder(fc$rates, fc$train_D, fc$train_E))
  expect_equal(sort(re$arr[3, 2, 5, ]), sort(fc$rates[3, 2, 5, ]))
  expect_length(re$years, 50)
})

test_that("reconciled models give coherent point and predictive draws and run in backtests", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  sp <- reconciled(lee_carter(), h3(), method = "wls")
  expect_s3_class(sp, "mt_spec_reconciled")
  fit <- mt_fit(sp, x)
  expect_equal(nrow(fit$strata), 2 * 3 + 3 * 3 * 6 - 2 * 3 * 1 - 0)
  fc <- mt_forecast(fit, h = 2, n_draws = 30, seed = 2)
  expect_lt(coherence_gap(fc, h3()), 1e-8)
  pr <- mt_forecast(fit, h = 2, n_draws = 30, seed = 2, predictive = TRUE)
  expect_lt(coherence_gap(pr, h3()), 1e-8)
  xa <- mt_aggregate_strata(x, h3())
  bt <- mt_backtest(reconciled(li_lee(common_by = "sex"), h3(), method = "bottom_up"), xa, origins = 2011,
                    h = 3, n_draws = 30)
  expect_true(all(c("both", "female") %in% bt$sex))
})

test_that("hierarchical Li-Lee pools bottom strata within each sex", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  xa <- mt_aggregate_strata(x, h3())
  fit <- mt_fit(li_lee(common_by = "sex"), xa)
  expect_equal(length(fit$fits$commons), 3)
  expect_error(mt_fit(li_lee(common_by = "nonsense"), xa), class = "mt_error_model_fit")
  expect_mt_model(li_lee(common_by = "sex"), data = x, n_draws = 10)
})

test_that("draws are aggregated over strata and ages", {
  x <- mt_example(deaths = "single", years = 2000:2012)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 2, n_draws = 10, seed = 1)
  fa <- mt_aggregate_draws(fc, hierarchy = mt_hierarchy(sex = "both"), age = age_groups(c(30, seq(40, 100, 10))))
  expect_equal(dim(fa$rates), c(8, 2, 3, 10))
  E <- fc$exposure_last[, 1, ]
  both <- (fc$rates[1:10, 1, 1, 1] %*% E[1:10, 1] + fc$rates[1:10, 1, 2, 1] %*% E[1:10, 2]) / sum(E[1:10, ])
  expect_equal(fa$rates[1, 1, "both", 1], as.numeric(both))
  expect_error(mt_aggregate_draws(fc, age = age_groups(c(30, 35.5))), class = "mt_error_argument")
})

test_that("rate-additive hierarchies (causes) reconcile on the rate scale, also in predictive draws", {
  x <- mt_example_causes(years = 2005:2019)
  h <- mt_hierarchy(cause = "all", additive = "rates")
  fit <- mt_fit(reconciled(lee_carter(), h, method = "wls"), x)
  base <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 50, seed = 1)
  st <- mt_strata(x)
  for (pred in c(FALSE, TRUE)) {
    fc <- mt_forecast(fit, h = 3, n_draws = 50, seed = 1, predictive = pred)
    for (sx in c("female", "male")) {
      i <- st$stratum[st$sex == sx & st$cause == "all"]
      j <- st$stratum[st$sex == sx & st$cause != "all"]
      parts <- apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum)
      expect_lt(max(abs(parts / fc$rates[, , i, ] - 1)), 1e-8)
      ratio <- stats::median(fc$rates[, , i, ] / base$rates[, , i, ])
      expect_gt(ratio, 0.8)
      expect_lt(ratio, 1.25)
    }
  }
})
