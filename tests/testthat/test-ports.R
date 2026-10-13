test_that("wrapped and multi-population models satisfy the contract", {
  expect_mt_model(li_lee())
  expect_mt_model(product_ratio(ratio_order = 2))
  skip_if_not_installed("demography")
  skip_if_not_installed("forecast")
  expect_mt_model(hyndman_ullah(order = 2))
  skip_if_not_installed("StMoMo")
  for (m in c("apc", "cbd", "lc", "m6")) expect_mt_model(stmomo(m), n_draws = 10)
})

test_that("StMoMo cohort models need single years of age", {
  skip_if_not_installed("StMoMo")
  x <- mortensor:::simulate_lc_dgp(n_t = 12)
  fit <- tryCatch(mt_fit(stmomo("apc"), x), error = function(e) e)
  expect_true(inherits(fit, "mt_error_model_fit") || all(fit$diagnostics$status == "failed"))
})

test_that("Li-Lee forecasts of strata stay close together", {
  x <- mt_example(deaths = "single")
  fc <- mt_forecast(mt_fit(li_lee(), x), h = 30, n_draws = 200, seed = 1)
  ratio <- apply(fc$rates[40, , 2, ] / fc$rates[40, , 1, ], 1, median)
  obs <- mean(fc$observed[40, 26:30, 2] / fc$observed[40, 26:30, 1])
  expect_lt(abs(log(ratio[30] / obs)), 0.3)
})

test_that("cause shares add up to all causes in every draw", {
  x <- mt_example_causes(years = 2000:2015)
  fit <- mt_fit(cause_shares(total = lee_carter(), ranks = c(2, 2, 2, 2)), x)
  fc <- mt_forecast(fit, h = 3, n_draws = 20, seed = 1)
  st <- fc$strata
  for (sx in c("female", "male")) {
    i <- st$stratum[st$sex == sx & st$cause == "all"]
    j <- st$stratum[st$sex == sx & st$cause != "all"]
    total_parts <- apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum)
    expect_lt(max(abs(total_parts / fc$rates[, , i, ] - 1)), 1e-3)
  }
  p <- mt_forecast(fit, h = 2, n_draws = 10, seed = 1, predictive = TRUE)
  expect_false(anyNA(p$rates))
  bt <- mt_backtest(cause_shares(total = rwd(), ranks = c(2, 2, 2, 2)), x, origins = 2012, h = 3, n_draws = 30)
  expect_setequal(unique(bt$cause), c("all", "cvd", "cancer", "other"))
  expect_error(mt_fit(cause_shares(var = "region"), x), class = "mt_error_argument")
  expect_error(cause_shares(ranks = 1:3), class = "mt_error_argument")
})

test_that("ILR transforms round trip", {
  V <- mortensor:::helmert_basis(4)
  expect_equal(crossprod(V), diag(3), tolerance = 1e-12)
  expect_equal(colSums(V), rep(0, 3), tolerance = 1e-12)
  S <- array(c(0.1, 0.2, 0.3, 0.4), c(1, 1, 1, 4))
  expect_equal(mortensor:::ilr_inverse_array(mortensor:::ilr_array(S, V), V), S, tolerance = 1e-12)
})

test_that("reconciliation makes parts add up to the total", {
  x <- mt_example_causes(years = 2000:2015)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 2, n_draws = 30, seed = 1)
  st <- fc$strata
  i <- st$stratum[st$sex == "male" & st$cause == "all"]
  j <- st$stratum[st$sex == "male" & st$cause != "all"]
  for (m in c("mint", "wls")) {
    r <- mt_reconcile(fc, var = "cause", total = "all", method = m)
    expect_lt(max(abs(apply(r$rates[, , j, ], c(1, 2, 4), sum) / r$rates[, , i, ] - 1)), 1e-8)
    expect_gt(r$reconciliation$max_relative_gap_before, 0)
  }
  expect_error(mt_reconcile(fc, var = "region", total = "all"), class = "mt_error_argument")
  expect_error(mt_reconcile(fc, var = "cause", total = "none"), class = "mt_error_argument")
})

test_that("coherence keeps the lower stratum below the upper one in every draw", {
  x <- mt_example(deaths = "single")
  fc <- mt_forecast(mt_fit(rwd(), x), h = 10, n_draws = 50, seed = 1)
  co <- mt_order_strata(fc, "sex", lower = "female", upper = "male")
  expect_true(all(co$rates[, , 1, ] <= co$rates[, , 2, ] + 1e-15))
  expect_equal((co$rates[, , 1, ] + co$rates[, , 2, ]), (fc$rates[, , 1, ] + fc$rates[, , 2, ]))
  expect_s3_class(co$coherence, "data.table")
  expect_error(mt_order_strata(fc, "region", "a", "b"), class = "mt_error_argument")
})

test_that("stacking weights lie on the simplex", {
  w <- mortensor:::simplex_least_squares(cbind(1:10, (1:10) * 2), (1:10) * 1.5)
  expect_equal(sum(w), 1)
  expect_equal(w, c(0.5, 0.5), tolerance = 1e-3)
  x <- mt_example(deaths = "single", years = 1990:2010)
  fit <- mt_fit(ensemble(lee_carter(), rwd(), weights = "stack", inner_origins = 2, inner_h = 3), x)
  expect_equal(sum(fit$weights), 1)
})

test_that("Whittaker smoothing, Kannisto closure and monotone deaths work", {
  x <- mt_example(years = 2001:2003, deaths = "single", size = 2e4)
  w <- mt_smooth(x, method = "whittaker")
  truth <- attr(x, "truth")
  sm <- mt_array(w, "deaths") / mt_array(w, "exposure")
  raw <- pmax(mt_array(x, "deaths"), 0.5) / mt_array(x, "exposure")
  expect_lt(median(abs(log(sm / truth))), median(abs(log(raw / truth))))
  skip_if_not_installed("mgcv")
  gs <- mt_smooth(x, method = "gam")
  smg <- mt_array(gs, "deaths") / mt_array(gs, "exposure")
  expect_lt(median(abs(log(smg / truth))), median(abs(log(raw / truth))))
  m <- truth[, 1, 1]
  ek <- tail(life_table(m, mt_grid(x)$age, open = "kannisto")$ex, 1)
  eg <- tail(life_table(m, mt_grid(x)$age, open = "gompertz")$ex, 1)
  expect_gt(ek, eg)
  expect_warning(u <- mt_ungroup(mt_example(years = 2001:2002), method = "monotone"),
                 class = "mt_warning_ungroup_open_counts")
  expect_equal(sum(mt_obs(u)$value), sum(mt_obs(mt_example(years = 2001:2002))$value))
})

test_that("bootstrap refits pool draws and widen intervals", {
  x <- mt_example(deaths = "single", years = 2000:2019)
  for (type in c("grid", "data")) {
    fit <- mt_fit(bootstrapped(lee_carter(), n_boot = 4, type = type), x)
    expect_s3_class(fit, "mt_fit_bootstrap")
    expect_length(fit$boot, 4)
    fc <- mt_forecast(fit, h = 3, n_draws = 40, seed = 1)
    expect_identical(dim(fc$rates), c(71L, 3L, 2L, 40L))
    expect_true(all(is.finite(fc$rates) & fc$rates > 0))
  }
  fc1 <- mt_forecast(fit, h = 3, n_draws = 40, seed = 2)
  fc2 <- mt_forecast(fit, h = 3, n_draws = 40, seed = 2)
  expect_identical(fc1$rates, fc2$rates)
  xg <- mt_example(years = 2000:2019)
  fit <- mt_fit(bootstrapped(rwd(), n_boot = 2, type = "data"), xg, prep = mt_prep(step_ungroup()))
  expect_match(fit$diagnostics$message[3], "2 of 2")
  expect_error(bootstrapped(lee_carter(), n_boot = 0), class = "mt_error_argument")
})

test_that("Tucker3 ranks can be chosen by an information criterion", {
  x <- mt_example(deaths = "single", years = 2000:2019)
  X <- log(mt_array(x, "deaths") / mt_array(x, "exposure"))
  sel <- tucker_select_ranks(sweep(X, c(1, 3), apply(X, c(1, 3), mean)), "bic")
  expect_length(sel$ranks, 3)
  expect_true(all(is.finite(sel$table$ic) | is.infinite(sel$table$ic)))
  fit <- mt_fit(tucker3(ranks = "bic"), x)
  expect_match(fit$spec$label, "BIC")
  expect_error(tucker3(ranks = "dic"), class = "mt_error_argument")
  expect_s3_class(mt_forecast(fit, h = 2, n_draws = 10, seed = 1), "mt_draws")
})
