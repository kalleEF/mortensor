# WLS and MinT with weights from the errors of inner forecasts
# (Wickramasuriya, Athanasopoulos and Hyndman 2019; Schaefer and Strimmer 2005)

# The shrinkage estimator of MinT written out from its formula: the sample
# covariance (1/n) sum e e' of errors of mean zero, its correlations shrunk
# towards the identity with lambda = sum_{i != j} Var(r_ij) / sum_{i != j} r_ij^2
shrink_errors_ref <- function(e) {
  n <- nrow(e)
  covm <- crossprod(e) / n
  corm <- stats::cov2cor(covm)
  xs <- sweep(e, 2, sqrt(diag(covm)), "/")
  v <- (crossprod(xs^2) - crossprod(xs)^2 / n) / (n * (n - 1))
  diag(v) <- 0
  r2 <- corm^2
  diag(r2) <- 0
  lambda <- max(min(sum(v) / sum(r2), 1), 0)
  list(W = lambda * diag(diag(covm)) + (1 - lambda) * covm, lambda = lambda)
}

# The same for draws: centred, divisor n - 1 (Schaefer and Strimmer 2005, target D on the correlations)
shrink_draws_ref <- function(x) {
  n <- nrow(x)
  V <- stats::cov(x)
  z <- scale(x)
  w <- array(NA_real_, c(n, ncol(x), ncol(x)))
  for (k in seq_len(n)) w[k, , ] <- outer(z[k, ], z[k, ])
  wbar <- apply(w, 2:3, mean)
  var_r <- n / (n - 1)^3 * apply(sweep(w, 2:3, wbar)^2, 2:3, sum)
  r <- stats::cor(x)
  off <- row(r) != col(r)
  lambda <- max(min(sum(var_r[off]) / sum(r[off]^2), 1), 0)
  list(W = lambda * diag(diag(V)) + (1 - lambda) * V, lambda = lambda)
}

test_that("the shrinkage estimators follow their formulas", {
  set.seed(1)
  e <- matrix(rnorm(4 * 9), 4) %*% diag(1:9) + rnorm(4)
  ref <- shrink_errors_ref(e)
  sh <- mortensor:::shrink_cov(e, center = FALSE)
  expect_equal(sh$W, ref$W, tolerance = 1e-12)
  expect_equal(sh$lambda, ref$lambda, tolerance = 1e-12)
  x <- matrix(rnorm(30 * 5), 30) %*% matrix(runif(25), 5)
  rd <- shrink_draws_ref(x)
  sd <- mortensor:::shrink_cov(x)
  expect_equal(sd$W, rd$W, tolerance = 1e-12)
  expect_equal(sd$lambda, rd$lambda, tolerance = 1e-12)
  expect_equal(mortensor:::shrink_cov(e, 1, center = FALSE)$W, diag(colMeans(e^2)), tolerance = 1e-12)
})

test_that("MinT with error covariance uses the shrunk covariance of the inner forecast errors", {
  x <- mt_subset(mt_example_regions(years = 2000:2014, regions = 2), ages = c(70, Inf))
  h <- mt_hierarchy(region = "all", quintile = "all")
  xa <- mt_aggregate_strata(x, h)
  sp <- reconciled(lee_carter(), h, method = "mint", covariance = "errors", years = 6, n_origins = 3, horizon = 2,
                   control = list(n_draws = 30), seed = 7)
  expect_match(sp$label, "error covariance")
  fit <- mt_fit(sp, xa)
  er <- fit$reconcile$errors
  expect_equal(er$origins, c(2008, 2010, 2012))
  expect_true(all(er$n == 3))
  # the inner forecasts, made independently: the base model refitted at each
  # origin, with the random number streams of the training
  seeds <- mortensor:::with_seed(7, mortensor:::rng_streams(3))
  for (i in 1:3) {
    o <- er$origins[i]
    ref <- mortensor:::with_stream(seeds[[i]], {
      f <- mt_fit(lee_carter(), mt_window(xa, to = o))
      test <- mt_window(xa, from = o + 1, to = o + 2)
      E <- mt_array(test, "exposure")
      fc <- suppressWarnings(mt_forecast(f, to = o + 2, n_draws = 30, exposure = E))
      mt_array(test, "deaths") - apply(fc$rates * as.vector(E), 1:3, mean)
    })
    expect_equal(unname(er$errors[, , , i]), unname(ref), tolerance = 1e-10)
  }
  # W of each age and inner horizon: the shrinkage estimate of its errors
  for (a in c(1, 4)) for (k in 1:2) {
    ref <- shrink_errors_ref(t(er$errors[a, k, , ]))
    w <- mortensor:::errors_weights(er, a, k, "mint")
    expect_equal(w$W, ref$W, tolerance = 1e-10, ignore_attr = TRUE)
    expect_equal(w$lambda, ref$lambda, tolerance = 1e-10)
    # forecast years beyond the inner horizons take the weights of the last one
    if (k == 2) expect_identical(mortensor:::errors_weights(er, a, 5, "mint"), w)
  }
  # the forecast is the projection of the base draws with these weights
  fc <- mt_forecast(fit, h = 4, n_draws = 30, seed = 9)
  expect_equal(fc$reconciliation$share_set_to_zero, 0)
  expect_lt(coherence_gap(fc, h), 1e-8)
  base <- mortensor:::with_seed(9, mt_forecast(mortensor:::unwrap_fit(fit, "mt_fit_reconciled"), h = 4, n_draws = 30))
  hs <- fit$reconcile$structure
  E <- mortensor:::expand_exposure(fit$exposure_last, dim(base$rates), dimnames(base$rates), hs)
  agg <- which(!hs$is_bottom)
  U <- matrix(0, length(agg), nrow(hs$S))
  U[cbind(seq_along(agg), agg)] <- 1
  U[, hs$bottom] <- -hs$S[agg, ]
  for (a in c(1, 4)) for (t in c(1, 2, 4)) {
    W <- shrink_errors_ref(t(er$errors[a, min(t, 2), , ]))$W
    Y <- t(base$rates[a, t, , ] * E[a, t, ])
    X <- Y - t(W %*% t(U) %*% solve(U %*% W %*% t(U), U %*% t(Y)))
    expect_equal(t(fc$rates[a, t, , ] * E[a, t, ]), X, tolerance = 1e-8, ignore_attr = TRUE)
  }
  pr <- mt_forecast(fit, h = 4, n_draws = 30, seed = 9, predictive = TRUE)
  expect_lt(coherence_gap(pr, h), 1e-8)
  # with diagonal weights (shrinkage intensity 1), MinT equals WLS with the mean squared errors
  fw <- mt_fit(reconciled(lee_carter(), h, method = "wls", covariance = "errors", years = 6, n_origins = 3,
                          horizon = 2, control = list(n_draws = 30), seed = 7), xa)
  fm <- mt_fit(reconciled(lee_carter(), h, method = "mint", shrink = 1, covariance = "errors", years = 6,
                          n_origins = 3, horizon = 2, control = list(n_draws = 30), seed = 7), xa)
  expect_identical(fw$reconcile$errors$errors, er$errors)
  expect_equal(mortensor:::errors_weights(er, 2, 1, "wls"), colMeans(t(er$errors[2, 1, , ])^2), ignore_attr = TRUE)
  rm1 <- mt_forecast(fm, h = 4, n_draws = 30, seed = 9)
  expect_equal(rm1$rates, mt_forecast(fw, h = 4, n_draws = 30, seed = 9)$rates, tolerance = 1e-10)
  expect_equal(rm1$reconciliation$mean_shrinkage, 1)
  expect_match(mortensor:::model_assumptions(fit$spec), "errors of inner forecasts")
})

test_that("error covariances work with fixed nodes and shared exposures, and are checked", {
  h <- mt_hierarchy(cause = "all", additive = "rates")
  x <- mt_example_causes(years = 2003:2016)
  fit <- mt_fit(reconciled(lee_carter(), h, method = "mint", covariance = "errors", years = 6, n_origins = 3,
                           horizon = 2, control = list(n_draws = 30)), x)
  # errors of rates: observed minus forecast rates
  expect_lt(max(abs(fit$reconcile$errors$errors), na.rm = TRUE), 1)
  fc <- mt_forecast(fit, h = 3, n_draws = 30, seed = 1)
  st <- mt_strata(x)
  for (sx in c("female", "male")) {
    i <- st$stratum[st$sex == sx & st$cause == "all"]
    j <- st$stratum[st$sex == sx & st$cause != "all"]
    expect_lt(max(abs(apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum) / fc$rates[, , i, ] - 1)), 1e-8)
  }
  xr <- mt_example_regions(years = 2003:2016, regions = 2)
  hr <- mt_hierarchy(region = "all")
  fx <- mt_fit(reconciled(lee_carter(), hr, method = "mint", covariance = "errors", fixed = "total", years = 6,
                          n_origins = 3, horizon = 2, control = list(n_draws = 30)), xr)
  ff <- mt_forecast(fx, h = 2, n_draws = 30, seed = 1)
  expect_lt(ff$reconciliation$max_relative_change_fixed, 1e-6)
  expect_lt(coherence_gap(ff, hr), 1e-8)
  expect_error(reconciled(lee_carter(), hr, method = "ols", covariance = "errors"), class = "mt_error_argument")
  expect_error(reconciled(lee_carter(), hr, method = "mint", covariance = "errors", n_origins = 1),
               class = "mt_error_argument")
  expect_identical(reconciled(lee_carter(), hr)$params$covariance, "draws")
})
