test_that("every model satisfies the contract", {
  expect_mt_model(rwd())
  expect_mt_model(rwd(innovations = "gaussian", drift_window = 10))
  expect_mt_model(suppressWarnings(rwd(observation_noise = TRUE)))
  expect_mt_model(suppressWarnings(rwd(observation_noise = TRUE, innovations = "gaussian")))
  expect_mt_model(lee_carter())
  expect_mt_model(lee_carter(method = "svd", jump_off = "observed"))
  expect_mt_model(tucker3())
  expect_mt_model(tucker3(centre = FALSE))
  expect_mt_model(suppressWarnings(tucker3(transform = "coda")))
  expect_mt_model(ensemble(lee_carter(), rwd()))
  skip_if_not_installed("forecast")
  expect_mt_model(ets_series(), h = 2, n_draws = 5)
  expect_mt_model(arima_series(), h = 2, n_draws = 5)
  expect_mt_model(arima_series(order = c(0, 1, 1), drift = TRUE), h = 2, n_draws = 5)
  expect_mt_model(tucker3(ts = "ets"), h = 2, n_draws = 5)
})

test_that("models need grid-resolution data and enough periods", {
  x <- mt_example(years = 2001:2010)
  expect_error(mt_fit(lee_carter(), x), class = "mt_error_model_data")
  fit <- mt_fit(lee_carter(), x, prep = mt_prep(step_ungroup()))
  expect_s3_class(fit, "mt_fit")
  expect_error(mt_fit(lee_carter(), mt_example(deaths = "single", years = 2001:2003)), class = "mt_error_model_data")
})

test_that("a failing stratum is recorded and the others are kept", {
  x <- mt_example(deaths = "single", years = 2001:2010)
  bad <- new_mt_spec("failing", min_periods = 2)
  local_mocked_bindings(fit_model = function(spec, D, E, ages, periods, ...) {
    if (D[1, 1] == mt_array(x, "deaths")[1, 1, 2]) stop("boom")
    list(last = log(D[, ncol(D)] / E[, ncol(E)]))
  })
  registerS3method("simulate_model", "mt_spec_failing", function(spec, fit, h, n_draws, ...) {
    array(exp(fit$last), c(length(fit$last), h, n_draws))
  }, envir = asNamespace("mortensor"))
  expect_warning(fit <- mt_fit(bad, x), class = "mt_warning_model_fit")
  expect_equal(fit$diagnostics$status, c("ok", "failed"))
  fc <- mt_forecast(fit, h = 2, n_draws = 3)
  expect_true(all(is.na(fc$rates[, , 2, ])))
  expect_false(anyNA(fc$rates[, , 1, ]))
})

test_that("Poisson Lee-Carter agrees with StMoMo", {
  skip_if_not_installed("StMoMo")
  ew <- StMoMo::EWMaleData
  ages <- 55:89
  yrs <- 1961:2011
  D <- ew$Dxt[as.character(ages), as.character(yrs)]
  E <- ew$Ext[as.character(ages), as.character(yrs)]
  mine <- mortensor:::lc_poisson(D, E)
  ref <- StMoMo::fit(StMoMo::lc(), data = ew, ages.fit = ages, years.fit = yrs, verbose = FALSE)
  expect_lt(max(abs(mine$a - ref$ax)), 1e-4)
  expect_lt(max(abs(mine$b - ref$bx)), 1e-5)
  expect_lt(max(abs(mine$k - ref$kt[1, ])), 1e-3)
})

test_that("SVD Lee-Carter reproduces the fitted rates of demography::lca", {
  skip_if_not_installed("StMoMo")
  skip_if_not_installed("demography")
  ew <- StMoMo::EWMaleData
  ages <- 55:89
  yrs <- 1961:2011
  D <- ew$Dxt[as.character(ages), as.character(yrs)]
  E <- ew$Ext[as.character(ages), as.character(yrs)]
  dd <- demography::demogdata(data = D / E, pop = E, ages = ages, years = yrs, type = "mortality",
                              label = "EW", name = "male")
  l <- demography::lca(dd, adjust = "dt", interpolate = FALSE)
  s <- mortensor:::lc_svd(D, E)
  mine <- s$a + outer(s$b, s$k)
  ref <- l$ax + outer(l$bx, as.numeric(l$kt))
  expect_lt(max(abs(mine - ref)), 1e-5)
})

test_that("Lee-Carter constraints hold and parameters can be extracted", {
  fit <- mt_fit(lee_carter(), mt_example(deaths = "single", years = 2001:2015))
  p <- mt_parameters(fit)
  for (s in c("female", "male")) {
    expect_equal(sum(p[stratum == s & parameter == "b"]$value), 1, tolerance = 1e-10)
    expect_equal(sum(p[stratum == s & parameter == "k"]$value), 0, tolerance = 1e-8)
  }
  boot <- mt_fit(lee_carter(bootstrap = TRUE, n_boot = 5), mt_example(deaths = "single", years = 2001:2015))
  expect_length(boot$fits[[1]]$boot, 5)
  expect_mt_model(lee_carter(bootstrap = TRUE, n_boot = 3), n_draws = 6)
})

test_that("HOOI recovers an exact low-rank array and CoDa round trips", {
  set.seed(1)
  U <- list(qr.Q(qr(matrix(rnorm(20), 10, 2))), qr.Q(qr(matrix(rnorm(16), 8, 2))), qr.Q(qr(matrix(rnorm(6), 3, 2))))
  G <- array(rnorm(8), c(2, 2, 2))
  X <- mortensor:::tucker_reconstruct(G, U)
  dec <- mortensor:::hooi(X, c(2, 2, 2))
  expect_lt(max(abs(mortensor:::tucker_reconstruct(dec$G, dec$U) - X)), 1e-8)
  ages <- single_ages(30, 100)
  m <- array(rep(3e-4 + 2e-5 * exp(0.1 * (30:100)), 6), c(71, 3, 2))
  z <- mortensor:::coda_forward(m, ages)
  back <- mortensor:::coda_backward(z, ages, open_rates = m[71, , ])
  expect_lt(max(abs(back[-71, , ] / m[-71, , ] - 1)), 0.02)
  closed <- single_ages(30, 99, open = FALSE)
  z2 <- mortensor:::coda_forward(m[1:70, , , drop = FALSE], closed)
  expect_equal(dim(z2)[1], 71)
  back2 <- mortensor:::coda_backward(z2, closed)
  expect_lt(max(abs(back2 / m[1:70, , ] - 1)), 0.02)
})

test_that("tucker3 works with one stratum", {
  x <- mt_example(deaths = "single", years = 2001:2015)
  x1 <- x
  x1$obs <- lapply(x1$obs, function(o) o[o$stratum == 1])
  x1$strata <- x1$strata[1]
  expect_mt_model(tucker3(), data = x1)
})

test_that("ensemble mixes component draws and drops failed components", {
  x <- mt_example(deaths = "single", years = 2001:2015)
  ens <- ensemble(a = lee_carter(), b = rwd(), weights = c(3, 1))
  fit <- mt_fit(ens, x)
  expect_equal(unname(fit$weights), c(0.75, 0.25))
  fc <- mt_forecast(fit, h = 2, n_draws = 400, seed = 1)
  la <- mt_forecast(fit$components$a, h = 2, n_draws = 400, seed = 2)
  lb <- mt_forecast(fit$components$b, h = 2, n_draws = 400, seed = 3)
  # the mixture is at least as wide as the narrower component
  w <- function(f) diff(stats::quantile(log(f$rates[40, 2, 1, ]), c(0.05, 0.95)))
  expect_gte(w(fc), min(w(la), w(lb)) * 0.9)
  expect_error(ensemble(1), class = "mt_error_argument")
  bad <- ensemble(lee_carter(), tucker3(ranks = c(3, 3, 2)), rwd())
  short <- mt_example(deaths = "single", years = 2001:2005)
  fit2 <- mt_fit(bad, short)
  expect_true(any(fit2$diagnostics$status == "failed") || length(fit2$components) == 3)
})

test_that("forecast draws are summarised and life expectancy is computed path by path", {
  x <- mt_example(deaths = "single", years = 2001:2015)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 50, seed = 1)
  s <- mt_summarise(fc, probs = c(0.1, 0.9))
  expect_equal(nrow(s), 71 * 3 * 2)
  expect_true(all(c("q0.1", "q0.9", "mean", "sex") %in% names(s)))
  e <- life_expectancy(fc, at = 65)
  expect_equal(nrow(e), 6)
  expect_true(all(e$q0.025 < e$q0.975))
  expect_true(all(e$mean > 15 & e$mean < 30))
  d <- mt_deaths(fc, seed = 1)
  expect_equal(dim(d), dim(fc$rates))
  p <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 50, seed = 1, predictive = TRUE)
  expect_true(p$predictive)
  expect_output(print(fc), "mt_draws")
})

test_that("the random walk plus noise separates noise from the trend", {
  set.seed(1)
  na <- 5
  nt <- 30
  level <- outer(log(seq(0.01, 0.05, length.out = na)), rep(1, nt)) + outer(rep(1, na), cumsum(rep(-0.02, nt)))
  E <- matrix(1e4, na, nt)
  D <- matrix(stats::rpois(na * nt, exp(level) * E), na, nt)
  expect_warning(sp <- rwd(observation_noise = TRUE), class = "mt_warning_deprecated")
  f <- fit_model(sp, D, E, seq_len(na), seq_len(nt))
  # the level is deterministic here, so most of the variation of the changes is Poisson noise
  expect_true(all(f$s2e > 0.8 * pmin(rowMeans(1 / D), 0.45 * f$sigma^2)))
  expect_true(all(f$s2eta < 0.5 * f$sigma^2))
  p_noise <- simulate_model(sp, f, h = 10, n_draws = 400)
  p_plain <- simulate_model(rwd(), fit_model(rwd(), D, E, seq_len(na), seq_len(nt)), h = 10, n_draws = 400)
  width <- function(p) mean(apply(log(p[, 10, ]), 1, stats::sd))
  expect_lt(width(p_noise), width(p_plain))
  expect_error(rwd(observation_noise = NA), class = "mt_error_argument")
})

test_that("ensemble draws come from one component for all strata", {
  x <- mt_example_causes(years = 2008:2019)
  h <- mt_hierarchy(cause = "all", additive = "rates")
  ens <- ensemble(a = reconciled(lee_carter(), h, method = "bottom_up"),
                     b = reconciled(rwd(), h, method = "bottom_up"))
  fc <- mt_forecast(mt_fit(ens, x), h = 2, n_draws = 30, seed = 1)
  st <- mt_strata(x)
  for (sx in c("female", "male")) {
    i <- st$stratum[st$sex == sx & st$cause == "all"]
    j <- st$stratum[st$sex == sx & st$cause != "all"]
    parts <- apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum)
    expect_lt(max(abs(parts / fc$rates[, , i, ] - 1)), 1e-8)
  }
})

test_that("every kind of model can be an ensemble component", {
  x <- mt_example(deaths = "single", years = 2000:2019)
  comps <- list(tucker3(ranks = c(2, 2, 2)), li_lee(common_by = "sex"), suppressWarnings(product_ratio()), ets_series(),
                bootstrapped(rwd(), n_boot = 2), ensemble(lee_carter(), rwd()),
                tuned(function(w) rwd(drift_window = w), data.frame(w = c(5, 10)), inner_origins = 2, inner_h = 2))
  if (requireNamespace("demography", quietly = TRUE)) comps <- c(comps, list(hyndman_ullah(order = 2)))
  if (requireNamespace("StMoMo", quietly = TRUE)) comps <- c(comps, list(stmomo("cbd")))
  for (cp in comps) {
    fc <- mt_forecast(mt_fit(ensemble(lc = lee_carter(), other = cp), x), h = 2, n_draws = 20, seed = 1,
                      predictive = TRUE)
    expect_true(all(is.finite(fc$rates)), label = cp$label)
  }
})

test_that("parameters, fitted values, residuals and summaries are available", {
  x <- mt_example(deaths = "single", years = 2001:2015)
  fit <- mt_fit(li_lee(common_by = NULL), x)
  p <- mt_parameters(fit)
  expect_true(all(c("B", "K", "a", "b1", "k1", "phi1") %in% p$parameter))
  pr <- mt_parameters(mt_fit(suppressWarnings(rwd(observation_noise = TRUE)), x))
  expect_true(all(c("drift", "sd", "var_level", "var_noise") %in% pr$parameter))
  f <- fitted(fit)
  expect_identical(dim(f), dim(fit$train_D))
  r <- residuals(fit)
  expect_true(all(c("stratum", "age", "period", "residual") %in% names(r)))
  expect_lt(abs(mean(r$residual)), 1)
  expect_output(print(summary(fit)), "Li-Lee")
  expect_output(print(summary(mt_fit(lee_carter(), x))), "Poisson deviance")
  expect_error(mt_parameters(mt_fit(tucker3(), x)), class = "mt_error_argument")
})
