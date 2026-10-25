# Observation model of predictive draws (Poisson given the rates, negative
# binomial for stanmomo(family = "nb")), the removed deprecated options,
# safeguards that warn, and the three-way Lee-Carter model

test_that("predictive draws are Poisson deaths given the rates of each draw", {
  x <- mt_example(deaths = "single", years = 2000:2019)
  fit <- mt_fit(lee_carter(), x)
  expect_null(fit$noise_sd)
  a <- mt_forecast(fit, h = 2, n_draws = 300, seed = 1)
  p <- mt_forecast(fit, h = 2, n_draws = 300, seed = 1, predictive = TRUE)
  expect_true(p$predictive)
  expect_null(p$nb_size)
  E <- mortensor:::expand_exposure(fit$exposure_last, dim(a$rates))
  mu <- a$rates * as.vector(E)
  D <- p$rates * as.vector(E)
  expect_equal(D, round(D), tolerance = 1e-8)
  # mean mu and variance mu: Pearson residuals have mean 0 and variance 1 (no
  # extra-Poisson noise, which would inflate the variance)
  z <- as.vector((D - mu) / sqrt(mu))
  expect_lt(abs(mean(z)), 0.02)
  expect_lt(abs(stats::var(z) - 1), 0.04)
  # models of observed rates add nothing
  rw <- mt_fit(rwd(), x)
  expect_equal(mt_forecast(rw, h = 2, n_draws = 20, seed = 1, predictive = TRUE)$rates,
               mt_forecast(rw, h = 2, n_draws = 20, seed = 1)$rates)
})

test_that("deaths are negative binomial with a finite size and Poisson otherwise", {
  set.seed(1)
  mu <- rep(50, 2e5)
  nb <- mortensor:::r_deaths(mu, 5)
  expect_equal(mean(nb), 50, tolerance = 0.01)
  expect_equal(stats::var(nb), 50 + 50^2 / 5, tolerance = 0.03)
  po <- mortensor:::r_deaths(mu, Inf)
  expect_equal(stats::var(po), 50, tolerance = 0.03)
  # sizes per stratum and draw of predictive_rates(): stratum 1 Poisson, stratum 2 size 2
  arr <- array(0.01, c(3, 2, 2, 2000))
  fit <- list(exposure_last = array(1e4, c(3, 1, 2)), smoothed_input = FALSE)
  size <- rbind(rep(Inf, 2000), rep(2, 2000))
  out <- mortensor:::predictive_rates(arr, fit, NULL, list(smooth = TRUE), size)
  D <- out * 1e4
  expect_equal(stats::var(as.vector(D[, , 1, ])), 100, tolerance = 0.1)
  expect_equal(stats::var(as.vector(D[, , 2, ])), 100 + 100^2 / 2, tolerance = 0.1)
  # wrappers that add the observation noise themselves use the sizes of the draws
  base <- list(spec = lee_carter(), fits = list(NULL, NULL), age = mt_axis(c(0, 1, 2), c(1, 2, 3)))
  fc <- list(nb_size = size)
  E <- array(1e4, c(3, 2, 2))
  noisy <- mortensor:::add_observation_noise(arr, base, fc, 1:2, E)
  expect_equal(stats::var(as.vector(noisy[, , 2, ] * 1e4)), 100 + 100^2 / 2, tolerance = 0.1)
  expect_equal(stats::var(as.vector(noisy[, , 1, ] * 1e4)), 100, tolerance = 0.1)
})

test_that("the calibration of recalibrated() scores negative binomial deaths for finite sizes", {
  n <- 4000
  mu <- c(100, 100, 50, 50)
  # observed cells: grid cell 1 (Poisson), grid cell 2 (size 4), grid cells 3 and 4 summed (size 4 each)
  P <- list(L = matrix(log(mu), 4, n), m = log(mu), h = rep(1, 4), E = rep(1, 4),
            C = Matrix::sparseMatrix(i = c(1, 2, 3, 3), j = 1:4, x = 1, dims = c(3, 4)), Ecell = rep(1, 3),
            pois = matrix(TRUE, 3, n), size = rbind(rep(Inf, n), rep(4, n), rep(4, n)))
  set.seed(1)
  P$U <- matrix(stats::runif(3 * n), 3)
  P$Z <- stats::qnorm(P$U)
  X <- exp(mortensor:::calibration_sample(0, P, "constant", FALSE))
  v <- apply(X, 1, stats::var)
  expect_equal(v, c(100, 100 + 100^2 / 4, 2 * (50 + 50^2 / 4)), tolerance = 0.1)
  expect_equal(rowMeans(X), c(100, 100, 100), tolerance = 0.02)
})

test_that("mt_deaths() does not add Poisson noise to predictive draws a second time", {
  x <- mt_example(deaths = "single", years = 2000:2015)
  fit <- mt_fit(lee_carter(), x)
  p <- mt_forecast(fit, h = 2, n_draws = 30, seed = 1, predictive = TRUE)
  d0 <- mt_deaths(p, noise = FALSE)
  expect_identical(mt_deaths(p, seed = 2), d0)
  expect_equal(d0, round(d0), tolerance = 1e-8)
  # draws of expected rates get Poisson noise
  a <- mt_forecast(fit, h = 2, n_draws = 30, seed = 1)
  e1 <- mt_deaths(a, seed = 2)
  expect_equal(e1, round(e1))
  expect_false(isTRUE(all.equal(e1, mt_deaths(a, noise = FALSE))))
})

test_that("stanmomo(family = \"nb\") draws deaths with the posterior dispersion of each draw", {
  skip_on_cran()
  skip_if_not_installed("StanMoMo")
  skip_if_not_installed("rstan")
  x <- mt_subset(mt_example(deaths = "single", years = 2006:2019), sex == "female")
  sa <- single_ages(60, 89, open = FALSE)
  x <- mt_regroup(mt_regroup(mt_set_grid(x, age = sa), age = sa),
                  age = age_groups(seq(60, 85, 5), open = FALSE, last_upper = 90))
  set.seed(1)
  fit <- suppressWarnings(mt_fit(stanmomo("lc", family = "nb", chains = 1, iter = 200, h_max = 2), x))
  f <- fit$fits[[1]]
  expect_true(length(f$nb_size) > 0 && all(f$nb_size > 0))
  p0 <- mt_forecast(fit, h = 1, n_draws = 50, seed = 1, predictive = TRUE)
  expect_true(all(p0$nb_size[1, ] %in% f$nb_size))
  # with strong overdispersion the deaths vary as the negative binomial
  fit$fits[[1]]$nb_size[] <- 2
  a <- mt_forecast(fit, h = 1, n_draws = 600, seed = 1)
  p <- mt_forecast(fit, h = 1, n_draws = 600, seed = 1, predictive = TRUE)
  expect_equal(p$nb_size[1, ], rep(2, 600))
  E <- as.vector(fit$exposure_last)
  mu <- a$rates[, 1, 1, ] * E
  D <- p$rates[, 1, 1, ] * E
  expect_equal(D, round(D), tolerance = 1e-8)
  z <- as.vector((D - mu) / sqrt(mu + mu^2 / 2))
  expect_lt(abs(stats::var(z) - 1), 0.25)
  expect_gt(stats::var(as.vector((D - mu) / sqrt(mu))), 5)
})

test_that("the deprecated functions and options removed in 0.3.0 are gone", {
  ours <- getNamespaceExports("mortensor")
  expect_false("product_ratio" %in% ours)
  expect_false("cause_shares" %in% ours)
  expect_false("observation_noise" %in% names(formals(rwd)))
  expect_false("transform" %in% names(formals(tucker3)))
  expect_error(rwd(observation_noise = TRUE))
  expect_error(tucker3(transform = "coda"))
  expect_no_warning(rwd())
  expect_no_warning(tucker3())
  expect_no_warning(li_lee())
})

test_that("arima_series() fits a given order: the local level model with drift", {
  skip_if_not_installed("forecast")
  x <- mt_example(deaths = "single", years = 2000:2015)
  fit <- mt_fit(arima_series(order = c(0, 1, 1), drift = TRUE), x)
  m <- fit$fits[[1]]$models[[1]]
  expect_equal(unname(m$arma[c(1, 6, 2)]), c(0, 1, 1))
  expect_true(all(c("ma1", "drift") %in% names(stats::coef(m))))
  expect_match(fit$spec$label, "ARIMA(0,1,1) with drift", fixed = TRUE)
  expect_false(fit$spec$smooth)
  expect_error(arima_series(order = c(0, 1)), class = "mt_error_argument")
  expect_error(arima_series(order = c(0, 1, 1), drift = NA), class = "mt_error_argument")
})

test_that("li_lee() warns when the bound of the AR(1) coefficient binds", {
  k <- 1.1^(1:20)
  expect_true(mortensor:::ar1_fit(k - mean(k))$bounded)
  set.seed(1)
  expect_false(mortensor:::ar1_fit(as.numeric(stats::arima.sim(list(ar = 0.5), 40)))$bounded)
  x <- mt_example(deaths = "single", years = 2000:2019)
  expect_no_warning(mt_fit(li_lee(), x))
  # the sexes drift apart at a growing pace: the deviations are explosive
  st <- mt_strata(x)
  x$obs$deaths <- data.table::copy(x$obs$deaths)
  male <- st$stratum[st$sex == "male"]
  x$obs$deaths[stratum == male, value := round(value * exp(0.01 * 1.25^(period_lo - 2000)))]
  expect_warning(mortensor:::fit_model(li_lee(), mt_array(x, "deaths"), mt_array(x, "exposure"), mt_grid(x)$age,
                                       2000:2019, strata = st),
                 class = "mt_warning_li_lee")
  # mt_fit() keeps the class in its warning about the fit
  expect_warning(fit <- mt_fit(li_lee(), x), class = "mt_warning_li_lee")
  expect_true(any(vapply(fit$fits$strata, function(s) s$ar[[1]]$bounded, logical(1))))
})

test_that("lee_carter() warns when it falls back to the ridge penalty or to SVD estimates", {
  x <- mt_subset(mt_example(deaths = "single", years = 2000:2015), sex == "female")
  D <- mt_array(x, "deaths")[, , 1]
  E <- mt_array(x, "exposure")[, , 1]
  ages <- mt_grid(x)$age
  real <- mortensor:::lc_poisson
  local_mocked_bindings(lc_poisson = function(D, E, init = NULL, maxit = 500, tol = 1e-8, ridge = 0) {
    if (ridge == 0) mt_abort("Lee-Carter fit diverged.", class = "model_fit")
    real(D, E, init, maxit, tol, ridge)
  })
  expect_warning(est <- mortensor:::fit_model(lee_carter(), D, E, ages, 2000:2015), "ridge",
                 class = "mt_warning_model_fit")
  expect_equal(est$method, "poisson")
  local_mocked_bindings(lc_poisson = function(...) mt_abort("Lee-Carter fit diverged.", class = "model_fit"))
  expect_warning(est <- mortensor:::fit_model(lee_carter(), D, E, ages, 2000:2015), "SVD",
                 class = "mt_warning_model_fit")
  expect_equal(est$method, "svd")
  expect_warning(mt_fit(lee_carter(), x), class = "mt_warning_model_fit")
})

test_that("tucker3() centres log rates as the three-way Lee-Carter model", {
  x <- mt_example(deaths = "single", years = 2000:2019)
  sp <- tucker3()
  expect_false(sp$params$residual_noise)
  expect_true(sp$params$centre)
  fit <- mt_fit(sp, x)
  y <- log(pmax(mt_array(x, "deaths"), 0.5) / mt_array(x, "exposure"))
  alpha <- apply(y, c(1, 3), mean)
  expect_equal(unname(fit$fits$alpha), unname(alpha))
  # the decomposition describes the deviations from a[x, i]: their mean over the years is zero
  dev <- sweep(log(fit$fits$fitted), c(1, 3), alpha)
  expect_lt(max(abs(apply(dev, c(1, 3), mean))), 1e-8)
  # ranks of the application of Russolillo et al. (2011); shares of
  # variation as in their Table 4 (squared core elements over the sum of
  # squares of the centred array), whose sum is the fit
  expect_equal(fit$fits$ranks, c(2, 2, 2))
  Xc <- sweep(y, c(1, 3), alpha)
  expect_equal(fit$fits$fit_pct, 100 * (1 - sum((log(fit$fits$fitted) - y)^2) / sum(Xc^2)))
  expect_equal(sum(fit$fits$core_pct), fit$fits$fit_pct)
  expect_equal(dim(fit$fits$core_pct), c(2, 2, 2))
  # one stratum, ranks 1: the Lee-Carter model (first SVD component of the centred log rates)
  x1 <- mt_subset(x, sex == "female")
  f1 <- mt_fit(tucker3(ranks = c(1, 1, 1)), x1)
  a1 <- rowMeans(y[, , 1])
  sv <- svd(y[, , 1] - a1, nu = 1, nv = 1)
  expect_equal(log(f1$fits$fitted[, , 1]), a1 + sv$d[1] * outer(sv$u[, 1], sv$v[, 1]), ignore_attr = TRUE,
               tolerance = 1e-8)
  # without centring, no means are removed
  f0 <- mt_fit(tucker3(centre = FALSE), x)
  expect_null(f0$fits$alpha)
  expect_error(tucker3(centre = NA), class = "mt_error_argument")
  fc <- mt_forecast(fit, h = 3, n_draws = 20, seed = 1)
  expect_true(all(is.finite(fc$rates) & fc$rates > 0))
})

test_that("the simulated Lee-Carter generator adds no noise beyond Poisson by default", {
  a <- mortensor:::simulate_lc_dgp(n_t = 5, seed = 2)
  b <- mortensor:::simulate_lc_dgp(n_t = 5, seed = 2, sd_obs = 0)
  expect_equal(mt_obs(a, "deaths"), mt_obs(b, "deaths"))
  c3 <- mortensor:::simulate_lc_dgp(n_t = 5, seed = 2, sd_obs = 0.03)
  expect_false(isTRUE(all.equal(mt_obs(a, "deaths")$value, mt_obs(c3, "deaths")$value)))
})
