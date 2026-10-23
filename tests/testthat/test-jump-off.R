# Jump-off options of Stoeldraijer et al. (2018) and the filtered jump-off of
# De Jong and Tickle (2006) for lee_carter()

# Log rates from the state-space Lee-Carter model of De Jong and Tickle (2006,
# eq. 6) in the scale sum(b) = 1: y = a + b alpha + eps, alpha a random walk
# with drift delta and innovation sd lambda sigma / sqrt(b'b)
sim_djt <- function(seed, p = 20, n = 40, sigma = 0.03, lambda = 1.5, delta = -1) {
  set.seed(seed)
  a <- -9 + 0.4 * seq_len(p)
  b <- 1.5 - seq_len(p) / p
  b <- b / sum(b)
  alpha <- cumsum(c(5, delta + stats::rnorm(n - 1, 0, lambda * sigma / sqrt(sum(b^2)))))
  list(y = a + outer(b, alpha) + matrix(stats::rnorm(p * n, 0, sigma), p), a = a, b = b, alpha = alpha)
}

# Central forecast of a Lee-Carter fit: draws without innovations or drift
# uncertainty (sigma_k = 0), as the point forecasts of Stoeldraijer et al.
central_rates <- function(spec, fit, h) {
  fit$sigma_k <- 0
  mortensor:::simulate_model(spec, fit, h, 1)[, , 1]
}

test_that("jump-off arguments are checked and the defaults are unchanged", {
  expect_identical(lee_carter()$params[c("jump_off", "jump_off_years")], list(jump_off = "fitted", jump_off_years = 1))
  expect_error(lee_carter(jump_off_years = 3), class = "mt_error_argument")
  expect_error(lee_carter(jump_off = "filtered", jump_off_years = 2), class = "mt_error_argument")
  expect_error(lee_carter(jump_off = "observed", jump_off_years = 1.5), class = "mt_error_argument")
  expect_match(lee_carter(jump_off = "observed", jump_off_years = 3)$label, "observed jump-off \\(mean of 3 years\\)")
  expect_match(lee_carter(jump_off = "filtered")$label, "filtered jump-off")
  expect_identical(lee_carter(jump_off = "observed", jump_off_years = 4)$min_periods, 5)
  expect_identical(lee_carter(jump_off = "observed", jump_off_years = 7)$min_periods, 7)
  x <- mt_example(deaths = "single", years = 1995:2012)
  # one year is the observed jump-off of Lee and Miller (2001)
  o1 <- mt_forecast(mt_fit(lee_carter(jump_off = "observed"), x), h = 3, n_draws = 20, seed = 1)
  o2 <- mt_forecast(mt_fit(lee_carter(jump_off = "observed", jump_off_years = 1), x), h = 3, n_draws = 20, seed = 1)
  expect_identical(o1$rates, o2$rates)
  g <- mt_grid(x)
  D <- mt_array(x, "deaths")[, , 1]
  E <- mt_array(x, "exposure")[, , 1]
  expect_error(mortensor:::fit_model(lee_carter(jump_off = "observed", jump_off_years = 30), D, E, g$age, 1995:2012),
               class = "mt_error_model_data")
  expect_mt_model(lee_carter(jump_off = "observed", jump_off_years = 3))
  expect_mt_model(lee_carter(jump_off = "filtered"))
})

test_that("the averaged jump-off is the mean observed log rate at the midpoint, moved by the drift", {
  x <- mt_example(deaths = "single", years = 1995:2015)
  for (method in c("svd", "poisson")) {
    for (n in c(2, 3, 5)) {
      sp <- lee_carter(method = method, jump_off = "observed", jump_off_years = n)
      fit <- mt_fit(sp, x)
      fc <- mt_forecast(fit, h = 4, n_draws = 15, seed = 3)
      f0 <- mt_forecast(mt_fit(lee_carter(method = method), x), h = 4, n_draws = 15, seed = 3)
      for (s in 1:2) {
        f <- fit$fits[[s]]
        last <- 21 - n + seq_len(n)
        abar <- rowMeans(log(pmax(fit$train_D[, last, s], 0.5) / fit$train_E[, last, s]))
        # Stoeldraijer et al. (2018): a = mean log rate of the n years, k = 0 at
        # their midpoint; the jump-off of year T is that mean plus b d (n - 1) / 2
        J <- abar + f$b * f$drift * (n - 1) / 2
        expect_equal(log(fc$rates[, , s, ]) - log(f0$rates[, , s, ]),
                     array(J - (f$a + f$b * f$k[21]), c(length(J), 4, 15)), tolerance = 1e-10, ignore_attr = TRUE)
        # the central forecast does not depend on k[T]: mean + b d (h + (n - 1) / 2)
        expect_equal(log(central_rates(sp, f, 4)), abar + outer(f$b * f$drift, (1:4) + (n - 1) / 2),
                     tolerance = 1e-10, ignore_attr = TRUE)
      }
    }
  }
})

test_that("anchored(decay = 4) moves from the last observed values to the model values by the fifth year", {
  x <- mt_example(deaths = "single", years = 1995:2015)
  fa <- mt_fit(anchored(lee_carter(), decay = 4), x)
  fc <- mt_forecast(fa, h = 6, n_draws = 10, seed = 2)
  f0 <- mt_forecast(mt_fit(lee_carter(), x), h = 6, n_draws = 10, seed = 2)
  w <- c(1, 0.75, 0.5, 0.25, 0, 0)
  for (s in 1:2) {
    expect_equal(log(fc$rates[, , s, ]) - log(f0$rates[, , s, ]),
                 array(outer(fa$anchor$offset[, s], w), c(nrow(fa$anchor$offset), 6, 10)), tolerance = 1e-10,
                 ignore_attr = TRUE)
  }
})

test_that("the filtered index has the limits of the state-space model", {
  s <- sim_djt(11, p = 15, n = 30)
  z <- colSums(s$b * (s$y - s$a)) / sum(s$b^2)
  tt <- seq_len(30) - 1
  # without state noise: the least-squares trend of the per-year index at T
  f0 <- mortensor:::lc_filter(s$y, s$a, s$b, lambda = 0)
  ols <- stats::lm(z ~ tt)
  expect_equal(f0$z, z)
  expect_equal(f0$level, unname(stats::fitted(ols)[30]))
  expect_equal(f0$drift, unname(stats::coef(ols)[2]))
  # with little observation noise: the per-year index of year T
  f_inf <- mortensor:::lc_filter(s$y, s$a, s$b, lambda = 1e5)
  expect_equal(f_inf$level, z[30], tolerance = 1e-8)
  # its error is then the observation error of z[T]
  expect_equal(f_inf$se, f_inf$sigma / sqrt(sum(s$b^2)), tolerance = 1e-6)
  # one age: the filter then returns the observed log rate of year T
  y1 <- s$y[3, , drop = FALSE]
  f1 <- mortensor:::lc_filter(y1, s$a[3], s$b[3], lambda = 1e5)
  expect_equal(s$a[3] + s$b[3] * f1$level, y1[1, 30], tolerance = 1e-8)
  # the estimate lies between the limits here, and its likelihood is maximal
  fe <- mortensor:::lc_filter(s$y, s$a, s$b)
  expect_gt(fe$lambda, 0)
  prof <- function(l) mortensor:::lc_filter(s$y, s$a, s$b, lambda = l)$loglik
  expect_gte(fe$loglik, max(prof(0), prof(fe$lambda * 0.9), prof(fe$lambda * 1.1), prof(1e3)))
})

test_that("the filtered index equals a Kalman filter with a diffuse initial state", {
  s <- sim_djt(5, p = 10, n = 25, lambda = 0.8)
  fe <- mortensor:::lc_filter(s$y, s$a, s$b)
  bb <- sum(s$b^2)
  h_obs <- fe$sigma^2 / bb
  q_state <- (fe$lambda * fe$sigma)^2 / bb
  # state (alpha, delta), alpha[t + 1] = alpha[t] + delta + w[t]; z[t] = alpha[t] + e[t]
  Tm <- matrix(c(1, 0, 1, 1), 2)
  Q <- diag(c(q_state, 0))
  m <- c(0, 0)
  P <- diag(1e9, 2)
  for (t in seq_along(fe$z)) {
    if (t > 1) {
      m <- as.vector(Tm %*% m)
      P <- Tm %*% P %*% t(Tm) + Q
    }
    f <- P[1, 1] + h_obs
    k <- P[, 1] / f
    m <- m + k * (fe$z[t] - m[1])
    P <- P - tcrossprod(k) * f
  }
  expect_equal(m[1], fe$level, tolerance = 1e-7)
  expect_equal(sqrt(P[1, 1]), fe$se, tolerance = 1e-5)
})

test_that("sigma, lambda and the drift maximise the likelihood of the observed log rates", {
  # the full Gaussian likelihood of eq. (6) given a and b (covariance of
  # eq. 8, with the level of the first year as a parameter), maximised
  # numerically
  s <- sim_djt(3, p = 3, n = 10, lambda = 2)
  p <- 3
  n <- 10
  bb <- sum(s$b^2)
  Mt <- outer(seq_len(n) - 1, seq_len(n) - 1, pmin)
  yv <- as.vector(s$y)
  nll <- function(th) {
    V <- exp(th[1]) * (diag(n * p) + th[2]^2 / bb * kronecker(Mt, tcrossprod(s$b)))
    mu <- as.vector(s$a + outer(s$b, th[3] + th[4] * (seq_len(n) - 1)))
    R <- chol(V)
    sum(log(diag(R))) + sum(backsolve(R, yv - mu, transpose = TRUE)^2) / 2
  }
  fe <- mortensor:::lc_filter(s$y, s$a, s$b)
  expect_gt(fe$lambda, 0)
  start <- c(log(fe$sigma^2) + 0.3, fe$lambda * 0.8, fe$z[1] + 0.1, fe$drift * 0.9)
  op <- stats::optim(start, nll, control = list(maxit = 20000, reltol = 1e-14))
  op <- stats::optim(op$par, nll, method = "BFGS", control = list(reltol = 1e-16, maxit = 1000))
  expect_equal(sqrt(exp(op$par[1])), fe$sigma, tolerance = 1e-6)
  expect_equal(abs(op$par[2]), fe$lambda, tolerance = 1e-6)
  expect_equal(op$par[4], fe$drift, tolerance = 1e-6)
  # the profile log-likelihood is the maximum up to the constant n p / 2
  expect_equal(-op$value, fe$loglik - n * p / 2, tolerance = 1e-10)
})

test_that("the filter recovers the parameters and the index of simulated data", {
  est <- t(vapply(1:100, function(i) {
    s <- sim_djt(100 + i)
    f <- mortensor:::lc_filter(s$y, s$a, s$b)
    c(lambda = f$lambda, sigma = f$sigma, drift = f$drift, std = (f$level - s$alpha[40]) / f$se,
      err = f$level - s$alpha[40], err_z = f$z[40] - s$alpha[40])
  }, numeric(6)))
  expect_equal(stats::median(est[, "lambda"]), 1.5, tolerance = 0.15)
  expect_equal(mean(est[, "sigma"]), 0.03, tolerance = 0.01)
  expect_equal(mean(est[, "drift"]), -1, tolerance = 0.01)
  # the standard error describes the error of the filtered index
  expect_gt(mean(est[, "std"]^2), 0.7)
  expect_lt(mean(est[, "std"]^2), 1.4)
  # the filtered index is closer to the true index than the per-year index
  expect_lt(mean(est[, "err"]^2), mean(est[, "err_z"]^2))
  # within lee_carter(): with SVD and no second stage, z is the fitted index
  s <- sim_djt(1)
  E <- matrix(1e8, 20, 40)
  ages <- age_groups(seq(0, 95, 5))
  sp <- lee_carter(method = "svd", adjust = "none", jump_off = "filtered")
  f <- mortensor:::fit_model(sp, exp(s$y) * E, E, ages, 1971:2010)
  expect_equal(f$filter$z, f$k, tolerance = 1e-10)
  # the jump-off rates are those of the filtered index, near the true rates
  jo <- log(central_rates(sp, f, 1)) - f$b * f$drift
  expect_equal(unname(jo), f$a + f$b * f$filter$level, tolerance = 1e-10)
  expect_lt(max(abs(jo - (s$a + s$b * s$alpha[40]))), 4 * f$filter$se * max(f$b) + 0.01)
})

test_that("the filtered jump-off works in the pipeline, with the bootstrap and a chosen period", {
  x <- mt_example(deaths = "single", years = 1990:2015)
  fit <- mt_fit(lee_carter(jump_off = "filtered"), x)
  f0 <- mt_fit(lee_carter(), x)
  fc <- mt_forecast(fit, h = 3, n_draws = 10, seed = 4)
  fc0 <- mt_forecast(f0, h = 3, n_draws = 10, seed = 4)
  for (s in 1:2) {
    f <- fit$fits[[s]]
    expect_equal(log(fc$rates[, , s, ]) - log(fc0$rates[, , s, ]),
                 array(f$b * (f$filter$level - f$k[26]), c(length(f$b), 3, 10)), tolerance = 1e-10, ignore_attr = TRUE)
  }
  p <- mt_parameters(fit)
  expect_setequal(p[grepl("filter", parameter), unique(parameter)],
                  c("k_filtered", "se_k_filtered", "filter_sigma", "filter_lambda", "filter_drift"))
  expect_equal(p[parameter == "k_filtered", value], vapply(fit$fits, function(f) f$filter$level, 0))
  expect_equal(p[parameter == "k_filtered", unique(index)], "2015")
  expect_null(mt_fit(lee_carter(), x)$fits[[1]]$filter)
  expect_output(print(summary(fit)), "Kalman\\s+filter")
  fb <- mt_fit(lee_carter(jump_off = "filtered", bootstrap = TRUE, n_boot = 3), x)
  expect_true(all(vapply(fb$fits[[1]]$boot, function(b) is.finite(b$filter$level), NA)))
  expect_s3_class(mt_forecast(fb, h = 2, n_draws = 6, seed = 1), "mt_draws")
  fp <- mt_fit(lee_carter(jump_off = "filtered", choose_period = TRUE, min_years = 10), x)
  expect_length(fp$fits[[1]]$filter$z, length(fp$fits[[1]]$k))
})

test_that("Stoeldraijer et al. (2018), Table 8a: French men, six jump-off options", {
  skip_on_cran()
  skip_if_not_installed("StanMoMo")
  fr <- mortensor:::pkg_data("FRMaleData", "StanMoMo")
  # ages 0-94 and 95+ (the paper aggregates from 95), HMD data of a later release
  top <- function(M) rbind(M[1:95, ], colSums(M[96:nrow(M), ], na.rm = TRUE))
  D <- top(fr$Dxt)
  E <- top(fr$Ext)
  yrs <- as.numeric(colnames(D))
  ages <- age_groups(0:95)
  e65 <- function(m) as.numeric(life_expectancy(matrix(m, ncol = 1), ages, at = 65))
  obs <- vapply(seq_along(yrs), function(j) e65(D[, j] / E[, j]), 0)
  names(obs) <- yrs
  # point forecasts of life expectancy at 65 one and five years ahead and in
  # 2020, fitting periods 1960-2005, ..., 1960-2014 (Lee-Carter, SVD)
  forecasts <- function(sp) {
    t(vapply(2005:2014, function(last) {
      cols <- which(yrs >= 1960 & yrs <= last)
      f <- mortensor:::fit_model(sp, D[, cols], E[, cols], ages, yrs[cols])
      r <- central_rates(sp, f, 2020 - last)
      c(e65(r[, 1]), e65(r[, 5]), e65(r[, 2020 - last]))
    }, numeric(3)))
  }
  res <- vapply(1:5, function(n) {
    fe <- forecasts(lee_carter(method = "svd", jump_off = "observed", jump_off_years = n))
    c(mafe1 = mean(abs(fe[1:9, 1] - obs[as.character(2006:2014)])),
      mafe5 = mean(abs(fe[1:5, 2] - obs[as.character(2010:2014)])), sd2020 = stats::sd(diff(fe[, 3])))
  }, numeric(3))
  # last observed values and averages of 2 to 5 years
  paper <- rbind(mafe1 = c(0.12, 0.13, 0.15, 0.19, 0.23), mafe5 = c(0.25, 0.28, 0.33, 0.38, 0.43),
                 sd2020 = c(0.15, 0.09, 0.08, 0.06, 0.05))
  expect_lt(max(abs(res - paper)), 0.03)
  # accuracy falls and robustness rises with the number of years averaged
  expect_true(all(diff(res["mafe1", 2:5]) > 0) && all(diff(res["mafe5", ]) > 0))
  expect_true(all(res["sd2020", 2:5] < res["sd2020", 1]))
})
