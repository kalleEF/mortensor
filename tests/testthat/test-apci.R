# APCI model (Richards, Currie, Kleinow and Ritchie 2019)

# Poisson deaths from the APCI model with parameters that satisfy the
# constraints of Appendix C; gamma = 0 for cohorts with four or fewer cells
sim_apci <- function(seed, ages = 60:69, years = 2000:2014, exposure = 1e5) {
  set.seed(seed)
  nx <- length(ages)
  ny <- length(years)
  cohorts <- seq(min(years) - max(ages), max(years) - min(ages))
  ix <- rep(seq_len(nx), ny)
  iy <- rep(seq_len(ny), each = nx)
  ic <- years[iy] - ages[ix] - cohorts[1] + 1
  w <- tabulate(ic, length(cohorts))
  est <- w > 4
  alpha <- -4.5 + 0.09 * (ages - 60)
  beta <- -0.02 + 0.0005 * (ages - 60)
  kappa <- as.vector(stats::residuals(stats::lm(cumsum(stats::rnorm(ny, 0, 0.02)) ~ years)))
  # a random walk without weighted level, linear and quadratic trend
  cc <- seq_along(cohorts)[est]
  Q <- cbind(1, cc, cc^2)
  g <- cumsum(stats::rnorm(sum(est), 0, 0.02))
  gamma <- numeric(length(cohorts))
  gamma[est] <- g - Q %*% solve(crossprod(Q, w[est] * Q), crossprod(Q, w[est] * g))
  eta <- alpha[ix] + beta[ix] * (years[iy] - mean(years)) + kappa[iy] + gamma[ic]
  E <- matrix(exposure * stats::runif(nx * ny, 0.5, 1.5), nx)
  D <- matrix(stats::rpois(nx * ny, as.vector(E) * exp(eta)), nx)
  list(D = D, E = E, alpha = alpha, beta = beta, kappa = kappa, gamma = gamma, estimated = est, w = w,
       cohorts = cohorts, ages = ages, years = years)
}

# mt_data of a simulated stratum
apci_data <- function(x) {
  long <- data.frame(age = rep(x$ages, length(x$years)), year = rep(x$years, each = length(x$ages)),
                     deaths = as.vector(x$D), exposure = as.vector(x$E))
  mt_data(deaths = long, exposure = long)
}

# The five constraints of Appendix C (eq. 17, 21, 19, 20, 22), each scaled
# by the size of its terms
apci_constraints <- function(kappa, gamma, years, cohorts, w) {
  cc <- cohorts - min(cohorts) + 1
  c(sum(kappa) / sum(abs(kappa)), sum((years - min(years)) * kappa) / sum(abs((years - min(years)) * kappa)),
    sum(w * gamma) / sum(abs(w * gamma)), sum(w * cc * gamma) / sum(abs(w * cc * gamma)),
    sum(w * cc^2 * gamma) / sum(abs(w * cc^2 * gamma)))
}

test_that("the constraints of Appendix C hold exactly and corner cohorts are 0", {
  x <- sim_apci(1)
  d <- apci_data(x)
  for (smooth in c(FALSE, TRUE)) {
    fit <- mt_fit(apci(smooth = smooth), d)
    f <- fit$fits[[1]]
    expect_identical(f$cohorts, x$cohorts)
    expect_identical(f$cohort_cells, x$w)
    expect_identical(f$estimated, x$w > 4)
    expect_true(all(f$gamma[!f$estimated] == 0))
    expect_lt(max(abs(apci_constraints(f$kappa, f$gamma, x$years, f$cohorts, f$cohort_cells))), 1e-12)
    p <- mt_parameters(fit)
    expect_equal(p[parameter == "kappa"]$value, f$kappa)
    expect_equal(p[parameter == "gamma"]$index, as.character(x$cohorts))
    expect_equal(p[parameter == "ybar"]$value, mean(x$years))
  }
  # the smoothed fit has lower effective dimension and its smoothing parameters
  ps <- mt_parameters(mt_fit(apci(), d))
  pu <- mt_parameters(mt_fit(apci(smooth = FALSE), d))
  expect_lt(ps[parameter == "ed"]$value, pu[parameter == "ed"]$value)
  expect_length(ps[parameter %in% c("lambda_alpha", "lambda_beta")]$value, 2)
  expect_length(pu[parameter %in% c("lambda_alpha", "lambda_beta")]$value, 0)
})

test_that("the unsmoothed fit equals a constrained glm() fit", {
  x <- sim_apci(2)
  f <- mortensor:::apci_estimate(x$D, x$E, x$ages, x$years, smooth = FALSE)
  # the model matrix from factors, the corner cohorts left out, and the
  # constraints imposed by a basis of their null space
  dat <- expand.grid(age = x$ages, year = x$years)
  dat$coh <- dat$year - dat$age
  Xa <- stats::model.matrix(~ 0 + factor(age), dat)
  Xk <- stats::model.matrix(~ 0 + factor(year), dat)
  Xg <- stats::model.matrix(~ 0 + factor(coh), dat)
  wc <- colSums(Xg)
  Xg <- Xg[, wc >= 5]
  cc <- (sort(unique(dat$coh)) - min(dat$coh) + 1)[wc >= 5]
  wc <- wc[wc >= 5]
  X <- cbind(Xa, Xa * (dat$year - mean(x$years)), Xk, Xg)
  na <- ncol(Xa)
  nk <- ncol(Xk)
  ng <- ncol(Xg)
  H <- rbind(c(rep(0, 2 * na), rep(1, nk), rep(0, ng)),
             c(rep(0, 2 * na), x$years - min(x$years), rep(0, ng)),
             c(rep(0, 2 * na + nk), wc),
             c(rep(0, 2 * na + nk), wc * cc),
             c(rep(0, 2 * na + nk), wc * cc^2))
  Z <- qr.Q(qr(t(H)), complete = TRUE)[, -(1:5)]
  XZ <- X %*% Z
  g <- stats::glm(as.vector(x$D) ~ 0 + XZ, offset = log(as.vector(x$E)), family = stats::poisson(),
                  control = stats::glm.control(epsilon = 1e-14, maxit = 100))
  th <- as.vector(Z %*% stats::coef(g))
  expect_equal(f$alpha, th[seq_len(na)], tolerance = 1e-8)
  expect_equal(f$beta, th[na + seq_len(na)], tolerance = 1e-8)
  expect_equal(f$kappa, th[2 * na + seq_len(nk)], tolerance = 1e-8)
  expect_equal(f$gamma[f$estimated], th[2 * na + nk + seq_len(ng)], tolerance = 1e-8)
  expect_equal(f$dev, stats::deviance(g), tolerance = 1e-10)
  expect_equal(as.vector(f$fitted), as.vector(stats::fitted(g)) / as.vector(x$E), tolerance = 1e-8)
  # the effective dimension of the unsmoothed model is its number of free parameters
  expect_equal(f$ed, ncol(XZ), tolerance = 1e-8)
})

test_that("the APCI fit recovers the parameters of simulated data", {
  x <- sim_apci(3, ages = 60:79, years = 1991:2015, exposure = 2e5)
  f <- mortensor:::apci_estimate(x$D, x$E, x$ages, x$years, smooth = FALSE)
  # standard errors from the information of the constrained Poisson likelihood
  des <- mortensor:::apci_design(x$ages, x$years, FALSE)
  mu <- as.vector(f$fitted * x$E)
  info <- as.matrix(Matrix::crossprod(des$X, des$X * mu))
  p <- ncol(info)
  cov <- solve(rbind(cbind(info, t(des$H)), cbind(des$H, matrix(0, 5, 5))))[seq_len(p), seq_len(p)]
  est <- c(f$alpha, f$beta, f$kappa, f$gamma[f$estimated])
  truth <- c(x$alpha, x$beta, x$kappa, x$gamma[x$estimated])
  z <- (est - truth) / sqrt(diag(cov))
  expect_lt(max(abs(z)), 4.5)
  expect_gt(stats::sd(z), 0.7)
  expect_lt(stats::sd(z), 1.4)
})

test_that("the unsmoothed model of the paper's data dimensions has 241 free parameters (Table 2)", {
  # ages 50-104, years 1971-2015: 55 alpha, 55 beta, 45 kappa and 91 gamma
  # (99 cohorts less 4 at each end) less 5 constraints
  set.seed(4)
  ages <- 50:104
  years <- 1971:2015
  E <- matrix(1e4, length(ages), length(years))
  m <- exp(-9 + 0.09 * (ages - 50))
  D <- matrix(stats::rpois(length(E), E * m), length(ages))
  f <- mortensor:::apci_estimate(D, E, ages, years, smooth = FALSE)
  expect_equal(f$ed, 241, tolerance = 1e-8)
  expect_equal(sum(f$estimated), 91)
})

test_that("forecasts continue beta (y - ybar) linearly and project kappa and gamma", {
  x <- sim_apci(5, ages = 60:74, years = 1996:2015)
  fit <- mt_fit(apci(smooth = FALSE), apci_data(x))
  f <- fit$fits[[1]]
  h <- 6
  fc <- mt_forecast(fit, h = h, n_draws = 50, seed = 1)
  lr <- log(fc$rates[, , 1, ])
  yrs <- 2015 + seq_len(h)
  # ages whose cohorts were estimated: the difference between two of them is
  # the same in every draw and linear in the year
  old <- c(72, 74) - 59
  coh <- outer(-x$ages, yrs, `+`)
  expect_true(all(coh[old, ] <= max(f$cohorts[f$estimated])))
  dif <- lr[old[2], , ] - lr[old[1], , ]
  expect_lt(max(apply(dif, 1, stats::sd)), 1e-12)
  det <- f$alpha[old[2]] - f$alpha[old[1]] + (f$beta[old[2]] - f$beta[old[1]]) * (yrs - f$ybar) +
    f$gamma[match(coh[old[2], ], f$cohorts)] - f$gamma[match(coh[old[1], ], f$cohorts)]
  expect_equal(unname(dif[, 1]), det, tolerance = 1e-12)
  # without innovations, the draws are the ARIMA point forecasts of kappa and
  # gamma added to alpha + beta (y - ybar)
  f0 <- f
  f0$kappa_ts$sigma2 <- 0
  f0$gamma_ts$sigma2 <- 0
  sim <- mortensor:::simulate_model(fit$spec, f0, h, 2)
  kf <- stats::predict(f$kappa_ts, n.ahead = h)$pred
  c_last <- max(f$cohorts[f$estimated])
  gf <- stats::predict(f$gamma_ts, n.ahead = max(coh) - c_last)$pred
  gam <- ifelse(coh <= c_last, f$gamma[match(coh, f$cohorts)], gf[pmax(coh - c_last, 1)])
  manual <- f$alpha + outer(f$beta, yrs - f$ybar) + matrix(kf, length(x$ages), h, byrow = TRUE) + gam
  expect_equal(log(sim[, , 1]), manual, tolerance = 1e-10)
  expect_equal(sim[, , 1], sim[, , 2])
  # the youngest cohorts with four or fewer cells (gamma = 0 in the fit) take
  # the projected gamma
  young <- coh > c_last & coh <= max(f$cohorts)
  expect_true(any(young))
  expect_true(all(abs(gam[young]) > 0))
})

test_that("ARIMA paths agree with predict() and orders follow Appendix E", {
  set.seed(6)
  y <- cumsum(cumsum(stats::rnorm(30, 0, 0.1)))
  for (ord in list(c(1, 2, 1), c(2, 1, 2), c(0, 1, 0))) {
    fa <- stats::arima(y, order = ord, include.mean = FALSE)
    pr <- stats::predict(fa, n.ahead = 5)
    paths <- with_seed_test(7, mortensor:::arima_paths(fa, 5, 20000))
    expect_lt(max(abs(rowMeans(paths) - pr$pred) / pr$se), 0.05)
    expect_lt(max(abs(apply(paths, 1, stats::sd) / pr$se - 1)), 0.03)
  }
  # order selection: the smallest AICc over p, q in 0:3 for each d, then the
  # smallest root mean squared error from the third residual
  sel <- mortensor:::apci_arima(y, d = 1:2)
  aicc <- function(o) {
    f <- tryCatch(suppressWarnings(stats::arima(y, order = o, include.mean = FALSE)), error = function(e) NULL)
    if (is.null(f)) return(Inf)
    k <- o[1] + o[3] + 1
    if (f$nobs - k - 1 > 0) f$aic + 2 * k * (k + 1) / (f$nobs - k - 1) else Inf
  }
  best <- lapply(1:2, function(d) {
    grid <- expand.grid(p = 0:3, q = 0:3)
    a <- mapply(function(p, q) aicc(c(p, d, q)), grid$p, grid$q)
    c(grid$p[which.min(a)], d, grid$q[which.min(a)])
  })
  rmse <- vapply(best, function(o) {
    sqrt(mean(stats::residuals(stats::arima(y, order = o, include.mean = FALSE))[-(1:2)]^2))
  }, 0)
  expect_equal(sel$order, best[[which.min(rmse)]])
  expect_false(any(c("intercept", "drift") %in% names(stats::coef(sel))))
  expect_equal(mortensor:::apci_arima(y, order = c(1, 1, 0))$order, c(1, 1, 0))
})

test_that("apci() checks its arguments and data and is described", {
  expect_error(apci(kappa_order = c(1, 1)), class = "mt_error_argument")
  expect_error(apci(gamma_order = c(-1, 1, 0)), class = "mt_error_argument")
  expect_identical(apci(kappa_order = c(1, 1, 2))$params$kappa_order, c(1L, 1L, 2L))
  x <- mt_example(deaths = "single", years = 2000:2015)
  # the open age interval of the example data
  expect_error(mt_fit(apci(), mt_subset(x, sex == "male", ages = c(80, Inf))), class = "mt_error_model_fit")
  d <- apci_data(sim_apci(8))
  fit <- mt_fit(apci(kappa_order = c(0, 1, 1), gamma_order = c(1, 1, 0)), d)
  p <- mt_parameters(fit)
  expect_equal(p[parameter == "order_kappa"]$value, c(0, 1, 1))
  expect_equal(p[parameter == "order_gamma"]$value, c(1, 1, 0))
  expect_match(summary(fit)$assumptions, "APCI")
  fc <- mt_forecast(fit, h = 4, n_draws = 30, seed = 1, predictive = TRUE)
  expect_true(all(is.finite(fc$rates)))
  expect_null(fc$nb_size)
})
