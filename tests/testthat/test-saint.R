# SAINT model (Jarner and Kryger 2011; Jarner and Jallbjorn 2022)

# Data simulated from the model: per sex a reference ("ref") and targets.
# Rates of the reference follow the trend (4.10) with M cumulated from the
# rates themselves (saint_trend_rates()); a target has the trend times
# exp(y_t' r(x)) with spreads `spread[[target]]` [year, 5] (zero if absent).
# noise = FALSE gives deaths equal to their expected values.
saint_world <- function(x = 40:100, years = 1996:2015, targets = c("A", "B"), spread = list(), noise = TRUE,
                        ref_size = 4e6, tgt_size = 4e5, sigma2 = c(female = 0.15, male = 0.1), seed = 1) {
  set.seed(seed)
  nt <- length(years)
  tt <- years - 2005
  R <- mortensor:::saint_regressors(x)
  rows <- list()
  truth <- list()
  for (sx in c("female", "male")) {
    theta <- rbind((if (sx == "female") -3.7 else -3.3) - 0.015 * tt, 0.105 + 0.0004 * tt, -0.01 + 0.0002 * tt)
    zeta <- log(if (sx == "female") 2e-4 else 5e-4) - 0.02 * tt
    mu <- mortensor:::saint_trend_rates(theta, zeta, sigma2[[sx]], x)
    truth[[sx]] <- list(theta = theta, zeta = zeta, sigma2 = sigma2[[sx]], mu = mu)
    for (r in c("ref", targets)) {
      y <- spread[[r]] %||% matrix(0, nt, 5)
      m <- mu * exp(R %*% t(y))
      E <- matrix((if (r == "ref") ref_size else tgt_size) * exp(-0.03 * (x - x[1])), length(x), nt)
      D <- if (noise) matrix(stats::rpois(length(m), m * E), length(x)) else m * E
      rows[[length(rows) + 1L]] <- data.table::data.table(region = r, sex = sx, year = rep(years, each = length(x)),
                                                          age = rep(x, nt), deaths = as.vector(D),
                                                          exposure = as.vector(E))
    }
  }
  d <- data.table::rbindlist(rows)
  x_data <- suppressWarnings(mt_data(deaths = d, exposure = d, by = c("region", "sex")))
  list(data = x_data, truth = truth)
}

test_that("spread regressors follow eq. 5.11 and are identified on the age range", {
  R <- unname(mortensor:::saint_regressors(20:100))
  expect_equal(dim(R), c(81, 5))
  at <- function(age, i) R[age - 19, i]
  expect_equal(c(at(20, 1), at(30, 1), at(40, 1), at(45, 1)), c(1, 0.5, 0, 0))
  expect_equal(c(at(60, 3), at(70, 3), at(80, 3)), c(1, 0.5, 0))
  expect_equal(c(at(80, 4), at(90, 4), at(100, 4)), c(1, 0.5, 0))
  expect_true(all(R[, 5] == 1))
  expect_identical(mortensor:::saint_spread_columns(20:100), 1:5)
  expect_identical(mortensor:::saint_spread_columns(40:100), 2:5)
  expect_identical(mortensor:::saint_spread_columns(60:100), 3:5)
  # ages 20-79: r4 is 1 at every age like the level r5 and is left out
  expect_identical(mortensor:::saint_spread_columns(20:79), c(1L, 2L, 3L, 5L))
  for (x in list(20:100, 40:100, 20:79, 75:110)) {
    cols <- mortensor:::saint_spread_columns(x)
    expect_equal(qr(mortensor:::saint_regressors(x)[, cols])$rank, length(cols))
  }
  # baseline (5.2): the excess slope needs ages on both sides of 75
  expect_identical(colnames(mortensor:::saint_design(20:100)), c("alpha", "beta", "kappa"))
  expect_identical(colnames(mortensor:::saint_design(75:100)), c("alpha", "beta"))
})

test_that("cumulated rates follow the cohort diagonals (eqs. 4.7-4.8)", {
  mt <- matrix(seq_len(12) / 100, 4, 3)
  M <- mortensor:::saint_cumulated(mt)
  expect_equal(M[1, ], c(0, 0, 0))
  # first year: earlier years take the rates of the first year
  expect_equal(M[, 1], c(0, cumsum(mt[1:3, 1])))
  # cell (age 3, year 3): rates of (age 1, year 1) and (age 2, year 2)
  expect_equal(M[3, 3], mt[1, 1] + mt[2, 2])
  # cell (age 4, year 2): (age 1, year 1 for year 0), (age 2, year 1), (age 3, year 2)... along the diagonal
  expect_equal(M[4, 2], mt[1, 1] + mt[2, 1] + mt[3, 1])
  expect_equal(M[4, 3], mt[1, 1] + mt[2, 1] + mt[3, 2])
})

test_that("the trend recovers its parameters from deaths equal to their expectations", {
  x <- 30:100
  tt <- seq(-7, 7)
  theta <- rbind(-3.6 - 0.015 * tt, 0.105 + 0.0005 * tt, -0.01 + 0.0002 * tt)
  zeta <- log(3e-4) - 0.02 * tt
  mu <- mortensor:::saint_trend_rates(theta, zeta, 0.2, x)
  E <- matrix(1e6 * exp(-0.03 * (x - 30)), length(x), length(tt))
  f <- mortensor:::saint_trend_fit(mu * E, E, x)
  expect_true(f$converged)
  expect_equal(f$sigma2, 0.2, tolerance = 1e-5)
  expect_equal(unname(f$theta), theta, tolerance = 1e-6)
  expect_equal(unname(f$zeta), zeta, tolerance = 1e-5)
  expect_equal(f$fitted, mu, tolerance = 1e-8)
  # mean frailty exp(-sigma2 M) equals 1 / (1 + sigma2 M0) with M0 = (exp(sigma2 M) - 1) / sigma2 (Example 4.2)
  M0 <- expm1(f$sigma2 * f$M) / f$sigma2
  expect_equal(1 / (1 + f$sigma2 * M0), exp(-f$sigma2 * f$M), tolerance = 1e-12)
})

test_that("background = FALSE fits the trend without background mortality", {
  x <- 30:100
  tt <- seq(-5, 5)
  theta <- rbind(-3.6 - 0.015 * tt, 0.105 + 0.0005 * tt, -0.01 + 0.0002 * tt)
  mu <- mortensor:::saint_trend_rates(theta, rep(-Inf, length(tt)), 0.2, x)
  E <- matrix(1e6 * exp(-0.03 * (x - 30)), length(x), length(tt))
  f <- mortensor:::saint_trend_fit(mu * E, E, x, background = FALSE)
  expect_equal(f$sigma2, 0.2, tolerance = 1e-5)
  expect_equal(unname(f$theta), theta, tolerance = 1e-6)
  expect_true(all(f$zeta == -Inf))
  w <- saint_world(noise = TRUE, targets = "A", years = 2001:2012)
  fit <- mt_fit(saint(reference = list(region = "ref"), by = "sex", background = FALSE), w$data)
  p <- mt_parameters(fit)
  expect_false(any(c("zeta", "sd_zeta") %in% p$parameter))
  fc <- mt_forecast(fit, h = 5, n_draws = 10, seed = 1)
  expect_true(all(is.finite(fc$rates) & fc$rates > 0))
  expect_error(saint(reference = "ref", background = NA), class = "mt_error_argument")
})

test_that("Newton-Raphson sweeps reach the estimates of the EM algorithm (Appendix B)", {
  x <- 30:100
  tt <- seq(-5, 5)
  theta <- rbind(-3.6 - 0.015 * tt, 0.105 + 0.0005 * tt, -0.01 + 0.0002 * tt)
  zeta <- log(3e-4) - 0.02 * tt
  mu <- mortensor:::saint_trend_rates(theta, zeta, 0.15, x)
  E <- matrix(5e5 * exp(-0.03 * (x - 30)), length(x), length(tt))
  set.seed(3)
  D <- matrix(stats::rpois(length(mu), mu * E), nrow(mu))
  for (s2 in c(0, 0.15, 0.4)) {
    em <- mortensor:::saint_trend_em(D, E, x, s2, tol = 1e-11)
    nr <- mortensor:::saint_trend_newton(D, E, x, s2)
    expect_true(em$converged && nr$converged)
    expect_equal(unname(nr$theta), unname(em$theta), tolerance = 1e-6)
    expect_equal(unname(nr$zeta), unname(em$zeta), tolerance = 1e-5)
    expect_equal(nr$loglik, em$loglik, tolerance = 1e-10)
  }
})

test_that("maximum likelihood of the trend recovers the parameters from Poisson deaths", {
  w <- saint_world(noise = TRUE, targets = character())
  D <- mt_array(w$data, "deaths")
  E <- mt_array(w$data, "exposure")
  for (sx in c("female", "male")) {
    f <- mortensor:::saint_trend_fit(D[, , paste0("ref.", sx)], E[, , paste0("ref.", sx)], 40:100)
    tr <- w$truth[[sx]]
    expect_lt(abs(f$sigma2 - tr$sigma2), 0.02)
    expect_lt(max(abs(f$theta[1, ] - tr$theta[1, ])), 0.02)
    expect_lt(max(abs(f$theta[2, ] - tr$theta[2, ])), 0.002)
    expect_lt(max(abs(log(f$fitted / tr$mu))), 0.03)
    # the estimate maximises the profile likelihood (4.14)
    for (d in c(-0.01, 0.01)) {
      expect_gt(f$loglik, mortensor:::saint_trend_newton(D[, , paste0("ref.", sx)], E[, , paste0("ref.", sx)],
                                                          40:100, f$sigma2 + d)$loglik)
    }
  }
})

test_that("a target identical to the reference has spread zero; known spreads are recovered", {
  nt <- 20
  y_b <- cbind(0, 0.05 * sin(seq_len(nt)), 0.1, -0.05, 0.02 * seq_len(nt) / nt)
  w <- saint_world(noise = FALSE, targets = c("same", "B"), spread = list(B = y_b), tgt_size = 4e6)
  # spreads without noise: their VAR(1) is not identified (warnings)
  fit <- suppressWarnings(mt_fit(saint(reference = list(region = "ref"), by = "sex"), w$data))
  f <- fit$fits
  lab <- fit$strata$label
  for (sx in c("female", "male")) {
    expect_null(f$spreads[[match(paste0("ref.", sx), lab)]])
    y_same <- f$spreads[[match(paste0("same.", sx), lab)]]$y
    expect_lt(max(abs(y_same)), 1e-7)
    yb <- f$spreads[[match(paste0("B.", sx), lab)]]$y
    expect_identical(colnames(yb), c("y2", "y3", "y4", "y5"))
    expect_equal(unname(yb), y_b[, 2:5], tolerance = 1e-6)
    # the reference stratum is forecast by the trend
    expect_equal(stats::fitted(fit)[, , paste0("ref.", sx)], f$trends[[sx]]$fitted, ignore_attr = TRUE)
  }
})

test_that("with Poisson deaths the spread of a target from the reference's intensity is within its error", {
  w <- saint_world(noise = TRUE, targets = "same", seed = 7)
  fit <- mt_fit(saint(reference = list(region = "ref"), by = "sex"), w$data)
  f <- fit$fits
  E <- mt_array(w$data, "exposure")
  for (sx in c("female", "male")) {
    s <- match(paste0("same.", sx), fit$strata$label)
    y <- f$spreads[[s]]$y
    mu <- f$trends[[sx]]$fitted
    z <- vapply(seq_len(nrow(y)), function(t) {
      info <- crossprod(f$R, f$R * (mu[, t] * E[, t, s]))
      y[t, ] / sqrt(diag(solve(info)))
    }, numeric(ncol(y)))
    expect_lt(max(abs(z)), 4.5)
  }
})

test_that("VAR(1) estimators of the spread follow their restrictions", {
  set.seed(11)
  n <- 4000
  A <- diag(c(0.5, 0.9, 0.2))
  Om <- matrix(c(1, 0.3, 0.1, 0.3, 1, -0.2, 0.1, -0.2, 1), 3) * 1e-3
  y <- matrix(0, n, 3)
  L <- t(chol(Om))
  for (t in 2:n) y[t, ] <- A %*% y[t - 1, ] + L %*% stats::rnorm(3)
  colnames(y) <- c("y3", "y4", "y5")
  d <- mortensor:::saint_var_fit(y, "diagonal")
  expect_equal(unname(diag(d$A)), diag(A), tolerance = 0.05)
  expect_true(all(d$A[row(d$A) != col(d$A)] == 0))
  expect_equal(d$Omega, Om, tolerance = 0.1, ignore_attr = TRUE)
  # unrestricted: conditional maximum likelihood, the least-squares regressions on the lagged values
  full <- mortensor:::saint_var_fit(y, "full")
  for (i in 1:3) {
    expect_equal(unname(full$A[i, ]), unname(stats::coef(stats::lm(y[-1, i] ~ y[-n, ] - 1))), tolerance = 1e-10)
  }
  e <- y[-1, ] - y[-n, ] %*% t(full$A)
  expect_equal(unname(full$Omega), unname(crossprod(e) / (n - 1)), tolerance = 1e-12)
  # fixed: A = a I
  expect_equal(unname(mortensor:::saint_var_fit(y, 0.7)$A), diag(0.7, 3))
  # diagonal coefficients bounded to [0, 0.99]
  rw <- apply(matrix(stats::rnorm(600), 200), 2, cumsum)
  b <- mortensor:::saint_var_fit(cbind(rw[, 1], (-1)^(1:200)), "diagonal")
  expect_equal(unname(diag(b$A)), c(0.99, 0))
  expect_identical(b$bound, c(TRUE, FALSE))
})

test_that("spread forecasts have mean A^h y_T and variance sum A^i Omega A^i' (2011, Section 2.3.3)", {
  A <- matrix(c(0.7, 0.1, -0.2, 0.5), 2)
  Om <- matrix(c(2, 0.5, 0.5, 1), 2) * 1e-3
  sp <- list(y = rbind(c(0, 0), c(0.2, -0.1)), A = A, Omega = Om)
  set.seed(5)
  p <- mortensor:::saint_var_paths(sp, 3, 40000)
  m3 <- A %*% A %*% A %*% c(0.2, -0.1)
  V3 <- Om + A %*% Om %*% t(A) + A %*% A %*% Om %*% t(A %*% A)
  expect_equal(rowMeans(p[, 3, ]), as.vector(m3), tolerance = 0.02)
  expect_equal(stats::cov(t(p[, 3, ])), V3, tolerance = 0.03)
})

test_that("A = 0.99 I gives a deviation half-life of about 69 years (2022, Section 5.4.1)", {
  expect_equal(round(log(0.5) / log(0.99)), 69)
  sp <- list(y = rbind(c(0.1, 0.2), c(0.2, -0.1)), A = diag(0.99, 2), Omega = matrix(0, 2, 2))
  p <- mortensor:::saint_var_paths(sp, 69, 1)
  expect_equal(p[, 69, 1] / c(0.2, -0.1), c(0.5, 0.5), tolerance = 1e-3)
})

test_that("the published VAR(1) of Danish spreads behaves as described (2011, Section 3.3)", {
  # A for Danish women and men and the spread of women in 2005 (level, slope, curvature)
  A_w <- matrix(c(0.6861, -0.1907, -0.2739, -0.1423, 0.8724, -0.1558, -0.2422, -0.1035, 0.5179), 3, byrow = TRUE)
  A_m <- matrix(c(0.7885, -0.1714, -0.1283, -0.2485, 0.6387, 0.0477, -0.0650, -0.0792, 0.9130), 3, byrow = TRUE)
  # "In both cases the A matrices give rise to stationary series"
  for (A in list(A_w, A_m)) expect_lt(max(Mod(eigen(A)$values)), 1)
  sp <- list(y = rbind(c(0.17, 0.06, -0.16)), A = A_w, Omega = matrix(0, 3, 3))
  p <- mortensor:::saint_var_paths(sp, 100, 1)[, , 1]
  # "The slope ... starts out positive but is forecasted to become negative before approaching zero"
  expect_gt(p[2, 1], 0)
  expect_true(any(p[2, 1:30] < 0))
  expect_lt(max(abs(p[, 100])), 0.01)
})

test_that("trend forecasts follow eqs. 4.15-4.16 with the indexes held fixed", {
  w <- saint_world(noise = TRUE, targets = "A")
  fit <- mt_fit(saint(reference = list(region = "ref"), by = "sex"), w$data)
  f <- fit$fits
  # no innovations, no drift: the indexes stay at their jump-off values
  for (g in seq_along(f$dynamics$groups)) {
    f$dynamics$groups[[g]]$drift[] <- 0
    f$dynamics$groups[[g]]$cov[] <- 0
    f$dynamics$groups[[g]]$sd_kappa <- 0
    f$dynamics$groups[[g]]$sd_zeta <- 0
  }
  h <- 4
  sim <- mortensor:::simulate_model(fit$spec, f, h, 2)
  tr <- f$trends[["female"]]
  nt <- ncol(tr$M)
  s2 <- tr$sigma2
  FT <- tr$baseline[, nt]
  GT <- exp(tr$zeta[nt])
  M0 <- (exp(s2 * tr$M[, nt]) - 1) / s2
  for (k in seq_len(h)) {
    M0 <- c(0, M0[-length(M0)] + FT[-length(FT)])
    expect_equal(sim[, k, match("ref.female", fit$strata$label), 1], FT / (1 + s2 * M0) + GT, tolerance = 1e-12)
  }
})

test_that("trend indexes follow a random walk with drift, or the error-correction model", {
  w <- saint_world(noise = TRUE, targets = "A")
  fit <- mt_fit(saint(reference = list(region = "ref"), by = "sex"), w$data)
  f <- fit$fits
  set.seed(2)
  p <- mortensor:::saint_index_paths(f, 10, 20000)
  dg <- f$dynamics$groups[["male"]]
  th <- f$trends[["male"]]$theta
  expect_equal(mean(p[[2]]$alpha[10, ]), unname(th[1, ncol(th)] + 10 * dg$drift[1]), tolerance = 1e-3)
  expect_equal(stats::var(p[[2]]$alpha[10, ]), 10 * dg$cov[1, 1], tolerance = 0.05)
  expect_equal(stats::sd(p[[2]]$zeta[10, ]), sqrt(10) * dg$sd_zeta, tolerance = 0.05)
  # error-correction model with known coefficients: lambda = 0.8
  ec <- list(A = matrix(c(-0.1, 0.1, 0, 0, 0, 0, -0.05, 0.15), 4), B = cbind(c(1, -1, 0, 0), c(0, 0, 1, -1)),
             CD = c(-0.4, 0.01), C_delta = c(-0.015, -0.015, 4e-4, 4e-4), Omega = diag(c(1e-4, 2e-4, 1e-7, 2e-7)),
             last = c(-3.7, -3.3, 0.11, 0.10), lead = 1L, other = 2L)
  expect_equal(sort(Re(eigen(diag(4) + ec$A %*% t(ec$B))$values)), c(0.8, 0.8, 1, 1))
  f$dynamics$ecm <- ec
  f0 <- f
  f0$dynamics$ecm$Omega <- matrix(0, 4, 4)
  p0 <- mortensor:::saint_index_paths(f0, 30, 1)
  # zero innovations: Y_T + C_delta h (eq. 5.8 with C_D the jump-off differences)
  expect_equal(p0[[1]]$alpha[, 1], -3.7 - 0.015 * (1:30), tolerance = 1e-12)
  expect_equal(p0[[2]]$beta[, 1], 0.10 + 4e-4 * (1:30), tolerance = 1e-12)
  set.seed(4)
  pe <- mortensor:::saint_index_paths(f, 60, 4000)
  gap <- pe[[1]]$alpha - pe[[2]]$alpha
  # the difference is stationary around C_D; each level is a random walk
  expect_equal(mean(gap[60, ]), -0.4, tolerance = 0.01)
  expect_equal(stats::var(gap[60, ]) / stats::var(gap[30, ]), 1, tolerance = 0.15)
  expect_gt(stats::var(pe[[1]]$alpha[60, ]) / stats::var(pe[[1]]$alpha[30, ]), 1.6)
})

test_that("the error-correction model is estimated by maximum likelihood given C_D and C_delta", {
  set.seed(8)
  n <- 40
  Y <- matrix(0, n, 4)
  Y[1, ] <- c(-3.7, -3.3, 0.11, 0.10)
  A <- matrix(c(-0.1, 0.1, 0, 0, 0, 0, -0.05, 0.15), 4)
  B <- cbind(c(1, -1, 0, 0), c(0, 0, 1, -1))
  for (t in 2:n) {
    Y[t, ] <- Y[t - 1, ] + A %*% (crossprod(B, Y[t - 1, ]) - c(-0.4, 0.01)) + c(-0.015, -0.015, 4e-4, 4e-4) +
      stats::rnorm(4, 0, c(0.01, 0.015, 3e-4, 4e-4))
  }
  ec <- mortensor:::saint_ecm_fit(Y[, 1], Y[, 2], Y[, 3], Y[, 4])
  expect_equal(ec$CD, c(Y[n, 1] - Y[n, 2], Y[n, 3] - Y[n, 4]))
  expect_equal(ec$C_delta, rep(c(mean(diff(Y[, 1])), mean(diff(Y[, 3]))), each = 2))
  expect_equal(unname(ec$lambda), c(1 + ec$A[1, 1] - ec$A[2, 1], 1 + ec$A[3, 2] - ec$A[4, 2]))
  # the coefficients maximise the concentrated Gaussian likelihood, -n/2 log det Omega(a)
  dY <- diff(Y) - matrix(ec$C_delta, n - 1, 4, byrow = TRUE)
  Z <- cbind(Y[-n, 1] - Y[-n, 2] - ec$CD[1], Y[-n, 3] - Y[-n, 4] - ec$CD[2])[, c(1, 1, 2, 2)]
  logdet <- function(a) as.numeric(determinant(crossprod(dY - Z * matrix(a, n - 1, 4, byrow = TRUE)))$modulus)
  a_hat <- ec$A[cbind(1:4, c(1, 1, 2, 2))]
  opt <- stats::optim(a_hat + 0.01, logdet, method = "BFGS", control = list(reltol = 1e-14))
  expect_lte(logdet(a_hat), opt$value + 1e-8)
  expect_equal(sort(Re(eigen(diag(4) + ec$A %*% t(ec$B))$values)), sort(unname(c(1, 1, ec$lambda))),
               tolerance = 1e-10)
})

test_that("saint() fits and forecasts: contract, pooled reference, parameters, error correction", {
  w <- saint_world(noise = TRUE, targets = "A", years = 2001:2012)
  r <- expect_mt_model(saint(reference = list(region = "ref"), by = "sex"), data = w$data)
  fc <- mt_forecast(r$fit, h = 3, n_draws = 10, seed = 1, predictive = TRUE)
  expect_true(all(is.finite(fc$rates)))
  p <- mt_parameters(r$fit)
  expect_setequal(unique(p$stratum), c("trend female", "trend male", "A.female", "A.male"))
  expect_true(all(c("sigma2", "alpha", "beta", "kappa", "zeta", "drift", "y2", "y5", "ar", "Omega") %in% p$parameter))
  expect_true(all(p[parameter == "ar"]$value >= 0 & p[parameter == "ar"]$value <= 0.99))
  expect_output(print(summary(r$fit)), "SAINT")
  # error correction between the sexes
  fe <- mt_fit(saint(reference = list(region = "ref"), by = "sex", ecm = "female"), w$data)
  pe <- mt_parameters(fe)
  expect_true("trend error correction" %in% pe$stratum)
  ec <- fe$fits$dynamics$ecm
  af <- fe$fits$trends$female$theta[1, ]
  expect_equal(ec$C_delta[1:2], rep((af[12] - af[1]) / 11, 2))
  expect_s3_class(mt_forecast(fe, h = 3, n_draws = 5, seed = 1), "mt_draws")
  # one trend from two pooled reference strata; both are targets of the pooled trend
  fp <- mt_fit(saint(reference = c("ref.female", "A.female"), spread_ar = 0.9), mt_subset(w$data, sex == "female"))
  D <- mt_array(mt_subset(w$data, sex == "female"), "deaths")
  E <- mt_array(mt_subset(w$data, sex == "female"), "exposure")
  tr <- mortensor:::saint_trend_fit(D[, , 1] + D[, , 2], E[, , 1] + E[, , 2], 40:100)
  expect_equal(fp$fits$trends[[1]]$theta, tr$theta)
  expect_false(any(vapply(fp$fits$spreads, is.null, logical(1))))
  expect_equal(unname(fp$fits$spreads[[1]]$A), diag(0.9, 4))
})

test_that("saint() checks ages and arguments", {
  w <- saint_world(noise = TRUE, targets = "A", years = 2001:2010)
  sp <- saint(reference = list(region = "ref"), by = "sex")
  # an open age interval or ages below 20 are refused
  x <- mt_example(deaths = "single", years = 2001:2010)
  expect_error(suppressWarnings(mt_fit(saint(reference = "female"), x)), class = "mt_error_model_fit")
  expect_error(mortensor:::saint_check_ages(mt_grid(x)$age), class = "mt_error_model_data")
  expect_error(mortensor:::saint_check_ages(single_ages(15, 100, open = FALSE)), class = "mt_error_model_data")
  expect_equal(mortensor:::saint_check_ages(single_ages(20, 100, open = FALSE)), 20:100)
  fx <- suppressWarnings(mt_fit(saint(reference = "female"), mt_subset(x, ages = c(30, 100))))
  expect_null(fx$fits$spreads[[1]])
  expect_error(saint(reference = 1), class = "mt_error_argument")
  expect_error(saint(reference = "ref", ecm = "female"), class = "mt_error_argument")
  expect_error(saint(reference = "ref", spread_ar = 1), class = "mt_error_argument")
  st <- mt_strata(w$data)
  expect_error(mortensor:::saint_strata(list(reference = "nope"), st, nrow(st)), class = "mt_error_argument")
  expect_error(mortensor:::saint_strata(list(reference = list(region = "ref"), by = "sex", ecm = "both"), st,
                                        nrow(st)), class = "mt_error_argument")
  expect_error(mortensor:::saint_strata(list(reference = c("ref.female"), by = "sex"), st, nrow(st)),
               class = "mt_error_argument")
  s <- mortensor:::saint_strata(list(reference = c(region = "ref"), by = "sex"), st, nrow(st))
  expect_identical(s$sole, st$region == "ref")
})
