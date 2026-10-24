# Common age effect model (Kleinow 2015)

# Rates from a CAE model with K populations, n ages, T years and order 2:
# beta orthonormal, kappa of each population with orthogonal rows (so that
# beta diagonalises every m m' exactly) and decreasing variation by order
cae_truth <- function(n = 8, nt = 30, K = 3, seed = 1) {
  set.seed(seed)
  beta <- qr.Q(qr(cbind(seq(1, 2, length.out = n), sin(seq_len(n)))))
  beta <- sweep(beta, 2, sign(colSums(beta)), "*")
  a <- sapply(seq_len(K), function(i) -9 + 0.09 * (seq_len(n) - 1) * 5 + 0.1 * i)
  kappa <- array(NA_real_, c(2, nt, K))
  for (i in seq_len(K)) {
    k1 <- cumsum(rnorm(nt, -0.3 - 0.05 * i, 0.15))
    k2 <- cumsum(rnorm(nt, 0, 0.1))
    k1 <- k1 - mean(k1)
    k2 <- k2 - mean(k2)
    k2 <- k2 - sum(k2 * k1) / sum(k1^2) * k1
    k2 <- k2 / sqrt(sum(k2^2)) * 0.2 * sqrt(sum(k1^2))
    kappa[, , i] <- rbind(k1, k2)
  }
  logm <- array(NA_real_, c(n, nt, K))
  for (i in seq_len(K)) logm[, , i] <- a[, i] + beta %*% kappa[, , i]
  list(beta = beta, a = a, kappa = kappa, logm = logm)
}

fit_arrays <- function(spec, D, E, strata = NULL, hierarchy = NULL) {
  mortensor:::fit_model(spec, D, E, NULL, seq_len(dim(D)[2]), strata = strata, hierarchy = hierarchy)
}

cae_data <- function() mt_subset(mt_example_regions(years = 2000:2019, regions = 2), ages = c(60, Inf))

test_that("cae() with one population is the Lee-Carter model by SVD", {
  x <- mt_subset(mt_regroup(mt_example(), age = "deaths"), sex == "female")
  f_cae <- mt_fit(cae(), x)
  f_lc <- mt_fit(lee_carter(method = "svd", adjust = "none"), x)
  expect_equal(fitted(f_cae), fitted(f_lc), tolerance = 1e-10)
  # the same age pattern and index up to the scale of the constraints
  b <- f_lc$fits[[1]]$b
  beta <- f_cae$fits$beta[[1]][, 1]
  expect_equal(beta, b / sqrt(sum(b^2)), tolerance = 1e-10)
  expect_equal(f_cae$fits$kappa[1, , 1] * sum(beta), f_lc$fits[[1]]$k, tolerance = 1e-8, ignore_attr = TRUE)
  # the multivariate random walk of one index is Lee-Carter's random walk with drift: identical draws
  fc_cae <- mt_forecast(f_cae, h = 6, n_draws = 50, seed = 3)
  fc_lc <- mt_forecast(f_lc, h = 6, n_draws = 50, seed = 3)
  expect_equal(fc_cae$rates, fc_lc$rates, tolerance = 1e-10)
})

test_that("cae() satisfies its identification constraints", {
  fit <- mt_fit(cae(order = 3), cae_data())
  f <- fit$fits
  beta <- f$beta[[1]]
  expect_equal(crossprod(beta), diag(3), tolerance = 1e-10)
  expect_true(all(colSums(beta) >= 0))
  # period effects of every population sum to zero over the years
  expect_lt(max(abs(apply(f$kappa, c(1, 3), sum))), 1e-9)
  # kappa = t(beta) m, with m the centred log rates
  y <- log(pmax(fit$train_D, 0.5) / fit$train_E)
  m <- sweep(y, c(1, 3), apply(y, c(1, 3), mean))
  expect_equal(f$kappa[, , 5], crossprod(beta, m[, , 5]), tolerance = 1e-10, ignore_attr = TRUE)
  # the columns minimise the off-diagonal criterion: rotations in any plane do not lower it
  Q <- lapply(seq_len(dim(m)[3]), function(s) tcrossprod(m[, , s]))
  crit <- function(B) {
    sum(vapply(Q, function(q) {
      a <- crossprod(B, q %*% B)
      sum(a^2) - sum(diag(a)^2)
    }, 0))
  }
  B <- mortensor:::cpca_ls(Q)$beta
  c0 <- crit(B)
  for (pq in list(c(1, 2), c(1, 3), c(2, 5), c(4, 7))) {
    for (th in c(-1e-3, 1e-3)) {
      R <- diag(ncol(B))
      R[pq, pq] <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2)
      expect_gt(crit(B %*% R), c0)
    }
  }
})

test_that("cae() recovers the parameters of data simulated from the model", {
  tr <- cae_truth()
  n <- dim(tr$logm)[1]
  nt <- dim(tr$logm)[2]
  # exact rates: the estimates are the true parameters
  E <- array(1e9, dim(tr$logm))
  f <- fit_arrays(cae(order = 2), exp(tr$logm) * E, E)
  expect_equal(f$beta[[1]], tr$beta, tolerance = 1e-8)
  expect_equal(f$a, tr$a, tolerance = 1e-8)
  expect_equal(f$kappa, tr$kappa, tolerance = 1e-8, ignore_attr = TRUE)
  # Poisson deaths with large exposures: estimates close to the truth
  set.seed(2)
  E <- array(2e6, dim(tr$logm))
  D <- array(rpois(length(E), exp(tr$logm) * E), dim(E))
  f <- fit_arrays(cae(order = 2), D, E)
  expect_lt(max(abs(f$beta[[1]] - tr$beta)), 0.05)
  expect_lt(max(abs(f$kappa[1, , ] - tr$kappa[1, , ])) / diff(range(tr$kappa[1, , ])), 0.02)
  # the random walk: drift and covariance of the year-to-year changes of the
  # period effects, stacked as order 1 and 2 of population 1, then of 2 and 3
  K <- rbind(f$kappa[, , 1], f$kappa[, , 2], f$kappa[, , 3])
  dK <- t(diff(t(K)))
  expect_equal(f$last, K[, nt], ignore_attr = TRUE)
  expect_equal(f$drift, rowMeans(dK))
  expect_equal(f$sigma, stats::cov(t(dK)), tolerance = 1e-10, ignore_attr = TRUE)
  # mean squared error and approximate BIC of eq. (4) and (5)
  res <- log(f$fitted) - log(pmax(D, 0.5) / E)
  N <- length(res)
  expect_equal(f$mse, mean(res^2))
  expect_equal(f$bic, N * log(mean(res^2)) + 2 * (3 * nt + n) * log(N))
})

test_that("cae() draws of the populations are correlated as estimated", {
  x <- cae_data()
  fit <- mt_fit(cae(order = 1), mt_subset(x, region == "R1", sex == "female"))
  f <- fit$fits
  fc <- mt_forecast(fit, h = 1, n_draws = 4000, seed = 1)
  # period effect of each population and draw in the first forecast year
  y <- log(fc$rates[, 1, , ]) - as.vector(f$a)
  k1 <- apply(y, 2:3, function(v) sum(v * f$beta[[1]][, 1]))
  sim_cor <- stats::cor(t(k1))
  est_cor <- stats::cov2cor(f$sigma)
  expect_lt(max(abs(sim_cor - est_cor)), 0.06)
  # covariance: innovations plus the uncertainty of the drift
  sim_cov <- stats::cov(t(k1))
  expect_equal(sim_cov, f$sigma * (1 + 1 / f$n_changes), tolerance = 0.08, ignore_attr = TRUE)
  expect_true(any(abs(est_cor[upper.tri(est_cor)]) > 0.3))
  # each population starts from its last period effect and moves by its drift
  nt <- dim(f$kappa)[2]
  dk <- f$kappa[1, nt, ] - f$kappa[1, 1, ]
  expect_equal(rowMeans(k1), f$kappa[1, nt, ] + dk / (nt - 1), tolerance = 0.02, ignore_attr = TRUE)
  # the correlations are those of the year-to-year changes of each population's period effect
  expect_equal(est_cor, stats::cor(diff(f$kappa[1, , ])), ignore_attr = TRUE)
})

test_that("cae(common_by = ) estimates the age effects of each group from its strata", {
  x <- cae_data()
  fit <- mt_fit(cae(common_by = "sex"), x)
  expect_length(fit$fits$beta, 2)
  women <- mt_fit(cae(), mt_subset(x, sex == "female"))
  is_f <- fit$strata$sex == "female"
  expect_equal(fit$fits$beta[[fit$fits$group[which(is_f)[1]]]], women$fits$beta[[1]], tolerance = 1e-10)
  expect_equal(fitted(fit)[, , is_f], fitted(women), tolerance = 1e-10, ignore_attr = TRUE)
  expect_error(mt_fit(cae(common_by = "nonsense"), x), class = "mt_error_model_fit")
  expect_mt_model(cae(order = 2, common_by = "sex"), data = x, n_draws = 10)
})

test_that("cae() estimates the age effects from the bottom strata of a hierarchy", {
  x <- cae_data()
  h <- mt_hierarchy(quintile = "all")
  xa <- mt_aggregate_strata(x, h)
  fa <- mt_fit(cae(order = 2), xa)
  fb <- mt_fit(cae(order = 2), x)
  expect_equal(fa$fits$beta, fb$fits$beta, tolerance = 1e-10)
  expect_equal(nrow(fa$strata), nrow(x$strata) + 4)
  rc <- mt_forecast(mt_fit(reconciled(cae(), h, method = "bottom_up"), x), h = 2, n_draws = 20, seed = 1)
  expect_lt(coherence_gap(rc, h), 1e-8)
})

test_that("cae() needs more years than age groups", {
  x <- mt_regroup(mt_example(years = 2005:2015), age = "deaths")
  expect_error(mt_fit(cae(), x), class = "mt_error_model_fit")
  D <- mt_array(x, "deaths")
  expect_error(fit_arrays(cae(), D, mt_array(x, "exposure")), class = "mt_error_model_data")
  expect_error(cae(order = 0), class = "mt_error_argument")
  expect_error(cae(common_by = 1), class = "mt_error_argument")
})

test_that("mt_parameters() and summary() describe cae() fits", {
  fit <- mt_fit(cae(order = 2), cae_data())
  p <- mt_parameters(fit)
  expect_setequal(unique(p$parameter), c("beta1", "beta2", "a", "kappa1", "kappa2", "drift_kappa1", "drift_kappa2",
                                         "sd_kappa1", "sd_kappa2", "cor_kappa1", "cor_kappa2", "mse", "bic"))
  expect_equal(p[parameter == "bic", value], fit$fits$bic)
  cr <- p[parameter == "cor_kappa1" & stratum == fit$strata$label[1]]
  expect_equal(cr[index == fit$strata$label[1], value], 1)
  expect_match(mortensor:::model_assumptions(cae(order = 2)), "Kleinow 2015")
  expect_output(print(summary(fit)), "Common age effect")
})

# index = "coherent": the forecasting model of Enchev, Kleinow and Cairns (2017, Table 10, Model 3)

# Exact rates of a CAE model of order 2 whose period effects follow that
# model: kappa1 of K populations a random walk with a common drift `d` and
# correlated innovations, kappa2 AR(1) processes with coefficients `phi`
cae_coherent_truth <- function(n = 8, nt = 1000, d = -0.03, phi = c(0.3, 0.6, 0.8), seed = 1) {
  set.seed(seed)
  K <- length(phi)
  beta <- qr.Q(qr(cbind(seq(1, 2, length.out = n), sin(seq_len(n)))))
  beta <- sweep(beta, 2, sign(colSums(beta)), "*")
  a <- sapply(seq_len(K), function(i) -6 + 0.3 * (seq_len(n) - 1) + 0.1 * i)
  S1 <- 0.15^2 * (0.5 + 0.5 * diag(K))
  S2 <- 0.1^2 * (0.3 + 0.7 * diag(K))
  k1 <- t(apply(d + t(chol(S1)) %*% matrix(rnorm(K * nt), K), 1, cumsum))
  e2 <- t(chol(S2)) %*% matrix(rnorm(K * nt), K)
  k2 <- e2
  for (t in 2:nt) k2[, t] <- phi * k2[, t - 1] + e2[, t]
  logm <- array(NA_real_, c(n, nt, K))
  for (i in seq_len(K)) logm[, , i] <- a[, i] + beta %*% rbind(k1[i, ], k2[i, ])
  list(beta = beta, d = d, phi = phi, S1 = S1, S2 = S2, logm = logm)
}

test_that("cae(index = \"coherent\") with one population and order 1 equals the random walk with drift", {
  x <- mt_subset(mt_regroup(mt_example(), age = "deaths"), sex == "female")
  f_rw <- mt_fit(cae(), x)
  f_co <- mt_fit(cae(index = "coherent"), x)
  expect_equal(f_co$fits$drift, f_rw$fits$drift)
  expect_equal(f_co$fits$common_drift, f_rw$fits$drift)
  expect_equal(f_co$fits$sigma, f_rw$fits$sigma)
  fc_rw <- mt_forecast(f_rw, h = 6, n_draws = 50, seed = 3)
  fc_co <- mt_forecast(f_co, h = 6, n_draws = 50, seed = 3)
  expect_equal(fc_co$rates, fc_rw$rates, tolerance = 1e-10)
  # order 2: the first period effect keeps its own drift, the second is an AR(1) through the origin
  f2 <- mt_fit(cae(order = 2, index = "coherent"), x)$fits
  k <- f2$kappa[, , 1]
  nt <- ncol(k)
  expect_equal(f2$drift, c((k[1, nt] - k[1, 1]) / (nt - 1), 0), ignore_attr = TRUE)
  expect_equal(f2$phi, c(1, sum(k[2, -1] * k[2, -nt]) / sum(k[2, -nt]^2)))
})

test_that("cae(index = \"coherent\") estimates a common drift, AR(1) coefficients and one covariance", {
  x <- cae_data()
  fit <- mt_fit(cae(order = 2, index = "coherent"), x)
  f <- fit$fits
  ns <- nrow(fit$strata)
  nt <- dim(f$kappa)[2]
  K <- matrix(aperm(f$kappa, c(1, 3, 2)), 2 * ns)
  # the common drift: mean year-to-year change of the first period effects, every population with the same weight
  d <- mean((f$kappa[1, nt, ] - f$kappa[1, 1, ]) / (nt - 1))
  expect_equal(f$common_drift, d)
  expect_equal(f$drift, rep(c(d, 0), ns))
  phi <- vapply(seq_len(ns), function(s) {
    k <- f$kappa[2, , s]
    sum(k[-1] * k[-nt]) / sum(k[-nt]^2)
  }, 0)
  expect_equal(f$phi, as.vector(rbind(1, phi)))
  expect_true(all(abs(phi) < 0.99))
  R <- K[, -1] - f$phi * K[, -nt] - f$drift
  expect_equal(f$sigma, tcrossprod(R) / (nt - 2))
  # the same age effects, period effects and fitted rates as index = "rwd"
  rw <- mt_fit(cae(order = 2), x)
  expect_equal(fitted(fit), fitted(rw))
  expect_equal(f$kappa, rw$fits$kappa)
  # first period effects with equal drifts: the common drift is that drift, the covariance that of "rwd"
  k1 <- f$kappa[1, , , drop = FALSE]
  own <- (k1[1, nt, ] - k1[1, 1, ]) / (nt - 1)
  k1[1, , ] <- k1[1, , ] - outer(seq_len(nt) - 1, own - d)
  eq <- mortensor:::cae_index_fit(k1, f$group, f$members, "coherent")
  eq_rw <- mortensor:::cae_index_fit(k1, f$group, f$members, "rwd")
  expect_equal(eq$drift, rep(d, ns))
  expect_equal(eq_rw$drift, rep(d, ns))
  expect_equal(eq$sigma, eq_rw$sigma)
  expect_equal(eq$sigma, rw$fits$sigma[seq(1, 2 * ns, by = 2), seq(1, 2 * ns, by = 2)])
})

test_that("cae(index = \"coherent\") recovers the parameters of data simulated from its forecasting model", {
  tr <- cae_coherent_truth()
  E <- array(1e20, dim(tr$logm))
  f <- fit_arrays(cae(order = 2, index = "coherent"), exp(tr$logm) * E, E)
  # age effects up to a small rotation: the sample period effects of the two orders are not exactly orthogonal
  expect_lt(max(abs(f$beta[[1]] - tr$beta)), 1e-3)
  n <- dim(tr$logm)[2] - 1
  expect_lt(abs(f$common_drift - tr$d), 4 * sqrt(sum(tr$S1) / 9 / n))
  expect_lt(max(abs(f$phi[c(2, 4, 6)] - tr$phi)), 4 * sqrt((1 - 0.3^2) / n))
  expect_equal(sqrt(diag(f$sigma)[c(1, 3, 5)]), sqrt(diag(tr$S1)), tolerance = 0.1)
  expect_equal(sqrt(diag(f$sigma)[c(2, 4, 6)]), sqrt(diag(tr$S2)), tolerance = 0.1)
  expect_lt(max(abs(stats::cov2cor(f$sigma[c(1, 3, 5), c(1, 3, 5)]) - stats::cov2cor(tr$S1))), 0.1)
  expect_lt(max(abs(stats::cov2cor(f$sigma[c(2, 4, 6), c(2, 4, 6)]) - stats::cov2cor(tr$S2))), 0.1)
  # innovations of the two orders are independent in the simulation
  expect_lt(max(abs(stats::cov2cor(f$sigma)[c(1, 3, 5), c(2, 4, 6)])), 0.1)
})

test_that("cae(index = \"coherent\") draws follow the forecasting model", {
  x <- cae_data()
  fit <- mt_fit(cae(order = 2, index = "coherent"), mt_subset(x, region == "R1", sex == "female"))
  f <- fit$fits
  h <- 8
  fc <- mt_forecast(fit, h = h, n_draws = 6000, seed = 1)
  # period effects of the draws [order and stratum, draw] in the last forecast year
  y <- log(fc$rates[, h, , ]) - as.vector(f$a)
  ns <- ncol(f$a)
  kd <- matrix(NA_real_, 2 * ns, dim(y)[3])
  for (s in seq_len(ns)) kd[(s - 1) * 2 + 1:2, ] <- crossprod(f$beta[[1]], y[, s, ])
  # mean phi^h last + drift (1 + phi + ... + phi^(h - 1)); covariance sigma[r, s] times the sum over k < h of
  # (phi_r phi_s)^k, plus h^2 times the variance of the common drift for the first period effects
  g <- function(z) ifelse(abs(z - 1) < 1e-12, h, (1 - z^h) / (1 - z))
  mu <- f$phi^h * f$last + f$drift * g(f$phi)
  V <- f$sigma * g(outer(f$phi, f$phi))
  first <- f$drift_group > 0
  V[first, first] <- V[first, first] + h^2 * drop(f$drift_weights %*% f$sigma %*% t(f$drift_weights)) / f$n_changes
  sdv <- sqrt(diag(V))
  expect_lt(max(abs(rowMeans(kd) - mu) / sdv), 0.06)
  expect_equal(sqrt(diag(stats::cov(t(kd)))), sdv, tolerance = 0.05)
  expect_lt(max(abs(stats::cor(t(kd)) - stats::cov2cor(V))), 0.05)
  # the common drift cancels in the differences of the first period effects: their expected value stays at
  # the last fitted year, and their variance is h times the variance of the differences of the innovations
  i <- which(first)[1:2]
  w <- c(1, -1)
  dk <- drop(w %*% kd[i, ])
  expect_lt(abs(mean(dk) - sum(w * f$last[i])) / stats::sd(dk), 0.06)
  expect_equal(stats::var(dk), h * drop(w %*% f$sigma[i, i] %*% w), tolerance = 0.06)
})

test_that("cae(index = \"coherent\") with common_by and hierarchies gives each group its own drift", {
  x <- cae_data()
  fit <- mt_fit(cae(common_by = "sex", index = "coherent"), x)
  expect_length(fit$fits$common_drift, 2)
  women <- mt_fit(cae(index = "coherent"), mt_subset(x, sex == "female"))
  g <- fit$fits$group[which(fit$strata$sex == "female")[1]]
  expect_equal(fit$fits$common_drift[g], women$fits$common_drift)
  expect_equal(unique(fit$fits$drift[fit$fits$group == g]), women$fits$common_drift)
  # aggregates take the common drift of the bottom strata
  h <- mt_hierarchy(quintile = "all")
  fa <- mt_fit(cae(order = 2, index = "coherent"), mt_aggregate_strata(x, h))
  fb <- mt_fit(cae(order = 2, index = "coherent"), x)
  expect_equal(fa$fits$common_drift, fb$fits$common_drift)
  expect_equal(fa$fits$drift[fa$fits$drift_group > 0], rep(fb$fits$common_drift, nrow(fa$strata)))
  rc <- mt_forecast(mt_fit(reconciled(cae(index = "coherent"), h, method = "bottom_up"), x), h = 2, n_draws = 20,
                    seed = 1)
  expect_lt(coherence_gap(rc, h), 1e-8)
  expect_mt_model(cae(order = 2, common_by = "sex", index = "coherent"), data = x, n_draws = 10)
})

test_that("cae(index = \"coherent\") bounds the AR(1) coefficients and reports its parameters", {
  # exact rates of order 2 whose second period effects are slow waves: phi above 0.99
  n <- 6
  nt <- 100
  beta <- qr.Q(qr(cbind(seq(1, 2, length.out = n), cos(seq_len(n)))))
  tt <- seq_len(nt)
  logm <- array(NA_real_, c(n, nt, 2))
  for (i in 1:2) logm[, , i] <- -5 + beta %*% rbind(-0.1 * i * tt, (0.5 + i) * cos(pi * tt / nt))
  E <- array(1e12, dim(logm))
  expect_warning(f <- fit_arrays(cae(order = 2, index = "coherent"), exp(logm) * E, E), class = "mt_warning_cae")
  expect_equal(f$phi[c(2, 4)], c(0.99, 0.99))
  expect_true(all(f$bounded[c(2, 4)]))
  expect_error(cae(index = "nonsense"), class = "mt_error_argument")
  expect_identical(cae()$params$index, "rwd")
  fit <- mt_fit(cae(order = 2, index = "coherent"), cae_data())
  p <- mt_parameters(fit)
  expect_true(all(c("phi_kappa2", "drift_kappa1", "drift_kappa2") %in% p$parameter))
  expect_false("phi_kappa1" %in% p$parameter)
  expect_equal(unique(p[parameter == "drift_kappa1", value]), fit$fits$common_drift)
  expect_equal(unique(p[parameter == "drift_kappa2", value]), 0)
  expect_equal(p[parameter == "phi_kappa2", value], fit$fits$phi[seq(2, length(fit$fits$phi), by = 2)])
  expect_match(mortensor:::model_assumptions(cae(order = 2, index = "coherent")), "Enchev")
  expect_match(fit$spec$label, "coherent")
})
