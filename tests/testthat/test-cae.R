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
