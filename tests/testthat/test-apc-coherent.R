# Simulated data of the coherent APC model: S strata, I age groups, J years
coh_sim_data <- function(S = 3, I = 8, J = 15, M = 1, E = 2e5, rho = 0.7, sd_d = 0.05, seed = 42) {
  set.seed(seed)
  K <- M * (I - 1) + J
  theta <- seq(-7, -2, length.out = I)
  phi <- cumsum(rnorm(J, -0.02, 0.02))
  psi <- cumsum(rnorm(K, 0, 0.03))
  dv <- matrix(0, S, J)
  for (s in seq_len(S)) {
    dv[s, 1] <- rnorm(1, 0, sd_d / sqrt(1 - rho^2))
    for (t in 2:J) dv[s, t] <- rho * dv[s, t - 1] + rnorm(1, 0, sd_d)
  }
  dv <- sweep(dv, 2, colMeans(dv))
  a <- seq(0.2, -0.2, length.out = S)
  P <- array(0, c(I, J, S))
  for (s in seq_len(S)) for (i in seq_len(I)) for (j in seq_len(J)) {
    P[i, j, s] <- plogis(theta[i] + a[s] + phi[j] + psi[M * (I - i) + j] + dv[s, j])
  }
  Ex <- array(E, c(I, J, S))
  list(D = array(rbinom(length(Ex), Ex, P), dim(Ex)), E = Ex, P = P, dv = dv, M = M)
}

test_that("the coherent APC sampler keeps the constraints in every draw", {
  sim <- coh_sim_data(J = 12)
  st <- coh_setup(sim$D, sim$E, sim$M, c(1L, 1L, 1L))
  set.seed(1)
  ch <- coh_chain(st, 150, 50, 1)
  S <- st$S
  J <- st$J
  dsum <- vapply(seq_len(J), function(j) max(abs(rowSums(ch$x[, st$idx$d[(seq_len(S) - 1L) * J + j]]))), 0)
  expect_lt(max(dsum), 1e-10)
  for (s in seq_len(S)) expect_lt(max(abs(rowSums(ch$x[, st$idx$th[(s - 1L) * st$I + seq_len(st$I)]]))), 1e-10)
  expect_lt(max(abs(rowSums(ch$x[, st$idx$ph]))), 1e-10)
  expect_lt(max(abs(rowSums(ch$x[, st$idx$ps]))), 1e-10)
  expect_gt(ch$accept, 0.5)
  # reference stratum: its deviation and level stay at zero
  st2 <- coh_setup(sim$D, sim$E, sim$M, c(1L, 1L, 1L), cohort_dev = TRUE, ref = 2L)
  ch2 <- coh_chain(st2, 100, 50, 1)
  expect_lt(max(abs(ch2$x[, st2$idx$d[J + seq_len(J)]])), 1e-10)
  expect_lt(max(abs(ch2$x[, st2$idx$c[st2$K + seq_len(st2$K)]])), 1e-10)
  expect_lt(max(abs(ch2$x[, st2$idx$a[2]])), 1e-10)
  expect_true(all(is.finite(ch2$h[, "rho_c"])))
})

test_that("with fixed hyperparameters the sampler agrees with importance sampling", {
  skip_on_cran()
  # small counts, so that the posterior is not Gaussian
  sim <- coh_sim_data(S = 2, I = 4, J = 6, E = 300, seed = 7)
  st <- coh_setup(sim$D, sim$E, sim$M, c(1L, 1L, 1L))
  h <- list(kappa = 2, lambda = 50, nu = 50, lambda_d = 100, rho = 0.6, lambda_c = 10, rho_c = 0.5)
  Qp <- coh_prior(st, h)
  set.seed(3)
  ch <- coh_chain(st, 12000, 1000, 1, init = list(h = h), fix_hyper = TRUE)
  eta_draws <- as.matrix(Matrix::tcrossprod(ch$x, st$X))
  eta_mcmc <- colMeans(eta_draws)
  batches <- split(seq_len(nrow(eta_draws)), cut(seq_len(nrow(eta_draws)), 40))
  bm <- vapply(batches, function(i) colMeans(eta_draws[i, , drop = FALSE]), numeric(ncol(eta_draws)))
  se_mcmc <- apply(bm, 1, sd) / sqrt(length(batches))
  # self-normalised importance sampling from the Gaussian approximation at
  # the mode, its covariance inflated by 1.1^2
  sym <- Matrix::Cholesky(Qp, perm = TRUE, LDL = FALSE, super = FALSE)
  ap <- coh_approx(st, coh_start(st, Qp, sym), Qp, sym)
  n <- 20000
  f <- 1.1
  Z <- vapply(seq_len(n), function(r) {
    z <- as.vector(Matrix::solve(ap$L, Matrix::solve(ap$L, rnorm(st$P), system = "Lt"), system = "Pt"))
    ap$m + f * ap$correct(z)
  }, numeric(st$P))
  lw <- vapply(seq_len(n), function(r) {
    v <- Z[, r] - ap$m
    coh_lp(st, Z[, r], Qp) + 0.5 * sum(v * as.vector(ap$H %*% v)) / f^2
  }, 0)
  w <- exp(lw - max(lw))
  w <- w / sum(w)
  ess <- 1 / sum(w^2)
  EZ <- as.matrix(st$X %*% Z)
  eta_is <- as.vector(EZ %*% w)
  se_is <- sqrt(as.vector((EZ - eta_is)^2 %*% w) / ess)
  expect_gt(ess, 5000)
  # the linear predictors of all 48 cells agree within Monte Carlo error
  expect_lt(max(abs(eta_mcmc - eta_is) / sqrt(se_mcmc^2 + se_is^2)), 4)
})

test_that("forecasts of the deviations and effects have the variances of their processes", {
  set.seed(11)
  S <- 4
  rho <- 0.8
  prec <- 25
  h <- 6
  sims <- replicate(20000, coh_extend_ar1(matrix(0, S, 3), h, rho, prec, NULL)[1, 3 + c(1, h)])
  # independent AR(1) series centred over the strata
  v <- function(k) (1 - 1 / S) * (1 - rho^(2 * k)) / ((1 - rho^2) * prec)
  expect_equal(apply(sims, 1, var), c(v(1), v(h)), tolerance = 0.04)
  ref <- replicate(20000, coh_extend_ar1(matrix(0, S, 3), h, rho, prec, 2L)[c(1, 2), 3 + h])
  expect_equal(var(ref[1, ]), (1 - rho^(2 * h)) / ((1 - rho^2) * prec), tolerance = 0.04)
  expect_identical(max(abs(ref[2, ])), 0)
  rw1 <- replicate(20000, coh_extend_rw(c(0, 0), h, prec, 1L)[2 + h])
  expect_equal(var(rw1), h / prec, tolerance = 0.04)
  rw2 <- replicate(20000, coh_extend_rw(c(0, 0), h, prec, 2L)[2 + h])
  expect_equal(var(rw2), sum((h - seq_len(h) + 1)^2) / prec, tolerance = 0.04)
})

test_that("the coherent APC model recovers simulated deviations", {
  skip_on_cran()
  sim <- coh_sim_data(I = 10, J = 20)
  st <- coh_setup(sim$D, sim$E, sim$M, c(1L, 1L, 1L))
  set.seed(2)
  ch <- coh_chain(st, 800, 300, 2)
  dhat <- matrix(colMeans(ch$x[, st$idx$d]), st$S, st$J, byrow = TRUE)
  expect_gt(cor(as.vector(sim$dv), as.vector(dhat)), 0.75)
  eta <- as.matrix(Matrix::tcrossprod(ch$x, coh_cell_design(st)))
  q <- apply(plogis(eta), 2, quantile, c(0.005, 0.995))
  expect_gt(mean(sim$P >= q[1, ] & sim$P <= q[2, ]), 0.9)
})

test_that("apc_bayes(coherent = TRUE) fits and forecasts all strata jointly", {
  skip_on_cran()
  x <- mt_example(years = 2000:2019)
  sp <- apc_bayes(coherent = TRUE, mcmc = mcmc_control("test"), seed = 1)
  expect_s3_class(sp, "mt_spec_apc_coherent")
  expect_true(sp$joint)
  fit <- mt_fit(sp, x, prep = mt_prep(step_ungroup()))
  expect_identical(fit$diagnostics$status, "ok")
  expect_lt(fit$diagnostics$max_rhat, 1.1)
  expect_equal(dim(fitted(fit)), c(length(fit$age), 20L, 2L))
  fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
  expect_false(anyNA(fc$rates))
  # the same seed gives the same fit
  fit2 <- mt_fit(sp, x, prep = mt_prep(step_ungroup()))
  expect_equal(fit$fits$groups$all$x, fit2$fits$groups$all$x)
  p <- mt_parameters(fit)
  expect_true(all(c("rho", "lambda_d", "rhat", "deviation", "period") %in% p$parameter))
  fr <- mt_fit(apc_bayes(coherent = TRUE, reference = list(sex = "male"), mcmc = mcmc_control("test"), seed = 1), x,
               prep = mt_prep(step_ungroup()))
  pr <- mt_parameters(fr)
  expect_lt(max(abs(pr[stratum == "male" & parameter == "deviation"]$value)), 1e-10)
})

test_that("apc_bayes(coherent = TRUE) checks its arguments", {
  expect_error(apc_bayes(coherent = TRUE, overdispersion = TRUE), class = "mt_error_argument")
  expect_error(apc_bayes(common_by = "sex"), class = "mt_error_argument")
  expect_error(apc_bayes(coherent = TRUE, reference = "male"), class = "mt_error_argument")
  x <- mt_example(years = 2000:2019)
  expect_error(mt_fit(apc_bayes(coherent = TRUE, common_by = "sex", mcmc = mcmc_control("test")), x,
                      prep = mt_prep(step_ungroup())), class = "mt_error_argument")
  expect_error(mt_fit(apc_bayes(coherent = TRUE, reference = list(sex = "none"), mcmc = mcmc_control("test")), x,
                      prep = mt_prep(step_ungroup())), class = "mt_error_argument")
})
