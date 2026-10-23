# Weights and joint and individual components of coda_causes() (two-step
# model 2S-CoDA of Kjærgaard, Ergemen, Kallestrup-Lamb, Oeppen and
# Lindahl-Jacobsen 2019, Section 2.4). The paper's data (France, INED;
# the Netherlands, Statistics Netherlands) are not available, so its numbers
# are not reproduced; the tests check limiting cases, identities and
# recovery on data simulated from the model.

closure_2s <- function(M) sweep(M, 2, colSums(M), "/")

# Data from the two-step model (eq. 12) for one group with three causes,
# single ages 50-70 (closed range, so the survivors at 71 are a part) and
# `nt` years. The components are chosen so that the weighted decomposition
# separates them exactly: individual loadings sum to zero on their cause and
# are orthogonal to the joint loadings in the metric of the part weights
# (eq. 10, found by fixed-point iteration, as they depend on the deaths);
# the indexes have mean zero and the individual indexes are orthogonal to the
# joint index in the metric of the year weights (eq. 11). `ind`: scale of
# the individual components. Deaths are the rates times exposures of 1e6, or
# Poisson counts with these means.
sim_2s <- function(rho = 0.2, weights = "deaths", nt = 30, ind = 0.25, poisson = FALSE, seed = 1) {
  ages <- single_ages(50, 70, open = FALSE)
  na <- length(ages)
  nc <- 3
  np <- na * nc + 1
  x <- (seq_len(na) - 1) / (na - 1)
  alpha <- c(log(0.5) + 3 * x, log(0.3) + 2.5 * x, log(0.2) + 2 * x, 3)
  alpha <- alpha - mean(alpha)
  bJ <- c(-1 + 2 * x, -0.5 + x, 0.3 - x, 0.5)
  bJ <- bJ - mean(bJ)
  bJ <- bJ / sqrt(sum(bJ^2))
  tt <- seq_len(nt)
  kJ <- 3 * (tt - mean(tt)) / nt
  wt <- if (is.null(rho)) rep(1, nt) else rho * (1 - rho)^(nt - tt)
  kI <- vapply(seq_len(nc), function(i) {
    v <- sin(tt * i / 3) + 0.2 * i * ((tt - mean(tt)) / nt)^2
    v <- v - mean(v)
    v <- v - sum(wt * v * kJ) / sum(wt * kJ^2) * kJ
    ind * v / max(abs(v))
  }, numeric(nt))
  rows <- lapply(seq_len(nc), function(i) (i - 1) * na + seq_len(na))
  wp <- rep(1, np)
  for (it in 1:100) {
    BI <- matrix(0, np, nc)
    for (i in seq_len(nc)) {
      A <- rbind(1, (wp * bJ)[rows[[i]]])
      v <- cos(pi * x * (i + 1)) + 0.3 * x^2
      v <- v - drop(t(A) %*% solve(A %*% t(A), A %*% v))
      BI[rows[[i]], i] <- 0.6 * v / sqrt(sum(v^2))
    }
    Z <- alpha + outer(bJ, kJ) + BI %*% t(kI)
    wp_new <- if (weights == "deaths") rowMeans(closure_2s(exp(Z))) else rep(1, np)
    wp_new <- wp_new / sum(wp_new)
    done <- max(abs(wp_new - wp)) < 1e-15
    wp <- wp_new
    if (done) break
  }
  R <- mortensor:::md_rates(Z, mortensor:::coda_layout(ages, "regression"), nc)
  Dc <- lapply(seq_len(nc), function(i) as.vector(R[, , i]) * 1e6)
  if (poisson) {
    set.seed(seed)
    Dc <- lapply(Dc, function(d) as.numeric(stats::rpois(length(d), d)))
  }
  grid <- expand.grid(age_lo = ages$lower, period = 2000 + tt - 1)
  long <- rbindlist(c(list(data.table(grid, cause = "all", deaths = Reduce(`+`, Dc))),
                      lapply(seq_len(nc), function(i) data.table(grid, cause = paste0("c", i), deaths = Dc[[i]]))))
  long$age_hi <- long$age_lo + 1
  ex <- data.table(long[, c("age_lo", "age_hi", "period", "cause")], exposure = 1e6)
  d <- withCallingHandlers(
    mt_data(deaths = long, exposure = ex, by = "cause", age = c("age_lo", "age_hi"), period = "period"),
    mt_warning_noninteger = function(w) invokeRestart("muffleWarning"))
  list(data = d, alpha = alpha, bJ = bJ, kJ = kJ, BI = BI, kI = kI, Z = Z, wp = wp, wt = wt)
}

# Largest relative gap between the sum of the causes and all causes, per sex
cause_gap <- function(arr, st) {
  max(vapply(unique(st$sex), function(sx) {
    i <- st$stratum[st$sex == sx & st$cause == "all"]
    j <- st$stratum[st$sex == sx & st$cause != "all"]
    sums <- if (length(dim(arr)) == 4) apply(arr[, , j, , drop = FALSE], c(1, 2, 4), sum) else
      apply(arr[, , j, drop = FALSE], 1:2, sum)
    tot <- if (length(dim(arr)) == 4) arr[, , i, ] else arr[, , i]
    max(abs(sums - tot) / pmax(tot, 1e-12))
  }, 0))
}

test_that("coda_causes() checks the options of the two-step model and keeps its defaults", {
  sp <- coda_causes()
  expect_identical(sp$params$weights, "equal")
  expect_null(sp$params$rho)
  expect_false(sp$params$individual)
  expect_identical(sp$label, "CoDa causes (multiple decrement)")
  s2 <- coda_causes(weights = "deaths", rho = 0.1, individual = TRUE)
  expect_match(s2$label, "parts weighted by deaths, years weighted \\(rho = 0.1\\), joint and individual components")
  expect_error(coda_causes(rho = 0), class = "mt_error_argument")
  expect_error(coda_causes(rho = 1), class = "mt_error_argument")
  expect_error(coda_causes(rho = c(0.1, 0.2)), class = "mt_error_argument")
  expect_error(coda_causes(rho = "0.1"), class = "mt_error_argument")
  expect_error(coda_causes(individual = NA), class = "mt_error_argument")
  expect_error(coda_causes(weights = "size"))
  x <- mt_example_causes(years = 1998:2012)
  expect_mt_model(s2, data = x)
  expect_mt_model(coda_causes(rank = 2, weights = "deaths", individual = TRUE, jump_off = "fitted"), data = x)
})

test_that("the generalized SVD minimises the weighted sum of squares and is the SVD with equal weights", {
  C <- with_seed_test(1, matrix(stats::rnorm(12 * 9), 12))
  C <- sweep(C, 2, colMeans(C))
  for (r in 1:2) {
    s0 <- mortensor:::coda_svd(C, r)
    expect_equal(mortensor:::coda_gsvd(C, r), s0, tolerance = 1e-12)
    expect_equal(mortensor:::coda_gsvd(C, r, rep(3, 12), rep(0.2, 9)), s0, tolerance = 1e-12)
  }
  wp <- with_seed_test(2, stats::runif(12, 0.1, 2))
  wt <- 0.1 * 0.9^(8:0)
  g <- mortensor:::coda_gsvd(C, 1, wp, wt)
  b <- g$b[, 1]
  k <- g$k[, 1]
  expect_equal(sum(b^2), 1)
  # first-order conditions: each of b and k is the weighted least-squares fit given the other
  expect_equal(b, drop(C %*% (wt * k)) / sum(wt * k^2))
  expect_equal(k, drop(crossprod(C, wp * b)) / sum(wp * b^2))
  # global minimum: alternating weighted least squares from random starts ends no lower
  W <- outer(wp, wt)
  sse <- function(b, k) sum(W * (C - outer(b, k))^2)
  als <- vapply(1:5, function(s) {
    kk <- with_seed_test(s, stats::rnorm(9))
    for (i in 1:300) {
      bb <- drop(C %*% (wt * kk)) / sum(wt * kk^2)
      kk <- drop(crossprod(C, wp * bb)) / sum(wp * bb^2)
    }
    sse(bb, kk)
  }, 0)
  expect_equal(als, rep(sse(b, k), 5), tolerance = 1e-8)
  expect_equal(g$explained, 1 - sse(b, k) / sum(W * C^2))
})

test_that("equal weights without individual components give the published model; rho near 0 gives equal weights", {
  x <- mt_example_causes(years = 1998:2012)
  f0 <- mt_fit(coda_causes(), x)
  f1 <- mt_fit(coda_causes(weights = "equal", rho = NULL, individual = FALSE), x)
  expect_identical(fitted(f1), fitted(f0))
  expect_identical(mt_forecast(f1, h = 3, n_draws = 10, seed = 1)$rates,
                   mt_forecast(f0, h = 3, n_draws = 10, seed = 1)$rates)
  expect_false(any(grepl("^w_|^explained_", mortensor:::coda_parameters(f0)$parameter)))
  f2 <- mt_fit(coda_causes(rho = 1e-9), x)
  for (g in 1:2) {
    a <- f0$fits$groups[[g]]
    b <- f2$fits$groups[[g]]
    expect_equal(b$b, a$b, tolerance = 1e-6)
    expect_equal(b$k, a$k, tolerance = 1e-6)
    expect_equal(mortensor:::coda_point(b, 5)$Z, mortensor:::coda_point(a, 5)$Z, tolerance = 1e-6)
  }
  expect_equal(fitted(f2), fitted(f0), tolerance = 1e-6)
})

test_that("individual components vanish when the causes share one component", {
  s <- sim_2s(ind = 0)
  base <- coda_causes(residual = FALSE, weights = "deaths", rho = 0.2, order = c(0, 1, 0))
  ind <- coda_causes(residual = FALSE, weights = "deaths", rho = 0.2, order = c(0, 1, 0), individual = TRUE)
  f0 <- mt_fit(base, s$data)
  f1 <- mt_fit(ind, s$data)
  g <- f1$fits$groups[[1]]
  expect_lt(max(abs(g$b[, -1] %*% t(g$k[, -1]))), 1e-10)
  expect_equal(fitted(f1), fitted(f0), tolerance = 1e-10)
  expect_equal(mortensor:::coda_point(g, 5)$Z, mortensor:::coda_point(f0$fits$groups[[1]], 5)$Z,
               tolerance = 1e-10)
})

test_that("the two-step model recovers its components from data simulated from it", {
  # weights of parts and years (2S-CoDA), then equal weights
  for (w in list(list(rho = 0.2, weights = "deaths"), list(rho = NULL, weights = "equal"))) {
    s <- sim_2s(rho = w$rho, weights = w$weights)
    sp <- coda_causes(residual = FALSE, weights = w$weights, rho = w$rho, individual = TRUE, order = c(0, 1, 0),
                      jump_off = "fitted")
    fit <- mt_fit(sp, s$data)
    g <- fit$fits$groups[[1]]
    if (w$weights == "deaths") expect_equal(g$decomp$w_part, s$wp, tolerance = 1e-12)
    if (!is.null(w$rho)) expect_equal(g$decomp$w_time, s$wt)
    expect_equal(g$alpha, s$alpha, tolerance = 1e-10)
    expect_equal(outer(g$b[, 1], g$k[, 1]), outer(s$bJ, s$kJ), tolerance = 1e-10)
    for (i in 1:3) {
      expect_equal(outer(g$b_individual[, i], g$k[, 1 + i]), outer(s$BI[, i], s$kI[, i]), tolerance = 1e-10)
    }
    expect_lt(max(abs(g$resid)), 1e-10)
    ex <- g$explained_causes
    expect_equal(unname(rowSums(ex)), rep(1, 3))
    expect_lt(max(ex[, "residual"]), 1e-12)
    # point forecasts: random walks with drift of the recovered indexes, eq. (7)
    h <- 4
    nt <- ncol(s$Z)
    ext <- function(k) k[nt] + seq_len(h) * (k[nt] - k[1]) / (nt - 1)
    truth <- s$alpha + outer(s$bJ, ext(s$kJ)) + s$BI %*% t(apply(s$kI, 2, ext))
    expect_equal(closure_2s(exp(mortensor:::coda_point(g, h)$Z)), closure_2s(exp(truth)), tolerance = 1e-10)
  }
  # Poisson deaths: the components are recovered up to the noise
  s <- sim_2s(poisson = TRUE, seed = 1)
  s0 <- sim_2s()
  g <- mt_fit(coda_causes(residual = FALSE, weights = "deaths", rho = 0.2, individual = TRUE), s$data)$fits$groups[[1]]
  expect_gt(stats::cor(g$b[, 1], s0$bJ), 0.99)
  expect_gt(stats::cor(g$k[, 1], s0$kJ), 0.99)
  expect_true(all(abs(diag(stats::cor(g$k[, 2:4], s0$kI))) > 0.98))
  expect_lt(max(abs(outer(g$b[, 1], g$k[, 1]) - outer(s0$bJ, s0$kJ))), 0.1 * max(abs(outer(s0$bJ, s0$kJ))))
})

test_that("two-step forecasts of the causes add up to all causes and their life-table deaths to one", {
  x <- mt_example_causes(years = 1998:2012)
  fit <- mt_fit(coda_causes(weights = "deaths", rho = 0.1, individual = TRUE), x)
  st <- fit$strata
  expect_lt(cause_gap(fitted(fit), st), 1e-12)
  fc <- mt_forecast(fit, h = 4, n_draws = 30, seed = 1)
  expect_lt(cause_gap(fc$rates, st), 1e-12)
  pr <- mt_forecast(fit, h = 4, n_draws = 30, seed = 1, predictive = TRUE)
  expect_lt(cause_gap(pr$rates, st), 1e-12)
  # multiple-decrement life table of each draw: deaths by age and cause sum to one
  ages <- mt_grid(x)$age
  i <- st$stratum[st$sex == "female" & st$cause == "all"]
  j <- st$stratum[st$sex == "female" & st$cause != "all"]
  tot <- matrix(fc$rates[, , i, ], length(ages))
  lt <- mortensor:::lt_columns(tot, ages)
  dxc <- vapply(j, function(c) lt$dx * matrix(fc$rates[, , c, ], length(ages)) / tot, tot)
  expect_equal(apply(dxc, 2, sum), rep(1, ncol(tot)), tolerance = 1e-10)
  # parameters: weights, individual components, explained variation
  g <- fit$fits$groups[[1]]
  p <- mt_parameters(fit)
  expect_true(all(c("w_part", "w_time", "b_cvd", "k_cancer", "arima_k_other", "explained_joint",
                    "explained_individual") %in% p$parameter))
  expect_equal(sum(p[stratum == "causes: female" & parameter == "w_part"]$value), 1)
  expect_equal(p[stratum == "causes: female" & parameter == "w_time"]$value, 0.1 * 0.9^(14:0))
  expect_length(p[stratum == "causes: female" & parameter == "b_open"]$value, 1)
  expect_equal(g$open$b[-1], rep(0, 3))
  expect_equal(unname(rowSums(g$explained_causes)), rep(1, 3))
  expect_match(summary(fit)$assumptions, "extensions of Kjaergaard et al. \\(2019")
  # a residual cause gets its own individual component
  fr <- mt_fit(coda_causes(individual = TRUE), mt_subset(x, cause != "other"))
  expect_true(fr$fits$residual)
  expect_identical(fr$fits$groups[[1]]$cause_labels, c("cvd", "cancer", "residual"))
  expect_true("b_residual" %in% mt_parameters(fr)$parameter)
})

test_that("bootstrap residuals of a weighted fit come from their own part and favour recent years", {
  R <- matrix(seq_len(40) + 100 * rep(1:4, 10), 4)
  set.seed(1)
  pooled <- mortensor:::coda_resample(R, list(w_part = NULL, w_time = NULL))
  expect_true(all(pooled %in% R))
  expect_false(all(vapply(1:4, function(p) all(pooled[p, ] %in% R[p, ]), TRUE)))
  by_part <- mortensor:::coda_resample(R, list(w_part = rep(0.25, 4), w_time = NULL))
  expect_true(all(vapply(1:4, function(p) all(by_part[p, ] %in% R[p, ]), TRUE)))
  recent <- mortensor:::coda_resample(R, list(w_part = rep(0.25, 4), w_time = 0.99 * 0.01^(9:0)))
  expect_gt(mean(recent == R[, 10]), 0.9)
})

test_that("coda_causes() fits the two-step model to the ONS causes of death", {
  skip_on_cran()
  causes <- names(ons_data)[-(1:7)]
  x <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile", "cause"),
               age = "age_group", value = list(deaths = c(all = "all_causes", causes), population = "population"))
  x <- mt_subset(mt_as_exposure(x, type = "midyear"), ages = c(50, Inf))
  x <- mt_group_levels(x, "imd_decile", function(d) rep("England", length(d)), name = "area")
  x <- mt_group_levels(x, "cause", ons_cause_groups, other = "rest", total = "all", additive = "rates")
  fit <- mt_fit(coda_causes(weights = "deaths", rho = 0.1, individual = TRUE), x)
  expect_true(all(fit$diagnostics$status == "ok"))
  fc <- mt_forecast(fit, h = 5, n_draws = 20, seed = 1)
  expect_lt(cause_gap(fc$rates, fit$strata), 1e-12)
  ex <- fit$fits$groups[[1]]$explained_causes
  expect_equal(nrow(ex), length(ons_cause_groups) + 1)
  expect_equal(unname(rowSums(ex)), rep(1, nrow(ex)))
})
