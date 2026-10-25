# Weights of ensemble(): stacking of predictive distributions with the CRPS
# (Yao et al. 2018), the bias weights of Kontis et al. (2017), and the error
# of the removed stacking of medians

# CRPS of a weighted empirical distribution (atoms x, probabilities p) by its
# definition, the integral of (F(z) - 1{z >= y})^2 (Gneiting and Raftery
# 2007, eq. 20), exact for step functions
crps_by_integral <- function(x, p, y) {
  z <- sort(unique(c(x, y)))
  cdf <- vapply(z, function(t) sum(p[x <= t]), 0)
  step <- as.numeric(z >= y)
  sum(((cdf - step)^2)[-length(z)] * diff(z))
}

# Mean CRPS over cells of the mixture with weights w of the draws X[[k]]
# (cells x draws), each draw of component k with probability w[k] / n_k
mixture_crps_direct <- function(X, y, w) {
  mean(vapply(seq_along(y), function(i) {
    x <- unlist(lapply(X, function(m) m[i, ]))
    p <- unlist(lapply(seq_along(X), function(k) rep(w[k] / ncol(X[[k]]), ncol(X[[k]]))))
    crps_by_integral(x, p, y[i])
  }, 0))
}

energy_form <- function(tm, w) sum(tm$a * w) - 0.5 * sum(w * (tm$B %*% w))

test_that("the kernel form of the mixture CRPS is the CRPS of the weighted pooled draws", {
  set.seed(1)
  X <- list(matrix(rnorm(3 * 7), 3), matrix(rnorm(3 * 5, 1, 2), 3), matrix(rexp(3 * 4), 3))
  y <- c(-0.5, 0.3, 2)
  tm <- mortensor:::crps_mixture_terms(X, y)
  for (w in list(c(0.2, 0.5, 0.3), c(1, 0, 0), c(0, 0.4, 0.6))) {
    expect_equal(energy_form(tm, w), mixture_crps_direct(X, y, w), tolerance = 1e-12)
  }
  # the expectations are means over all draws and pairs of draws
  B <- outer(1:3, 1:3, Vectorize(function(k, l) {
    mean(vapply(1:3, function(i) mean(abs(outer(X[[k]][i, ], X[[l]][i, ], "-"))), 0))
  }))
  expect_equal(tm$B, B, tolerance = 1e-12)
  # rational weights: the pooled sample with replicated draws, scored like the backtest
  Z <- list(X[[1]][, 1:5], X[[2]])
  tz <- mortensor:::crps_mixture_terms(Z, y)
  pooled <- cbind(Z[[1]], Z[[2]], Z[[2]])
  expect_equal(energy_form(tz, c(1, 2) / 3), mean(mortensor:::score_samples(pooled, y)$crps), tolerance = 1e-12)
  # a vertex gives the CRPS of that component
  expect_equal(energy_form(tm, c(0, 1, 0)), mean(mortensor:::score_samples(X[[2]], y)$crps), tolerance = 1e-12)
})

test_that("stacking weights minimise the mean CRPS of the mixture (2- and 3-simplex)", {
  set.seed(2)
  n_cells <- 40
  y <- rnorm(n_cells)
  X2 <- list(matrix(rnorm(n_cells * 30, 0.4, 0.6), n_cells), matrix(rnorm(n_cells * 30, -0.5, 1.5), n_cells))
  tm2 <- mortensor:::crps_mixture_terms(X2, y)
  w2 <- mortensor:::simplex_qp(-tm2$B, tm2$a)
  expect_equal(sum(w2), 1)
  # fine grid on the 2-simplex with the CRPS of the pooled draws by its definition
  grid <- seq(0, 1, by = 0.001)
  f_grid <- vapply(grid, function(v) mixture_crps_direct(X2, y, c(v, 1 - v)), 0)
  expect_equal(w2[1], grid[which.min(f_grid)], tolerance = 1.5e-3)
  expect_lte(mixture_crps_direct(X2, y, w2), min(f_grid) + 1e-12)
  # 3 components: optim over the softmax of the weights and a grid of the simplex
  X3 <- c(X2, list(matrix(rnorm(n_cells * 25, 0, 1), n_cells)))
  tm3 <- mortensor:::crps_mixture_terms(X3, y)
  w3 <- mortensor:::simplex_qp(-tm3$B, tm3$a)
  sm <- function(z) exp(z - max(z)) / sum(exp(z - max(z)))
  op <- stats::optim(c(0, 0, 0), function(z) mixture_crps_direct(X3, y, sm(z)), method = "BFGS",
                     control = list(reltol = 1e-14))
  expect_equal(w3, sm(op$par), tolerance = 1e-3)
  g <- expand.grid(a = seq(0, 1, 0.02), b = seq(0, 1, 0.02))
  g <- g[g$a + g$b <= 1 + 1e-9, ]
  f3 <- apply(g, 1, function(v) energy_form(tm3, c(v[1], v[2], 1 - v[1] - v[2])))
  expect_lte(energy_form(tm3, w3), min(f3) + 1e-12)
  best <- unlist(g[which.min(f3), ])
  expect_equal(w3[1:2], unname(best), tolerance = 0.03)
  # optimality conditions: equal gradients on the support, larger outside it
  grad <- as.vector(tm3$a - tm3$B %*% w3)
  on <- w3 > 1e-10
  expect_lt(diff(range(grad[on])), 1e-10)
  if (any(!on)) expect_true(all(grad[!on] >= max(grad[on]) - 1e-10))
})

test_that("the quadratic program on the simplex is solved exactly", {
  # every support of the weights by its optimality conditions
  by_supports <- function(Q, q) {
    K <- length(q)
    best <- Inf
    for (m in seq_len(2^K - 1)) {
      S <- which(bitwAnd(m, 2^(seq_len(K) - 1)) > 0)
      A <- rbind(cbind(Q[S, S, drop = FALSE], 1), c(rep(1, length(S)), 0))
      sol <- tryCatch(solve(A, c(-q[S], 1)), error = function(e) NULL)
      if (is.null(sol) || any(sol[seq_along(S)] < -1e-12)) next
      w <- numeric(K)
      w[S] <- pmax(sol[seq_along(S)], 0)
      f <- 0.5 * sum(w * (Q %*% w)) + sum(q * w)
      if (f < best) {
        best <- f
        out <- w
      }
    }
    out
  }
  set.seed(3)
  for (r in 1:40) {
    K <- sample(2:5, 1)
    X <- lapply(seq_len(K), function(k) matrix(rnorm(15 * 20, rnorm(1, 0, 0.5), runif(1, 0.3, 2)), 15))
    tm <- mortensor:::crps_mixture_terms(X, rnorm(15))
    expect_equal(mortensor:::simplex_qp(-tm$B, tm$a), by_supports(-tm$B, tm$a), tolerance = 1e-8)
  }
  skip_if_not_installed("quadprog")
  # strictly convex problems against quadprog
  for (r in 1:20) {
    K <- sample(2:6, 1)
    L <- matrix(rnorm(K * K), K)
    Q <- crossprod(L) + diag(0.1, K)
    q <- rnorm(K)
    ref <- quadprog::solve.QP(Q, -q, cbind(1, diag(K)), c(1, rep(0, K)), meq = 1)$solution
    expect_equal(mortensor:::simplex_qp(Q, q), pmax(ref, 0), tolerance = 1e-8)
  }
})

test_that("a component from the true distribution gets nearly all the weight", {
  set.seed(4)
  n_cells <- 400
  y <- rnorm(n_cells)
  truth <- matrix(rnorm(n_cells * 200), n_cells)
  biased <- matrix(rnorm(n_cells * 200, 1.5), n_cells)
  tm <- mortensor:::crps_mixture_terms(list(biased, truth), y)
  w <- mortensor:::simplex_qp(-tm$B, tm$a)
  expect_gt(w[2], 0.95)
})

test_that("identical components give a valid minimum without error", {
  set.seed(5)
  y <- rnorm(20)
  A <- matrix(rnorm(20 * 50), 20)
  C <- matrix(rnorm(20 * 50, 2), 20)
  # the minimum is that of the problem with one copy of A
  crps_a <- mean(mortensor:::score_samples(A, y)$crps)
  tm1 <- mortensor:::crps_mixture_terms(list(A, C), y)
  min_ac <- energy_form(tm1, mortensor:::simplex_qp(-tm1$B, tm1$a))
  cases <- list(list(X = list(A, A), min = crps_a), list(X = list(A, A, A), min = crps_a),
                list(X = list(A, A, C), min = min_ac), list(X = list(C, A, A), min = min_ac))
  for (cs in cases) {
    tm <- mortensor:::crps_mixture_terms(cs$X, y)
    w <- expect_silent(mortensor:::simplex_qp(-tm$B, tm$a))
    expect_equal(sum(w), 1)
    expect_true(all(w >= 0))
    expect_equal(energy_form(tm, w), cs$min, tolerance = 1e-12)
  }
})

test_that("stacking weights come from the draws of the inner backtest", {
  x <- mt_example(deaths = "single", years = 1990:2010)
  ens <- ensemble(lc = lee_carter(), rw = rwd(), weights = "crps_stack", inner_origins = 2, inner_h = 3)
  fit <- mt_fit(ens, x)
  expect_equal(sum(fit$weights), 1)
  bts <- mortensor:::ensemble_inner_backtests(ens$params, x, NULL)
  common <- mortensor:::ensemble_inner_common(bts)
  X <- lapply(1:2, function(i) do.call(rbind, common[[paste0("draws", i)]]))
  # the kept draws are those the backtest scored
  expect_equal(mortensor:::score_samples(X[[1]], common$observed)$crps, common$crps1, tolerance = 1e-12)
  tm <- mortensor:::crps_mixture_terms(X, common$observed)
  expect_equal(unname(fit$weights), mortensor:::simplex_qp(-tm$B, tm$a), tolerance = 1e-12)
  expect_equal(energy_form(tm, c(0, 1)), mean(common$crps2), tolerance = 1e-12)
  # the scores of the backtest are those of mt_backtest()
  bt <- mt_backtest(rwd(), x, origins = 2005:2006, h = 3, n_draws = 100)
  kd <- mortensor:::backtest_impl(rwd(), x, origins = 2005:2006, h = 3, n_draws = 100, keep_draws = TRUE)
  expect_equal(kd$crps, bt$crps)
  # stacking through the inner weights on toy tables with draws
  cells <- function(dr) {
    data.table::data.table(origin = 2010, stratum = 1L, age_lo = 1:3, age_hi = 2:4, period = 2011, target = "rate",
                           observed = c(0, 0.5, 1), crps = 1, median = 0, draws = dr)
  }
  set.seed(6)
  d1 <- lapply(1:3, function(i) rnorm(40))
  d2 <- lapply(1:3, function(i) rnorm(40, 1))
  w <- mortensor:::ensemble_inner_weights(list(weights = "crps_stack"), list(cells(d1), NULL, cells(d2)))
  tm <- mortensor:::crps_mixture_terms(list(do.call(rbind, d1), do.call(rbind, d2)), c(0, 0.5, 1))
  expect_equal(w[c(1, 3)], mortensor:::simplex_qp(-tm$B, tm$a))
  expect_equal(w[2], 0)
})

test_that("bias weights follow Kontis et al. (2017)", {
  # observed and projected life expectancy (at the mean forecast rates) in
  # the withheld years of an inner backtest
  le <- function(origin, stratum, observed, plugin) {
    data.table::data.table(origin = origin, stratum = stratum, age_lo = NA_real_, age_hi = NA_real_,
                           period = origin + seq_along(observed), target = "life_expectancy",
                           observed = observed, crps = 0.1, median = plugin, plugin = plugin)
  }
  obs <- c(80, 80.5, 81, 81.2)
  a <- le(2000, 1L, obs, c(79.5, 80.2, 80.4, 80.5))
  b <- le(2000, 1L, obs, c(80.3, 81, 81.4, 81.8))
  cc <- le(2000, 1L, obs, c(78, 78.2, 78.1, 78))
  w <- mortensor:::ensemble_inner_weights(list(weights = "bias"), list(a, b, cc))
  # hand computation: mean bias 0.525, -0.45 and 2.6 years
  expect_equal(w, exp(-c(0.525, 0.45, 2.6)))
  expect_equal(w / sum(w), exp(-abs(c(mean(obs - a$plugin), mean(obs - b$plugin), mean(obs - cc$plugin)))) /
                 sum(exp(-c(0.525, 0.45, 2.6))))
  # positive and negative errors cancel within the withheld years, as in the paper
  zig <- le(2000, 1L, obs, obs + c(1, -1, 1, -1))
  expect_equal(mortensor:::ensemble_inner_weights(list(weights = "bias"), list(zig, a))[1], 1)
  # several strata (package-specific): the mean over strata of |bias|
  a2 <- rbind(a, le(2000, 2L, obs, obs - 1))
  b2 <- rbind(b, le(2000, 2L, obs, obs + 0.2))
  w2 <- mortensor:::ensemble_inner_weights(list(weights = "bias"), list(a2, b2, NULL))
  expect_equal(w2, c(exp(-mean(c(0.525, 1))), exp(-mean(c(0.45, 0.2))), 0))
  # rows of rates are not used
  rate_rows <- data.table::data.table(origin = 2000, stratum = 1L, age_lo = 1, age_hi = 2, period = 2001,
                                      target = "rate", observed = 0, crps = 0.1, median = 0, plugin = NA_real_)
  expect_equal(mortensor:::ensemble_inner_weights(list(weights = "bias"), list(rbind(a, rate_rows), b)),
               exp(-c(0.525, 0.45)))
})

test_that("the projected life expectancy is that of the mean forecast rates", {
  o <- data.table::data.table(period_lo = 2001, age_lo = c(60, 70, 80), age_hi = c(70, 80, 90),
                              value = c(20, 50, 90))
  E <- c(1000, 800, 400)
  D <- matrix(c(18, 45, 95, 22, 52, 85, 21, 49, 99, 19, 55, 88), 3)
  qs <- mortensor:::quantity_scores(D, o, E, "life_expectancy", mortensor:::wis_probs(), c(60, 90), plugin = TRUE)
  ax <- mt_axis(c(60, 70, 80), c(70, 80, 90), kind = "age")
  expect_equal(qs$plugin, as.vector(life_expectancy(rowMeans(D) / E, ax)))
  expect_equal(qs$observed, as.vector(life_expectancy(o$value / E, ax)))
  expect_null(mortensor:::quantity_scores(D, o, E, "life_expectancy", mortensor:::wis_probs(), c(60, 90))$plugin)
})

test_that("bias weights use one origin 13 years back by default", {
  ens <- ensemble(lc = lee_carter(), rw = rwd(), weights = "bias")
  expect_equal(ens$params[c("inner_origins", "inner_h")], list(inner_origins = 1, inner_h = 13))
  expect_equal(ensemble(lee_carter(), weights = "crps")$params[c("inner_origins", "inner_h")],
               list(inner_origins = 3, inner_h = 5))
  x <- mt_example(deaths = "single", years = 1990:2012)
  ens <- ensemble(lc = lee_carter(), rw = rwd(), weights = "bias", inner_h = 6)
  fit <- mt_fit(ens, x)
  bts <- mortensor:::ensemble_inner_backtests(ens$params, x, NULL)
  b <- bts[[2]][bts[[2]]$target == "life_expectancy"]
  expect_equal(sort(unique(b$period)), 2007:2012)
  expect_equal(unique(b$origin), 2006)
  common <- mortensor:::ensemble_inner_common(bts, "life_expectancy")
  bias <- vapply(1:2, function(i) {
    mean(abs(tapply(common$observed - common[[paste0("plugin", i)]], common$stratum, mean)))
  }, 0)
  expect_equal(unname(fit$weights), exp(-bias) / sum(exp(-bias)), tolerance = 1e-12)
})

test_that("top_k works with every weight option", {
  x <- mt_example(deaths = "single", years = 1990:2012)
  for (wt in c("equal", "crps", "crps_stack", "bias")) {
    ens <- ensemble(lc = lee_carter(), rw = rwd(), tu = tucker3(), weights = wt, top_k = 2, inner_origins = 2,
                    inner_h = 3)
    fit <- mt_fit(ens, x)
    expect_lte(sum(fit$weights > 0), 2)
    expect_equal(sum(fit$weights), 1)
  }
  # stacking among the kept components
  ens <- ensemble(lc = lee_carter(), rw = rwd(), tu = tucker3(), weights = "crps_stack", top_k = 2,
                  inner_origins = 2, inner_h = 3)
  bts <- mortensor:::ensemble_inner_backtests(ens$params, x, NULL)
  common <- mortensor:::ensemble_inner_common(bts)
  keep <- rank(vapply(1:3, function(i) mean(common[[paste0("crps", i)]]), 0)) <= 2
  common_k <- mortensor:::ensemble_inner_common(replace(bts, !keep, list(NULL)))
  X <- lapply(which(keep), function(i) do.call(rbind, common_k[[paste0("draws", i)]]))
  tm <- mortensor:::crps_mixture_terms(X, common_k$observed)
  w <- numeric(3)
  w[keep] <- mortensor:::simplex_qp(-tm$B, tm$a)
  expect_equal(unname(mt_fit(ens, x)$weights), w, tolerance = 1e-12)
})

test_that("weights = \"stack\" is removed and failed inner backtests give equal weights", {
  expect_error(ensemble(lee_carter(), rwd(), weights = "stack"), class = "mt_error_argument")
  expect_warning(w <- mortensor:::ensemble_inner_weights(list(weights = "crps_stack", inner_origins = 3, inner_h = 5),
                                                         list(NULL, NULL)),
                 class = "mt_warning_ensemble_weights")
  expect_equal(w, c(1, 1))
  expect_error(ensemble(lee_carter(), weights = "mean"), class = "mt_error_argument")
  expect_match(model_assumptions(ensemble(lee_carter(), rwd(), weights = "crps_stack")), "Yao et al. 2018")
  expect_match(model_assumptions(ensemble(lee_carter(), rwd(), weights = "bias")), "Kontis et al. 2017")
})

test_that("ensemble() warns when inner backtests score fewer origins than requested", {
  # 11 years: the first inner origin (2003) has 4 training years, fewer than Lee-Carter needs
  x <- mt_example(deaths = "single", years = 2000:2010)
  ens <- ensemble(lc = lee_carter(), rw = rwd(), weights = "crps", inner_origins = 3, inner_h = 5)
  expect_warning(mt_fit(ens, x), class = "mt_warning_ensemble_weights")
  x2 <- mt_example(deaths = "single", years = 1990:2012)
  expect_no_warning(mt_fit(ens, x2), class = "mt_warning_ensemble_weights")
})
