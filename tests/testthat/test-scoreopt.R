# Score-optimal reconciliation (reconciled(method = "energy_score"))

h3 <- function() mt_hierarchy(sex = "both", region = "all", quintile = "all")

# Toy hierarchy: total = A + B, A = a1 + a2, B = b1 + b2
toy_S <- function() rbind(c(1, 1, 1, 1), c(1, 1, 0, 0), c(0, 0, 1, 1), diag(4))

# A quickly trained energy-score fit on 2 regions (54 nodes) at ages 75+;
# fits are cached within this file
es_cache <- new.env()
es_small <- function(control = list(max_iter = 30, n_draws = 30), seed = 1) {
  key <- paste(deparse(list(control, seed)), collapse = "")
  if (is.null(es_cache[[key]])) {
    x <- mt_subset(mt_example_regions(years = 2000:2014, regions = 2), ages = c(75, Inf))
    xa <- mt_aggregate_strata(x, h3())
    sp <- reconciled(lee_carter(), h3(), method = "energy_score", years = 5, n_origins = 2, horizon = 2,
                     control = control, seed = seed)
    es_cache[[key]] <- list(xa = xa, fit = mt_fit(sp, xa))
  }
  es_cache[[key]]
}

test_that("the energy score and its gradient agree with ProbReco::total_score()", {
  skip_if_not_installed("ProbReco")
  S <- toy_S()
  n <- nrow(S)
  m <- ncol(S)
  nt <- 4
  q <- 25
  set.seed(1)
  Y <- vapply(seq_len(nt), function(t) as.vector(S %*% (c(5, 3, 2, 4) + stats::rnorm(4))), numeric(n))
  mu <- as.vector(S %*% c(5, 3, 2, 4))
  X1 <- lapply(seq_len(nt), function(t) matrix(stats::rnorm(n * q, mu), n))
  X2 <- lapply(seq_len(nt), function(t) matrix(stats::rnorm(n * q, mu), n))
  Gvec <- stats::runif(m + m * n)
  # each function returns the first set of draws, then the second
  mk <- function(t) {
    calls <- new.env()
    calls$k <- 0
    function() {
      calls$k <- calls$k + 1
      if (calls$k %% 2 == 1) X1[[t]] else X2[[t]]
    }
  }
  ref <- ProbReco::total_score(lapply(seq_len(nt), function(t) Y[, t]), lapply(seq_len(nt), mk), S, Gvec)
  z <- list(x1 = do.call(cbind, X1), x2 = do.call(cbind, X2), per = rep(seq_len(nt), each = q))
  mine <- mortensor:::scoreopt_score(Y, z, S, Gvec)
  expect_equal(mine$value, ref$value, tolerance = 1e-10)
  expect_equal(mine$grad, as.vector(ref$grad), tolerance = 1e-9)
})

test_that("the optimisation reproduces ProbReco::scoreopt() and its d and G reconcile as S (d + G y)", {
  skip_if_not_installed("ProbReco")
  S <- toy_S()
  n <- nrow(S)
  m <- ncol(S)
  nt <- 6
  q <- 20
  set.seed(2)
  Y <- vapply(seq_len(nt), function(t) as.vector(S %*% (c(5, 3, 2, 4) + stats::rnorm(4))), numeric(n))
  # biased base forecasts with independent noise, so that d and G move away from OLS
  mu <- as.vector(S %*% c(5, 3, 2, 4)) * c(1.2, 0.9, 1.1, 1, 1, 1, 1)
  prob <- lapply(seq_len(nt), function(t) function() matrix(stats::rnorm(n * q, mu, 0.7), n))
  data <- lapply(seq_len(nt), function(t) Y[, t])
  ctl <- list(max_iter = 60, eta = 0.01, beta1 = 0.9, beta2 = 0.999, tol = 1e-4, epsilon = 1e-8)
  set.seed(3)
  ref <- suppressWarnings(ProbReco::scoreopt(data, prob, S, control = list(maxIter = ctl$max_iter, eta = ctl$eta)))
  # ProbReco checks its inputs by calling every function three times before the first step
  set.seed(3)
  for (k in 1:3) for (f in prob) f()
  draw <- function() {
    x1 <- lapply(prob, function(f) f())
    x2 <- lapply(prob, function(f) f())
    list(x1 = do.call(cbind, x1), x2 = do.call(cbind, x2), per = rep(seq_len(nt), each = q))
  }
  G0 <- solve(crossprod(S), t(S))
  mine <- mortensor:::scoreopt_sgd(Y, draw, S, c(rep(0, m), as.vector(G0)), ctl)
  expect_equal(mine$d, ref$d, tolerance = 1e-6)
  expect_equal(mine$G, ref$G, tolerance = 1e-6)
  expect_equal(mine$value, ref$val, tolerance = 1e-6)
  expect_gt(max(abs(ref$G - G0)), 0.05)
  # mortensor's forecast step with ProbReco's d and G: S (d + G y) draw by draw
  yhat <- aperm(array(stats::rnorm(n * 2 * 3 * 15, mu, 0.7), c(n, 2, 3, 15)), c(2, 3, 1, 4))
  hs <- list(S = S, bottom = 4:7, is_bottom = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE))
  so <- list(d = matrix(ref$d, m, 2), G = array(ref$G, c(m, n, 2)))
  rec <- mortensor:::scoreopt_apply(yhat, hs, "rates", so, array(1, c(2, 3, n)), dim(yhat))
  for (a in 1:2) for (t in 1:3) {
    expect_equal(rec$rates[a, t, , ], pmax(S %*% pmax(ref$d + ref$G %*% yhat[a, t, , ], 0), 1e-12), tolerance = 1e-12)
  }
})

test_that("score-optimal reconciled draws are coherent in every draw", {
  r <- es_small()
  fit <- r$fit
  so <- fit$reconcile$scoreopt
  expect_equal(dim(so$G), c(20, 54, 4))
  expect_equal(nrow(so$table), 4)
  expect_true(all(so$table$forecasts >= 2))
  expect_true(all(is.finite(so$table$energy_opt)))
  fc <- mt_forecast(fit, h = 3, n_draws = 40, seed = 1)
  expect_lt(coherence_gap(fc, h3()), 1e-8)
  expect_true(all(fc$rates > 0))
  expect_identical(fc$reconciliation$method, "energy_score")
  pr <- mt_forecast(fit, h = 3, n_draws = 40, seed = 1, predictive = TRUE)
  expect_lt(coherence_gap(pr, h3()), 1e-8)
  # deaths: aggregates equal the sums of their bottom strata in every draw
  hs <- mortensor:::hierarchy_structure(fc$strata, h3())
  D <- mt_deaths(fc, noise = FALSE)
  for (j in c(1, 17, 40)) {
    full <- matrix(D[, 3, hs$bottom, j], dim(D)[1]) %*% t(hs$S)
    expect_lt(max(abs(full - matrix(D[, 3, , j], dim(D)[1])) / matrix(D[, 3, , j], dim(D)[1])), 1e-8)
  }
})

test_that("without optimisation steps the method is OLS reconciliation", {
  r <- es_small(control = list(max_iter = 0, n_draws = 20))
  expect_true(all(r$fit$reconcile$scoreopt$table$iterations == 0))
  expect_equal(r$fit$reconcile$scoreopt$table$energy_opt, r$fit$reconcile$scoreopt$table$energy_ols)
  ols <- mt_fit(reconciled(lee_carter(), h3(), method = "ols"), r$xa)
  fe <- mt_forecast(r$fit, h = 2, n_draws = 30, seed = 4)
  fo <- mt_forecast(ols, h = 2, n_draws = 30, seed = 4)
  expect_equal(fo$reconciliation$share_set_to_zero, 0)
  expect_equal(fe$rates, fo$rates, tolerance = 1e-8)
  # a zero learning rate leaves OLS as well
  z <- es_small(control = list(max_iter = 5, eta = 0, n_draws = 20))
  expect_equal(z$fit$reconcile$scoreopt$G, r$fit$reconcile$scoreopt$G, tolerance = 1e-12)
})

test_that("training is reproducible with a seed and moves d and G away from OLS", {
  a <- es_small()$fit$reconcile$scoreopt
  r <- es_small()
  b <- mt_fit(reconciled(lee_carter(), h3(), method = "energy_score", years = 5, n_origins = 2, horizon = 2,
                         control = list(max_iter = 30, n_draws = 30), seed = 1), r$xa)$reconcile$scoreopt
  c2 <- es_small(seed = 2)$fit$reconcile$scoreopt
  expect_identical(a$G, b$G)
  expect_identical(a$d, b$d)
  expect_false(isTRUE(all.equal(a$G, c2$G)))
  ols <- es_small(control = list(max_iter = 0, n_draws = 20))$fit$reconcile$scoreopt
  expect_gt(max(abs(a$G - ols$G)), 1e-3)
  expect_true(all(a$table$iterations > 0))
})

test_that("rate-additive hierarchies (causes) are reconciled on the rate scale", {
  x <- mt_subset(mt_example_causes(years = 2005:2019), ages = c(80, 90))
  h <- mt_hierarchy(cause = "all", additive = "rates")
  sp <- reconciled(lee_carter(), h, method = "energy_score", years = 5, n_origins = 2, horizon = 2,
                   control = list(max_iter = 20, n_draws = 20))
  fit <- mt_fit(sp, x)
  fc <- mt_forecast(fit, h = 2, n_draws = 30, seed = 1)
  st <- mt_strata(x)
  for (sx in c("female", "male")) {
    i <- st$stratum[st$sex == sx & st$cause == "all"]
    j <- st$stratum[st$sex == sx & st$cause != "all"]
    parts <- apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum)
    expect_lt(max(abs(parts / fc$rates[, , i, ] - 1)), 1e-8)
  }
})

test_that("score-optimal reconciliation refuses what it cannot do", {
  xa <- es_small()$xa
  fc <- mt_forecast(mt_fit(lee_carter(), xa), h = 2, n_draws = 10, seed = 1)
  expect_error(mt_reconcile(fc, h3(), method = "energy_score"), class = "mt_error_argument")
  expect_error(reconciled(lee_carter(), h3(), method = "energy_score", fixed = "total"), class = "mt_error_argument")
  expect_error(reconciled(lee_carter(), h3(), method = "energy_score", control = list(maxIter = 10)),
               class = "mt_error_argument")
  expect_error(reconciled(lee_carter(), h3(), method = "energy_score", horizon = 0), class = "mt_error_argument")
  sp <- reconciled(lee_carter(), h3(), method = "energy_score", years = 6)
  expect_equal(sp$min_periods, lee_carter()$min_periods + 6)
  # 6 training years cannot hold 6 years of training forecasts after 5 years of base fit
  expect_error(mt_fit(sp, mt_window(xa, from = 2005, to = 2010)), class = "mt_error_model_data")
})
