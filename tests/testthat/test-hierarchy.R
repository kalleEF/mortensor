h3 <- function() mt_hierarchy(sex = "both", region = "all", quintile = "all")

test_that("hierarchies are defined and printed", {
  h <- mt_hierarchy(sex = "both", region = list(total = "DE", groups = list(East = c("BE", "SN"))))
  expect_s3_class(h, "mt_hierarchy")
  expect_output(print(h), "groups East")
  expect_error(mt_hierarchy("both"), class = "mt_error_hierarchy")
  expect_error(mt_hierarchy(sex = 1), class = "mt_error_hierarchy")
})

test_that("aggregate strata are sums of the bottom strata", {
  x <- mt_example_regions(years = 2010:2013, regions = 3)
  xa <- mt_aggregate_strata(x, h3())
  st <- mt_strata(xa)
  expect_equal(nrow(st), 3 * 4 * 6)
  expect_equal(st$stratum, seq_len(nrow(st)))
  D <- mt_array(xa, "deaths")
  expect_equal(sum(D[, , "both.all.all"]), sum(mt_obs(x, "deaths")$value))
  expect_equal(D[, , "female.all.Q1"], apply(D[, , paste0("female.R", 1:3, ".Q1")], 1:2, sum))
  v <- mt_validate(xa, hierarchy = h3())
  expect_false(any(v$severity == "error"))
  # an inconsistent observed total is flagged
  bad <- xa
  bad$obs$deaths <- data.table::copy(bad$obs$deaths)
  tot <- st$stratum[st$label == "both.all.all"]
  bad$obs$deaths[stratum == tot, value := value * 2]
  vb <- mt_validate(bad, hierarchy = h3())
  expect_true(any(vb$severity == "error" & grepl("median 2.000", vb$detail)))
  # intermediate groups
  hg <- mt_hierarchy(region = list(total = "all", groups = list(north = c("R1", "R2"))))
  xg <- mt_aggregate_strata(x, hg)
  expect_true(all(c("female.north.Q3", "male.all.Q5") %in% mt_strata(xg)$label))
  expect_identical(mt_aggregate_strata(xa, h3())$strata, xa$strata)
})

test_that("reconciliation makes all nodes add up for every method", {
  x <- mt_example_regions(years = 2005:2016, regions = 2)
  xa <- mt_aggregate_strata(x, h3())
  fc <- mt_forecast(mt_fit(lee_carter(), xa), h = 2, n_draws = 40, seed = 1)
  expect_gt(coherence_gap(fc, h3()), 0.01)
  for (m in c("wls", "mint", "structural", "ols", "bottom_up")) {
    r <- mt_reconcile(fc, h3(), method = m, seed = 1)
    expect_lt(coherence_gap(r, h3()), 1e-8)
    expect_true(all(r$rates > 0))
    expect_equal(r$reconciliation$method, m)
  }
  r0 <- mt_reconcile(fc, h3(), dependence = "none")
  expect_lt(coherence_gap(r0, h3()), 1e-8)
  expect_error(mt_reconcile(fc, mt_hierarchy(cause = "all")), class = "mt_error_argument")
})

test_that("reconciliation keeps the draws as they are by default and never reorders joint models", {
  x <- mt_example_regions(years = 2005:2016, regions = 2)
  xa <- mt_aggregate_strata(x, h3())
  fc <- mt_forecast(mt_fit(lee_carter(), xa), h = 2, n_draws = 30, seed = 1)
  b <- mortensor:::hierarchy_structure(fc$strata, h3())$bottom
  # bottom-up without reordering: the bottom draws are those of the base model
  rb <- mt_reconcile(fc, h3(), method = "bottom_up")
  expect_equal(rb$rates[, , b, ], fc$rates[, , b, ], tolerance = 1e-10)
  expect_identical(mt_reconcile(fc, h3())$rates, mt_reconcile(fc, h3(), dependence = "none")$rates)
  re <- mt_reconcile(fc, h3(), method = "bottom_up", dependence = "empirical", seed = 1)
  expect_false(isTRUE(all.equal(re$rates[, , b, ], fc$rates[, , b, ])))
  expect_lt(coherence_gap(re, h3()), 1e-8)
  expect_identical(reconciled(lee_carter(), h3())$params$dependence, "none")
  expect_identical(reconciled(lee_carter(), h3(), dependence = "auto")$params$dependence, "none")
  expect_identical(reconciled(lee_carter(), h3(), dependence = "empirical")$params$dependence, "empirical")
  # joint models carry the dependence between strata: their paths are not reordered
  expect_warning(sp <- reconciled(li_lee(common_by = "sex"), h3(), dependence = "empirical"),
                 class = "mt_warning_reconcile")
  expect_identical(sp$params$dependence, "none")
  fj <- mt_forecast(mt_fit(li_lee(common_by = "sex"), xa), h = 2, n_draws = 30, seed = 1)
  expect_warning(rj <- mt_reconcile(fj, h3(), method = "bottom_up", dependence = "empirical", seed = 1),
                 class = "mt_warning_reconcile")
  expect_equal(rj$rates[, , b, ], fj$rates[, , b, ], tolerance = 1e-10)
  expect_true(mortensor:::spec_is_joint(anchored(li_lee())))
  expect_true(mortensor:::spec_is_joint(ensemble(li_lee(), tucker3())))
  expect_false(mortensor:::spec_is_joint(ensemble(li_lee(), lee_carter())))
  expect_false(mortensor:::spec_is_joint(stratified("sex", female = li_lee(), male = lee_carter())))
})

test_that("the empirical copula keeps each stratum's marginal distribution", {
  x <- mt_example_regions(years = 2005:2016, regions = 2)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 2, n_draws = 50, seed = 1)
  re <- with_seed_test(1, mortensor:::copula_reorder(fc$rates, fc$train_D, fc$train_E))
  expect_equal(sort(re$arr[3, 2, 5, ]), sort(fc$rates[3, 2, 5, ]))
  expect_length(re$years, 50)
})

test_that("reconciled models give coherent point and predictive draws and run in backtests", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  sp <- reconciled(lee_carter(), h3(), method = "wls")
  expect_s3_class(sp, "mt_spec_reconciled")
  fit <- mt_fit(sp, x)
  expect_equal(nrow(fit$strata), 2 * 3 + 3 * 3 * 6 - 2 * 3 * 1 - 0)
  fc <- mt_forecast(fit, h = 2, n_draws = 30, seed = 2)
  expect_lt(coherence_gap(fc, h3()), 1e-8)
  pr <- mt_forecast(fit, h = 2, n_draws = 30, seed = 2, predictive = TRUE)
  expect_lt(coherence_gap(pr, h3()), 1e-8)
  xa <- mt_aggregate_strata(x, h3())
  bt <- mt_backtest(reconciled(li_lee(common_by = "sex"), h3(), method = "bottom_up"), xa, origins = 2011,
                    h = 3, n_draws = 30)
  expect_true(all(c("both", "female") %in% bt$sex))
})

test_that("hierarchical Li-Lee pools bottom strata within each sex", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  xa <- mt_aggregate_strata(x, h3())
  fit <- mt_fit(li_lee(common_by = "sex"), xa)
  expect_equal(length(fit$fits$commons), 3)
  expect_error(mt_fit(li_lee(common_by = "nonsense"), xa), class = "mt_error_model_fit")
  expect_mt_model(li_lee(common_by = "sex"), data = x, n_draws = 10)
})

test_that("draws are aggregated over strata and ages", {
  x <- mt_example(deaths = "single", years = 2000:2012)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 2, n_draws = 10, seed = 1)
  fa <- mt_aggregate_draws(fc, hierarchy = mt_hierarchy(sex = "both"), age = age_groups(c(30, seq(40, 100, 10))))
  expect_equal(dim(fa$rates), c(8, 2, 3, 10))
  E <- fc$exposure_last[, 1, ]
  both <- (fc$rates[1:10, 1, 1, 1] %*% E[1:10, 1] + fc$rates[1:10, 1, 2, 1] %*% E[1:10, 2]) / sum(E[1:10, ])
  expect_equal(fa$rates[1, 1, "both", 1], as.numeric(both))
  expect_error(mt_aggregate_draws(fc, age = age_groups(c(30, 35.5))), class = "mt_error_argument")
})

test_that("rate-additive hierarchies (causes) reconcile on the rate scale, also in predictive draws", {
  x <- mt_example_causes(years = 2005:2019)
  h <- mt_hierarchy(cause = "all", additive = "rates")
  fit <- mt_fit(reconciled(lee_carter(), h, method = "wls"), x)
  base <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 50, seed = 1)
  st <- mt_strata(x)
  for (pred in c(FALSE, TRUE)) {
    fc <- mt_forecast(fit, h = 3, n_draws = 50, seed = 1, predictive = pred)
    for (sx in c("female", "male")) {
      i <- st$stratum[st$sex == sx & st$cause == "all"]
      j <- st$stratum[st$sex == sx & st$cause != "all"]
      parts <- apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum)
      expect_lt(max(abs(parts / fc$rates[, , i, ] - 1)), 1e-8)
      ratio <- stats::median(fc$rates[, , i, ] / base$rates[, , i, ])
      expect_gt(ratio, 0.8)
      expect_lt(ratio, 1.25)
    }
  }
})

test_that("reconciled() warns when the exposures of the totals do not match the hierarchy", {
  x <- mt_example_causes(years = 2005:2019)
  expect_warning(mt_fit(reconciled(lee_carter(), mt_hierarchy(cause = "all")), x),
                 class = "mt_warning_hierarchy_exposure")
  expect_no_warning(mt_fit(reconciled(lee_carter(), mt_hierarchy(cause = "all", additive = "rates")), x))
  xr <- mt_aggregate_strata(mt_example_regions(years = 2005:2019, regions = 2),
                            mt_hierarchy(sex = "both", region = "all"))
  expect_no_warning(mt_fit(reconciled(lee_carter(), mt_hierarchy(sex = "both", region = "all")), xr))
})

# A small hierarchy (total, three groups, eight bottom series) for the
# quadratic programs of the reconciliation
small_hierarchy <- function() {
  groups <- list(1:3, 4:5, 6:8)
  sagg <- rbind(rep(1, 8), t(vapply(groups, function(g) replace(numeric(8), g, 1), numeric(8))))
  list(S = rbind(sagg, diag(8)), Ut = cbind(diag(4), -sagg), b = 5:12, agg = 1:4)
}

# Bottom values of the quadratic program of Wickramasuriya, Turlach and
# Hyndman (2020), Eq. (3), by quadprog; nodes in `keep` as equality constraints
qp_bottom <- function(y, Lam, S, keep = integer()) {
  li <- solve(Lam)
  nb <- ncol(S)
  A <- diag(nb)
  b0 <- numeric(nb)
  if (length(keep)) {
    A <- cbind(t(S[keep, , drop = FALSE]), A)
    b0 <- c(y[keep], b0)
  }
  quadprog::solve.QP(t(S) %*% li %*% S, as.vector(t(S) %*% li %*% y), A, b0, meq = length(keep))$solution
}

test_that("non-negative reconciliation solves the quadratic program of Wickramasuriya et al. (2020)", {
  skip_if_not_installed("quadprog")
  hs <- small_hierarchy()
  set.seed(3)
  n_neg <- 0
  for (i in 1:60) {
    y <- as.vector(hs$S %*% (stats::rgamma(8, 2) * 10)) * exp(stats::rnorm(12, 0, 0.6))
    y[hs$agg] <- y[hs$agg] * stats::runif(4, 0.2, 1)
    full <- i %% 2 == 0
    Lam <- if (full) {
      crossprod(matrix(stats::rnorm(144), 12)) / 12 + diag(stats::runif(12, 0.5, 2))
    } else {
      diag(stats::runif(12, 0.5, 3))
    }
    W <- if (full) Lam else diag(Lam)
    keep <- list(integer(), 1L, c(2L, 9L))[[i %% 3 + 1]]
    Wk <- mortensor:::keep_weights(W, keep)
    x0 <- mortensor:::project_rows(matrix(y, 1), Wk, hs$Ut)
    if (all(x0[, hs$b] >= 0)) {
      # without negative values the projection is the optimum of the program
      expect_equal(as.vector(x0[, hs$b]), qp_bottom(y, Lam, hs$S, keep), tolerance = 1e-7)
      next
    }
    n_neg <- n_neg + 1
    nn <- mortensor:::nnls_reconcile(matrix(y, 1), Wk, hs$Ut, hs$b, x0)
    expect_true(nn$converged)
    ref <- pmax(qp_bottom(y, Lam, hs$S, keep), 0)
    expect_lt(max(abs(nn$X[, hs$b] - ref)) / max(abs(y)), 1e-8)
    expect_true(all(nn$X[, hs$b] >= -1e-10 * max(abs(y))))
    expect_equal(as.vector(nn$X), as.vector(hs$S %*% nn$X[1, hs$b]), tolerance = 1e-10)
    expect_equal(nn$held[1, ], nn$X[1, hs$b] == 0)
    # fixed nodes are kept as equality constraints
    if (length(keep)) expect_lt(max(abs(nn$X[1, keep] / y[keep] - 1)), 1e-7)
  }
  expect_gt(n_neg, 20)
})

test_that("non-negative reconciliation of forecast draws is the optimum of every draw", {
  skip_if_not_installed("quadprog")
  x <- mt_example_regions(years = 2005:2019, regions = 2)
  h <- h3()
  xa <- mt_aggregate_strata(x, h)
  fc <- mt_forecast(mt_fit(lee_carter(), xa), h = 2, n_draws = 30, seed = 1)
  tot <- which(fc$strata$label == "both.all.all")
  low <- fc
  # aggregates over regions or quintiles far below the sums of their parts
  st <- fc$strata
  agg <- which(st$region == "all" | st$quintile == "all")
  low$rates[, , agg, ] <- 0.3 * low$rates[, , agg, ]
  hs <- mortensor:::hierarchy_structure(fc$strata, h)
  E <- mortensor:::expand_exposure(fc$exposure_last, dim(fc$rates), dimnames(fc$rates), hs)
  for (fx in list(NULL, "total")) {
    rc <- suppressWarnings(mt_reconcile(low, h, method = "wls", fixed = fx))
    expect_gt(rc$reconciliation$share_set_to_zero, 0)
    # rates are floored at 1e-12, which matters only for aggregates of zeros
    expect_lt(coherence_gap(rc, h), 1e-4)
    # the age with most zeros: the bottom deaths of each draw solve the
    # quadratic program with the variances of the base draws as weights and
    # the fixed total as an equality constraint
    a <- which.max(apply(rc$rates[, 1, hs$bottom, ] <= 1e-12, 1, sum))
    Y <- t(low$rates[a, 1, , ] * E[a, 1, ])
    Lam <- diag(pmax(apply(Y, 2, stats::var), 1e-20))
    Xb <- t(rc$rates[a, 1, hs$bottom, ] * E[a, 1, hs$bottom])
    keep <- if (is.null(fx)) integer() else tot
    ref <- t(apply(Y, 1, function(y) pmax(qp_bottom(y, Lam, hs$S, keep), 0)))
    expect_lt(max(abs(Xb - ref)) / max(abs(Y)), 1e-7)
  }
})
