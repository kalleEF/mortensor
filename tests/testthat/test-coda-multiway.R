# Multiway CoDa model of Bergeron-Boucher, Simonacci, Oeppen and Gallo (2018):
# coda_multiway(). No code or numbers of the paper can be reproduced here
# (Canadian Human Mortality Database data); the tests check parameter
# recovery on simulated compositions of known structure, the least-squares
# optimum of the Tucker3 fit, the reduction to coda() for one population, the
# shared time index of the forecasts and the draws.

closure_cols_mw <- function(M) sweep(M, 2, colSums(M), "/")
unit <- function(v) v / sqrt(sum(v^2))

# Life-table deaths of `ns` populations with a known rank-1 Tucker3
# structure in clr space, eq. (5): clr(d[x, t, i]) = alpha[x, i] + g b[x]
# k[t] gamma[i] (+ noise); the log rate of the open interval is linear in
# k. Returns the deaths and exposures that give these compositions and the
# true parameters.
sim_multiway <- function(noise = 0, nt = 20, gamma = c(1, 1.3, 0.7, 1.1), g = 8, seed = 1) {
  set.seed(seed)
  ages <- age_groups(c(0, 1, seq(5, 85, 5)))
  na <- length(ages)
  ns <- length(gamma)
  layout <- mortensor:::coda_layout(ages)
  b <- unit(scale(-cos(seq(0, pi, length.out = na)), scale = FALSE)[, 1])
  k <- cumsum(c(0, 0.1 + stats::rnorm(nt - 1, 0, 0.02)))
  k <- unit(k - mean(k))
  gamma <- unit(gamma)
  x <- ages$lower
  alpha <- vapply(seq_len(ns), function(i) {
    m0 <- exp(-9.5 + 0.2 * i + 0.09 * x)
    as.vector(mortensor:::coda_clr(matrix(m0), layout))
  }, numeric(na))
  Z <- array(NA_real_, c(na, nt, ns))
  open_log <- matrix(NA_real_, nt, ns)
  m <- array(NA_real_, c(na, nt, ns))
  for (i in seq_len(ns)) {
    # noise in clr space: centred over the parts of each composition
    e <- matrix(stats::rnorm(na * nt, 0, noise), na)
    Z[, , i] <- alpha[, i] + g * gamma[i] * outer(b, k) + sweep(e, 2, colMeans(e))
    open_log[, i] <- log(0.2) + 0.05 * i - 0.3 * k
    m[, , i] <- mortensor:::coda_layout_rates(Z[, , i], layout, open_log[, i])
  }
  E <- array(1e5, dim(m))
  list(D = m * E, E = E, m = m, ages = ages, periods = 2001:(2000 + nt), layout = layout, Z = Z, b = b, k = k,
       gamma = gamma, g = g, alpha = alpha, open_log = open_log,
       strata = data.table::data.table(stratum = seq_len(ns), label = paste0("pop", seq_len(ns))))
}

fit_mw <- function(s, spec = coda_multiway()) {
  fit_model(spec, s$D, s$E, s$ages, s$periods, strata = s$strata)
}

test_that("coda_multiway() satisfies the model contract", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  expect_mt_model(coda_multiway(), data = x)
  expect_mt_model(coda_multiway(ranks = c(2, 2, 2), jump_off = "observed"), data = x)
  expect_mt_model(coda_multiway(order = c(0, 1, 1)), data = x)
  # single ages from 80: the Kannisto extension, as in coda()
  expect_mt_model(coda_multiway())
  expect_mt_model(coda_multiway(open = "regression"))
  # a hierarchy of strata and separate models per sex, as the paper fits women and men
  expect_mt_model(coda_multiway(), data = mt_aggregate_strata(x, mt_hierarchy(sex = "both", region = "all",
                                                                                quintile = "all")))
  expect_mt_model(stratified("sex", female = coda_multiway(), male = coda_multiway()), data = x)
})

test_that("coda_multiway() checks its arguments and replaces tucker3(transform = \"coda\")", {
  sp <- coda_multiway()
  expect_s3_class(sp, "mt_spec_coda_multiway")
  expect_true(sp$joint)
  expect_identical(sp$label, "Multiway CoDa")
  expect_equal(sp$params$ranks, c(1, 1, 1))
  expect_equal(sp$params$order, c(0, 1, 0))
  expect_identical(sp$params$jump_off, "fitted")
  expect_match(coda_multiway(ranks = c(2, 2, 2), order = c(0, 1, 1), jump_off = "observed")$label,
               "ranks 2/2/2.*ARIMA\\(0,1,1\\) with drift.*observed jump-off")
  expect_error(coda_multiway(ranks = c(1, 1)), class = "mt_error_argument")
  expect_error(coda_multiway(ranks = c(1, 0, 1)), class = "mt_error_argument")
  expect_error(coda_multiway(ranks = c(1, 1.5, 1)), class = "mt_error_argument")
  expect_error(coda_multiway(order = c(0, 3, 0)), class = "mt_error_argument")
  expect_error(coda_multiway(drift = NA), class = "mt_error_argument")
  expect_error(coda_multiway(jump_off = "actual"))
  expect_warning(tucker3(transform = "coda"), "coda_multiway\\(\\)", class = "mt_warning_deprecated")
})

test_that("the multiway fit recovers a known rank-1 structure exactly", {
  s <- sim_multiway()
  f <- fit_mw(s)
  # the compositions of the rates are those of the simulation
  expect_equal(f$alpha, s$alpha, tolerance = 1e-9)
  expect_equal(drop(f$U[[1]]), s$b, tolerance = 1e-9)
  expect_equal(drop(f$U[[2]]), s$k, tolerance = 1e-9)
  expect_equal(drop(f$U[[3]]), s$gamma, tolerance = 1e-9)
  expect_equal(drop(f$G), s$g, tolerance = 1e-9)
  expect_equal(f$fit_pct, 100, tolerance = 1e-9)
  # b and k sum to zero (clr and centring over the years)
  expect_equal(sum(f$U[[1]]), 0, tolerance = 1e-12)
  expect_equal(sum(f$U[[2]]), 0, tolerance = 1e-12)
  # fitted rates and the regression of the open interval on the index
  expect_equal(f$fitted, s$m, tolerance = 1e-9, ignore_attr = TRUE)
  for (i in 1:4) {
    expect_equal(f$open[[i]]$a, log(0.2) + 0.05 * i, tolerance = 1e-9)
    expect_equal(f$open[[i]]$b, -0.3, tolerance = 1e-9)
  }
  # the index follows a random walk with drift: drift = mean change
  expect_equal(f$ts[[1]]$order, c(0, 1, 0))
  expect_equal(f$ts[[1]]$const, mean(diff(s$k)))
})

test_that("the multiway fit is the least-squares Tucker3 solution with noise", {
  s <- sim_multiway(noise = 0.03, seed = 2)
  f <- fit_mw(s)
  X <- s$Z - array(f$alpha[, rep(1:4, each = 20)], dim(s$Z))
  b <- drop(f$U[[1]])
  k <- drop(f$U[[2]])
  gm <- drop(f$U[[3]])
  # the parameters are close to the truth
  expect_gt(abs(sum(b * s$b)), 0.99)
  expect_gt(abs(sum(k * s$k)), 0.99)
  expect_gt(abs(sum(gm * s$gamma)), 0.99)
  # stationarity of the best rank-1 approximation: each factor is the
  # normalised contraction of the centred array with the other two, and the
  # core is the contraction with all three
  contract <- function(n) {
    switch(n, apply(X, 1, function(M) sum(M * outer(k, gm))), apply(X, 2, function(M) sum(M * outer(b, gm))),
           apply(X, 3, function(M) sum(M * outer(b, k))))
  }
  g <- drop(f$G)
  expect_equal(b, contract(1) / g, tolerance = 1e-8)
  expect_equal(k, contract(2) / g, tolerance = 1e-8)
  expect_equal(gm, contract(3) / g, tolerance = 1e-8)
  expect_equal(f$fit_pct, 100 * g^2 / sum(X^2), tolerance = 1e-10)
  # the global optimum: no higher-order power iteration from random starts
  # reaches a larger core
  set.seed(5)
  best <- 0
  for (start in 1:10) {
    u <- unit(stats::rnorm(dim(X)[1]))
    v <- unit(stats::rnorm(dim(X)[2]))
    w <- unit(stats::rnorm(dim(X)[3]))
    for (it in 1:100) {
      u <- unit(apply(X, 1, function(M) sum(M * outer(v, w))))
      v <- unit(apply(X, 2, function(M) sum(M * outer(u, w))))
      w <- unit(apply(X, 3, function(M) sum(M * outer(u, v))))
    }
    best <- max(best, abs(sum(X * outer(outer(u, v), w))))
  }
  expect_lte(best, g + 1e-8)
  # ranks 2/2/2: the core is the projection of the centred array on the factors
  f2 <- fit_mw(s, coda_multiway(ranks = c(2, 2, 2)))
  G <- X
  for (n in 1:3) G <- mortensor:::mode_product(G, t(f2$U[[n]]), n)
  expect_equal(f2$G, G, tolerance = 1e-10)
  for (n in 1:3) expect_equal(crossprod(f2$U[[n]]), diag(2), tolerance = 1e-10)
  expect_gt(f2$fit_pct, f$fit_pct)
  expect_equal(f2$fit_pct, sum(f2$core_pct))
  # signs: increasing time indexes, markers with a positive sum
  expect_true(all(f2$U[[2]][20, ] > f2$U[[2]][1, ]))
  expect_true(all(colSums(f2$U[[3]]) > 0))
})

test_that("with one population the multiway model is the CoDa model of coda()", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  D <- mt_array(x, "deaths")[, , 1, drop = FALSE]
  E <- mt_array(x, "exposure")[, , 1, drop = FALSE]
  ages <- mt_grid(x)$age
  for (jo in c("fitted", "observed")) {
    for (open in c("kannisto", "regression")) {
      cd <- fit_model(coda(order = c(0, 1, 0), jump_off = jo, open = open), D[, , 1], E[, , 1], ages, 1995:2012)
      mw <- fit_model(coda_multiway(jump_off = jo, open = open), D, E, ages, 1995:2012)
      expect_equal(mw$alpha[, 1], cd$alpha, tolerance = 1e-10)
      expect_equal(mw$U[[1]], cd$b, tolerance = 1e-8)
      expect_equal(drop(mw$U[[2]] * drop(mw$G)), drop(cd$k), tolerance = 1e-8)
      expect_equal(drop(mw$U[[3]]), 1)
      expect_equal(mw$fitted[, , 1], cd$fitted, tolerance = 1e-8, ignore_attr = TRUE)
      p_cd <- mortensor:::coda_point(cd, 10)
      p_mw <- mortensor:::coda_multiway_point(mw, 1, 10)
      expect_equal(p_mw$Z, p_cd$Z, tolerance = 1e-8, ignore_attr = TRUE)
      if (open == "regression") expect_equal(p_mw$open, p_cd$open, tolerance = 1e-8)
    }
  }
  # also with the ARIMA(0,1,1) index of coda() and its first-step correction
  cd <- fit_model(coda(order = c(0, 1, 1)), D[, , 1], E[, , 1], ages, 1995:2012)
  mw <- fit_model(coda_multiway(order = c(0, 1, 1), jump_off = "observed"), D, E, ages, 1995:2012)
  expect_equal(mortensor:::coda_multiway_point(mw, 1, 10)$Z, mortensor:::coda_point(cd, 10)$Z, tolerance = 1e-5,
               ignore_attr = TRUE)
})

test_that("multiway forecasts share the time index and draws simulate it with resampled innovations", {
  s <- sim_multiway(noise = 0.02, seed = 3)
  spec <- coda_multiway()
  f <- fit_mw(s, spec)
  h <- 6
  n <- 200
  # point forecasts: every population moves along g b gamma[i] (k[T + j] - k[T])
  pts <- lapply(1:4, function(i) mortensor:::coda_multiway_point(f, i, h))
  gm <- drop(f$U[[3]])
  for (i in 2:4) {
    expect_equal(pts[[i]]$Z - f$alpha[, i], gm[i] / gm[1] * (pts[[1]]$Z - f$alpha[, 1]), tolerance = 1e-10)
  }
  expect_equal(drop(pts[[1]]$k), f$k[20] + (1:h) * mean(diff(drop(f$k))))
  # draws: alpha[, i] + loadings of i times paths of the index, which come
  # from the random walk with drift with resampled innovations
  K <- with_seed_test(7, mortensor:::index_arima_paths(f$ts[[1]], h, n))
  sim <- with_seed_test(7, mortensor:::simulate_model(spec, f, h, n))
  expect_equal(dim(sim), c(length(s$ages), h, 4, n))
  for (i in 1:4) {
    Zi <- mortensor:::coda_clr(matrix(sim[, , i, ], length(s$ages)), f$layout)
    expect_equal(Zi, f$alpha[, i] + f$load[[i]] %*% t(as.vector(K)), tolerance = 1e-8, ignore_attr = TRUE)
    # the compositions sum to one and are those of the clr
    dx <- mortensor:::lt_columns(matrix(sim[, , i, ], length(s$ages)), s$ages)$dx
    expect_equal(colSums(dx), rep(1, h * n), tolerance = 1e-12)
    expect_equal(dx, closure_cols_mw(exp(Zi)), tolerance = 1e-10)
    # the open interval follows the index paths
    expect_equal(log(matrix(sim[length(s$ages), , i, ], h)), f$open[[i]]$a + f$open[[i]]$b * K, tolerance = 1e-8)
  }
  # the median of the draws is the point forecast; intervals widen
  med <- apply(K, 1, stats::median)
  expect_lt(max(abs(med - drop(pts[[1]]$k))), 0.25 * stats::sd(K[h, ]))
  expect_gt(stats::sd(K[h, ]), stats::sd(K[1, ]))
})

test_that("the observed jump-off starts every population at its observed deaths", {
  s <- sim_multiway(noise = 0.02, seed = 4)
  fo <- fit_mw(s, coda_multiway(jump_off = "observed"))
  ff <- fit_mw(s, coda_multiway())
  for (i in 1:4) {
    po <- mortensor:::coda_multiway_point(fo, i, 3)
    pf <- mortensor:::coda_multiway_point(ff, i, 3)
    # the same change from the observed deaths as from the fitted deaths of the last year
    zfit_T <- ff$alpha[, i] + ff$load[[i]] %*% ff$k[20, ]
    expect_equal(po$Z - fo$Z_last[, i], pf$Z - drop(zfit_T), tolerance = 1e-10, ignore_attr = TRUE)
  }
  # without drift the forecast stays at the observed rates of the last year
  fo$ts[[1]]$const <- 0
  for (i in 1:4) {
    po <- mortensor:::coda_multiway_point(fo, i, 2)
    m <- mortensor:::coda_layout_rates(po$Z, fo$layout, po$open)
    expect_equal(m[, 1], s$m[, 20, i], tolerance = 1e-9)
  }
})

test_that("ranks above what the data allow are reduced with a warning, in the parameters and the label", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  # the population rank is capped by the product of the other two
  expect_warning(fit <- mt_fit(coda_multiway(ranks = c(2, 2, 9)), x), "c\\(2, 2, 9\\).*c\\(2, 2, 4\\) are used",
                 class = "mt_warning_argument")
  expect_equal(fit$fits$ranks, c(2, 2, 4))
  p <- mt_parameters(fit)
  expect_equal(p[parameter == "rank"]$value, c(2, 2, 4))
  expect_identical(p[parameter == "rank"]$index, c("age", "time", "population"))
  expect_true("gamma4" %in% p$parameter)
  expect_false("gamma5" %in% p$parameter)
  expect_identical(fit$spec$label, "Multiway CoDa (ranks 2/2/4)")
  expect_identical(mt_forecast(fit, h = 2, n_draws = 5, seed = 1)$model, "Multiway CoDa (ranks 2/2/4)")
  # within the limits, no warning and the label of the specification
  expect_no_warning(f2 <- mt_fit(coda_multiway(ranks = c(2, 2, 4)), x))
  expect_identical(f2$spec$label, "Multiway CoDa (ranks 2/2/4)")
  expect_equal(mt_parameters(f2)[parameter == "rank"]$value, c(2, 2, 4))
  # the years: 6 years allow 5 time indexes
  expect_warning(f5 <- mt_fit(coda_multiway(ranks = c(2, 8, 4)), mt_subset(x, periods = c(2009, 2014))),
                 "c\\(2, 5, 4\\) are used", class = "mt_warning_argument")
  expect_equal(f5$fits$ranks, c(2, 5, 4))
})

test_that("mt_parameters() and summary() describe the multiway fit", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  fit <- mt_fit(coda_multiway(), x)
  p <- mt_parameters(fit)
  expect_true(all(c("b", "k", "gamma", "g", "order_k", "sd_k", "explained", "a", "a_open", "b_open") %in%
                    p$parameter))
  expect_equal(nrow(p[parameter == "gamma"]), nrow(mt_strata(x)))
  expect_equal(p[parameter == "explained"]$value, fit$fits$fit_pct / 100)
  expect_true(all(abs(p[parameter == "a", sum(value), by = stratum]$V1 - 1) < 1e-12))
  p2 <- mt_parameters(mt_fit(coda_multiway(ranks = c(2, 2, 2)), x))
  expect_true(all(c("b1", "b2", "k1", "k2", "gamma1", "gamma2", "arima_k1", "b_open") %in% p2$parameter))
  expect_equal(nrow(p2[parameter == "g"]), 8)
  expect_match(mortensor:::coda_assumptions(coda_multiway()), "Tucker3")
})
