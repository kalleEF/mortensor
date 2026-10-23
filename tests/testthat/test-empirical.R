# The Schaake shuffle of stratum paths (dependence = "empirical"): the
# shuffle of Clark et al. (2004) with the stratum as margin, ranked by the
# log age-standardised rate

h3e <- function() mt_hierarchy(sex = "both", region = "all", quintile = "all")

# Paths of each stratum as a set: the columns [age x year] of each stratum,
# sorted, so that two forecasts compare equal when every stratum holds the
# same paths in any order
path_sets <- function(a) {
  lapply(seq_len(dim(a)[3]), function(s) {
    A <- matrix(a[, , s, ], ncol = dim(a)[4])
    A[, do.call(order, as.data.frame(t(A))), drop = FALSE]
  })
}

# Log age-standardised rates [period, stratum] of rates R [age, period,
# stratum] with standard weights w
log_asr <- function(R, w) {
  apply(R, c(2, 3), function(r) log(sum(w * r)))
}

# The same [stratum, draw] of draws [age, period, stratum, draw] in period t
log_asr_draws <- function(A, w, t) {
  apply(A[, t, , , drop = FALSE], c(3, 4), function(r) log(sum(w * r)))
}

test_that("the shuffle of stratum paths reproduces the example of Clark et al. (2004, Section 3b)", {
  # one age and one forecast year: each path is one value, ranked by its log rate
  X <- c(15.3, 11.2, 8.8, 11.9, 7.5, 9.7, 8.3, 12.5, 10.3, 10.1)
  Y <- c(10.7, 9.3, 6.8, 11.3, 12.2, 13.6, 8.9, 9.9, 11.8, 12.9)
  Xss <- c(10.1, 8.8, 7.5, 10.3, 11.9, 15.3, 8.3, 9.7, 11.2, 12.5)
  # log rates whose annual changes in the years 2001-2010 are Y
  y <- array(cumsum(c(0, Y)), c(1, 11, 1))
  E <- array(1, c(1, 11, 1))
  er <- mortensor:::with_seed(1, mortensor:::empirical_reorder(array(X, c(1, 1, 1, 10)), y, E, 2000:2010))
  expect_identical(er$info$dates, 10L)
  expect_identical(er$info$years, 2001:2010)
  expect_setequal(er$info$draw_year, 2001:2010)
  expect_equal(unname(er$info$template[1, ]), Y)
  # the draw with the date 2000 + i holds the i-th value of Xss
  expect_identical(as.vector(er$arr)[order(er$info$draw_year)], Xss)
  expect_identical(er$info$blocks, 1L)
})

test_that("every stratum keeps its paths, and their summaries take the ranks of the template", {
  set.seed(11)
  na <- 4
  H <- 3
  ns <- 5
  M <- 23
  # rates [age, year, stratum, draw] with the age pattern of a stratum and paths that differ by draw
  arr <- array(NA_real_, c(na, H, ns, M), dimnames = list(as.character(1:na), NULL, letters[1:ns], NULL))
  for (s in 1:ns) for (j in 1:M) {
    arr[, , s, j] <- exp(-6 + 0.8 * (1:na) + rnorm(1, 0, 0.1) + outer(rnorm(na, 0, 0.05), cumsum(rnorm(H, 0, 0.1)),
                                                                       "+"))
  }
  # 9 training years with exposures that differ by age and stratum
  E <- array(runif(na * 9 * ns, 1e3, 1e4), c(na, 9, ns), dimnames = list(as.character(1:na), NULL, letters[1:ns]))
  D <- array(rpois(na * 9 * ns, 0.01 * E), dim(E), dimnames(E))
  D[1, 2, 3] <- 0
  y <- mortensor:::schaake_log_rates(NULL, E, D)
  er <- mortensor:::with_seed(4, mortensor:::empirical_reorder(arr, y, E, 2001:2009))
  # each stratum holds its own paths, whole
  expect_identical(path_sets(er$arr), path_sets(arr))
  for (s in 1:ns) for (j in 1:M) expect_identical(er$arr[, , s, j], arr[, , s, er$perm[s, j]])
  # the template: annual changes of the log age-standardised rate, standard the mean exposure by age
  w <- apply(E, 1, mean)
  w <- w / sum(w)
  expect_equal(unname(er$info$standard), unname(w))
  z <- log_asr(pmax(D, 0.5) / E, w)
  tpl <- t(z[-1, ] - z[-9, ])
  expect_equal(unname(er$info$template), unname(tpl), tolerance = 1e-12)
  expect_identical(er$info$years, 2002:2009)
  expect_identical(er$info$dates, 8L)
  expect_identical(er$info$blocks, 3L)
  expect_equal(as.vector(table(er$info$block)), c(8, 8, 7))
  # the draws are ranked by their mean annual change from the last observed year
  chg <- (log_asr_draws(er$arr, w, H) - z[9, ]) / H
  for (b in 1:3) {
    j <- which(er$info$block == b)
    i <- er$info$draw_year[j] - 2001
    expect_false(anyDuplicated(i) > 0)
    for (s in 1:ns) expect_identical(rank(chg[s, j]), rank(tpl[s, i]))
  }
})

test_that("a comonotone template gives comonotone strata, and short templates leave the draws as they are", {
  set.seed(5)
  arr <- array(rexp(2 * 2 * 3 * 7), c(2, 2, 3, 7))
  # log rates alpha + k t^2, k > 0: the annual change rises with the year in every stratum
  y <- array(NA_real_, c(2, 8, 3))
  for (a in 1:2) for (s in 1:3) y[a, , s] <- -5 + a + runif(1, 0.001, 0.01) * (1:8)^2
  E <- array(1000, c(2, 8, 3))
  er <- mortensor:::with_seed(2, mortensor:::empirical_reorder(arr, y, E, 1:8))
  w <- c(0.5, 0.5)
  chg <- log_asr_draws(er$arr, w, 2)
  r <- t(apply(chg, 1, rank))
  expect_true(all(r[1, ] == r[2, ] & r[2, ] == r[3, ]))
  expect_identical(r[1, ], rank(er$info$draw_year))
  # three training years give two dates; two give a warning and keep the order
  e3 <- mortensor:::with_seed(1, mortensor:::empirical_reorder(arr, y[, 6:8, , drop = FALSE], E[, 6:8, ], 6:8))
  expect_identical(e3$info$dates, 2L)
  expect_identical(e3$info$blocks, 4L)
  expect_identical(path_sets(e3$arr), path_sets(arr))
  expect_warning(e2 <- mortensor:::empirical_reorder(arr, y[, 7:8, , drop = FALSE], E[, 7:8, ], 7:8),
                 class = "mt_warning_reconcile")
  expect_identical(e2$arr, arr)
  expect_null(e2$perm)
  # a stratum with a missing template value keeps its order, with a warning
  ym <- y
  ym[2, 4, 2] <- NA
  em <- mortensor:::with_seed(1, mortensor:::empirical_reorder(arr, ym, E, 1:8))
  expect_identical(em$arr[, , 2, ], arr[, , 2, ])
  expect_identical(em$info$uncoupled, "2")
  expect_warning(mortensor:::warn_schaake(em$info), class = "mt_warning_reconcile")
})

test_that("tied template values are ranked at random", {
  set.seed(8)
  arr <- array(rexp(1 * 1 * 2 * 8), c(1, 1, 2, 8))
  # annual changes 1, 0, 1, 0, ... in the eight years after the first
  y <- array(c(0, 1, 1, 2, 2, 3, 3, 4, 4), c(1, 9, 2))
  E <- array(1, c(1, 9, 2))
  e1 <- mortensor:::with_seed(1, mortensor:::empirical_reorder(arr, y, E, 1:9))
  e2 <- mortensor:::with_seed(2, mortensor:::empirical_reorder(arr, y, E, 1:9))
  expect_equal(e1$info$tied, 1)
  expect_false(identical(e1$arr, e2$arr))
  for (er in list(e1, e2)) {
    # template values are differences of log(exp(y)), whose rounding depends on the maths library: 1 up to rounding
    up <- (abs(er$info$template[1, as.character(er$info$draw_year)] - 1) < 1e-12) + 0
    for (s in 1:2) {
      v <- er$arr[1, 1, s, ]
      expect_true(all(outer(v, v, ">=")[outer(up, up, ">")]))
    }
  }
})

test_that("mt_reconcile(dependence = \"empirical\") takes its template from the training years only", {
  x <- mt_example_regions(years = 2000:2019, regions = 2)
  h <- h3e()
  xa <- mt_aggregate_strata(x, h)
  train <- mt_window(xa, to = 2012)
  fc <- mt_forecast(mt_fit(lee_carter(), train), h = 3, n_draws = 30, seed = 1)
  re <- mt_reconcile(fc, h, method = "bottom_up", dependence = "empirical", seed = 2)
  info <- re$empirical
  expect_identical(info$dates, 12L)
  expect_equal(info$years, 2001:2012)
  expect_true(all(info$draw_year >= 2001 & info$draw_year <= 2012))
  for (b in unique(info$block)) expect_false(anyDuplicated(info$draw_year[info$block == b]) > 0)
  expect_lt(coherence_gap(re, h), 1e-8)
  b <- mortensor:::hierarchy_structure(fc$strata, h)$bottom
  expect_equal(path_sets(re$rates[, , b, ]), path_sets(fc$rates[, , b, ]), tolerance = 1e-12)
  # the template from the observed deaths and exposures of the training years
  E <- mt_array(train, "exposure")
  w <- apply(E, 1, mean)
  w <- w / sum(w)
  z <- log_asr(pmax(mt_array(train, "deaths"), 0.5) / E, w)
  expect_equal(unname(info$template), unname(t(z[-1, ] - z[-13, ])), tolerance = 1e-10)
  # later observations do not enter the template
  xm <- xa
  xm$obs$deaths <- data.table::copy(xm$obs$deaths)
  xm$obs$deaths[period_lo > 2012, value := 3 * value]
  fm <- mt_forecast(mt_fit(lee_carter(), mt_window(xm, to = 2012)), h = 3, n_draws = 30, seed = 1)
  expect_identical(mt_reconcile(fm, h, method = "bottom_up", dependence = "empirical", seed = 2)$rates, re$rates)
  # the summary of the draws ranks as the template within each block
  zh <- log_asr_draws(re$rates, w, 3)
  for (blk in unique(info$block)) {
    j <- which(info$block == blk)
    for (s in b[c(1, 9, 20)]) {
      expect_identical(rank(zh[s, j]), rank(unname(info$template[s, as.character(info$draw_year[j])])))
    }
  }
  # the draws need the training data of all their strata
  fb <- mt_forecast(mt_fit(lee_carter(), mt_window(x, to = 2012)), h = 3, n_draws = 30, seed = 1)
  expect_error(mt_reconcile(mt_aggregate_draws(fb, hierarchy = h), h, dependence = "empirical"),
               class = "mt_error_reconcile")
})

test_that("bottom-up with \"empirical\" keeps every stratum's draws and scores of dependence = \"none\"", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  h <- h3e()
  xa <- mt_aggregate_strata(x, h)
  none <- mt_fit(reconciled(lee_carter(), h, method = "bottom_up"), xa)
  emp <- mt_fit(reconciled(lee_carter(), h, method = "bottom_up", dependence = "empirical"), xa)
  b <- emp$reconcile$structure$bottom
  for (pred in c(FALSE, TRUE)) {
    f0 <- mt_forecast(none, h = 4, n_draws = 60, seed = 2, predictive = pred)
    f1 <- mt_forecast(emp, h = 4, n_draws = 60, seed = 2, predictive = pred)
    # whole paths of each bottom stratum, also with the observation noise of predictive draws
    expect_identical(path_sets(f1$rates[, , b, ]), path_sets(f0$rates[, , b, ]))
    expect_false(isTRUE(all.equal(f1$rates[, , b, ], f0$rates[, , b, ])))
    expect_lt(coherence_gap(f1, h), 1e-8)
  }
  expect_identical(f1$empirical$dates, 14L)
  # backtests: the scores of every bottom stratum are those of "none", other scores of the totals
  bt0 <- mt_backtest(reconciled(lee_carter(), h, method = "bottom_up"), xa, origins = 2010, h = 3, n_draws = 40)
  bt1 <- mt_backtest(reconciled(lee_carter(), h, method = "bottom_up", dependence = "empirical"), xa, origins = 2010,
                     h = 3, n_draws = 40)
  lb <- xa$strata$label[b]
  key <- c("label", "age_lo", "period")
  s0 <- bt0[label %in% lb][order(label, age_lo, period)]
  s1 <- bt1[label %in% lb][order(label, age_lo, period)]
  expect_identical(s1[, c(key, "cov50", "cov95"), with = FALSE], s0[, c(key, "cov50", "cov95"), with = FALSE])
  # the CRPS sums the draws of a cell in their order, which the shuffle changes: equal up to rounding
  expect_equal(s1$crps, s0$crps, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(bt1[label == "both.all.all"]$crps, bt0[label == "both.all.all"]$crps)))
  expect_match(mortensor:::model_assumptions(emp$spec), "stratum")
})

test_that("\"empirical\" works with single ages, strata with zero deaths and projection methods", {
  x <- mt_subset(mt_example_regions(years = 2003:2014, regions = 2, ages = "single"), ages = c(80, 90))
  h <- mt_hierarchy(region = "all", quintile = "all")
  xa <- mt_aggregate_strata(x, h)
  fit <- mt_fit(reconciled(lee_carter(), h, method = "wls", dependence = "empirical"), xa)
  fc <- mt_forecast(fit, h = 2, n_draws = 30, seed = 1, predictive = TRUE)
  expect_true(all(is.finite(fc$rates)))
  expect_lt(coherence_gap(fc, h), 1e-8)
  expect_identical(fc$empirical$margins, nrow(xa$strata))
  expect_length(fc$empirical$standard, length(fit$age))
  # with a projection, the paths are reordered and then reconciled, as by
  # mt_reconcile() after the base forecast
  u <- mt_forecast(fit, h = 2, n_draws = 30, seed = 1)
  re <- mortensor:::with_seed(1, {
    base <- mt_forecast(mortensor:::unwrap_fit(fit, "mt_fit_reconciled"), h = 2, n_draws = 30)
    mt_reconcile(base, h, method = "wls", dependence = "empirical")
  })
  expect_lt(coherence_gap(u, h), 1e-8)
  expect_identical(u$rates, re$rates)
  expect_identical(u$empirical, re$empirical)
  # a small bottom stratum with about one death per age and year: deaths of
  # zero count as 0.5 in the template
  set.seed(3)
  xs <- x
  s0 <- xs$strata$stratum[xs$strata$label == "female.R1.Q5"]
  for (m in c("deaths", "exposure")) {
    xs$obs[[m]] <- data.table::copy(xs$obs[[m]])
    xs$obs[[m]][stratum == s0, value := value / 600]
  }
  xs$obs$deaths[stratum == s0, value := as.numeric(stats::rpois(.N, value))]
  expect_gt(sum(xs$obs$deaths$value == 0), 20)
  fs <- mt_fit(reconciled(rwd(), h, method = "bottom_up", dependence = "empirical"), mt_aggregate_strata(xs, h))
  f0 <- mt_fit(reconciled(rwd(), h, method = "bottom_up"), mt_aggregate_strata(xs, h))
  # random walks of so few deaths reach rates above 1 (a warning)
  ps <- suppressWarnings(mt_forecast(fs, h = 2, n_draws = 30, seed = 1), classes = "mt_warning_implausible")
  p0 <- suppressWarnings(mt_forecast(f0, h = 2, n_draws = 30, seed = 1), classes = "mt_warning_implausible")
  expect_true(all(is.finite(ps$empirical$template)))
  expect_identical(ps$empirical$margins, 20L)
  b <- fs$reconcile$structure$bottom
  expect_identical(path_sets(ps$rates[, , b, ]), path_sets(p0$rates[, , b, ]))
  expect_false(identical(ps$rates[, , "female.R1.Q5", ], p0$rates[, , "female.R1.Q5", ]))
  expect_lt(coherence_gap(ps, h), 1e-8)
})
