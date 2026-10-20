# The Schaake shuffle (Clark, Gangopadhyay, Hay, Rajagopalan and Wilby 2004)

h3s <- function() mt_hierarchy(sex = "both", region = "all", quintile = "all")
sorted_margins_s <- function(a) apply(a, 1:3, sort)

# Template [age, horizon, stratum, date] computed directly from log rates y
# [age, year, stratum]: the change from year i to year i + h
template_of <- function(y, H) {
  n <- dim(y)[2] - H
  tpl <- array(NA_real_, c(dim(y)[1], H, dim(y)[3], n))
  for (a in seq_len(dim(y)[1])) for (h in seq_len(H)) for (s in seq_len(dim(y)[3])) for (i in seq_len(n)) {
    tpl[a, h, s, i] <- y[a, i + h, s] - y[a, i, s]
  }
  tpl
}

test_that("the Schaake shuffle reproduces the example of Clark et al. (2004, Section 3b)", {
  # ten ensemble members X and ten historical observations Y of one station
  X <- c(15.3, 11.2, 8.8, 11.9, 7.5, 9.7, 8.3, 12.5, 10.3, 10.1)
  Y <- c(10.7, 9.3, 6.8, 11.3, 12.2, 13.6, 8.9, 9.9, 11.8, 12.9)
  Xss <- c(10.1, 8.8, 7.5, 10.3, 11.9, 15.3, 8.3, 9.7, 11.2, 12.5)
  ec <- mortensor:::ecc_reorder(array(X, c(1, 1, 1, 10)), array(Y, c(1, 1, 1, 10)))
  expect_identical(as.vector(ec$rates), Xss)
  # the same through the template of observed changes: log rates whose
  # one-year changes from the ten start years are Y, a forecast of one year
  y <- array(cumsum(c(0, Y)), c(1, 11, 1))
  sr <- mortensor:::with_seed(1, mortensor:::schaake_reorder(array(X, c(1, 1, 1, 10)), y, 2000:2010))
  expect_identical(sr$info$dates, 10L)
  expect_identical(sr$info$years, 2000:2009)
  expect_setequal(sr$info$draw_year, 2000:2009)
  # the draw with the date of start year 2000 + i - 1 holds the i-th value of Xss
  expect_equal(as.vector(sr$arr)[order(sr$info$draw_year)], Xss)
  expect_identical(sr$info$blocks, 1L)
})

test_that("the Schaake shuffle keeps every margin and takes the ranks of the template in each block", {
  set.seed(1)
  d <- c(3, 4, 5, 50)
  arr <- array(rgamma(prod(d), 2), d, dimnames = list(NULL, NULL, letters[1:5], NULL))
  # log rates of 12 training years: n = 12 - 4 = 8 dates
  y <- aperm(apply(array(rnorm(3 * 12 * 5, -0.02, 0.05), c(3, 12, 5)), c(1, 3), cumsum), c(2, 1, 3))
  sr <- mortensor:::with_seed(2, mortensor:::schaake_reorder(arr, y, 2001:2012))
  expect_identical(sorted_margins_s(sr$arr), sorted_margins_s(arr))
  expect_identical(sr$info$dates, 8L)
  expect_identical(sr$info$blocks, 7L)
  expect_equal(sr$info$margins, prod(d[1:3]))
  expect_equal(as.vector(table(sr$info$block)), c(rep(8, 6), 2))
  # `draw`: the base draw each value comes from
  idx <- cbind(as.vector(slice.index(arr, 1)), as.vector(slice.index(arr, 2)), as.vector(slice.index(arr, 3)),
               as.vector(sr$draw))
  expect_identical(sr$arr, array(arr[idx], d, dimnames = dimnames(arr)))
  tpl <- template_of(y, 4)
  for (b in seq_len(sr$info$blocks)) {
    j <- which(sr$info$block == b)
    dates <- sr$info$draw_year[j] - 2000
    # distinct dates within a block, from the training years whose next 4 years are observed
    expect_false(anyDuplicated(dates) > 0)
    expect_true(all(sr$info$draw_year[j] >= 2001 & sr$info$draw_year[j] + 4 <= 2012))
    for (a in 1:3) for (h in 1:4) for (s in 1:5) {
      expect_identical(rank(sr$arr[a, h, s, j]), rank(tpl[a, h, s, dates]))
      # each block keeps the values of its own draws
      expect_identical(sort(sr$arr[a, h, s, j]), sort(arr[a, h, s, j]))
    }
  }
  # fewer draws than dates: one block of distinct dates chosen at random
  s2 <- mortensor:::with_seed(3, mortensor:::schaake_reorder(arr[, , , 1:5, drop = FALSE], y, 2001:2012))
  expect_identical(s2$info$blocks, 1L)
  expect_length(unique(s2$info$draw_year), 5)
  expect_identical(sorted_margins_s(s2$arr), sorted_margins_s(arr[, , , 1:5, drop = FALSE]))
})

test_that("a comonotone template gives comonotone draws", {
  set.seed(4)
  d <- c(4, 3, 2, 9)
  arr <- array(rexp(prod(d)), d)
  # log rates alpha + k t^2 with k > 0: the change over h years from year i
  # rises with i for every margin and horizon
  tt <- (1:12)^2
  y <- array(NA_real_, c(4, 12, 2))
  for (a in 1:4) for (s in 1:2) y[a, , s] <- rnorm(1) + runif(1, 0.001, 0.01) * tt
  sr <- mortensor:::with_seed(5, mortensor:::schaake_reorder(arr, y, 1:12))
  expect_identical(sr$info$dates, 9L)
  r <- matrix(apply(sr$arr, 1:3, rank), nrow = d[4])
  expect_true(all(r == rank(sr$info$draw_year)))
  expect_identical(sorted_margins_s(sr$arr), sorted_margins_s(arr))
})

test_that("tied template values are ranked at random (Clark et al. 2004, Section 4c)", {
  set.seed(7)
  d <- c(1, 1, 2, 8)
  arr <- array(rnorm(prod(d)), d)
  # one-year changes 1, 0, 1, 0, ... from the eight start years
  y <- array(c(0, 1, 1, 2, 2, 3, 3, 4, 4), c(1, 9, 2))
  s1 <- mortensor:::with_seed(1, mortensor:::schaake_reorder(arr, y, 1:9))
  s2 <- mortensor:::with_seed(1, mortensor:::schaake_reorder(arr, y, 1:9))
  s3 <- mortensor:::with_seed(2, mortensor:::schaake_reorder(arr, y, 1:9))
  expect_identical(s1$arr, s2$arr)
  expect_false(identical(s1$arr, s3$arr))
  expect_equal(s1$info$tied, 1)
  for (sr in list(s1, s3)) {
    expect_identical(sorted_margins_s(sr$arr), sorted_margins_s(arr))
    chg <- (sr$info$draw_year %% 2 == 1) + 0
    for (s in 1:2) {
      v <- sr$arr[1, 1, s, ]
      # a draw with a larger template value never gets a smaller value
      expect_true(all(outer(v, v, ">=")[outer(chg, chg, ">")]))
    }
  }
})

test_that("margins with missing template values keep their order, and short templates are refused", {
  set.seed(6)
  d <- c(2, 2, 3, 10)
  arr <- array(rnorm(prod(d)), d, dimnames = list(NULL, NULL, c("a", "b", "c"), NULL))
  y <- array(rnorm(2 * 10 * 3), c(2, 10, 3))
  y[1, 4, 2] <- NA
  sr <- mortensor:::with_seed(1, mortensor:::schaake_reorder(arr, y, 1:10))
  expect_identical(sr$arr[1, , 2, ], arr[1, , 2, ])
  expect_identical(sr$info$uncoupled, "b")
  expect_warning(mortensor:::warn_schaake(sr$info), class = "mt_warning_reconcile")
  expect_false(identical(sr$arr[2, , 2, ], arr[2, , 2, ]))
  expect_error(mortensor:::schaake_reorder(arr, y[, 1:3, ], 1:3), class = "mt_error_reconcile")
})

test_that("mt_reconcile(dependence = \"schaake\") takes its template from the training years only", {
  x <- mt_example_regions(years = 2000:2019, regions = 2)
  h <- h3s()
  xa <- mt_aggregate_strata(x, h)
  train <- mt_window(xa, to = 2012)
  fc <- mt_forecast(mt_fit(lee_carter(), train), h = 3, n_draws = 30, seed = 1)
  rs <- mt_reconcile(fc, h, method = "bottom_up", dependence = "schaake", seed = 2)
  info <- rs$schaake
  expect_identical(info$dates, 10L)
  expect_equal(info$years, 2000:2009)
  expect_true(all(info$draw_year >= 2000 & info$draw_year + 3 <= 2012))
  for (b in unique(info$block)) expect_false(anyDuplicated(info$draw_year[info$block == b]) > 0)
  expect_lt(coherence_gap(rs, h), 1e-8)
  expect_null(rs$component)
  # the ranks of the bottom strata are those of the observed changes of the training data
  y <- log(pmax(mt_array(train, "deaths"), 0.5) / mt_array(train, "exposure"))
  tpl <- template_of(y, 3)
  b <- mortensor:::hierarchy_structure(fc$strata, h)$bottom
  for (blk in unique(info$block)) {
    j <- which(info$block == blk)
    i <- info$draw_year[j] - 1999
    for (s in b[c(1, 7)]) for (a in c(1, 6)) for (t in 1:3) {
      expect_identical(rank(rs$rates[a, t, s, j]), rank(tpl[a, t, s, i]))
    }
  }
  expect_equal(sorted_margins_s(rs$rates[, , b, ]), sorted_margins_s(fc$rates[, , b, ]), tolerance = 1e-10)
  # later observations do not enter the template
  xm <- xa
  xm$obs$deaths <- data.table::copy(xm$obs$deaths)
  xm$obs$deaths[period_lo > 2012, value := 3 * value]
  fm <- mt_forecast(mt_fit(lee_carter(), mt_window(xm, to = 2012)), h = 3, n_draws = 30, seed = 1)
  expect_identical(mt_reconcile(fm, h, method = "bottom_up", dependence = "schaake", seed = 2)$rates, rs$rates)
  # the Schaake shuffle needs the training data of all strata of the draws:
  # draws of the bottom strata summed to the aggregates have none for these
  fb <- mt_forecast(mt_fit(lee_carter(), mt_window(x, to = 2012)), h = 3, n_draws = 30, seed = 1)
  fa <- mt_aggregate_draws(fb, hierarchy = h)
  expect_error(mt_reconcile(fa, h, dependence = "schaake"), class = "mt_error_reconcile")
})

test_that("reconciled(dependence = \"schaake\") shuffles the base draws before reconciling", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  h <- h3s()
  xa <- mt_aggregate_strata(x, h)
  fit <- mt_fit(reconciled(lee_carter(), h, method = "bottom_up", dependence = "schaake"), xa)
  fc <- mt_forecast(fit, h = 3, n_draws = 40, seed = 3)
  base <- mt_forecast(mt_fit(lee_carter(), xa), h = 3, n_draws = 40, seed = 3)
  b <- mortensor:::hierarchy_structure(fc$strata, h)$bottom
  expect_equal(sorted_margins_s(fc$rates[, , b, ]), sorted_margins_s(base$rates[, , b, ]), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(fc$rates[, , b, ], base$rates[, , b, ])))
  expect_lt(coherence_gap(fc, h), 1e-8)
  expect_identical(fc$schaake$dates, 12L)
  expect_true(all(fc$schaake$draw_year + 3 <= 2014))
  # the same as the shuffle of the base draws with mt_reconcile() after the base forecast
  rs <- mortensor:::with_seed(3, {
    f0 <- mt_forecast(mortensor:::unwrap_fit(fit, "mt_fit_reconciled"), h = 3, n_draws = 40)
    mt_reconcile(f0, h, method = "bottom_up", dependence = "schaake")
  })
  expect_identical(fc$rates, rs$rates)
  # predictive draws: Poisson deaths of the bottom strata given the shuffled rates, coherent
  pr <- mt_forecast(fit, h = 3, n_draws = 40, seed = 3, predictive = TRUE)
  expect_true(all(is.finite(pr$rates)))
  expect_lt(coherence_gap(pr, h), 1e-8)
  # a backtest refits the template at each origin
  bt <- mt_backtest(reconciled(lee_carter(), h, method = "wls", dependence = "schaake"), xa, origins = 2011,
                    h = 3, n_draws = 30)
  expect_true(all(is.finite(bt$crps)))
  expect_match(mortensor:::model_assumptions(fit$spec), "Schaake shuffle")
})

test_that("joint models are not shuffled, and \"empirical\" is deprecated", {
  h <- h3s()
  expect_warning(sp <- reconciled(li_lee(common_by = "sex"), h, dependence = "schaake"),
                 class = "mt_warning_reconcile")
  expect_identical(sp$params$dependence, "none")
  expect_warning(se <- reconciled(lee_carter(), h, dependence = "empirical"), class = "mt_warning_deprecated")
  expect_identical(se$params$dependence, "empirical")
  x <- mt_example_regions(years = 2005:2016, regions = 2)
  fc <- mt_forecast(mt_fit(lee_carter(), mt_aggregate_strata(x, h)), h = 2, n_draws = 20, seed = 1)
  expect_warning(mt_reconcile(fc, h, dependence = "empirical", seed = 1), class = "mt_warning_deprecated")
  expect_no_warning(mt_reconcile(fc, h, dependence = "schaake", seed = 1))
})

test_that("predictive draws after the Schaake shuffle take the noise of the base draw of each value", {
  x <- mt_example_regions(years = 2005:2016, regions = 2)
  h <- h3s()
  xa <- mt_aggregate_strata(x, h)
  # Poisson deaths for values from Lee-Carter draws, none for those of rwd()
  fit <- mt_fit(reconciled(ensemble(lc = lee_carter(), rw = rwd()), h, method = "bottom_up", dependence = "schaake"),
                xa)
  u <- mt_forecast(fit, h = 2, n_draws = 40, seed = 1)
  p <- mt_forecast(fit, h = 2, n_draws = 40, seed = 1, predictive = TRUE)
  expect_null(u$component)
  base <- mortensor:::with_seed(1, mt_forecast(mortensor:::unwrap_fit(fit, "mt_fit_reconciled"), h = 2, n_draws = 40))
  expect_setequal(unique(base$component), 1:2)
  s <- mortensor:::hierarchy_structure(fit$strata, h)$bottom[3]
  E <- fit$exposure_last[, 1, s]
  for (a in c(2, 8)) for (t in 1:2) {
    src <- vapply(u$rates[a, t, s, ], function(v) which.min(abs(base$rates[a, t, s, ] - v)), 1L)
    expect_equal(sort(src), seq_len(40))
    deaths <- p$rates[a, t, s, ] * E[a]
    expect_identical(abs(deaths - round(deaths)) < 1e-6, base$component[src] == 1L)
  }
  expect_lt(coherence_gap(p, h), 1e-8)
})
