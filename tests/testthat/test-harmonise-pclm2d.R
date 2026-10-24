# The PCLM over age and period of mt_ungroup(dims = 2) and mt_smooth(dims =
# 2) in array form (Currie, Durban and Eilers 2006) gives the results of the
# dense implementation of mortensor 0.2.0 (helper-pclm-dense.R)

test_that("the banded Cholesky factor, solution and inverse equal the dense ones", {
  set.seed(1)
  m <- 4L
  nb <- 6L
  h <- 2L
  X <- matrix(stats::rnorm(60 * m * nb), 60)
  blk <- (seq_len(m * nb) - 1L) %/% m
  band <- abs(outer(blk, blk, "-")) <= h
  D <- crossprod(X) * band + diag(m * nb)
  bd <- mortensor:::band_layout(nb, h, m)
  A <- mortensor:::band_from_dense(D, bd)
  expect_identical(mortensor:::band_dense(A, bd), D)
  U <- mortensor:::band_chol(A, bd)
  R <- chol(D)
  expect_equal(mortensor:::band_dense(U, bd)[upper.tri(D, diag = TRUE)], R[upper.tri(R, diag = TRUE)],
               tolerance = 1e-12)
  b <- stats::rnorm(m * nb)
  expect_equal(mortensor:::band_solve(U, bd, b), as.numeric(solve(D, b)), tolerance = 1e-12)
  expect_equal(mortensor:::band_mv(A, bd, b), as.numeric(D %*% b), tolerance = 1e-12)
  expect_equal(mortensor:::band_quad(A, bd, b), sum(b * (D %*% b)), tolerance = 1e-12)
  S <- mortensor:::band_dense(mortensor:::band_inverse(U, bd), bd)
  expect_equal(S[band], chol2inv(R)[band], tolerance = 1e-12)
  # not positive definite
  D[1, 1] <- -1
  expect_null(mortensor:::band_chol(mortensor:::band_from_dense(D, bd), bd))
})

test_that("the starting values for a tensor basis equal the dense least-squares solution", {
  # six periods and six B-splines over periods: B'B is not singular
  b <- block_inputs(mt_example(years = 2001:2006))
  B <- kronecker(b$bp$B, b$ba$B)
  z <- stats::runif(nrow(B))
  tb <- mortensor:::pclm_tensor(b$ba$B, b$bp$B)
  expect_lt(norm_diff(mortensor:::solve_lsq(tb, z), ref_solve_lsq(B, z)), 1e-8)
  # some cells left out: B'B in array form
  keep <- rep(TRUE, nrow(B))
  keep[seq(3, nrow(B), 7)] <- FALSE
  expect_lt(norm_diff(mortensor:::solve_lsq(mortensor:::basis_rows(tb, keep), z[keep]),
                      ref_solve_lsq(B[keep, ], z[keep])), 1e-8)
  expect_equal(mortensor:::basis_eta(tb, z[seq_len(ncol(B))]), as.numeric(B %*% z[seq_len(ncol(B))]),
               tolerance = 1e-12)
})

test_that("the 2D PCLM equals the dense reference with fixed smoothing parameters", {
  x <- mt_example(years = 2001:2006)
  # coefficients, standard errors and fitted values of one stratum
  b <- block_inputs(x)
  compare_block(b, pclm_ctl(c(100, 10)))
  compare_block(b, pclm_ctl(c(1, 1e4), preserve_totals = TRUE))
  compare_block(b, pclm_ctl(c(100, 10), open_split = "stationary"))
  # mt_ungroup() with the package code and with the reference
  ref <- function(...) with_mocked_bindings(mt_ungroup(...), ungroup_rates = ref_ungroup_rates)
  expect_same_2d(mt_ungroup(x, dims = 2, lambda = c(100, 10)), ref(x, dims = 2, lambda = c(100, 10)), "ungroup")
  # deaths in two-year periods: observations over several periods of the grid
  x2 <- mt_regroup(x, period = c(2001, 2003, 2005), measures = "deaths")
  expect_true(any(mt_obs(x2, "deaths")$period_hi - mt_obs(x2, "deaths")$period_lo == 2))
  compare_block(block_inputs(x2), pclm_ctl(c(100, 10)))
  expect_same_2d(mt_ungroup(x2, dims = 2, lambda = c(100, 10)), ref(x2, dims = 2, lambda = c(100, 10)), "ungroup")
})

test_that("the 2D PCLM equals the dense reference on a grid with an open interval", {
  # deaths and exposures in 5-year groups to 85+, grid in single years to
  # 85+: the exposures are ungrouped first, the open interval of the grid
  # is expanded to omega in the PCLM and its exposure split by the PCLM of
  # the exposures
  x <- mt_example(years = 2001:2005)
  x <- mt_set_grid(mt_regroup(x, age = seq(30, 85, 5)), age = single_ages(30, 85))
  ref <- function(...) with_mocked_bindings(mt_ungroup(...), ungroup_rates = ref_ungroup_rates)
  expect_same_2d(mt_ungroup(x, dims = 2, lambda = c(10, 100)), ref(x, dims = 2, lambda = c(10, 100)), "ungroup")
  u <- mt_ungroup(x, dims = 1)
  x$obs$exposure <- u$obs$exposure
  b <- block_inputs(x, 2L)
  expect_gt(length(b$fa$lower), length(mt_grid(x)$age))
  compare_block(b, pclm_ctl(c(10, 100)))
})

test_that("the 2D P-spline smoothing equals the dense reference with fixed smoothing parameters", {
  x <- mt_example(years = 2001:2006, deaths = "single")
  ref <- function(...) with_mocked_bindings(mt_smooth(...), ungroup_rates = ref_ungroup_rates)
  expect_same_2d(mt_smooth(x, dims = 2, lambda = c(100, 10)), ref(x, dims = 2, lambda = c(100, 10)), "smooth")
  # one cell per observation, the open interval 100+ expanded to omega
  b <- block_inputs(x)
  compare_block(b, pclm_ctl(c(100, 10), keep_observed = FALSE))
})

test_that("the 2D PCLM and smoothing equal the dense reference with chosen smoothing parameters", {
  # The default range goes up to 1e8. Above about 1e5 the system is badly
  # conditioned and the criterion nearly flat, and the fits the search makes
  # there carry rounding errors into the next fits and into the choice: the
  # values chosen by the reference itself change by up to several per cent
  # when only the order of the observations changes. In a range up to 1e4,
  # the search and the chosen values are reproducible to rounding error.
  x <- mt_subset(mt_example(years = 2001:2006), sex == "female")
  ref <- function(f, ...) with_mocked_bindings(f(...), ungroup_rates = ref_ungroup_rates)
  u <- mt_ungroup(x, dims = 2, criterion = "aic", lambda_range = c(1e-2, 1e4))
  expect_false(any(attr(u, "ungroup")$at_bound))
  expect_same_2d(u, ref(mt_ungroup, x, dims = 2, criterion = "aic", lambda_range = c(1e-2, 1e4)), "ungroup")
  skip_on_cran()
  # coefficients and standard errors of the chosen fit
  compare_block(block_inputs(x), pclm_ctl(NULL, criterion = "aic", lambda_range = c(1e-2, 1e4)))
  xs <- mt_subset(mt_example(years = 2001:2006, deaths = "single"), sex == "female")
  expect_same_2d(mt_smooth(xs, dims = 2, criterion = "aic", lambda_range = c(1e-2, 1e4)),
                 ref(mt_smooth, xs, dims = 2, criterion = "aic", lambda_range = c(1e-2, 1e4)), "smooth")
  # the grid open at 85+ of the previous tests
  xo <- mt_subset(mt_example(years = 2001:2005), sex == "female")
  xo <- mt_set_grid(mt_regroup(xo, age = seq(30, 85, 5)), age = single_ages(30, 85))
  expect_same_2d(mt_ungroup(xo, dims = 2, criterion = "aic", lambda_range = c(1e-2, 1e4)),
                 ref(mt_ungroup, xo, dims = 2, criterion = "aic", lambda_range = c(1e-2, 1e4)), "ungroup")
})

test_that("the 2D PCLM and smoothing equal the dense reference on mt_example()", {
  skip_on_cran()
  # the package example (30 years, 71 grid ages, omega 110), one stratum
  x <- mt_subset(mt_example(), sex == "male")
  compare_block(block_inputs(x), pclm_ctl(c(1e4, 1e4)))
  xs <- mt_subset(mt_example(deaths = "single"), sex == "male")
  compare_block(block_inputs(xs), pclm_ctl(c(1e2, 1e4), keep_observed = FALSE))
})
