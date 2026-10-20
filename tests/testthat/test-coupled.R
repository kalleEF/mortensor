# Ensemble copula coupling (Schefzik, Thorarinsdottir and Gneiting 2013)

coupled_data <- function() mt_subset(mt_example_regions(years = 2000:2019, regions = 2), ages = c(60, Inf))

sorted_margins <- function(a) apply(a, 1:3, sort)
margin_ranks <- function(a) apply(a, 1:3, rank)

test_that("ECC keeps the values of each margin and takes the ranks of the template", {
  set.seed(1)
  d <- c(3, 4, 5, 40)
  base <- array(rgamma(prod(d), 2), d)
  tpl <- array(rnorm(prod(d)), d)
  ec <- mortensor:::ecc_reorder(base, tpl)
  expect_identical(sorted_margins(ec$rates), sorted_margins(base))
  expect_identical(margin_ranks(ec$rates), margin_ranks(tpl))
  # Spearman correlations between margins are those of the template
  sp <- function(a) stats::cor(t(matrix(a, ncol = d[4])), method = "spearman")
  expect_equal(sp(ec$rates), sp(tpl))
  # `draw`: the base draw each value comes from
  expect_identical(ec$rates, array(base[cbind(as.vector(slice.index(base, 1)), as.vector(slice.index(base, 2)),
                                              as.vector(slice.index(base, 3)), as.vector(ec$draw))], d))
  expect_equal(ec$info$margins, prod(d[1:3]))
  expect_equal(ec$info$tied, 0)
  # only the target strata are reordered
  ec2 <- mortensor:::ecc_reorder(base, tpl, targets = 2)
  expect_identical(ec2$rates[, , -2, ], base[, , -2, ])
  expect_identical(margin_ranks(ec2$rates[, , 2, , drop = FALSE]), margin_ranks(tpl[, , 2, , drop = FALSE]))
})

test_that("ECC with a comonotone template gives comonotone draws", {
  set.seed(2)
  d <- c(4, 3, 2, 30)
  u <- rnorm(d[4])
  tpl <- array(rep(u, each = prod(d[1:3])), d) * array(runif(prod(d[1:3]), 1, 2), d) + array(rnorm(prod(d[1:3])), d)
  base <- array(rexp(prod(d)), d)
  out <- mortensor:::ecc_reorder(base, tpl)$rates
  r <- matrix(margin_ranks(out), ncol = prod(d[1:3]))
  expect_true(all(r == rank(u)))
  expect_identical(sorted_margins(out), sorted_margins(base))
})

test_that("ECC ranks tied template values at random", {
  set.seed(3)
  d <- c(2, 2, 1, 50)
  tpl <- array(sample(1:5, prod(d), replace = TRUE), d)
  base <- array(rnorm(prod(d)), d)
  e1 <- mortensor:::with_seed(1, mortensor:::ecc_reorder(base, tpl))
  e2 <- mortensor:::with_seed(1, mortensor:::ecc_reorder(base, tpl))
  e3 <- mortensor:::with_seed(2, mortensor:::ecc_reorder(base, tpl))
  expect_identical(e1$rates, e2$rates)
  expect_false(identical(e1$rates, e3$rates))
  expect_equal(e1$info$tied, 1)
  expect_identical(sorted_margins(e1$rates), sorted_margins(base))
  # a higher template value always gets a value at least as high
  for (i in 1:2) for (j in 1:2) {
    t1 <- tpl[i, j, 1, ]
    v <- e1$rates[i, j, 1, ]
    hi <- outer(t1, t1, ">")
    expect_true(all(outer(v, v, ">=")[hi]))
  }
})

test_that("ECC leaves margins with missing values in their order", {
  set.seed(4)
  d <- c(2, 2, 3, 10)
  base <- array(rnorm(prod(d)), d, dimnames = list(NULL, NULL, c("a", "b", "c"), NULL))
  tpl <- array(rnorm(prod(d)), d)
  tpl[, , 2, ] <- NA
  ec <- mortensor:::ecc_reorder(base, tpl)
  expect_identical(ec$rates[, , 2, ], base[, , 2, ])
  expect_identical(ec$info$uncoupled, "b")
  expect_identical(unname(margin_ranks(ec$rates[, , 3, , drop = FALSE])), margin_ranks(tpl[, , 3, , drop = FALSE]))
})

test_that("coupled() forecasts have the base model's margins and the template's ranks", {
  x <- coupled_data()
  fit <- mt_fit(coupled(lee_carter(), template = cae()), x)
  fc <- mt_forecast(fit, h = 4, n_draws = 60, seed = 11)
  # the draws of the separate models with the same random numbers
  fb_fit <- mt_fit(lee_carter(), x)
  ft_fit <- mt_fit(cae(), x)
  sep <- mortensor:::with_seed(11, list(
    base = mt_forecast(fb_fit, h = 4, n_draws = 60)$rates,
    tpl = mt_forecast(ft_fit, h = 4, n_draws = 60)$rates))
  expect_identical(sorted_margins(fc$rates), sorted_margins(sep$base))
  expect_identical(margin_ranks(fc$rates), margin_ranks(sep$tpl))
  expect_identical(fc$coupling$method, "ECC-R")
  expect_identical(fc$coupling$template, "CAE (order 1)")
  # predictive draws: Poisson deaths given the coupled rates
  fp <- mt_forecast(fit, h = 4, n_draws = 60, seed = 11, predictive = TRUE)
  expect_true(all(is.finite(fp$rates) & fp$rates >= 0))
  E <- mortensor:::expand_exposure(fit$exposure_last, dim(fp$rates), dimnames(fp$rates))
  deaths <- fp$rates * as.vector(E)
  expect_equal(deaths, round(deaths), tolerance = 1e-8)
  # in a cell with many deaths, the observed rate of a draw follows its coupled rate
  big <- which(apply(E * fc$rates[, , , 1], 3, max) > 300)[1]
  expect_gt(stats::cor(fp$rates[7, 4, big, ], fc$rates[7, 4, big, ]), 0.8)
  # the template gives the strata the dependence of the joint model
  st <- c(1, 6)
  cor_of <- function(a) stats::cor(log(a[3, 4, st[1], ]), log(a[3, 4, st[2], ]), method = "spearman")
  expect_equal(cor_of(fc$rates), cor_of(sep$tpl))
})

test_that("reconciled(dependence = \"ecc\") couples before reconciling", {
  x <- coupled_data()
  h <- mt_hierarchy(quintile = "all")
  a <- mt_forecast(mt_fit(reconciled(lee_carter(), h, method = "bottom_up", dependence = "ecc", template = cae()), x),
                   h = 3, n_draws = 40, seed = 5)
  b <- mt_forecast(mt_fit(reconciled(coupled(lee_carter(), cae()), h, method = "bottom_up"), x), h = 3,
                   n_draws = 40, seed = 5)
  expect_identical(a$rates, b$rates)
  expect_lt(coherence_gap(a, h), 1e-8)
  ap <- mt_forecast(mt_fit(reconciled(lee_carter(), h, method = "bottom_up", dependence = "ecc", template = cae()), x),
                    h = 3, n_draws = 40, seed = 5, predictive = TRUE)
  expect_true(all(is.finite(ap$rates)))
  expect_error(reconciled(lee_carter(), h, dependence = "ecc"), class = "mt_error_argument")
  expect_error(reconciled(lee_carter(), h, template = cae()), class = "mt_error_argument")
  # coupled draws carry dependence: the Schaake shuffle is not applied
  expect_true(mortensor:::spec_is_joint(coupled(lee_carter(), cae())))
  expect_warning(reconciled(coupled(lee_carter(), cae()), h, dependence = "schaake"),
                 class = "mt_warning_reconcile")
})

test_that("coupled() around a reconciled model keeps the draws coherent", {
  x <- coupled_data()
  h <- mt_hierarchy(quintile = "all")
  fit <- mt_fit(coupled(reconciled(lee_carter(), h, method = "wls"), template = reconciled(cae(), h, "bottom_up")), x)
  fc <- mt_forecast(fit, h = 2, n_draws = 30, seed = 1)
  expect_lt(coherence_gap(fc, h), 1e-8)
  expect_error(mt_fit(coupled(lee_carter(), template = reconciled(cae(), h, "bottom_up")), x),
               class = "mt_error_coupled")
})

test_that("mt_couple() couples draws and checks them", {
  x <- coupled_data()
  fb <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 50, seed = 1)
  ft <- mt_forecast(mt_fit(cae(), x), h = 3, n_draws = 50, seed = 2)
  fc <- mt_couple(fb, ft, seed = 3)
  expect_identical(sorted_margins(fc$rates), sorted_margins(fb$rates))
  expect_identical(margin_ranks(fc$rates), margin_ranks(ft$rates))
  expect_identical(dimnames(fc$rates), dimnames(fb$rates))
  expect_error(mt_couple(fb, mt_forecast(mt_fit(cae(), x), h = 2, n_draws = 50, seed = 2)), class = "mt_error_coupled")
  expect_error(mt_couple(fb, mt_forecast(mt_fit(cae(), x), h = 3, n_draws = 40, seed = 2)), class = "mt_error_coupled")
  ftp <- mt_forecast(mt_fit(cae(), x), h = 3, n_draws = 50, seed = 2, predictive = TRUE)
  expect_warning(mt_couple(fb, ftp), class = "mt_warning_coupled")
  h <- mt_hierarchy(quintile = "all")
  rc <- mt_reconcile(mt_forecast(mt_fit(lee_carter(), mt_aggregate_strata(x, h)), h = 3, n_draws = 50, seed = 1), h)
  tp <- mt_forecast(mt_fit(cae(), mt_aggregate_strata(x, h)), h = 3, n_draws = 50, seed = 2)
  expect_warning(mt_couple(rc, tp), class = "mt_warning_coupled")
})

test_that("coupled() warns about templates without dependence between strata and reports parameters", {
  expect_warning(coupled(lee_carter(), template = rwd()), class = "mt_warning_coupled")
  x <- coupled_data()
  fit <- mt_fit(coupled(lee_carter(), cae()), x)
  p <- mt_parameters(fit)
  expect_setequal(unique(p$component), c("base", "template"))
  expect_true(all(c("b", "beta1") %in% p$parameter))
  expect_match(mortensor:::model_assumptions(fit$spec), "ECC-R")
})
