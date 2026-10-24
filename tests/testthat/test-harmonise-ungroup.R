# Largest relative gap between observed group totals of deaths and the sums
# of the ungrouped cells inside each group
group_gap <- function(x, u) {
  D0 <- mt_obs(x, "deaths")
  Du <- mt_obs(u, "deaths")
  est <- vapply(seq_len(nrow(D0)), function(i) {
    sum(Du$value[Du$stratum == D0$stratum[i] & Du$period_lo == D0$period_lo[i] &
                   Du$age_lo >= D0$age_lo[i] - 1e-9 & Du$age_hi <= D0$age_hi[i] + 1e-9])
  }, 0)
  max(abs(est / D0$value - 1))
}

noise_free_fixture <- function(omega = 101) {
  ages <- 30:(omega - 1)
  m <- 2e-5 * exp(0.1 * ages) + 3e-4
  e <- 1e5 * exp(-0.03 * (ages - 30)) * (1 + 0.2 * sin(ages / 3))
  d <- m * e
  lo <- c(seq(30, 85, 5), 90)
  hi <- c(seq(35, 90, 5), omega)
  dg <- data.frame(year = 2000, age_lo = lo, age_hi = hi, deaths = tapply(d, findInterval(ages, lo), sum))
  ex <- data.frame(year = 2000, age_lo = ages, age_hi = ages + 1, exposure = e)
  list(m = m, x = suppressWarnings(mt_data(deaths = dg, exposure = ex, age = c("age_lo", "age_hi"))))
}

test_that("rate-based PCLM recovers a smooth schedule without a level shift", {
  f <- noise_free_fixture()
  # fixed, moderate smoothing: no level shift and group totals fitted without rescaling
  u <- mt_ungroup(f$x, lambda = 1, preserve_totals = FALSE)
  est <- mt_array(u, "deaths")[, 1, 1] / mt_array(u, "exposure")[, 1, 1]
  rel <- est / f$m - 1
  expect_lt(max(abs(rel)), 0.01)
  expect_lt(abs(mean(rel)), 0.002)
  dg <- attr(u, "ungroup")
  expect_equal(dg$method, "pclm_rate_1d")
  expect_lt(dg$max_rel_diff, 0.01)
  # smoothing selected by BIC; the published PCLM (default) reproduces the
  # total over all groups, not each group
  u2 <- mt_ungroup(f$x)
  est2 <- mt_array(u2, "deaths")[, 1, 1] / mt_array(u2, "exposure")[, 1, 1]
  # BIC on noise-free data prefers strong smoothing (the deviance is close to
  # 0): the curve misses the lowest group (where the rate bends) by a few
  # per cent and is close to the truth elsewhere
  rel2 <- abs(est2 / f$m - 1)
  expect_lt(max(rel2), 0.15)
  expect_lt(max(rel2[f$x$obs$exposure$age_lo >= 40]), 0.03)
  expect_equal(sum(mt_obs(u2, "deaths")$value), sum(mt_obs(f$x, "deaths")$value), tolerance = 1e-6)
  expect_gt(group_gap(f$x, u2), 1e-6)
  expect_equal(group_gap(f$x, u2), attr(u2, "ungroup")$max_rel_diff, tolerance = 1e-6)
  # rescaled to the group totals on request
  u3 <- mt_ungroup(f$x, preserve_totals = TRUE)
  expect_lt(group_gap(f$x, u3), 1e-10)
  est3 <- mt_array(u3, "deaths")[, 1, 1] / mt_array(u3, "exposure")[, 1, 1]
  expect_lt(max(abs(est3 / f$m - 1)), 0.05)
  expect_error(mt_ungroup(f$x, preserve_totals = NA), class = "mt_error_argument")
})

test_that("counts-only PCLM without offset matches group totals and warns about the open interval", {
  f <- noise_free_fixture()
  x <- f$x
  x$obs$exposure <- NULL
  expect_warning(u <- mt_ungroup(x), class = "mt_warning_ungroup_counts")
  expect_equal(sum(mt_obs(u, "deaths")$value), sum(mt_obs(x, "deaths")$value), tolerance = 1e-6)
  expect_warning(ut <- mt_ungroup(x, preserve_totals = TRUE), class = "mt_warning_ungroup_counts")
  expect_lt(group_gap(x, ut), 1e-10)
  expect_error(suppressWarnings(mt_ungroup(x, dims = 2)), class = "mt_error_ungroup")
})

test_that("open intervals are expanded to omega and exposures can be ungrouped first", {
  x <- mt_example(years = 2001:2002)
  xr <- suppressWarnings(mt_regroup(x, age = c(30, seq(35, 90, 5))))
  xr <- mt_set_grid(xr, age = single_ages(30, 100))
  expect_warning(u <- mt_ungroup(xr, exposure_method = "monotone"), class = "mt_warning_ungroup_open_counts")
  e_true <- mt_array(x, "exposure")
  e_est <- mt_array(u, "exposure")
  expect_equal(sum(e_est), sum(e_true), tolerance = 1e-8)
  expect_lt(median(abs(log(e_est[1:60, , ] / e_true[1:60, , ]))), 0.05)
  expect_false(anyNA(mt_array(u, "deaths")))
  # the PCLM of the exposures chooses its own smoothing parameter: a fixed
  # `lambda` for the deaths leaves the exposures unchanged
  e_bic <- mt_array(suppressWarnings(mt_ungroup(xr)), "exposure")
  for (lam in list(1e6, c(100, 10))) {
    u_lam <- suppressWarnings(mt_ungroup(xr, lambda = lam, dims = length(lam)))
    expect_identical(mt_array(u_lam, "exposure"), e_bic)
  }
})

test_that("ungrouped rates are close to the simulated truth", {
  x <- mt_example(years = 2001:2003)
  u <- mt_ungroup(x)
  est <- mt_array(u, "deaths") / mt_array(u, "exposure")
  truth <- attr(x, "truth")
  le <- abs(log(est / truth))
  # deaths in 5-year groups and 90+, exposures in single years and 100+
  ages <- mt_grid(u)$age$lower
  expect_lt(median(le[ages < 90, , ]), 0.03)
  expect_lt(median(le[ages >= 90, , ]), 0.03)
  prov <- unique(mt_obs(u, "deaths")$provenance)
  expect_true("pclm" %in% prov)
})

test_that("2D ungrouping runs and keeps totals", {
  x <- mt_example(years = 2001:2004)
  u <- mt_ungroup(x, dims = 2, lambda = c(100, 10))
  expect_equal(sum(mt_obs(u, "deaths")$value), sum(mt_obs(x, "deaths")$value), tolerance = 1e-6)
  expect_equal(attr(u, "ungroup")$method, rep("pclm_rate_2d", 2))
  ut <- mt_ungroup(x, dims = 2, lambda = c(100, 10), preserve_totals = TRUE)
  expect_lt(group_gap(x, ut), 1e-10)
  # rescaling puts the year-to-year noise of the groups back on the surface
  rough <- function(v) {
    r <- log(mt_array(v, "deaths") / mt_array(v, "exposure"))[1:60, , ]
    mean(abs(apply(r, c(1, 3), diff, differences = 2)))
  }
  expect_lt(rough(u), rough(ut))
})

test_that("the exposure of the open interval is split with stationary weights or by the PCLM of the exposures", {
  x <- mt_example(years = 2001:2002)
  D <- mt_obs(x, "deaths")[stratum == 1 & period_lo == 2001]
  E <- mt_obs(x, "exposure")[stratum == 1 & period_lo == 2001]
  # grid open at 100+, expanded to omega = 110 inside the PCLM
  u <- mt_ungroup(x, open_split = "stationary")
  dg <- attr(u, "ungroup")
  expect_true(all(dg$open_converged))
  expect_true(all(dg$open_iterations >= 1 & dg$open_iterations < 50))
  g <- mt_grid(u)
  e <- E$value[match(g$age$lower, E$age_lo)]
  idx <- data.table::data.table(a1 = findInterval(D$age_lo, g$age$lower),
                                a2 = findInterval(pmin(D$age_hi, 110) - 1e-9, g$age$lower), p1 = 1L, p2 = 1L)
  o <- cbind(D, idx)
  fa <- mortensor:::fine_axis(g$age, 110)
  bas <- mortensor:::pclm_basis((fa$lower + fa$upper) / 2, min(fa$lower), max(fa$upper))
  open <- which(fa$grid == length(g$age))
  # "stationary": at the fixed point the weights are the stationary weights of the fitted rates
  ctl <- list(lambda = dg$lambda_age[1], criterion = "bic", kr = 2, deg = 3, pord = 2, lambda_range = c(1e-2, 1e8),
              preserve_totals = FALSE, open_split = "stationary")
  r <- mortensor:::ungroup_rate_block(o, list(e), fa, bas$B, list(bas$P), ctl, 110, NULL)
  rate <- exp(as.numeric(bas$B %*% r$fit$beta))[open]
  w <- mortensor:::open_weights(rate)
  dopen <- r$deaths[r$a == length(g$age)]
  expect_equal(dopen, sum(e[length(e)] * w * rate), tolerance = 1e-4)
  # "pclm": the weights are the PCLM estimates of the exposures in the open
  # interval (exposures of all grid cells as counts, then 15 years with 0
  # counts beyond omega, Rizzi et al. 2015), with a smoothing parameter
  # chosen by BIC for the exposures; the fixed `lambda` of the deaths in
  # `ctl` does not enter
  up <- mt_ungroup(x, open_split = "pclm")
  expect_true(all(is.na(attr(up, "ungroup")$open_iterations)))
  ctl$open_split <- "pclm"
  r <- mortensor:::ungroup_rate_block(o, list(e), fa, bas$B, list(bas$P), ctl, 110, NULL)
  lo <- c(fa$lower, 110:124)
  hi <- c(fa$upper, 111:125)
  Ce <- 1 * outer(seq_len(length(e) + 1), c(fa$grid, rep(length(e) + 1, 15)), "==")
  be <- mortensor:::pclm_basis((lo + hi) / 2, min(lo), max(hi))
  fe <- mortensor:::pclm_fit(c(e, 0), Ce, hi - lo, be$B, list(be$P), lambda = NULL, criterion = "bic")
  we <- fe$gamma[open] / sum(fe$gamma[open])
  expect_equal(mortensor:::open_weights_pclm(e, fa, ctl), we)
  for (lam in list(NULL, 1e6, c(100, 10))) {
    expect_identical(mortensor:::open_weights_pclm(e, fa, modifyList(ctl, list(lambda = lam))),
                     mortensor:::open_weights_pclm(e, fa, ctl))
  }
  expect_false(isTRUE(all.equal(fe$lambda, ctl$lambda)))
  rate <- exp(as.numeric(bas$B %*% r$fit$beta))[open]
  expect_equal(r$deaths[r$a == length(g$age)], sum(e[length(e)] * we * rate), tolerance = 1e-8)
  # the exposures of the open interval decline with age
  expect_true(all(diff(we) < 0))
  expect_error(mt_ungroup(x, open_split = "none"))
})

test_that("an open deaths interval at 90+ is extrapolated as documented (French women, HMD)", {
  skip_if_not_installed("demography")
  fr <- suppressWarnings(as_mt_data(mortensor:::pkg_data("fr.mort", "demography"), series = "female", years = 2006))
  lo <- c(0, 1, seq(5, 90, 5))
  rate <- function(v) (mt_array(v, "deaths") / mt_array(v, "exposure"))[, 1, 1]
  # deaths in 5-year groups and 90+, exposures in single years to 110+
  x <- fr
  x$obs$deaths <- mt_regroup(fr, age = lo)$obs$deaths
  x$grid <- NULL
  u <- mt_ungroup(x, age = mt_grid(fr)$age)
  ages <- mt_grid(fr)$age$lower
  le <- log(rate(u) / rate(fr))
  old <- ages >= 90 & ages <= 109
  # ?mt_ungroup: median absolute log error 0.09 at 90-109, rates too high,
  # most at the highest ages
  expect_lt(median(abs(le[old])), 0.12)
  expect_gt(median(le[old]), 0)
  expect_gt(mean(le[ages >= 104 & ages <= 109]), 0.2)
  expect_lt(median(abs(le[ages >= 30 & ages < 90])), 0.03)
  # grid open at 100+: the PCLM expands it to omega and splits its exposure
  # with stationary weights (default) or by the PCLM of the exposures
  f100 <- mt_regroup(fr, age = 0:100)
  x100 <- f100
  x100$obs$deaths <- x$obs$deaths
  x100$grid <- NULL
  for (os in c("stationary", "pclm")) {
    u100 <- mt_ungroup(x100, age = single_ages(0, 100), open_split = os)
    if (os == "stationary") expect_true(attr(u100, "ungroup")$open_converged)
    e100 <- log(rate(u100) / rate(f100))
    expect_lt(median(abs(e100[90:99 + 1])), 0.06)
    # ?mt_ungroup: the rate of 100+ is too high by about 0.15 on the log scale
    expect_gt(e100[101], 0.05)
    expect_lt(e100[101], 0.25)
  }
})

test_that("the spline methods removed in 0.3.0 stop with an argument error", {
  x <- mt_example(years = 2018:2019)
  for (m in c("hyman_kannisto", "natural_gompertz", "gam")) {
    expect_error(mt_ungroup(x, method = m), class = "simpleError")
  }
})

test_that("forecast draws disaggregate to single ages with the group rates preserved", {
  x <- mt_example(years = 2000:2019, deaths = "single")
  lo <- c(seq(30, 85, 5), 90)
  xg <- mt_regroup(x, age = lo)
  fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 2, n_draws = 5, seed = 1)
  e <- mt_array(x, "exposure")[, 20, ]
  grp <- findInterval(mt_grid(x)$age$lower, lo)
  group_rates <- function(r) {
    out <- array(NA_real_, dim(fc$rates))
    for (s in 1:2) for (t in 1:2) for (j in 1:5) {
      out[, t, s, j] <- rowsum(r[, t, s, j] * e[, s], grp) / rowsum(e[, s], grp)
    }
    out
  }
  # proportional allocation by a standard schedule
  sh <- mt_disaggregate(fc, mt_grid(x)$age, exposure = e, shape_rates = mt_array(x, "deaths")[, 20, ] / e)
  expect_identical(dim(sh$rates), c(71L, 2L, 2L, 5L))
  expect_equal(as.vector(group_rates(sh$rates)), as.vector(fc$rates), tolerance = 1e-6)
  expect_error(mt_disaggregate(fc, mt_grid(x)$age, exposure = e), class = "mt_error_argument")
  # the spline methods were removed in 0.3.0
  for (m in c("hyman_kannisto", "natural_gompertz")) {
    expect_error(mt_disaggregate(fc, mt_grid(x)$age, exposure = e, method = m), class = "simpleError")
  }
})
