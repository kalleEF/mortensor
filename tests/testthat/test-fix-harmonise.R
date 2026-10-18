# Regression tests for fixes in the harmonisation code (ungrouping, smoothing,
# disaggregation of draws, shock periods)

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

test_that("spline ungrouping works with grouped exposures", {
  x <- mt_example(years = 2010:2011)
  E <- mt_obs(x, "exposure")
  lo <- ifelse(E$age_lo >= 100, 100, 30 + 5 * floor((E$age_lo - 30) / 5))
  Eg <- data.table::data.table(stratum = E$stratum, age_lo = lo, age_hi = ifelse(lo >= 100, Inf, lo + 5),
                               period_lo = E$period_lo, period_hi = E$period_hi, value = E$value)
  Eg <- Eg[, list(value = sum(value), provenance = "observed"),
           by = c("stratum", "age_lo", "age_hi", "period_lo", "period_hi")]
  xg <- x
  xg$obs$exposure <- Eg
  xg$grid <- NULL
  for (m in c("hyman_kannisto", "natural_gompertz", "gam")) {
    if (m == "gam") skip_if_not_installed("mgcv")
    u <- suppressWarnings(mt_ungroup(xg, method = m), classes = "mt_warning_deprecated")
    dg <- attr(u, "ungroup")
    expect_false(anyNA(dg$at_bound))
    expect_equal(sum(mt_obs(u, "exposure")$value), sum(Eg$value), tolerance = 1e-6)
    expect_lt(group_gap(xg, u), 1e-6)
  }
})

test_that("spline ungrouping keeps group totals next to a missing group and above a closed top group", {
  x <- mt_example(years = 2010:2011)
  xm <- x
  xm$obs$deaths <- mt_obs(x, "deaths")[!(age_lo == 40 & period_lo == 2011 & stratum == 1)]
  xc <- x
  xc$obs$deaths <- mt_obs(x, "deaths")[is.finite(age_hi)]
  for (m in c("hyman_kannisto", "natural_gompertz", "gam")) {
    if (m == "gam") skip_if_not_installed("mgcv")
    u <- suppressWarnings(mt_ungroup(xm, method = m), classes = "mt_warning_deprecated")
    expect_lt(group_gap(xm, u), 1e-6)
    expect_equal(nrow(mt_obs(u)[stratum == 1 & period_lo == 2011 & age_lo >= 40 & age_lo < 45]), 0L)
    uc <- suppressWarnings(mt_ungroup(xc, method = m), classes = "mt_warning_deprecated")
    expect_lt(group_gap(xc, uc), 1e-6)
    expect_equal(max(mt_obs(uc)$age_hi), 90)
  }
})

test_that("a missing deaths cell does not break smoothing of the other cells", {
  x <- mt_example(years = 2006:2011, deaths = "single")
  xm <- x
  xm$obs$deaths <- mt_obs(x)[!(age_lo == 40 & period_lo == 2010 & stratum == 1)]
  args <- list(list(method = "whittaker"), list(method = "whittaker", along = "period"), list(method = "gam"))
  for (a in args) {
    if (a$method == "gam") skip_if_not_installed("mgcv")
    s <- do.call(mt_smooth, c(list(xm), a))
    D <- mt_obs(s)
    expect_equal(nrow(D), nrow(mt_obs(xm)))
    expect_false(anyNA(D$value))
    expect_equal(unique(D$provenance), "smoothed")
    # the other cells are smoothed much as without the missing cell
    s0 <- do.call(mt_smooth, c(list(x), a))
    v <- merge(D, mt_obs(s0), by = c("stratum", "age_lo", "period_lo"))
    expect_lt(stats::median(abs(log(v$value.x / v$value.y))), 0.02)
  }
})

test_that("smoothing over time works with a single age group", {
  x <- mt_regroup(mt_example(deaths = "single", years = 2000:2010), age = 30)
  s <- mt_smooth(x, along = "period")
  expect_equal(nrow(mt_obs(s)), nrow(mt_obs(x)))
  expect_false(anyNA(mt_obs(s)$value))
  expect_equal(unique(mt_obs(s)$provenance), "smoothed")
})

test_that("an observed age 0 is kept out of the rate PCLM and does not distort ages 1-6", {
  ages <- 0:99
  m <- ifelse(ages == 0, 5e-3, 4e-4 * exp(-0.6 * ages) + 1e-4 + 2e-5 * exp(0.09 * ages))
  e <- rep(1e5, 100)
  lo <- c(0, 1, seq(5, 90, 5))
  dg <- data.frame(year = rep(2000:2001, each = length(lo)), age_lo = lo, age_hi = c(lo[-1], Inf),
                   deaths = as.numeric(tapply(m * e, findInterval(ages, lo), sum)))
  ex <- data.frame(year = rep(2000:2001, each = 100), age_lo = ages, age_hi = c(ages[-100] + 1, Inf),
                   exposure = e)
  x <- suppressWarnings(mt_data(deaths = dg, exposure = ex, age = c("age_lo", "age_hi")))
  for (d in 1:2) {
    # rescaled to the group totals, so that the check is about age 0 and not
    # about how well the smooth curve follows the steep decline after it
    u <- mt_ungroup(x, dims = d, lambda = if (d == 2) c(10, 100), preserve_totals = TRUE)
    r <- mt_array(u, "deaths")[, , 1] / mt_array(u, "exposure")[, , 1]
    expect_lt(max(abs(log(r[2:7, ] / m[2:7]))), 0.5)
    expect_equal(as.vector(r[1, ]), rep(5e-3, 2), tolerance = 1e-10)
    expect_equal(unique(mt_obs(u)[age_lo == 0]$provenance), "observed")
    expect_lt(group_gap(x, u), 1e-8)
    # the published PCLM (default) keeps the observed age 0 as well
    u0 <- mt_ungroup(x, dims = d, lambda = if (d == 2) c(10, 100))
    r0 <- mt_array(u0, "deaths")[, , 1] / mt_array(u0, "exposure")[, , 1]
    expect_equal(as.vector(r0[1, ]), rep(5e-3, 2), tolerance = 1e-10)
    expect_lt(max(abs(log(r0[2:7, ] / m[2:7]))), 1)
  }
})

test_that("mt_smooth_check uses the smoothed cells only", {
  x <- mt_example(years = 2009:2011, deaths = "single")
  s <- mt_smooth(x, ages = c(50, 90))
  ck <- mt_smooth_check(x, s)
  D <- mt_array(x)
  Ds <- mt_array(s)
  sm <- mt_grid(x)$age$lower >= 50 & mt_grid(x)$age$upper <= 90
  pear <- (D - Ds) / sqrt(Ds)
  # sum of squared Pearson residuals over n - edf (McCullagh and Nelder 1989)
  edf <- tapply(attr(s, "smooth")$ed, attr(s, "smooth")$stratum, sum)
  n <- sum(sm) * 3
  expect_equal(ck$edf, as.vector(edf), tolerance = 1e-10)
  expect_equal(ck$dispersion, c(sum(pear[sm, , 1]^2), sum(pear[sm, , 2]^2)) / (n - edf), tolerance = 1e-10,
               ignore_attr = TRUE)
  # without the diagnostics of the fits, divided by n
  s2 <- s
  attr(s2, "smooth") <- NULL
  ck2 <- mt_smooth_check(x, s2)
  expect_equal(ck2$dispersion, c(mean(pear[sm, , 1]^2), mean(pear[sm, , 2]^2)), tolerance = 1e-10)
  expect_true(all(is.na(ck2$edf)))
})

test_that("disaggregated draws handle zero shape rates, missing draws and keep exposure names", {
  x <- mt_example(years = 2010:2019, deaths = "single")
  x <- mt_aggregate_strata(x, mt_hierarchy(sex = "total"))
  lo <- c(seq(30, 85, 5), 90)
  xg <- mt_regroup(x, age = lo)
  fc <- mt_forecast(mt_fit(lee_carter(), xg), h = 2, n_draws = 5, seed = 1)
  e <- mt_array(x, "exposure")[, 10, ]
  age <- mt_grid(x)$age
  # zero shape rates in a group: flat shape (factor 1) instead of NaN
  xz <- x
  xz$obs$deaths <- data.table::copy(xz$obs$deaths)
  xz$obs$deaths[age_lo < 35 & period_lo >= 2015, value := 0]
  sh <- mt_shape_rates(xz, periods = 2020:2021)
  fs <- mt_disaggregate(fc, age, exposure = e, method = "shape", shape_rates = sh)
  expect_false(anyNA(fs$rates))
  expect_equal(fs$rates[age$lower < 35, , , ], fc$rates[rep(1, 5), , , ], ignore_attr = TRUE)
  # missing draws of a stratum stay missing; a missing group leaves the others
  fc2 <- fc
  fc2$rates[, , 1, ] <- NA
  fc2$rates[3, 1, 2, 1] <- NA
  # default (proportional allocation by the shape): missing groups stay missing
  fm <- mt_disaggregate(fc2, age, exposure = e, shape_rates = sh)
  expect_true(all(is.na(fm$rates[, , 1, ])))
  expect_false(anyNA(fm$rates[, , 2:3, 2:5]))
  # deprecated spline: the spline goes through the other groups
  fn <- suppressWarnings(mt_disaggregate(fc2, age, exposure = e, method = "hyman_kannisto"),
                         classes = "mt_warning_deprecated")
  expect_true(all(is.na(fn$rates[, , 1, ])))
  grp <- findInterval(age$lower, lo)
  expect_true(all(is.na(fn$rates[grp == 3, 1, 2, 1])))
  expect_false(anyNA(fn$rates[grp != 3, 1, 2, 1]))
  expect_false(anyNA(fn$rates[, , 2:3, 2:5]))
  ok <- grp != 3
  wm <- rowsum(fn$rates[ok, 1, 2, 1] * e[ok, 2], grp[ok]) / rowsum(e[ok, 2], grp[ok])
  expect_equal(as.vector(wm), as.vector(fc2$rates[-3, 1, 2, 1]), tolerance = 1e-6)
  # exposures stored with names: no warning about matching by position
  fs2 <- mt_disaggregate(fc, age, exposure = e, shape_rates = sh)
  expect_identical(dimnames(fs2$exposure_last)$age, age$labels)
  expect_identical(dimnames(fs2$exposure_last)$stratum, dimnames(fc$rates)$stratum)
  expect_no_warning(suppressWarnings(mt_reconcile(fs2, var = "sex", total = "total"), classes = "mt_warning_reconcile"),
                    class = "mt_warning_exposure_order")
})

test_that("interpolated shocks at the edges of the data are removed with a warning", {
  x <- mt_example(years = 2001:2008, deaths = "single")
  expect_warning(y <- mt_prep_apply(mt_prep(step_shocks(c(2005, 2008))), x), class = "mt_warning_prep")
  expect_equal(max(mt_grid(y)$period$lower), 2007)
  expect_equal(attr(y, "excluded_periods"), 2008)
  expect_equal(unique(mt_obs(y)[period_lo == 2005]$provenance), "interpolated")
  expect_warning(y2 <- mt_prep_apply(mt_prep(step_shocks(2001:2002)), x), class = "mt_warning_prep")
  expect_equal(min(mt_grid(y2)$period$lower), 2003)
  expect_error(mt_prep_apply(mt_prep(step_shocks(2001:2008)), x), class = "mt_error_prep")
  # interior shocks are interpolated without a warning
  expect_no_warning(mt_prep_apply(mt_prep(step_shocks(2005)), x))
  # backtest origins inside the shock years do not fail
  expect_warning(bt <- mt_backtest(rwd(), x, origins = c(2005, 2006), h = 2, n_draws = 20,
                                   prep = mt_prep(step_shocks(2006:2007))), class = "mt_warning_prep")
  expect_setequal(unique(bt$origin), c(2005, 2006))
})

test_that("dims = 2 warns for the one-dimensional ungrouping methods", {
  x <- mt_example(years = 2010:2011)
  expect_warning(suppressWarnings(mt_ungroup(x, method = "hyman_kannisto", dims = 2),
                                  classes = "mt_warning_deprecated"),
                 class = "mt_warning_ungroup_dims")
})


test_that("the 1D GAM over age is quasi-Poisson with the exposures as offset", {
  skip_if_not_installed("mgcv")
  x <- mt_example(years = 2010, deaths = "single")
  s <- mt_smooth(x, method = "gam")
  g <- mt_grid(x)
  df <- data.frame(d = mt_array(x)[, 1, 1], e = mt_array(x, "exposure")[, 1, 1],
                   age = mortensor:::cell_ages(g$age$lower, g$age$upper))
  fit <- mgcv::gam(d ~ s(age, k = 20) + offset(log(e)), family = stats::quasipoisson(), data = df, method = "REML")
  expect_equal(as.vector(mt_array(s)[, 1, 1]), as.vector(stats::fitted(fit)), tolerance = 1e-6)
  expect_equal(attr(s, "smooth")$dispersion[1], fit$scale, tolerance = 1e-6)
})

test_that("an open interval is placed at its lower bound plus half the width of the interval below", {
  expect_equal(mortensor:::cell_ages(c(80, 85, 90), c(85, 90, Inf)), c(82.5, 87.5, 92.5))
  expect_equal(mortensor:::cell_ages(c(99, 100), c(100, Inf)), c(99.5, 100.5))
  expect_equal(mortensor:::cell_ages(90, Inf), 90.5)
  # the Whittaker smoother over age uses these ages for its divided differences
  x <- mt_regroup(mt_example(years = 2010, deaths = "single"), age = c(seq(30, 85, 5), 90))
  s <- mt_smooth(x, method = "whittaker", lambda = 10)
  D <- mt_array(x)[, 1, 1]
  E <- mt_array(x, "exposure")[, 1, 1]
  w <- D
  P <- crossprod(mortensor:::divided_differences(c(seq(32.5, 87.5, 5), 92.5), 2))
  z <- solve(diag(w) + 10 * P, w * log(D / E))
  expect_equal(as.vector(log(mt_array(s)[, 1, 1] / E)), as.vector(z), tolerance = 1e-8)
})
