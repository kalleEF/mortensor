# Pessimism of forecasts by cause of death (Wilmoth 1995)

# Data with one age group and the causes c1, ... of each example in two years
two_year_data <- function(first, second, years = c(1951, 1990)) {
  one <- function(k, j, v) {
    data.table::data.table(example = k, year = years[j], cause = paste0("c", seq_along(v)), deaths = 10 * v)
  }
  rows <- data.table::rbindlist(lapply(names(first), function(k) rbind(one(k, 1, first[[k]]), one(k, 2, second[[k]]))))
  rows[, `:=`(age_lo = 60, age_hi = 65, exposure = 1e5)]
  mt_data(deaths = rows, exposure = rows, by = c("example", "cause"), age = c("age_lo", "age_hi"), period = "year",
          value = list(deaths = "deaths", exposure = "exposure"))
}

# First-to-last forecasts of every stratum (Wilmoth 1995, equations 8 and 9),
# one path: m(t) = m(t0) (m(t0) / m(tau))^((t - t0) / (t0 - tau))
first_to_last_draws <- function(x, years) {
  obs <- mt_observed(x, n_draws = 1)
  m <- obs$observed
  nt <- dim(m)[2]
  tau <- obs$periods[1]
  t0 <- obs$periods[nt]
  r <- m[, nt, ] / m[, 1, ]
  arr <- array(NA_real_, c(dim(m)[1], length(years), dim(m)[3], 1))
  for (j in seq_along(years)) arr[, j, , 1] <- m[, nt, ] * r^((years[j] - t0) / (t0 - tau))
  mortensor:::new_mt_draws(arr, obs$age, years, obs$strata, model = "first-to-last", origin = t0,
                           fit = list(observed = m, periods = obs$periods))
}

# Draws of causes with the given strata and rates [age, period, stratum, draw]
cause_draws <- function(fc, strata, rates, observed) {
  out <- fc
  out$strata <- strata
  out$rates <- rates
  out$observed <- observed
  out
}

test_that("the m-ratios of Wilmoth's Table 1 are reproduced", {
  # four hypothetical distributions of deaths by cause (percent) in 1951 and
  # 1990, m-ratios of the first-to-last forecasts in 2010 and 2040
  q <- list(es = c(1, 30, 30, 39), ew = c(10, 30, 30, 30), rs = c(20, 25, 25, 30), rw = c(25, 25, 25, 25))
  p <- list(es = c(20, 25, 25, 30), ew = c(25, 25, 25, 25), rs = c(1, 30, 30, 39), rw = c(10, 30, 30, 30))
  x <- two_year_data(q, p)
  r <- mt_pessimism(x, by = "age", to = 2040)[measure == "m_ratio" & period %in% c(2010, 2040)]
  m <- r$mean
  expect_equal(round(m[c(1, 2, 3, 5, 7, 8)], 2), c(1.65, 9.92, 1.08, 1.11, 1.05, 1.17))
  expect_equal(round(m[c(4, 6)], 1), c(1.4, 1.3))
  # equation (11)
  eq11 <- unlist(lapply(names(q), function(k) {
    vapply(c(2010, 2040), function(t) sum(p[[k]] / 100 * (p[[k]] / q[[k]])^((t - 1990) / 39)), 1)
  }))
  expect_equal(m, eq11, tolerance = 1e-12)
  # the contributions of the causes add up to the m-ratio minus 1 (equation 14)
  cc <- mt_pessimism(x, by = "cause", to = 2040)[period %in% c(2010, 2040)]
  expect_equal(cc[, sum(mean), by = .(example, period)]$V1, m - 1, tolerance = 1e-12)
  expect_equal(cc[example == "es" & period == 2010]$mean,
               p$es / 100 * ((p$es / q$es)^(20 / 39) - 1), tolerance = 1e-12)
})

test_that("the approximation is equation (12) and the first-to-last forecast", {
  x <- mt_subset(mt_example_causes(years = 2000:2019), ages = c(50, Inf))
  # first and second halves, midpoints 2004.5 and 2014.5
  a <- mt_pessimism(x, by = "age", to = 2030)[measure == "m_ratio" & period == 2030]
  D <- mt_array(x, "deaths")
  st <- mt_strata(x)
  for (s in c("female", "male")) {
    ids <- st$stratum[st$sex == s & st$cause != "all"]
    q <- apply(D[, 1:10, ids], c(1, 3), sum)
    q <- q / rowSums(q)
    p <- apply(D[, 11:20, ids], c(1, 3), sum)
    p <- p / rowSums(p)
    expect_equal(a[sex == s]$mean, unname(rowSums(p * (p / q)^((2030 - 2014.5) / 10))), tolerance = 1e-12)
  }
  # the default periods are the halves of the observed years
  expect_equal(mt_pessimism(x, to = 2030), mt_pessimism(x, periods = list(2000:2009, 2010:2019), to = 2030))
  # with an odd number of years the middle year is left out
  x19 <- mt_subset(x, periods = c(2001, 2019))
  expect_equal(mt_pessimism(x19, h = 3), mt_pessimism(x19, periods = list(2001:2009, 2011:2019), h = 3))
  # the measures of the approximation equal those of first-to-last forecasts
  # of the same data from the first and the last year
  years <- 2020:2040
  fc <- first_to_last_draws(x, years)
  for (b in c("total", "age", "cause")) {
    expect_equal(mt_pessimism(fc, by = b, at = 60), mt_pessimism(x, by = b, at = 60, periods = list(2000, 2019),
                                                                 to = 2040), tolerance = 1e-10)
  }
})

test_that("forecasts by cause are more pessimistic under the first-to-last approximation (Appendix B)", {
  skip_on_cran()
  causes <- names(mortensor::ons_data)[-(1:7)]
  x <- mt_data(deaths = mortensor::ons_data, population = mortensor::ons_data, by = c("sex", "imd_decile", "cause"),
               age = "age_group", value = list(deaths = c(all = "all_causes", causes), population = "population"))
  x <- mt_as_exposure(mt_subset(x, ages = c(50, Inf), imd_decile %in% c("1", "10")), type = "midyear")
  x <- mt_group_levels(x, "cause", mortensor::ons_cause_groups, other = "rest", total = "all", additive = "rates")
  # the twelve groups add up to all causes: no warning
  expect_no_warning(a <- mt_pessimism(x, by = "age", to = 2035))
  m <- a[measure == "m_ratio"]
  expect_true(all(m$mean > 1))
  # the m-ratio rises with the horizon
  expect_true(all(m[, all(diff(mean) > 0), by = .(sex, imd_decile, age)]$V1))
  e <- mt_pessimism(x, to = 2035, at = 50)
  expect_true(all(e[measure == "difference"]$mean > 0))
  expect_equal(e[measure == "difference"]$mean, e[measure == "total"]$mean - e[measure == "causes"]$mean)
  # contributions add up to the m-ratio minus 1
  cc <- mt_pessimism(x, by = "cause", to = 2035)
  s <- cc[, list(s = sum(mean)), by = .(sex, imd_decile, period, age)]
  expect_equal(s$s, m$mean - 1, tolerance = 1e-12)
})

test_that("one cause, or causes with the same proportional trend, give no pessimism", {
  base <- mt_example(deaths = "single", years = 2000:2019)
  long <- as.data.table(base)
  d <- long[measure == "deaths"]
  ex <- long[measure == "exposure"]
  # one cause equal to all causes
  one <- rbind(cbind(d, cause = "all"), cbind(d, cause = "only"))
  xe <- rbind(cbind(ex, cause = "all"), cbind(ex, cause = "only"))
  x1 <- mt_data(deaths = one[, .(sex, cause, period_lo, age_lo, age_hi, deaths = value)],
                exposure = xe[, .(sex, cause, period_lo, age_lo, age_hi, exposure = value)],
                by = c("sex", "cause"), age = c("age_lo", "age_hi"), period = "period_lo")
  a <- mt_pessimism(x1, by = "age", to = 2030)
  expect_equal(a[measure == "m_ratio"]$mean, rep(1, nrow(a) / 3), tolerance = 1e-14)
  expect_equal(mt_pessimism(x1, to = 2030)[measure == "difference"]$mean, rep(0, 22), tolerance = 1e-12)
  expect_equal(max(abs(mt_pessimism(x1, by = "cause", to = 2030)$mean)), 0, tolerance = 1e-14)
  # three causes with constant shares: the same proportional change
  sh <- c(cvd = 2, cancer = 3, other = 5)
  sh_all <- c(all = 10, sh)
  cells <- d[, .(sex, period_lo, age_lo, age_hi, value)]
  parts <- rbindlist(lapply(names(sh_all), function(k) {
    cells[, .(sex, period_lo, age_lo, age_hi, cause = k, deaths = sh_all[[k]] * value)]
  }))
  xe3 <- rbindlist(lapply(names(sh_all), function(k) {
    ex[, .(sex, period_lo, age_lo, age_hi, cause = k, exposure = 10 * value)]
  }))
  x3 <- mt_data(deaths = parts, exposure = xe3, by = c("sex", "cause"), age = c("age_lo", "age_hi"),
                period = "period_lo")
  expect_equal(mt_pessimism(x3, by = "age", to = 2030)[measure == "m_ratio"]$mean, rep(1, 71 * 11 * 2),
               tolerance = 1e-12)
  expect_equal(mt_pessimism(x3, to = 2030)[measure == "difference"]$mean, rep(0, 22), tolerance = 1e-10)

  # forecasts: the cause forecasts are the total forecast times fixed shares
  fc <- mt_forecast(mt_fit(lee_carter(), base), h = 5, n_draws = 20, seed = 1)
  st1 <- data.table::data.table(stratum = 1:2, sex = c("female", "male"), cause = "only",
                                label = c("female.only", "male.only"))
  only <- cause_draws(fc, st1, fc$rates, fc$observed)
  p1 <- mt_pessimism(fc, causes = only, pairing = "draw", by = "age", paths = TRUE)
  expect_identical(p1[measure == "m_ratio"]$value, rep(1, nrow(p1) / 3))
  st3 <- data.table::data.table(stratum = 1:6, sex = rep(c("female", "male"), 3),
                                cause = rep(names(sh), each = 2))
  st3[, label := paste(sex, cause, sep = ".")]
  w <- rep(sh / 10, each = 2)
  r3 <- fc$rates[, , rep(1:2, 3), , drop = FALSE] * rep(w, each = prod(dim(fc$rates)[1:2]))
  o3 <- fc$observed[, , rep(1:2, 3), drop = FALSE] * rep(w, each = prod(dim(fc$observed)[1:2]))
  prop <- cause_draws(fc, st3, r3, o3)
  pc <- mt_pessimism(fc, causes = prop, pairing = "draw", by = "cause", paths = TRUE)
  expect_lt(max(abs(pc$value)), 1e-12)
  pt <- mt_pessimism(fc, causes = prop, pairing = "draw", paths = TRUE)
  expect_lt(max(abs(pt[measure == "difference"]$value)), 1e-10)
})

test_that("signs follow Wilmoth: m-ratio above 1 and positive differences mean pessimism (equation 7)", {
  # one open age group (life expectancy 1 / m); total multiplier s, causes
  # s * c(0.9, 1.05) with shares 1/2 at the jump-off: sum p log(s_i / s) < 0,
  # so the m-ratio falls below 1 at first and exceeds it later (Appendix B)
  age <- mt_axis(0, Inf, "0+")
  yrs <- 2001:2030
  s <- 0.98
  m0 <- 0.02
  tot <- array(m0 * s^seq_along(yrs), c(1, length(yrs), 1, 1))
  cau <- array(NA_real_, c(1, length(yrs), 2, 1))
  cau[1, , 1, 1] <- m0 / 2 * (0.9 * s)^seq_along(yrs)
  cau[1, , 2, 1] <- m0 / 2 * (1.05 * s)^seq_along(yrs)
  obs <- array(m0 / 2, c(1, 1, 2), dimnames = list("0+", "2000", c("c1", "c2")))
  fx <- mortensor:::new_mt_draws(tot, age, yrs, data.table::data.table(stratum = 1L, label = "all"), model = "total",
                                 origin = 2000, fit = list(observed = array(m0, c(1, 1, 1)), periods = 2000))
  fcs <- mortensor:::new_mt_draws(cau, age, yrs, data.table::data.table(stratum = 1:2, cause = c("c1", "c2"),
                                                                        label = c("c1", "c2")),
                                  model = "causes", origin = 2000, fit = list(observed = obs, periods = 2000))
  m <- mt_pessimism(fx, causes = fcs, pairing = "draw", by = "age")[measure == "m_ratio"]$mean
  expect_equal(m, 0.5 * (0.9^seq_along(yrs) + 1.05^seq_along(yrs)), tolerance = 1e-12)
  expect_lt(m[1], 1)
  expect_gt(m[length(m)], 1)
  e <- mt_pessimism(fx, causes = fcs, pairing = "draw")
  dif <- e[measure == "difference"]$mean
  expect_equal(dif, 1 / tot[1, , 1, 1] - 1 / (m * tot[1, , 1, 1]), tolerance = 1e-10)
  expect_identical(sign(dif), sign(m - 1))
})

test_that("forecasts held in one object or two give the same measures", {
  x <- mt_example_causes(years = 2000:2019)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 30, seed = 1)
  p <- mt_pessimism(fc, at = 50)
  expect_equal(nrow(p), 2 * 5 * 3)
  expect_named(p, c("sex", "period", "measure", "at", "mean", "q0.025", "q0.5", "q0.975"))
  expect_equal(mt_pessimism(fc, causes = fc, pairing = "draw", at = 50), p)
  # independent pairing shuffles the draws of the causes: means are unchanged
  pind <- mt_pessimism(fc, causes = fc, at = 50)
  expect_equal(pind$mean, p$mean, tolerance = 1e-12)
  expect_equal(pind[measure != "difference"], p[measure != "difference"], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(pind[measure == "difference"]$q0.025, p[measure == "difference"]$q0.025)))
  # paths and their summary
  pp <- mt_pessimism(fc, at = 50, paths = TRUE)
  expect_equal(nrow(pp), 30 * nrow(p))
  expect_equal(pp[, mean(value), by = .(sex, period, measure)]$V1, p$mean, tolerance = 1e-12)
  e <- life_expectancy(fc, at = 50, paths = TRUE)
  expect_equal(pp[measure == "total" & sex == "female" & period == 2020]$value, unname(e["2020", "female.all", ]))
  # m-ratio and contributions path by path
  ma <- mt_pessimism(fc, by = "age", paths = TRUE)
  r <- fc$rates
  expect_equal(ma[measure == "m_ratio" & sex == "male" & period == 2022 & age == "80"]$value,
               unname((r["80", "2022", "male.cvd", ] + r["80", "2022", "male.cancer", ] +
                         r["80", "2022", "male.other", ]) / r["80", "2022", "male.all", ]))
  mc <- mt_pessimism(fc, by = "cause", paths = TRUE)
  s <- mc[, list(s = sum(value)), by = .(sex, period, age, draw)]
  expect_equal(s$s, ma[measure == "m_ratio"][order(sex, period, draw, match(age, fc$age$labels))]$value - 1,
               tolerance = 1e-12)
  obs <- fc$observed[, "2019", ]
  share <- obs["80", "male.cancer"] / sum(obs["80", c("male.cvd", "male.cancer", "male.other")])
  expect_equal(mc[sex == "male" & period == 2022 & age == "80" & cause == "cancer"]$value,
               unname(r["80", "2022", "male.cancer", ] / r["80", "2022", "male.all", ] - share))
})

test_that("the causes of coda_causes() add up to its all-cause strata", {
  skip_on_cran()
  x <- mt_example_causes(years = 2005:2019)
  total <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 20, seed = 1)
  fc <- mt_forecast(mt_fit(coda_causes(), x), h = 3, n_draws = 20, seed = 2)
  p <- mt_pessimism(total, causes = fc, pairing = "draw", paths = TRUE)
  e <- life_expectancy(fc, paths = TRUE)
  expect_equal(p[measure == "causes" & sex == "male" & period == 2021]$value, unname(e["2021", "male.all", ]),
               tolerance = 1e-10)
  et <- life_expectancy(total, paths = TRUE)
  expect_equal(p[measure == "total" & sex == "male" & period == 2021]$value, unname(et["2021", "male.all", ]))
})

test_that("Poisson draws give intervals of the approximation", {
  x <- mt_subset(mt_example_causes(years = 2000:2019), ages = c(60, Inf))
  a1 <- mt_pessimism(x, to = 2030)
  a <- mt_pessimism(x, to = 2030, n_draws = 50, seed = 3)
  expect_identical(a, mt_pessimism(x, to = 2030, n_draws = 50, seed = 3))
  expect_true(all(a$q0.025 <= a$q0.5 & a$q0.5 <= a$q0.975))
  expect_true(all(a[measure == "difference"]$q0.975 > a[measure == "difference"]$q0.025))
  expect_true(all(a1$mean >= a$q0.025 - 1e-8 & a1$mean <= a$q0.975 + 1e-8))
})

test_that("missing causes, zero deaths and wrong inputs are reported", {
  x <- mt_subset(mt_example_causes(years = 2000:2019), ages = c(50, Inf))
  # causes that do not add up to all causes
  xs <- mt_subset(x, cause != "other")
  expect_warning(mt_pessimism(xs, h = 2), class = "mt_warning_pessimism")
  fc <- mt_forecast(mt_fit(rwd(), xs), h = 2, n_draws = 10, seed = 1)
  expect_warning(mt_pessimism(fc), class = "mt_warning_pessimism")
  # a cause without deaths in the first period
  q <- list(a = c(0, 50, 50))
  p <- list(a = c(10, 40, 50))
  expect_warning(r <- mt_pessimism(two_year_data(q, p), by = "age", h = 5), class = "mt_warning_pessimism")
  expect_true(all(is.finite(r$mean)))
  # a cause without deaths in both periods contributes nothing
  z <- mt_pessimism(two_year_data(list(a = c(0, 50, 50)), list(a = c(0, 40, 60))), by = "cause", h = 5)
  expect_identical(z[cause == "c1"]$mean, rep(0, 5))
  # inputs
  fx <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 10, seed = 1)
  expect_error(mt_pessimism(fx, var = "region"), class = "mt_error_argument")
  expect_error(mt_pessimism(fx, total_level = "total"), class = "mt_error_argument")
  expect_error(mt_pessimism(mt_subset(x, cause == "all")), class = "mt_error_argument")
  expect_error(mt_pessimism(x, causes = fx), class = "mt_error_argument")
  expect_error(mt_pessimism(x, periods = list(2000:2010, 2010:2019)), class = "mt_error_argument")
  expect_error(mt_pessimism(x, periods = list(1990:1999, 2010:2019)), class = "mt_error_argument")
  expect_error(mt_pessimism(x, periods = 2000:2019), class = "mt_error_argument")
  fy <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 5, seed = 1)
  expect_error(mt_pessimism(fx, causes = fy), class = "mt_error_argument")
  f60 <- mt_forecast(mt_fit(rwd(), mt_subset(x, ages = c(60, Inf))), h = 2, n_draws = 10, seed = 1)
  expect_error(mt_pessimism(fx, causes = f60), class = "mt_error_argument")
  agg <- fx
  agg$observed <- NULL
  expect_error(mt_pessimism(agg, by = "cause"), class = "mt_error_pessimism")
  expect_error(mt_pessimism(fx, at = 52.5), class = "mt_error_lifetable")
})
