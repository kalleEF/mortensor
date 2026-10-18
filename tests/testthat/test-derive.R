test_that("life table identities hold", {
  m <- 3e-4 + 2e-5 * exp(0.1 * (30:100))
  lt <- life_table(m, single_ages(30, 100))
  expect_equal(lt$lx[-1], lt$lx[-nrow(lt)] - lt$dx[-nrow(lt)])
  expect_equal(sum(lt$dx), lt$lx[1])
  expect_equal(lt$Tx[1], sum(lt$Lx))
  expect_equal(lt$ex, lt$Tx / lt$lx)
  expect_equal(tail(lt$qx, 1), 1)
})

test_that("D10: life expectancy uses interval widths (single and abridged ages agree)", {
  ages <- 0:109
  m <- 2e-5 * exp(0.095 * ages) + 2e-4
  e_single <- life_table(c(m[1:100], m[101]), single_ages(0, 100))$ex[1]
  grp <- c(0, 1, seq(5, 100, 5))
  ax <- age_groups(grp)
  # abridged rates from the single-year schedule with stationary weights
  lt <- life_table(c(m[1:100], m[101]), single_ages(0, 100))
  gid <- findInterval(0:100, grp)
  mg <- tapply(lt$dx, gid, sum) / tapply(lt$Lx, gid, sum)
  e_abr <- life_table(as.numeric(mg), ax)$ex[1]
  expect_lt(abs(e_single - e_abr), 0.1)
})

test_that("Gompertz closure of the open interval matches the exact integral", {
  b <- 0.1
  a <- 1e-5
  m <- a * exp(b * (30:99 + 0.5))
  open_rate <- 0.5
  lt_g <- life_table(c(m, open_rate), single_ages(30, 100), open = "gompertz")
  # exact remaining life expectancy at 100 under mu(x) = a exp(b x)
  surv <- function(t) exp(-a / b * (exp(b * (100 + t)) - exp(b * 100)))
  e100 <- stats::integrate(surv, 0, 40)$value
  expect_lt(abs(tail(lt_g$ex, 1) - e100), 0.05)
  expect_error(life_table(1:3 / 10, single_ages(0, 5)), class = "mt_error_lifetable")
})

test_that("life_expectancy works on arrays and checks `at`", {
  x <- mt_example(years = 2001:2003, deaths = "single")
  m <- mt_array(x, "deaths") / mt_array(x, "exposure")
  e <- life_expectancy(m, mt_grid(x)$age, at = 65)
  expect_equal(dim(e), c(3, 2))
  expect_true(all(e > 15 & e < 30))
  expect_error(life_expectancy(m, mt_grid(x)$age, at = 65.5), class = "mt_error_lifetable")
  # a missing rate gives a missing life expectancy, not certain death
  m[1, 3, 1] <- NA
  e0 <- life_expectancy(m, mt_grid(x)$age)
  expect_true(is.na(e0[3, 1]))
  expect_false(anyNA(e0[-3, ]))
})

test_that("expected deaths use future exposures", {
  x <- mt_example(years = 2001:2012, deaths = "single")
  fc <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 10, seed = 1)
  mu <- mt_deaths(fc, noise = FALSE)
  expect_equal(mu[, 1, , 1], fc$rates[, 1, , 1] * fc$exposure_last[, 1, ])
  E2 <- fc$exposure_last[, 1, ] * 2
  expect_equal(mt_deaths(fc, exposure = E2, noise = FALSE), 2 * mu)
  expect_error(suppressWarnings(mt_deaths(fc, exposure = matrix(1, 3, 2))), class = "mt_error_argument")
})

test_that("HMD tables are converted", {
  d <- data.frame(Year = 2000, Age = 0:2, Female = c(10, 1, 1), Male = c(12, 2, 1), Total = c(22, 3, 2),
                  OpenInterval = c(FALSE, FALSE, TRUE))
  e <- transform(d, Female = 1000, Male = 1000, Total = 2000)
  x <- mt_from_hmd(d, e)
  expect_equal(mt_strata(x)$sex, c("female", "male"))
  expect_equal(tail(mt_grid(x)$age$labels, 1), "2+")
  expect_error(mt_read_hmd("SWE", username = "", password = ""), class = "mt_error_adapter")
})

test_that("probability of dying, age-standardised rates and YLL are computed per path", {
  x <- mt_example(deaths = "single", years = 2001:2015)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 40, seed = 1)
  p <- prob_dying(fc, from = 30, to = 70)
  expect_equal(nrow(p), 6)
  expect_true(all(p$q0.5 > 0 & p$q0.5 < 1))
  # equals 1 - l70 / l30 of the median path's life table within sampling error
  m <- apply(fc$rates[, 1, 1, ], 1, median)
  lt <- life_table(m, fc$age)
  expect_equal(p$q0.5[1], 1 - lt$lx[41] / lt$lx[1], tolerance = 0.05)
  a <- asr(fc)
  expect_true(all(a$mean > 0))
  esp <- standard_population()
  expect_equal(sum(esp$weight), 1e5)
  w <- mortensor:::standard_weights(fc$age, esp)
  expect_equal(sum(w), sum(esp$weight[esp$age_lo >= 30]))
  y <- yll(fc)
  expect_true(all(y$mean > 0))
  expect_error(prob_dying(fc, from = 70, to = 30), class = "mt_error_argument")
})

test_that("forecasts are exported as summaries and draws", {
  x <- mt_example(deaths = "single", years = 2005:2015)
  fc <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 5, seed = 1)
  f1 <- mt_export(fc, tempdir())
  f2 <- mt_export(fc, tempdir(), what = "draws")
  on.exit(unlink(c(f1, f2)))
  s <- data.table::fread(f1)
  d <- data.table::fread(f2)
  expect_equal(nrow(s), 71 * 2 * 2)
  expect_equal(nrow(d), 71 * 2 * 2 * 5)
  expect_match(basename(f1), "origin2015_summary.csv")
  expect_error(mt_export(fc, file.path(tempdir(), "missing-dir")), class = "mt_error_argument")
})

test_that("standard weights cover open grid intervals and paths can be returned", {
  esp <- standard_population("esp2013")
  w <- standard_weights(age_groups(c(seq(50, 85, 5), 90)), esp)
  expect_equal(w[length(w)], sum(esp$weight[esp$age_lo >= 90]))
  expect_equal(sum(w), sum(esp$weight[esp$age_lo >= 50]))
  w1 <- standard_weights(single_ages(0, 100), esp)
  expect_equal(sum(w1), sum(esp$weight))
  x <- mt_example(deaths = "single", years = 2005:2019)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 20, seed = 1)
  a <- asr(fc, paths = TRUE)
  expect_identical(dim(a), c(3L, 2L, 20L))
  expect_equal(apply(a, 1:2, stats::median)[, 1], asr(fc)[stratum == 1]$q0.5, ignore_attr = TRUE)
  expect_identical(dim(prob_dying(fc, 30, 70, paths = TRUE)), c(3L, 2L, 20L))
})

# Draws of one stratum and period with given rates and exposures of the
# last training year
hand_draws <- function(age, m, E) {
  n <- length(age)
  st <- data.table::data.table(stratum = 1L, label = "all")
  fit <- list(exposure_last = array(E, c(n, 1, 1), dimnames = list(age$labels, "2019", "all")))
  mortensor:::new_mt_draws(array(m, c(n, 1, 1, 1)), age, 2020, st, model = "test", origin = 2019, fit = fit)
}

test_that("asr() standardises directly on the age groups of the standard", {
  age <- single_ages(80, 92)
  m <- exp(-4 + 0.1 * (0:12))
  E <- seq(1000, 400, length.out = 13)
  dr <- hand_draws(age, m, E)
  # ESP 2013: 80-84 (2500), 85-89 (1500), and the open cell 92+ joins 90-94
  # and 95+ (800 + 200) with the cells 90 and 91
  grp <- rep(1:3, c(5, 5, 3))
  W <- c(2500, 1500, 1000)
  rate_g <- tapply(m * E, grp, sum) / tapply(E, grp, sum)
  expect_equal(asr(dr)$mean, sum(W * rate_g) / sum(W) * 1e5, tolerance = 1e-12)
  # the weights of the standard shared among the grid cells
  w <- c(rep(500, 5), rep(300, 5), 160, 160, 480 + 200)
  expect_equal(mortensor:::standard_weights(age, standard_population()), w)
  expect_equal(asr(dr, age_groups = "grid")$mean, sum(w * m) / sum(w) * 1e5, tolerance = 1e-12)
  # given exposures replace those of the last training year
  E2 <- matrix(rev(E), 13, 1, dimnames = list(age$labels, "all"))
  rate_2 <- tapply(m * rev(E), grp, sum) / tapply(rev(E), grp, sum)
  expect_equal(asr(dr, exposure = E2)$mean, sum(W * rate_2) / sum(W) * 1e5, tolerance = 1e-12)
  # a group without exposures: the plain mean of its rates
  E0 <- E
  E0[1:5] <- 0
  rate_0 <- c(mean(m[1:5]), rate_g[2:3])
  expect_equal(asr(hand_draws(age, m, E0))$mean, sum(W * rate_0) / sum(W) * 1e5, tolerance = 1e-12)
  # coarser grid cells collect the standard groups they cover, cells that
  # straddle a bound join the neighbouring groups, and a standard group the
  # grid covers in part counts with the covered share of its weight
  esp <- standard_population("esp2013")
  g <- mortensor:::asr_groups(age_groups(c(0, 5, 10, 90)), esp)
  expect_identical(g$group, 1:4)
  expect_equal(g$weight, c(5000, 5500, sum(esp$weight[esp$age_lo >= 10 & esp$age_lo < 90]), 1000))
  g2 <- mortensor:::asr_groups(mt_axis(c(0, 3, 7), c(3, 7, 12), kind = "age"), esp)
  expect_identical(g2$group, rep(1L, 3))
  expect_equal(g2$weight, 1000 + 4000 + 5500 + 5500 * 2 / 5)
  # on the standard's own age groups both ways agree
  ag <- age_groups(c(0, 1, seq(5, 95, 5)))
  mm <- exp(seq(-7, -1, length.out = 21))
  d5 <- hand_draws(ag, mm, seq(9e4, 1e3, length.out = 21))
  expect_equal(asr(d5)$mean, asr(d5, age_groups = "grid")$mean, tolerance = 1e-12)
})

test_that("yll() uses the life expectancy at the average age at death", {
  age <- age_groups(seq(50, 90, 5))
  m <- exp(seq(-5.5, -1.5, length.out = 9))
  E <- seq(5e4, 1e4, length.out = 9)
  dr <- hand_draws(age, m, E)
  lt <- life_table(m, age)
  # a_x = n / 2: the life expectancy halfway through each closed interval
  e_death <- c((lt$ex[-9] + lt$ex[-1]) / 2, lt$ex[9])
  expect_equal(yll(dr)$mean, sum(m * E * e_death), tolerance = 1e-10)
  # a reference life table, interpolated in the same way
  ref <- seq(30, 5, length.out = 9)
  expect_equal(yll(dr, reference_ex = ref)$mean, sum(m * E * c((ref[-9] + ref[-1]) / 2, ref[9])),
               tolerance = 1e-10)
  expect_error(yll(dr, reference_ex = 1:3), class = "mt_error_argument")
})

test_that("mt_observed() keeps the observed rates out of the bootstrap draws", {
  x <- mt_example(deaths = "single", years = 2015:2016)
  o1 <- mt_observed(x, n_draws = 1)
  expect_equal(o1$rates[, , , 1], o1$observed, ignore_attr = TRUE)
  o5 <- mt_observed(x, n_draws = 5)
  expect_equal(o5$observed, o1$observed)
  expect_false(isTRUE(all.equal(o5$rates[, , , 1], o5$observed, check.attributes = FALSE)))
  # cells without deaths get no variation
  zero <- which(o5$observed == 0)
  if (length(zero)) expect_true(all(matrix(o5$rates, ncol = 5)[zero, ] == 0))
})

test_that("observed inequality indices use the population weights of each year", {
  x <- mt_example_regions(years = 2010:2012, regions = 1)
  x <- mt_subset(x, sex == "female", region == unique(mt_strata(x)$region)[1])
  iq <- mt_inequality(x, "quintile", quantity = "rate", age = mt_grid(x)$age$labels[10])[index == "sii"]
  ob <- mt_observed(x)
  r <- ob$observed[10, , ]
  E <- apply(ob$exposure, c(2, 3), sum)
  manual <- vapply(seq_len(nrow(E)), function(t) {
    w <- E[t, ] / sum(E[t, ])
    rk <- cumsum(w) - w / 2
    stats::coef(stats::lm(r[t, ] ~ rk, weights = w))[2]
  }, 0)
  expect_equal(iq$mean, manual, tolerance = 0.05)
  expect_equal(nrow(iq), 3)
})

test_that("mt_inequality() leaves out levels not in `order`", {
  x <- mt_example_regions(years = 2010:2012, regions = 2)
  xa <- mt_aggregate_strata(x, mt_hierarchy(quintile = list(total = "all")))
  a <- mt_inequality(x, "quintile", order = paste0("Q", 5:1), quantity = "asr")
  b <- mt_inequality(xa, "quintile", order = paste0("Q", 5:1), quantity = "asr")
  expect_equal(b$q0.5, a$q0.5, tolerance = 0.01)  # Poisson draws of the observed rates
  expect_error(mt_inequality(x, "quintile", order = c("Q1", "Q9")), class = "mt_error_argument")
})


# Integrated hazard over [x, x + 1] of a Gompertz or Kannisto hazard
tail_H <- function(a, b, x, law) {
  if (law == "kannisto") (log1p(a * exp(b * (x + 1))) - log1p(a * exp(b * x))) / b else a / b * exp(b * x) * expm1(b)
}

test_that("Gompertz and Kannisto closures are fitted by binomial maximum likelihood and run to omega", {
  set.seed(2)
  m <- 1e-5 * exp(0.1 * (30:99 + 0.5)) * exp(stats::rnorm(70, 0, 0.05))
  age <- single_ages(30, 100)
  for (law in c("gompertz", "kannisto")) {
    lt <- life_table(c(m, 0.5), age, open = law, k = 10, omega = 105)
    tl <- attr(lt, "tail")
    expect_true(tl$converged)
    expect_equal(c(tl$age_from, tl$age_to), c(90, 100))
    # the same by hand: binomial likelihood of the life-table survivors and
    # deaths at ages 90-99, maximised by optim()
    x <- 90:99
    q <- m[61:70] / (1 + m[61:70] / 2)
    l <- cumprod(c(1, 1 - m[1:69] / (1 + m[1:69] / 2)))[61:70]
    # parameters log mu(95) and b, which are nearly uncorrelated
    nll <- function(p) {
      H <- tail_H(exp(p[1] - 95 * p[2]), p[2], x, law)
      -sum(l * q * log(-expm1(-H)) - l * (1 - q) * H)
    }
    p_fit <- c(log(tl$a) + 95 * tl$b, tl$b)
    o <- stats::optim(p_fit + c(0.2, -0.01), nll, method = "BFGS", control = list(reltol = 1e-15, maxit = 500))
    expect_equal(p_fit, o$par, tolerance = 1e-5)
    expect_lte(nll(p_fit), o$value + 1e-12)
    # open interval: years 100-104 with the fitted integrated hazard, closed at 105
    H <- tail_H(tl$a, tl$b, 100:105, law)
    lo <- c(1, cumprod(exp(-H[1:5])))
    e_open <- sum(lo[1:5] * (1 - exp(-H[1:5])) / H[1:5]) + lo[6] / H[6]
    expect_equal(tail(lt$ex, 1), e_open, tolerance = 1e-10)
  }
  # the closure at omega makes the result insensitive to omega
  e110 <- tail(life_table(c(m, 0.5), age, open = "gompertz")$ex, 1)
  e130 <- tail(life_table(c(m, 0.5), age, open = "gompertz", omega = 130)$ex, 1)
  expect_lt(abs(e110 - e130), 0.01)
  # observed exposures weight the fit; draws are fitted column by column
  E <- 1e5 * exp(-0.08 * (0:70))
  lte <- life_table(c(m, 0.5), age, open = "kannisto", exposure = E)
  expect_false(isTRUE(all.equal(attr(lte, "tail")$b, attr(life_table(c(m, 0.5), age, open = "kannisto"), "tail")$b)))
  M <- cbind(c(m, 0.5), c(m * 1.1, 0.55))
  e <- life_expectancy(M, age, at = 100, open = "kannisto")
  expect_equal(e[2], tail(life_table(M[, 2], age, open = "kannisto")$ex, 1))
  expect_error(life_table(c(m, 0.5), age, open = "kannisto", exposure = 1:3), class = "mt_error_lifetable")
})

test_that("the closures reproduce the parameters of Thatcher, Kannisto and Vaupel (1998)", {
  # probabilities of dying at ages 80-98 of the cohorts born in 1871-1880 (11
  # countries; their Table 6.1); for a cohort the life-table survivors are the
  # numbers reaching each age, which their maximum-likelihood fit uses
  q <- list(male = c(0.1224, 0.1326, 0.1452, 0.1565, 0.1702, 0.1822, 0.1965, 0.2116, 0.2276, 0.2424, 0.2580, 0.2749,
                     0.2907, 0.3073, 0.3253, 0.3434, 0.3571, 0.3691, 0.3817),
            female = c(0.0987, 0.1074, 0.1180, 0.1291, 0.1409, 0.1520, 0.1653, 0.1787, 0.1926, 0.2070, 0.2225, 0.2354,
                       0.2524, 0.2687, 0.2867, 0.3050, 0.3228, 0.3358, 0.3504))
  # their Appendix D, Tables 5 (Kannisto) and 1 (Gompertz): a and b
  pub <- list(kannisto = list(male = c(3.39306e-05, 1.04287e-01), female = c(2.31232e-05, 1.05869e-01)),
              gompertz = list(male = c(2.17126e-04, 7.97629e-02), female = c(1.22095e-04, 8.41400e-02)))
  for (law in names(pub)) {
    for (sx in names(q)) {
      m <- c(q[[sx]] / (1 - q[[sx]] / 2), 0.6)
      tl <- attr(life_table(m, single_ages(80, 99), open = law, k = 19), "tail")
      # the published probabilities have four decimals
      expect_equal(tl$a, pub[[law]][[sx]][1], tolerance = 3e-3)
      expect_equal(tl$b, pub[[law]][[sx]][2], tolerance = 3e-4)
    }
  }
})
