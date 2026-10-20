# Hazard laws of Feehan (2018) for the open age interval, and the
# out-of-sample check of Thatcher, Kannisto and Vaupel (1998)

# Integrated hazard over [x, x + 1] from the closed forms of Feehan (2018,
# Appendix E), natural parameters p = (alpha, beta, gamma, delta)
feehan_H <- function(law, p, x) {
  a <- p[1]
  b <- p[2]
  g <- p[3]
  dl <- p[4]
  erf <- function(z) 2 * stats::pnorm(z * sqrt(2)) - 1
  switch(law,
         makeham = a / b * (exp(b * (x + 1)) - exp(b * x)) + g,                                  # eq. 18
         weibull = a / b * ((x + 1)^b - x^b),                                                    # eq. 53
         beard = a / (b * dl) * log((1 + dl * exp(b * (x + 1))) / (1 + dl * exp(b * x))),        # eq. 46
         logistic = g + a / (b * dl) * log((1 + dl * exp(b * (x + 1))) / (1 + dl * exp(b * x))), # eq. 35
         perks = g - g / b * log((dl * exp(b * (x + 1)) + 1) / (dl * exp(b * x) + 1)) +          # eq. 49
           a / (b * dl) * log((dl * exp(b * (x + 1)) + 1) / (dl * exp(b * x) + 1)),
         # eq. 58, for gamma < 0
         logquad = sqrt(pi) * exp(a - b^2 / (4 * g)) / (2 * sqrt(-g)) *
           (erf((b + 2 * g * x) / (2 * sqrt(-g))) - erf((b + 2 * g * (x + 1)) / (2 * sqrt(-g)))),
         lynch_brown = {                                                                          # eq. 61
           k0 <- g * (x - dl)
           k1 <- g * (x - dl + 1)
           a + b / (2 * g) * (2 * k1 * atan(k1) - 2 * k0 * atan(k0) + log(k0^2 + 1) - log(k1^2 + 1))
         })
}

# parameters of the order of those at ages 80-100
truth <- list(makeham = c(5e-5, 0.1, 0.02, NA), weibull = c(1e-14, 7.5, NA, NA), beard = c(4e-5, 0.11, NA, 1e-4),
              logquad = c(-17, 0.27, -0.001, NA), logistic = c(4e-5, 0.11, 0.01, 1e-4),
              perks = c(4e-5, 0.11, 0.01, 1e-4), lynch_brown = c(0.3, 0.2, 0.08, 95))

test_that("integrated hazards and their derivatives agree with Feehan (2018, Appendix E)", {
  x <- 80:110
  centre <- 90
  for (law in names(truth)) {
    th <- mortensor:::oldage_internal(law, truth[[law]], centre)
    expect_equal(mortensor:::oldage_natural(law, th, centre)[, 1], truth[[law]], tolerance = 1e-12,
                 ignore_attr = TRUE)
    hz <- mortensor:::oldage_H(law, th, x - centre, centre)
    expect_equal(as.vector(hz$H), feehan_H(law, truth[[law]], x), tolerance = 1e-12)
    # derivatives with respect to the internal parameters: central differences
    for (k in seq_len(nrow(th))) {
      e <- 1e-6 * max(1, abs(th[k, 1]))
      tp <- tm <- th
      tp[k, 1] <- th[k, 1] + e
      tm[k, 1] <- th[k, 1] - e
      Hk <- function(t) mortensor:::oldage_H(law, t, x - centre, centre)$H
      num <- (Hk(tp) - Hk(tm)) / (2 * e)
      expect_equal(as.vector(hz$dH[[k]]), as.vector(num), tolerance = 1e-6)
    }
  }
  # the log-quadratic hazard by quadrature also for gamma > 0
  th <- mortensor:::oldage_internal("logquad", c(-12, 0.05, 2e-4, NA), centre)
  ref <- vapply(x, function(z) {
    stats::integrate(function(t) exp(-12 + 0.05 * t + 2e-4 * t^2), z, z + 1, rel.tol = 1e-12)$value
  }, 0)
  expect_equal(as.vector(mortensor:::oldage_H("logquad", th, x - centre, centre)$H), ref, tolerance = 1e-12)
})

test_that("each law recovers its parameters by maximum likelihood", {
  x <- 80:104
  for (law in names(truth)) {
    H <- feehan_H(law, truth[[law]], x)
    N <- 1e6 * exp(-cumsum(c(0, H[-length(H)])))
    # expected deaths: the maximum-likelihood estimate is the truth
    f <- mortensor:::oldage_fit(law, x, matrix(N), matrix(N * -expm1(-H)))
    expect_true(f$converged)
    expect_equal(f$par[, 1], truth[[law]], tolerance = 1e-5, ignore_attr = TRUE)
  }
  # binomial deaths of a large cohort: estimates close to the truth, and the
  # maximum of the likelihood as found by optim()
  set.seed(11)
  for (law in c("makeham", "beard", "logquad", "weibull")) {
    H <- feehan_H(law, truth[[law]], x)
    N <- numeric(length(x))
    d <- numeric(length(x))
    alive <- 1e8
    for (i in seq_along(x)) {
      N[i] <- alive
      d[i] <- stats::rbinom(1, alive, -expm1(-H[i]))
      alive <- alive - d[i]
    }
    f <- mortensor:::oldage_fit(law, x, matrix(N), matrix(d))
    expect_true(f$converged)
    p <- truth[[law]]
    ok <- !is.na(p)
    expect_equal(f$par[ok, 1], p[ok], tolerance = 0.05, ignore_attr = TRUE)
    # the fitted hazard over the ages of the fit is within 1% of the truth
    expect_equal(as.vector(mortensor:::oldage_fit_H(f, x)), H, tolerance = 0.01)
    centre <- f$centre
    nll <- function(th) {
      -mortensor:::oldage_loglik(mortensor:::oldage_H(law, matrix(th), x - centre, centre)$H, matrix(N), matrix(d))
    }
    o <- stats::optim(f$theta[, 1] + 0.01, nll, method = "BFGS", control = list(reltol = 1e-15, maxit = 2000))
    expect_lte(nll(f$theta[, 1]), o$value + 1e-6)
  }
})

# Probabilities of dying at ages 80-108 of the cohorts born in 1871-1880 (11
# countries) and the numbers reaching age 80 (Thatcher, Kannisto and Vaupel
# 1998, Tables 6.1, 7.3 and 7.4)
thatcher_q <- list(
  male = c(0.1224, 0.1326, 0.1452, 0.1565, 0.1702, 0.1822, 0.1965, 0.2116, 0.2276, 0.2424, 0.2580, 0.2749, 0.2907,
           0.3073, 0.3253, 0.3434, 0.3571, 0.3691, 0.3817, 0.4288, 0.4568, 0.4622, 0.4547, 0.4941, 0.4506, 0.5199,
           0.4204, 0.5055, 0.5556),
  female = c(0.0987, 0.1074, 0.1180, 0.1291, 0.1409, 0.1520, 0.1653, 0.1787, 0.1926, 0.2070, 0.2225, 0.2354, 0.2524,
             0.2687, 0.2867, 0.3050, 0.3228, 0.3358, 0.3504, 0.3738, 0.3971, 0.407, 0.4191, 0.4332, 0.4392, 0.4878,
             0.4804, 0.4908, 0.5600))
thatcher_n80 <- c(male = 3199159, female = 4799863)

test_that("the closures reproduce the Weibull, quadratic and logistic parameters of Thatcher et al. (1998)", {
  # their Appendix D: Weibull a x^b (Table 2), quadratic exp(a + b x + c x^2)
  # fitted at 85-98 (Table 6) and the logistic model (Table 4), whose
  # constant c is about 1e-5 for the cohorts, so that the Beard law (c = 0)
  # gives their a, b and sigma (the gamma-Makeham form c + a e^(bx) / (1 +
  # sigma^2 a / b (e^(bx) - 1)); the female sigma is printed as 3.47119e-03
  # instead of 3.47119e-01)
  pub <- list(weibull = list(male = c(8.06558e-15, 6.93343), female = c(9.88729e-16, 7.36009)),
              logquad = list(male = c(-1.92367e+01, 3.25880e-01, -1.39941e-03),
                             female = c(-1.73745e+01, 2.75321e-01, -1.08971e-03)),
              logistic = list(male = c(2.99457e-05, 1.05917e-01, 3.33775e-01),
                              female = c(1.84435e-05, 1.08791e-01, 3.47119e-01)))
  for (sx in names(thatcher_q)) {
    q <- thatcher_q[[sx]][1:19]
    m <- c(q / (1 - q / 2), 0.6)
    tw <- attr(life_table(m, single_ages(80, 99), open = "weibull", k = 19), "tail")
    # the published probabilities have four decimals
    expect_equal(tw$alpha, pub$weibull[[sx]][1], tolerance = 5e-3)
    expect_equal(tw$beta - 1, pub$weibull[[sx]][2], tolerance = 3e-4)
    tq <- attr(life_table(m, single_ages(80, 99), open = "logquad", k = 14), "tail")
    expect_equal(tq$age_from, 85)
    # the three parameters are strongly correlated (their standard errors
    # exceed the values); the hazard they give is close
    expect_equal(c(tq$alpha, tq$beta, tq$gamma), pub$logquad[[sx]], tolerance = 5e-3)
    hz <- function(p, x) exp(p[1] + p[2] * x + p[3] * x^2)
    expect_equal(hz(c(tq$alpha, tq$beta, tq$gamma), 85:98), hz(pub$logquad[[sx]], 85:98), tolerance = 2e-3)
    tb <- attr(life_table(m, single_ages(80, 99), open = "beard", k = 19), "tail")
    thatcher <- c(a = tb$alpha / (1 + tb$delta), b = tb$beta, sigma = sqrt(tb$delta * tb$beta / tb$alpha))
    expect_equal(thatcher, pub$logistic[[sx]], tolerance = 6e-3, ignore_attr = TRUE)
  }
})

# The pooled cohort data of Thatcher et al. as an mt_data object: N reaching
# each age from the numbers at 80 and the probabilities of dying, d = N q,
# person-years N - d / 2
thatcher_data <- function(sx) {
  q <- thatcher_q[[sx]]
  N <- thatcher_n80[[sx]] * cumprod(c(1, 1 - q))[seq_along(q)]
  d <- N * q
  suppressWarnings(mt_data(deaths = data.frame(age = 80:108, year = 1875, deaths = d),
                           exposure = data.frame(age = 80:108, year = 1875, exposure = N - d / 2)))
}

test_that("mt_open_check() reproduces the out-of-sample comparisons of Thatcher et al. (1998, Tables 6.3-6.8)", {
  x <- thatcher_data("male")
  # the reconstructed numbers reaching ages 90 and 100 match their Table 7.3
  # (441 690 and 7 408)
  chk <- mt_open_check(x, cut = 99, laws = c("gompertz", "weibull", "kannisto", "beard"), by_age = TRUE)
  expect_equal(chk[law == "kannisto" & age == 100, N], 7408, tolerance = 1e-3)
  quad <- mt_open_check(x, cut = 99, from = 85, laws = "logquad", by_age = TRUE)
  chk <- rbind(chk, quad)[age >= 101]
  # chi-squared at ages 101-108, male cohorts (their Tables 6.4 and 6.8; the
  # logistic model with c of about 1e-5 is the Beard law)
  pub <- list(gompertz = c(36.187, 60.107, 25.205, 52.092, 15.370, 37.932, 13.304, 5.386),
              weibull = c(6.462, 22.483, 5.893, 25.417, 4.383, 21.352, 5.523, 1.701),
              kannisto = c(11.170, 0.178, 4.376, 1.607, 1.280, 3.981, 0.010, 0.218),
              beard = c(14.883, 0.744, 6.081, 0.950, 1.879, 3.285, 0.003, 0.339),
              logquad = c(17.373, 1.245, 7.226, 0.623, 2.354, 2.831, 0.029, 0.466))
  for (l in names(pub)) expect_equal(chk[law == l, chisq], pub[[l]], tolerance = 0.01)
  # minus log-likelihood (their Tables 6.3 and 6.6)
  expect_equal(chk[law == "kannisto", nll], c(2783.275, 1491.170, 820.016, 411.693, 227.040, 108.826, 63.076, 31.022),
               tolerance = 2e-3)
  expect_equal(chk[law == "gompertz", nll], c(2795.829, 1521.022, 830.314, 436.332, 233.868, 124.635, 69.263, 33.394),
               tolerance = 2e-3)
  # their conclusion: at ages 101-108 the logistic, Kannisto and quadratic
  # models fit better than Gompertz and Weibull, for both sexes
  for (sx in names(thatcher_q)) {
    # fitted at 80-98 (85-98), scored at 101-108 as in their Tables 6.6-6.9
    all <- rbind(mt_open_check(thatcher_data(sx), cut = 99, laws = c("gompertz", "weibull", "kannisto", "beard"),
                               by_age = TRUE),
                 mt_open_check(thatcher_data(sx), cut = 99, from = 85, laws = "logquad", by_age = TRUE))[age >= 101]
    tot <- all[, list(nll = sum(nll)), by = law]
    expect_lt(max(tot[law %in% c("kannisto", "beard", "logquad"), nll]),
              min(tot[law %in% c("gompertz", "weibull"), nll]))
  }
})

test_that("mt_open_check() ranks the true law of simulated deaths first, with errors at the noise level", {
  set.seed(3)
  x <- 80:109
  # a Kannisto hazard with old-age levels of recent decades; cohorts of 10
  # million at age 80, one replication per "year"
  H <- (log1p(2e-5 * exp(0.11 * (x + 1))) - log1p(2e-5 * exp(0.11 * x))) / 0.11
  q <- -expm1(-H)
  reps <- 30
  sims <- rbindlist(lapply(seq_len(reps), function(r) {
    N <- numeric(length(x))
    d <- numeric(length(x))
    alive <- 1e7
    for (i in seq_along(x)) {
      N[i] <- alive
      d[i] <- stats::rbinom(1, alive, q[i])
      alive <- alive - d[i]
    }
    data.table(age = x, year = 2000 + r, deaths = d, exposure = N - d / 2, N = N)
  }))
  sim <- mt_data(deaths = sims[, list(age, year, deaths)], exposure = sims[, list(age, year, exposure)])
  laws <- c("gompertz", "kannisto", "weibull", "beard", "logquad")
  res <- mt_open_check(sim, cut = 99, laws = laws)
  expect_equal(nrow(res), reps * length(laws))
  expect_true(all(res$converged))
  # the noise level: chi-squared per test age with the true probabilities
  test <- sims[age >= 99]
  noise <- test[, mean((deaths - N * q[age - 79])^2 / (N * q[age - 79] * (1 - q[age - 79])))]
  # the Kannisto law extrapolates best on average and ranks first more often
  # than any law that does not contain it (Beard contains it: delta = alpha)
  mean_nll <- res[, list(m = mean(nll)), by = law]
  expect_equal(mean_nll[which.min(m), law], "kannisto")
  first <- res[, list(f = mean(rank == 1)), by = law]
  expect_gt(first[law == "kannisto", f], max(first[law %in% c("gompertz", "weibull", "logquad"), f]))
  # its errors are at the noise level; the Gompertz and Weibull laws miss
  ks <- res[law == "kannisto", mean(chisq / n_ages)]
  expect_lt(abs(ks - noise), 0.3)
  expect_gt(res[law == "gompertz", mean(chisq / n_ages)], 20 * noise)
  expect_gt(res[law == "weibull", mean(chisq / n_ages)], 10 * noise)
  # AIC of the fits below the cut
  expect_equal(res[, aic], -2 * res[, loglik] + 2 * res[, n_par])
  expect_true(all(res[, delta_aic] >= 0))
  expect_equal(res[, min(delta_aic), by = period]$V1, rep(0, reps))
  # one row per test age with by_age = TRUE, summing to the totals
  ba <- mt_open_check(sim, cut = 99, laws = "kannisto", by_age = TRUE)
  expect_equal(ba[, sum(chisq), by = period]$V1, res[law == "kannisto", chisq], tolerance = 1e-10)
  expect_equal(ba[period == 2001, q_obs], test[year == 2001, deaths / N])
  expect_error(mt_open_check(sim, cut = 99, laws = "makeham_x"), class = "mt_error_argument")
  expect_error(mt_open_check(sim, cut = 84, laws = "perks"), class = "mt_error_argument")
})

test_that("life_table() closes the open interval with every law and selects one by AIC", {
  x <- 80:99
  p <- truth$beard
  m <- feehan_H("beard", p, x)
  age <- single_ages(80, 100)
  for (law in names(truth)) {
    lt <- life_table(c(m, 0.5), age, open = law, k = 15, omega = 105)
    tl <- attr(lt, "tail")
    expect_true(tl$converged)
    expect_equal(c(tl$age_from, tl$age_to), c(85, 100))
    # the open interval: years 100-104 with the fitted integrated hazard,
    # closed at 105, from the closed form with the reported parameters
    Hf <- feehan_H(law, c(tl$alpha, tl$beta, tl$gamma, tl$delta), 100:105)
    if (law == "logquad" && tl$gamma > 0) {
      Hf <- vapply(100:105, function(z) {
        stats::integrate(function(t) exp(tl$alpha + tl$beta * t + tl$gamma * t^2), z, z + 1, rel.tol = 1e-12)$value
      }, 0)
    }
    lo <- c(1, cumprod(exp(-Hf[1:5])))
    e_open <- sum(lo[1:5] * (1 - exp(-Hf[1:5])) / Hf[1:5]) + lo[6] / Hf[6]
    expect_equal(tail(lt$ex, 1), e_open, tolerance = 1e-8)
  }
  # schedules fitted column by column
  M <- cbind(c(m, 0.5), c(m * 1.1, 0.55))
  e <- life_expectancy(M, age, at = 100, open = "logquad", k = 15)
  expect_equal(e[2], tail(life_table(M[, 2], age, open = "logquad", k = 15)$ex, 1))
  # AIC needs observed numbers
  expect_error(life_table(c(m, 0.5), age, open = "aic"), class = "mt_error_lifetable")
  expect_error(life_expectancy(M, age, open = "unknown"), class = "mt_error_lifetable")
  E <- 2e5 * exp(-cumsum(c(0, m)))
  set.seed(5)
  D <- stats::rpois(21, c(m, 0.5) * E)
  la <- life_table(D / E, age, open = "aic", exposure = E, k = 15)
  ta <- attr(la, "tail")
  expect_setequal(ta$law, mortensor:::oldage_laws)
  expect_equal(ta$aic, -2 * ta$loglik + 2 * mortensor:::oldage_npar[ta$law], ignore_attr = TRUE)
  expect_equal(ta[selected == TRUE, law], ta[which.min(aic), law])
  best <- ta[selected == TRUE, law]
  expect_equal(la$ex, life_table(D / E, age, open = best, exposure = E, k = 15)$ex)
  # the log-likelihood by hand: binomial deaths among those reaching each age
  lb <- attr(life_table(D / E, age, open = "beard", exposure = E, k = 15), "tail")
  xs <- 85:99
  d <- D[6:20]
  N <- E[6:20] + d / 2
  Hb <- feehan_H("beard", c(lb$alpha, lb$beta, NA, lb$delta), xs)
  expect_equal(lb$loglik, sum(d * log(-expm1(-Hb)) - (N - d) * Hb), tolerance = 1e-10)
})
