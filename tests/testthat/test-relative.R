# Relative model of subpopulations and a reference (Villegas and Haberman 2014)

# Data simulated from the model. The reference has single ages `ref_ages` in
# years `ref_years` (eq. 16; cohort effects without linear trend with
# `cohort = TRUE`); the subpopulations g1-g3 have 5-year age groups from
# `sub_lower` in years `sub_years` and rates relative to the geometric mean of
# the reference rates over their age group (eqs. 17-18). `same_grid = TRUE`:
# the reference has the age groups and years of the subpopulations and is a
# stratum "ref" of the same data, and `twin = TRUE` adds a stratum "twin"
# with its deaths and exposures. noise = FALSE gives expected deaths;
# `by_sex = TRUE` gives the same data for women and men (stratum variable
# "sex").
rel_world <- function(ref_ages = 40:89, ref_years = 1986:2015, sub_lower = seq(50, 80, 5), sub_years = 2001:2015,
                      cohort = FALSE, noise = FALSE, same_grid = FALSE, twin = FALSE, ref_size = 1e7, sub_size = 2e6,
                      by_sex = FALSE, seed = 1) {
  set.seed(seed)
  if (same_grid) {
    ref_ages <- sub_lower
    ref_years <- sub_years
  }
  na_r <- length(ref_ages)
  nt_r <- length(ref_years)
  b_r <- 1 + 0.5 * cos((ref_ages - 40) / 15)
  b_r <- b_r / sum(b_r)
  tt <- ref_years - mean(ref_years)
  k_r <- -1.2 * tt + 0.01 * tt^2
  k_r <- k_r - mean(k_r)
  a_r <- -9.6 + 0.09 * (ref_ages - 40)
  cohorts <- (ref_years[1] - ref_ages[na_r]):(ref_years[nt_r] - ref_ages[1])
  g_r <- NULL
  eta_r <- a_r + outer(b_r, k_r)
  if (cohort) {
    # without mean and linear trend (eqs. 21 and 26)
    g_r <- 0.06 * sin((cohorts - 1900) / 7) + 0.03 * cos((cohorts - 1900) / 13)
    cc <- cohorts - mean(cohorts)
    g_r <- g_r - mean(g_r) - sum(cc * g_r) / sum(cc^2) * cc
    eta_r <- eta_r + g_r[col(eta_r) - row(eta_r) + na_r]
  }
  Er <- matrix(ref_size * exp(-0.02 * (ref_ages - 40)), na_r, nt_r)
  Dr <- if (noise) matrix(stats::rpois(length(Er), exp(eta_r) * Er), na_r) else exp(eta_r) * Er
  na <- length(sub_lower)
  nt <- length(sub_years)
  lab <- paste0(sub_lower, "-", sub_lower + 4)
  W <- diag(na)
  if (!same_grid) {
    W <- t(vapply(sub_lower, function(lo) as.numeric(ref_ages >= lo & ref_ages < lo + 5) / 5, numeric(na_r)))
  }
  lmb <- W %*% eta_r[, match(sub_years, ref_years)]
  beta <- c(0.22, 0.2, 0.17, 0.14, 0.12, 0.09, 0.06)[seq_len(na)]
  beta <- beta / sum(beta)
  ts <- sub_years - mean(sub_years)
  alpha <- cbind(0.5 - 0.04 * seq_len(na), 0.05 + 0 * seq_len(na), -0.4 + 0.03 * seq_len(na))
  kappa <- cbind(0.03 * ts + 0.05 * sin(ts), -0.01 * ts, -0.04 * ts + 0.03 * cos(ts))
  kappa <- sweep(kappa, 2, colMeans(kappa))
  sub <- data.table::rbindlist(lapply(1:3, function(g) {
    m <- exp(lmb + alpha[, g] + outer(beta, kappa[, g]))
    E <- matrix(sub_size * exp(-0.03 * (sub_lower - 50)), na, nt)
    D <- if (noise) matrix(stats::rpois(length(m), m * E), na) else m * E
    data.table::data.table(pop = paste0("g", g), year = rep(sub_years, each = na), age = rep(lab, nt),
                           deaths = as.vector(D), exposure = as.vector(E))
  }))
  ref <- data.table::data.table(pop = "ref", year = rep(ref_years, each = na_r),
                                age = rep(if (same_grid) lab else ref_ages, nt_r), deaths = as.vector(Dr),
                                exposure = as.vector(Er))
  truth <- list(a_r = a_r, b_r = b_r, k_r = k_r, g_r = g_r, alpha = alpha, beta = beta, kappa = kappa)
  if (same_grid) {
    d <- rbind(ref, sub, if (twin) transform(ref, pop = "twin"))
    return(list(data = suppressWarnings(mt_data(deaths = d, exposure = d, by = "pop")), truth = truth))
  }
  by <- "pop"
  if (by_sex) {
    two <- function(d) rbind(transform(d, sex = "female"), transform(d, sex = "male"))
    ref <- two(ref)
    sub <- two(sub)
    by <- c("pop", "sex")
  }
  list(ref = suppressWarnings(mt_data(deaths = ref, exposure = ref, by = by)),
       sub = suppressWarnings(mt_data(deaths = sub, exposure = sub, by = by)), truth = truth)
}

# Largest deviation of the estimates from the parameters of rel_world()
rel_errors <- function(m, truth) {
  c(a_r = max(abs(m$ref$alpha - truth$a_r)), b_r = max(abs(m$ref$beta - truth$b_r)),
    k_r = max(abs(m$ref$kappa - truth$k_r)),
    g_r = if (!is.null(truth$g_r)) max(abs(m$ref$gamma - truth$g_r)) else 0,
    alpha = max(abs(m$sub$alpha - truth$alpha)), beta = max(abs(m$sub$beta - truth$beta)),
    kappa = max(abs(m$sub$kappa - truth$kappa)))
}

test_that("relative_model() checks its arguments", {
  expect_error(relative_model(reference = 1), class = "mt_error_argument")
  expect_error(relative_model(reference = list("England")), class = "mt_error_argument")
  expect_error(relative_model(reference = "ref", by = 1), class = "mt_error_argument")
  expect_error(relative_model(reference = "ref", cohort = NA), class = "mt_error_argument")
  expect_error(relative_model(reference = "ref", deviation = "ar1"))
  expect_error(relative_model(reference = "ref", trend_from = 1990.5), class = "mt_error_argument")
  sp <- relative_model(reference = list(quintile = "England"), by = "sex", deviation = "rwd", trend_from = 1980)
  expect_s3_class(sp, "mt_spec_relative")
  expect_match(sp$label, "Lee-Carter reference, by sex, independent random walks, reference trend from 1980")
  expect_match(relative_model("ref", cohort = TRUE, deviation = "constant")$label,
               "age-period-cohort reference, constant deviations")
})

test_that("the effective numbers of parameters are those of Table 1", {
  # k = 7 age groups, n = 27 years (1981-2007), m = 5 quintiles; the
  # reference with ages 10-99 (k' = 90) and years 1961-2009 (n' = 49)
  np <- mortensor:::relative_npar(90, 49, FALSE, 7, 27, 5)
  expect_equal(unname(np), c(227, 171))
  # with cohort effects, all 138 cohorts counted (Table 1 reports 353)
  expect_equal(unname(mortensor:::relative_npar(90, 49, TRUE, 7, 27, 5)[1]), 363)
})

test_that("expected deaths give the parameters back (Lee-Carter reference, wider ages and longer period)", {
  w <- rel_world()
  fit <- mt_fit(relative_model(reference = w$ref), w$sub)
  m <- fit$fits$models[[1]]
  expect_true(m$converged)
  expect_lt(max(rel_errors(m, w$truth)), 1e-7)
  expect_equal(fitted(fit), fit$observed, tolerance = 1e-8, ignore_attr = TRUE)
  # the geometric mean of the reference rates over the single ages of each
  # group (eq. 18)
  expect_equal(m$W[1, ], c(rep(0, 10), rep(0.2, 5), rep(0, 35)))
  expect_equal(m$ref$years, 1986:2015)
  expect_identical(m$jj, 16:30)
})

test_that("expected deaths give the parameters back with cohort effects, and the constraints hold", {
  w <- rel_world(cohort = TRUE)
  fit <- mt_fit(relative_model(reference = w$ref, cohort = TRUE), w$sub)
  m <- fit$fits$models[[1]]
  expect_true(m$converged)
  err <- rel_errors(m, w$truth)
  expect_lt(max(err), 1e-5)
  g <- m$ref$gamma
  cc <- seq_along(g) - mean(seq_along(g))
  expect_equal(c(sum(m$ref$beta), sum(m$ref$kappa), sum(g), sum(cc * g), sum(m$sub$beta)), c(1, 0, 0, 0, 1),
               tolerance = 1e-10)
  expect_equal(colSums(m$sub$kappa), rep(0, 3), tolerance = 1e-10)
})

test_that("Poisson deaths give estimates close to the parameters", {
  for (coh in c(FALSE, TRUE)) {
    w <- rel_world(cohort = coh, noise = TRUE, ref_size = 1e10, sub_size = 1e9)
    m <- mt_fit(relative_model(reference = w$ref, cohort = coh), w$sub)$fits$models[[1]]
    err <- rel_errors(m, w$truth)
    # parameters span about 0.5 (alpha, kappa), 35 (kappa') and 0.15 (gamma')
    expect_lt(err[["alpha"]], 0.01)
    expect_lt(err[["kappa"]], 0.01)
    expect_lt(err[["beta"]], 0.01)
    expect_lt(err[["k_r"]], 0.05)
    expect_lt(err[["g_r"]], 0.01)
  }
})

test_that("a reference stratum of the data: recovery, its forecast, and a subpopulation equal to it", {
  w <- rel_world(same_grid = TRUE, twin = TRUE)
  fit <- mt_fit(relative_model(reference = "ref"), w$data)
  f <- fit$fits
  lab <- fit$strata$label
  expect_identical(f$sole, lab == "ref")
  m <- f$models[[1]]
  expect_identical(lab[m$targets], c("g1", "g2", "g3", "twin"))
  m3 <- m
  m3$sub <- list(alpha = m$sub$alpha[, 1:3], beta = m$sub$beta, kappa = m$sub$kappa[, 1:3])
  expect_lt(max(rel_errors(m3, w$truth)), 1e-7)
  expect_equal(fitted(fit), fit$observed, tolerance = 1e-8, ignore_attr = TRUE)
  # the subpopulation with the deaths and exposures of the reference
  # deviates by nothing
  expect_lt(max(abs(m$sub$alpha[, 4]), abs(m$sub$kappa[, 4])), 1e-8)
  # the reference stratum is forecast by the reference model; with constant
  # deviations every subpopulation keeps its last ratio to it
  fc <- mt_forecast(mt_fit(relative_model(reference = "ref", deviation = "constant"), w$data), h = 4, n_draws = 5,
                    seed = 1)
  lr <- log(fc$rates[, , "g1", ]) - log(fc$rates[, , "ref", ])
  last <- m$sub$alpha[, 1] + m$sub$beta * m$sub$kappa[nrow(m$sub$kappa), 1]
  expect_equal(lr, array(last, dim(lr)), tolerance = 1e-10, ignore_attr = TRUE)
  k <- fc$rates[, , "ref", ]
  expect_true(all(abs(log(k[, 2, ]) - log(k[, 1, ]) - (m$ref$beta %o% rep(1, 5)) *
                        (log(k[1, 2, ]) - log(k[1, 1, ]))[col(matrix(0, 7, 5))] / m$ref$beta[1]) < 1e-10))
})

test_that("without subpopulations, the reference fit is the Poisson Lee-Carter model", {
  w <- rel_world(same_grid = TRUE)
  x <- mt_subset(w$data, pop == "ref")
  expect_identical(mt_strata(x)$label, "ref")
  fit <- mt_fit(relative_model(reference = "ref"), x)
  lc <- mt_fit(lee_carter(), x)
  expect_equal(fitted(fit), fitted(lc), tolerance = 1e-7)
  m <- fit$fits$models[[1]]
  expect_equal(c(m$ref$drift, m$ref$sd), c(lc$fits[[1]]$drift, lc$fits[[1]]$sigma_k), tolerance = 1e-6)
  expect_length(m$targets, 0)
})

test_that("without subpopulations and with cohort effects, the reference fit is StMoMo's rh(approxConst = TRUE)", {
  skip_on_cran()
  skip_if_not_installed("StMoMo")
  ew <- mortensor:::pkg_data("EWMaleData", "StMoMo")
  ages <- 75:89
  yrs <- 1996:2011
  Dr <- ew$Dxt[as.character(ages), as.character(yrs)]
  Er <- ew$Ext[as.character(ages), as.character(yrs)]
  st <- StMoMo::fit(StMoMo::rh(approxConst = TRUE), Dxt = Dr, Ext = Er, ages = ages, years = yrs, verbose = FALSE,
                    tolerance = 1e-10, iterMax = 1e5)
  empty <- array(0, c(length(ages), length(yrs), 0))
  fit_ref <- function(...) {
    mortensor:::relative_fit_group(Dr, Er, empty, empty, diag(length(ages)), seq_along(yrs), TRUE, ...)
  }
  m <- fit_ref(tol = 1e-11)
  expect_true(m$converged)
  expect_equal(m$ref$fitted, unname(stats::fitted(st, type = "rates")), tolerance = 1e-5, ignore_attr = TRUE)
  # StMoMo's last transformation leaves sum(bx) slightly away from 1
  sb <- sum(st$bx)
  expect_equal(m$ref$beta, unname(st$bx[, 1]) / sb, tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(m$ref$kappa, unname(st$kt[1, ]) * sb, tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(m$ref$gamma, unname(st$gc), tolerance = 1e-3, ignore_attr = TRUE)
  expect_equal(m$loglik_ref, st$loglik, tolerance = 1e-8)
  expect_equal(m$npar_ref, st$npar)
  # the accelerated iterations reach the fixed point of the plain ones
  plain <- fit_ref(tol = 1e-11, accelerate = FALSE, maxit = 1e5)
  expect_true(plain$iterations > 10 * m$iterations)
  expect_equal(m$ref$eta, plain$ref$eta, tolerance = 1e-8)
})

test_that("estimates maximise the total log-likelihood (eq. 25): the subpopulations inform the reference", {
  w <- rel_world(noise = TRUE)
  m <- mt_fit(relative_model(reference = w$ref), w$sub)$fits$models[[1]]
  md_r <- mortensor:::model_inputs(w$ref, list(min_periods = 2, label = "r"))
  md <- mortensor:::model_inputs(w$sub, list(min_periods = 2, label = "s"))
  Dr <- md_r$D[, , 1]
  Er <- md_r$E[, , 1]
  loglik <- function(k_r, b_r, sub = TRUE) {
    er <- m$ref$alpha + outer(b_r, k_r)
    l <- sum(Dr * er - Er * exp(er))
    if (!sub) return(l)
    lmb <- m$W %*% er[, m$jj]
    for (s in 1:3) {
      es <- lmb + m$sub$alpha[, s] + outer(m$sub$beta, m$sub$kappa[, s])
      l <- l + sum(md$D[, , s] * es - md$E[, , s] * exp(es))
    }
    l
  }
  score <- function(f, v, i, eps = 1e-5) {
    up <- v
    dn <- v
    up[i] <- up[i] + eps
    dn[i] <- dn[i] - eps
    (f(up) - f(dn)) / (2 * eps)
  }
  k_r <- m$ref$kappa
  b_r <- m$ref$beta
  # scores of kappa' and beta' of the total log-likelihood vanish; the
  # information per parameter is of the order of the deaths (about 1e5)
  sk <- vapply(c(1, 10, 20, 30), function(i) score(function(v) loglik(v, b_r), k_r, i), 0)
  sb <- vapply(c(1, 20, 40), function(i) score(function(v) loglik(k_r, v), b_r, i), 0)
  expect_lt(max(abs(sk)), 0.05)
  expect_lt(max(abs(sb)), 50)
  # the reference part alone is not at its maximum there: the subpopulations
  # moved kappa'
  sk_ref <- vapply(16:30, function(i) score(function(v) loglik(v, b_r, sub = FALSE), k_r, i), 0)
  expect_gt(max(abs(sk_ref)), 1)
})

test_that("forecasts: random walks of the indexes, and the processes compared in the paper", {
  w <- rel_world(same_grid = TRUE, noise = TRUE)
  fits <- lapply(c(mrwd = "mrwd", rwd = "rwd", constant = "constant"), function(d) {
    mt_fit(relative_model(reference = "ref", deviation = d), w$data)
  })
  out <- expect_mt_model(relative_model(reference = "ref"), w$data, h = 3, n_draws = 20)
  fc <- lapply(fits, function(f) mt_forecast(f, h = 6, n_draws = 50, seed = 3))
  # the reference stratum has the same draws whatever the subpopulation process
  expect_identical(fc$mrwd$rates[, , "ref", ], fc$rwd$rates[, , "ref", ])
  expect_identical(fc$mrwd$rates[, , "ref", ], fc$constant$rates[, , "ref", ])
  # the first subpopulation has the same draws with correlated and with
  # independent innovations (Section 2.1); the others do not
  expect_equal(fc$mrwd$rates[, , "g1", ], fc$rwd$rates[, , "g1", ], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(fc$mrwd$rates[, , "g2", ], fc$rwd$rates[, , "g2", ])))
  # drifts and covariance of the year-to-year changes of the indexes (eq. 4)
  m <- fits$mrwd$fits$models[[1]]
  k <- m$sub$kappa
  nt <- nrow(k)
  expect_equal(m$sub$drift, (k[nt, ] - k[1, ]) / (nt - 1))
  expect_equal(m$sub$sigma, stats::cov(diff(k)))
  expect_equal(m$ref$drift, (m$ref$kappa[nt] - m$ref$kappa[1]) / (nt - 1))
  # log rates of a subpopulation minus those of the reference follow its
  # index: alpha + beta kappa in every draw
  lr <- log(fc$mrwd$rates[, , "g2", ]) - log(fc$mrwd$rates[, , "ref", ])
  kap <- (lr[1, , ] - m$sub$alpha[1, 2]) / m$sub$beta[1]
  expect_equal(lr[3, , ], m$sub$alpha[3, 2] + m$sub$beta[3] * kap, tolerance = 1e-10, ignore_attr = TRUE)
  # with the drift only (no innovations) the index moves linearly
  m0 <- m
  m0$sub$sigma[] <- 0
  m0$sub$root[] <- 0
  p0 <- mortensor:::relative_target_paths(m0$sub, 4, 2, "mrwd")
  expect_equal(p0[, , 1], m$sub$last + outer(m$sub$drift, 1:4))
})

test_that("trend_from estimates the random walk of the reference index from the later years", {
  w <- rel_world()
  fit <- mt_fit(relative_model(reference = w$ref, trend_from = 2001), w$sub)
  r <- fit$fits$models[[1]]$ref
  k <- r$kappa[r$years >= 2001]
  expect_equal(r$drift, (k[length(k)] - k[1]) / (length(k) - 1))
  expect_equal(r$sd, stats::sd(diff(k)))
  expect_equal(r$trend_years, c(2001, 2015))
  expect_error(mt_fit(relative_model(reference = w$ref, trend_from = 2014), w$sub), class = "mt_error_model_fit")
})

test_that("cohorts after the last fitted one are needed only beyond the age distance", {
  w <- rel_world(cohort = TRUE, noise = TRUE)
  fit <- mt_fit(relative_model(reference = w$ref, cohort = TRUE), w$sub)
  r <- fit$fits$models[[1]]$ref
  # youngest subpopulation age 50, youngest reference age 40
  expect_identical(r$first_used, 11L)
  expect_true(all(c("phi", "mu", "sd", "last") %in% names(r$cohort_arima)))
  fit0 <- fit
  fit0$fits$models[[1]]$ref$cohort_arima <- NULL
  fc <- mt_forecast(fit0, h = 10, n_draws = 10, seed = 1)
  expect_true(all(is.finite(fc$rates)))
  expect_error(mt_forecast(fit0, h = 11, n_draws = 10, seed = 1), class = "mt_error_model_fit")
  fc <- mt_forecast(fit, h = 25, n_draws = 10, seed = 1)
  expect_true(all(is.finite(fc$rates) & fc$rates > 0))
  # the first ten years do not depend on the extrapolated cohorts
  fc10 <- mt_forecast(fit0, h = 10, n_draws = 10, seed = 1)
  expect_equal(mt_forecast(fit, h = 10, n_draws = 10, seed = 1)$rates, fc10$rates)
})

test_that("groups of `by` have a reference each, from the data or from an mt_data object", {
  w <- rel_world()
  w2 <- rel_world(by_sex = TRUE)
  ref2 <- w2$ref
  sub2 <- w2$sub
  fit <- mt_fit(relative_model(reference = ref2, by = "sex"), sub2)
  expect_identical(fit$fits$groups, c("female", "male"))
  m <- fit$fits$models
  expect_equal(m[[1]]$sub$alpha, m[[2]]$sub$alpha)
  expect_lt(max(rel_errors(m[[2]], w$truth)), 1e-7)
  # errors inside the joint fit are reported as failures of the fit
  expect_error(mt_fit(relative_model(reference = w$ref, by = "sex"), sub2), "does not have",
               class = "mt_error_model_fit")
  expect_error(mt_fit(relative_model(reference = mt_subset(ref2, sex == "male"), by = "sex"), sub2),
               "no stratum for group", class = "mt_error_model_fit")
  expect_error(mt_fit(relative_model(reference = list(pop = "nobody")), w$sub), class = "mt_error_model_fit")
  err <- tryCatch(mortensor:::relative_external(list(reference = w$ref, by = "sex"), mt_strata(sub2),
                                                mt_grid(sub2)$age, 2001:2015), error = function(e) e)
  expect_s3_class(err, "mt_error_argument")
})

test_that("the reference must cover the ages and years of the data", {
  w <- rel_world()
  expect_error(mt_fit(relative_model(reference = mt_subset(w$ref, ages = c(55, 90))), w$sub),
               class = "mt_error_model_fit")
  expect_error(mt_fit(relative_model(reference = mt_subset(w$ref, periods = c(2005, 2015))), w$sub),
               class = "mt_error_model_fit")
  # cohort effects need single ages: here the reference has the 5-year groups
  # of the data
  g <- rel_world(same_grid = TRUE)
  expect_error(mt_fit(relative_model(reference = "ref", cohort = TRUE), g$data), class = "mt_error_model_fit")
  err <- tryCatch(mortensor:::relative_check_single(c(50, 55), c(55, 60)), error = function(e) e)
  expect_s3_class(err, "mt_error_model_data")
  err <- tryCatch(mortensor:::relative_age_weights(mt_grid(w$sub)$age, mt_grid(mt_subset(w$ref, ages = c(52, 90)))$age),
                  error = function(e) e)
  expect_s3_class(err, "mt_error_model_data")
})

test_that("later years of the reference are not used, so backtests use past data only", {
  w <- rel_world()
  x <- mt_subset(w$sub, periods = c(2001, 2010))
  fit <- mt_fit(relative_model(reference = w$ref), x)
  expect_equal(range(fit$fits$models[[1]]$ref$years), c(1986, 2010))
  bt <- mt_backtest(relative_model(reference = w$ref), w$sub, origins = 2012, h = 3, n_draws = 20)
  expect_true(nrow(bt) > 0)
  expect_false(mortensor:::fits_bottom_alone(relative_model(reference = "ref")))
  expect_true(mortensor:::fits_bottom_alone(relative_model(reference = w$ref)))
})

test_that("parameters, BIC and the description of the model", {
  w <- rel_world(cohort = TRUE, noise = TRUE)
  fit <- mt_fit(relative_model(reference = w$ref, cohort = TRUE), w$sub)
  p <- mt_parameters(fit)
  expect_setequal(unique(p$stratum), c("reference all", "subpopulations all", "g1", "g2", "g3"))
  ref <- p[stratum == "reference all"]
  expect_equal(nrow(ref[parameter == "gamma"]), 50 + 30 - 1)
  expect_identical(ref[parameter == "gamma"]$index[1], "1897")
  stat <- function(part, nm) p[stratum == part & parameter == nm]$value
  for (part in c("reference all", "subpopulations all")) {
    expect_equal(stat(part, "bic"), stat(part, "loglik") - 0.5 * stat(part, "npar") * log(stat(part, "nobs")))
  }
  expect_equal(stat("subpopulations all", "npar"), 7 * 3 + 7 + 15 * 3 - 1 - 3)
  expect_equal(p[stratum == "g1" & parameter == "cor" & index == "g1"]$value, 1)
  # log-likelihood with the terms log d! of eq. (24)
  m <- fit$fits$models[[1]]
  mu <- fitted(fit) * fit$train_E
  expect_equal(m$loglik_sub, sum(stats::dpois(round(fit$train_D), mu, log = TRUE)), tolerance = 1e-8)
  s <- summary(fit)
  expect_match(s$assumptions, "Villegas and Haberman 2014")
  expect_match(mortensor:::model_assumptions(relative_model("ref", deviation = "constant")), "level differences")
})
