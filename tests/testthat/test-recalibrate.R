test_that("the calibration estimator recovers a known spread factor", {
  set.seed(1)
  n <- 200
  G <- 400
  L <- matrix(stats::rnorm(G * n, 0, 0.05), G)
  pieces <- list(list(label = "a", L = L, m = rep(0, G), h = rep(1:4, G / 4), E = rep(1, G),
                      C = Matrix::Diagonal(G), Ecell = rep(1, G), y = stats::rnorm(G, 0, 0.1),
                      noise = rep(FALSE, n), size = rep(Inf, n)))
  P <- with_seed_test(2, mortensor:::stack_pieces(pieces, n))
  r <- mortensor:::calibration_fit(P, "constant", FALSE)
  expect_equal(r$table$scale_h1, 2, tolerance = 0.15)
  expect_lt(r$table$crps_after, r$table$crps_before)
  expect_gt(r$table$cov95_after, 0.9)
  ri <- mortensor:::calibration_fit(P, "constant", FALSE, target = "interval", level = 0.9)
  expect_equal(ri$table$scale_h1, 2, tolerance = 0.2)
  rb <- mortensor:::calibration_fit(P, "constant", TRUE)
  expect_equal(rb$table$scale_h1, 2, tolerance = 0.15)
  expect_lt(abs(rb$table$drift), 0.02)
})

test_that("the transformation keeps the median and the ranks of the draws", {
  sh <- mortensor:::calibration_shift(log(2), 1:3, "constant", FALSE)
  expect_equal(sh$s, 2)
  expect_equal(sh$b, 0)
  sh0 <- mortensor:::calibration_shift(0, 1:3, "constant", FALSE)
  expect_equal(sh0$s, 1)
  expect_equal(sh0$b, 0)
})

test_that("recalibrated() wraps any model and forecasts on the grid", {
  x <- mt_example(deaths = "single", years = 1995:2016)
  sp <- recalibrated(lee_carter(), years = 6, n_origins = 3, n_draws = 30)
  expect_s3_class(sp, "mt_spec_recalibrated")
  expect_equal(sp$min_periods, lee_carter()$min_periods + 6)
  # one spread factor for all horizons by default
  expect_identical(sp$params$scale, "constant")
  r <- expect_mt_model(sp, data = x, n_draws = 20)
  cal <- r$fit$calibration
  expect_true(all(c("group", "scale_h1", "drift", "crps_before", "crps_after", "cov95_after") %in% names(cal)))
  expect_equal(cal$drift, 0)
  expect_lte(cal$crps_after, cal$crps_before + 1e-12)
  expect_equal(r$fit$calibration_origins, c(2010, 2012, 2014))
  expect_identical(r$fc$model, "Lee-Carter, recalibrated")
  pr <- mt_forecast(r$fit, h = 3, n_draws = 20, seed = 1, predictive = TRUE)
  expect_true(all(is.finite(pr$rates) & pr$rates >= 0))
  # summary and parameters of the base model
  expect_output(print(summary(r$fit)), "calibration \\(all\\): spread factor")
  pa <- mt_parameters(r$fit)
  expect_true(all(c("a", "b", "k", "scale_h1") %in% pa$parameter))
  # separate parameters by sex, pooled ones as well
  fs <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 2, by = "sex", n_draws = 30), x)
  expect_setequal(fs$calibration$group, c("all", "female", "male"))
  expect_identical(fs$spec$label, "Lee-Carter, recalibrated by sex")
  expect_identical(recalibrated(rwd(), target = "interval", level = 0.9)$label,
                   "Random walk with drift, recalibrated (90% interval)")
  expect_error(recalibrated(lee_carter(), years = 1), class = "mt_error_argument")
  # the power law of the factor over the horizon was removed in 0.3.0
  expect_error(recalibrated(lee_carter(), scale = "horizon"), class = "simpleError")
  expect_error(mt_fit(recalibrated(lee_carter(), years = 6, by = "region", n_draws = 30), x),
               class = "mt_error_argument")
  expect_warning(mt_fit(recalibrated(lee_carter(), years = 20, n_draws = 20), x), class = "mt_warning_calibration")
  short <- mt_subset(x, periods = c(1995, 1999))
  expect_error(mt_fit(recalibrated(lee_carter(), years = 3), short), class = "mt_error_model_data")
})

test_that("recalibrated ensembles, reconciled models and backtests work", {
  x <- mt_example(deaths = "single", years = 1995:2014)
  ens <- recalibrated(ensemble(lc = lee_carter(), rw = rwd()), years = 5, n_origins = 2, n_draws = 30)
  r <- expect_mt_model(ens, data = x, n_draws = 20)
  expect_equal(nrow(r$fit$calibration), 1)
  bt <- mt_backtest(recalibrated(rwd(), years = 5, n_origins = 2, n_draws = 30), x, origins = 2011, h = 3,
                    n_draws = 30)
  expect_equal(unique(bt$model), "Random walk with drift, recalibrated")
  xr <- mt_example_regions(years = 2000:2014, regions = 2)
  h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
  fit <- mt_fit(recalibrated(reconciled(lee_carter(), h, method = "wls"), years = 5, n_origins = 2, n_draws = 30),
                mt_aggregate_strata(xr, h))
  fc <- mt_forecast(fit, h = 2, n_draws = 20, seed = 1)
  expect_lt(coherence_gap(fc, h), 1e-8)
  pr <- mt_forecast(fit, h = 2, n_draws = 20, seed = 1, predictive = TRUE)
  expect_lt(coherence_gap(pr, h), 1e-8)
})

test_that("recalibrated() estimates factors by age band and on life expectancy", {
  x <- mt_example(deaths = "single", years = 1985:2015)
  fb <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 3, n_draws = 30, ages = c(60, 80)), x)
  expect_identical(fb$calibration$band, c("<60", "60-80", "80+"))
  fa <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 3, n_draws = 30), x)
  expect_equal(sum(fb$calibration$cells), fa$calibration$cells)
  fc <- mt_forecast(fb, h = 12, n_draws = 40, seed = 1)
  expect_true(all(is.finite(fc$rates)))
  # each band gets its own factor in the forecast
  base <- mt_forecast(mt_fit(lee_carter(), x), h = 12, n_draws = 40, seed = 1)
  ratio <- function(a) {
    w <- function(r) stats::IQR(log(r))
    w(fc$rates[a, 1, 1, ]) / w(base$rates[a, 1, 1, ])
  }
  expect_equal(ratio(10), fb$calibration$scale_h1[1], tolerance = 0.05)
  expect_equal(ratio(60), fb$calibration$scale_h1[3], tolerance = 0.05)
  expect_error(recalibrated(rwd(), ages = c(80, 60)), class = "mt_error_argument")
  expect_error(recalibrated(rwd(), ages = 60, quantity = "asr"), class = "mt_error_argument")
  fe <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 3, n_draws = 30, quantity = "life_expectancy"), x)
  # one life expectancy per stratum and test year of each origin (6 + 4 + 2 years)
  expect_equal(fe$calibration$scored, 2 * 12)
  expect_match(fe$spec$label, "on life expectancy")
})

test_that("recalibrated() reports intervals of the factors and damps the drift", {
  x <- mt_example(deaths = "single", years = 1985:2015)
  spec <- recalibrated(lee_carter(), years = 6, n_origins = 3, n_draws = 30, n_boot = 5, bias = TRUE)
  fit <- mt_fit(spec, x)
  cal <- fit$calibration
  expect_true(all(c("scale_h1_lo", "scale_h1_hi", "drift_lo", "drift_hi") %in% names(cal)))
  expect_lte(cal$scale_h1_lo, cal$scale_h1_hi)
  expect_output(print(summary(fit)), "spread factor .* \\(")
  sh <- mortensor:::calibration_shift(c(0, 0.02), 1:20, "constant", TRUE, h_max = 8)
  expect_equal(sh$b[8:20], rep(0.16, 13))
  expect_warning(mt_fit(recalibrated(rwd(), years = 3, n_origins = 1, n_draws = 20), x),
                 class = "mt_warning_calibration")
})

test_that("recalibrated() warns when the forecast is much longer than the calibrated horizons", {
  x <- mt_example(deaths = "single", years = 1995:2016)
  fit <- mt_fit(recalibrated(rwd(), years = 4, n_origins = 2, n_draws = 20), x)
  # origins 2012 and 2014: horizons 1 to 4 are calibrated
  expect_equal(fit$calibration_hmax, 4)
  expect_no_warning(mt_forecast(fit, h = 8, n_draws = 20, seed = 1), class = "mt_warning_recalibrate")
  expect_warning(mt_forecast(fit, h = 9, n_draws = 20, seed = 1), class = "mt_warning_recalibrate")
  expect_warning(mt_forecast(fit, to = 2026, n_draws = 20, seed = 1), class = "mt_warning_recalibrate")
  # beyond the calibrated horizons the factor and the drift correction of horizon 4 are kept
  sh <- mortensor:::calibration_shift(c(log(1.5), 0.02), 1:9, "constant", TRUE, h_max = 4)
  expect_equal(sh$s, 1.5)
  expect_equal(sh$b[4:9], rep(0.08, 6))
})

# Inner forecasts of `origins` with cells at horizons 1 to length(factors):
# draws N(0, sd) per cell, outcomes N(0, sd * factor of the horizon), no
# observation noise; large exposures keep the floor of half a death away
per_horizon_pieces <- function(factors, origins = 1:3, cells = 300, n = 200, sd = 0.05) {
  lapply(origins, function(o) {
    h <- rep(seq_along(factors), each = cells)
    k <- length(h)
    list(label = "a", origin = o, L = matrix(stats::rnorm(k * n, 0, sd), k), m = rep(0, k), h = h,
         E = rep(1e4, k), C = Matrix::Diagonal(k), Ecell = rep(1e4, k), y = stats::rnorm(k, 0, sd * factors[h]),
         noise = rep(FALSE, n), size = rep(Inf, n))
  })
}

calibration_params <- function(scale, h_reach = NULL, bias = FALSE) {
  list(scale = scale, bias = bias, target = "crps", level = 0.95, by = NULL, ages = NULL, seed = 2, n_boot = 0,
       h_reach = h_reach)
}

test_that("factors per horizon recover known spread factors of each horizon", {
  set.seed(1)
  true <- c(1, 1.5, 2, 3)
  pieces <- per_horizon_pieces(true)
  strata <- data.table::data.table(stratum = 1L, label = "a")
  est <- mortensor:::estimate_calibration(pieces, strata, calibration_params("per_horizon", 4))
  expect_equal(est$table$horizon, 1:4)
  expect_equal(est$table$origins, rep(3, 4))
  expect_equal(est$table$forecasts, rep(3, 4))
  expect_equal(est$table$scale, true, tolerance = 0.08)
  expect_true(all(est$table$crps_after <= est$table$crps_before + 1e-12))
  theta <- est$theta$all$`all ages`
  expect_equal(exp(theta[, 1]), est$table$scale)
  # the forecast transformation takes the factor of each horizon and none beyond
  sh <- mortensor:::calibration_shift(theta, c(1, 3, 4, 4), "per_horizon", FALSE)
  expect_equal(sh$s, est$table$scale[c(1, 3, 4, 4)])
  expect_equal(sh$b, 0)
  expect_error(mortensor:::calibration_shift(theta, 5, "per_horizon", FALSE), class = "mt_error_recalibrate")
  # a drift per horizon: the shift of horizon h is b_h * h
  shb <- mortensor:::calibration_shift(cbind(log(c(1, 2)), c(0.01, -0.02)), 1:2, "per_horizon", TRUE)
  expect_equal(shb$s, c(1, 2))
  expect_equal(shb$b, c(0.01, -0.04))
  # a horizon reached by one origin has no factor
  one <- per_horizon_pieces(1:2, origins = 1)
  expect_error(mortensor:::estimate_calibration(one, strata, calibration_params("per_horizon", 2)),
               class = "mt_error_calibration")
})

test_that("with the same factor at every horizon, factors per horizon agree with the constant factor", {
  set.seed(2)
  pieces <- per_horizon_pieces(rep(2, 4), cells = 200, n = 100)
  strata <- data.table::data.table(stratum = 1L, label = "a")
  const <- mortensor:::estimate_calibration(pieces, strata, calibration_params("constant"))
  per <- mortensor:::estimate_calibration(pieces, strata, calibration_params("per_horizon", 4))
  expect_equal(const$table$scale_h1, 2, tolerance = 0.05)
  expect_equal(per$table$scale, rep(const$table$scale_h1, 4), tolerance = 0.08)
  # one horizon: the factor per horizon is the constant factor, exactly
  single <- lapply(pieces, function(pc) {
    pc$h[] <- 1
    pc
  })
  c1 <- mortensor:::estimate_calibration(single, strata, calibration_params("constant"))
  p1 <- mortensor:::estimate_calibration(single, strata, calibration_params("per_horizon", 1))
  expect_identical(p1$table$scale, c1$table$scale_h1)
  expect_identical(p1$table$crps_after, c1$table$crps_after)
  # the same with a drift (Nelder-Mead on both parameters)
  small <- per_horizon_pieces(2, cells = 100, n = 50)
  cb <- mortensor:::estimate_calibration(small, strata, calibration_params("constant", bias = TRUE))
  pb <- mortensor:::estimate_calibration(small, strata, calibration_params("per_horizon", 1, bias = TRUE))
  expect_identical(c(pb$table$scale, pb$table$drift), c(cb$table$scale_h1, cb$table$drift))
})

test_that("recalibrated(scale = \"per_horizon\") needs two origins per horizon and applies the factor of each", {
  x <- mt_example(deaths = "single", years = 1995:2016)
  sp <- recalibrated(rwd(), years = 6, n_origins = 3, n_draws = 20, scale = "per_horizon")
  expect_identical(sp$label, "Random walk with drift, recalibrated per horizon")
  fit <- mt_fit(sp, x)
  # origins 2010, 2012 and 2014: horizons 1 to 4 are reached by two or more
  expect_equal(fit$calibration_origins, c(2010, 2012, 2014))
  expect_equal(fit$calibration$horizon, 1:4)
  expect_equal(fit$calibration$origins, c(3, 3, 2, 2))
  expect_equal(fit$calibration_hmax, 6 - floor(6 / 3))
  fc <- mt_forecast(fit, h = 4, n_draws = 30, seed = 1)
  base <- mt_forecast(mt_fit(rwd(), x), h = 4, n_draws = 30, seed = 1)
  ratio <- vapply(1:4, function(t) stats::IQR(log(fc$rates[5, t, 1, ])) / stats::IQR(log(base$rates[5, t, 1, ])), 0)
  expect_equal(ratio, fit$calibration$scale, tolerance = 1e-8)
  expect_error(mt_forecast(fit, h = 5, n_draws = 20, seed = 1), class = "mt_error_recalibrate")
  expect_error(mt_forecast(fit, to = 2021, n_draws = 20, seed = 1), class = "mt_error_recalibrate")
  expect_output(print(summary(fit)), "spread factors at horizons 1-4")
  pa <- mt_parameters(fit)
  expect_equal(pa[pa$parameter == "scale", ]$value, fit$calibration$scale)
  expect_error(recalibrated(rwd(), years = 3, n_origins = 1, scale = "per_horizon"), class = "mt_error_argument")
  # one origin leaves enough training years
  short <- mt_subset(x, periods = c(2012, 2016))
  expect_error(suppressWarnings(mt_fit(recalibrated(rwd(), years = 4, n_origins = 3, n_draws = 20,
                                                    scale = "per_horizon"), short)), class = "mt_error_calibration")
  # life expectancy as the scored quantity, per horizon: one value per
  # stratum and origin at each horizon
  fe <- mt_fit(recalibrated(rwd(), years = 6, n_origins = 3, n_draws = 20, scale = "per_horizon",
                            quantity = "life_expectancy"), x)
  expect_equal(fe$calibration$scored, 2 * c(3, 3, 2, 2))
})

# A toy model whose draws all equal the last observed rates; it records the
# last year of every fit
toy_last_model <- function(seen) {
  registerS3method("fit_model", "mt_spec_toylast", function(spec, D, E, ages, periods, ...) {
    seen$last <- c(seen$last, max(periods))
    list(last = log(pmax(D[, ncol(D)], 0.5) / E[, ncol(E)]))
  }, envir = asNamespace("mortensor"))
  registerS3method("simulate_model", "mt_spec_toylast", function(spec, fit, h, n_draws, ...) {
    array(exp(fit$last), c(length(fit$last), h, n_draws))
  }, envir = asNamespace("mortensor"))
  mortensor:::new_mt_spec("toylast", min_periods = 2, label = "last rates", smooth = FALSE)
}

test_that("the error method adds resampled error vectors of past forecasts to the median", {
  x <- mt_example(deaths = "single", years = 1995:2016)
  seen <- new.env()
  seen$last <- numeric()
  toy <- toy_last_model(seen)
  sp <- recalibrated(toy, years = 6, n_origins = 5, n_draws = 10, method = "errors")
  expect_identical(sp$label, "last rates, recalibrated by past errors")
  expect_false(sp$smooth)
  fit <- mt_fit(sp, x)
  # each inner fit saw the data up to its origin only
  expect_setequal(unique(seen$last), c(2010:2014, 2016))
  expect_equal(fit$calibration_origins, 2010:2014)
  expect_equal(fit$calibration$origins, c(5, 5, 4, 3, 2))
  expect_equal(fit$calibration_hmax, 5)
  # fewer than 20 origins: the bounds are the smallest and largest errors
  expect_equal(fit$calibration$tuning, rep(1, 5))
  # errors of the toy: observed log rate of year o + h minus that of year o
  y <- log(pmax(mt_array(x, "deaths"), 0.5) / mt_array(x, "exposure"))
  er <- fit$calibration_errors
  for (hv in 1:5) {
    for (o in er$origins[[hv]]) {
      expect_equal(unname(er$V[[hv]][, , as.character(o)]),
                   unname(y[, as.character(o + hv), ] - y[, as.character(o), ]), tolerance = 1e-10)
    }
  }
  # every draw at every horizon is the median plus one whole error vector (all
  # ages and strata of one origin)
  fc <- mt_forecast(fit, h = 5, n_draws = 40, seed = 3)
  m <- y[, "2016", ]
  used <- matrix(NA_integer_, 5, 40)
  for (t in 1:5) {
    for (j in 1:40) {
      d <- log(fc$rates[, t, , j]) - m
      hit <- which(vapply(seq_along(er$origins[[t]]), function(k) {
        max(abs(d - er$tuning[t] * er$V[[t]][, , k])) < 1e-10
      }, logical(1)))
      used[t, j] <- if (length(hit)) hit[1] else NA_integer_
    }
  }
  expect_false(anyNA(used))
  # resampled with replacement: several origins at each horizon
  expect_true(all(apply(used, 1, function(u) length(unique(u))) > 1))
  # the draws describe observed rates: predictive draws add nothing
  expect_identical(mt_forecast(fit, h = 5, n_draws = 40, seed = 3, predictive = TRUE)$rates, fc$rates)
  expect_false(mortensor:::fit_adds_noise(fit))
  expect_error(mt_forecast(fit, h = 6, n_draws = 10, seed = 1), class = "mt_error_recalibrate")
  expect_output(print(summary(fit)), "resampled errors of inner forecasts at horizons 1-5")
})

test_that("the inner forecasts of the error method use no data after their origin", {
  x <- mt_example(deaths = "single", years = 1995:2016)
  x2 <- x
  d <- data.table::copy(x$obs$deaths)
  d$value <- ifelse(d$period_lo > 2010, round(d$value * 1.2), d$value)
  x2$obs$deaths <- d
  sp <- recalibrated(lee_carter(), years = 6, n_origins = 5, n_draws = 20, method = "errors")
  f1 <- mt_fit(sp, x)
  f2 <- mt_fit(sp, x2)
  y1 <- log(pmax(mt_array(x, "deaths"), 0.5) / mt_array(x, "exposure"))
  y2 <- log(pmax(mt_array(x2, "deaths"), 0.5) / mt_array(x2, "exposure"))
  for (hv in 1:5) {
    yr <- as.character(2010 + hv)
    e1 <- f1$calibration_errors$V[[hv]][, , "2010"]
    e2 <- f2$calibration_errors$V[[hv]][, , "2010"]
    # the forecast of origin 2010 (observed minus error) is unchanged, and
    # its errors follow the changed observations
    expect_equal(y1[, yr, ] - e1, y2[, yr, ] - e2, tolerance = 1e-10)
    expect_false(isTRUE(all.equal(e1, e2)))
  }
  pa <- mt_parameters(f1)
  expect_equal(pa[pa$parameter == "tuning", ]$value, f1$calibration$tuning)
})

test_that("the tuning factor of the error method is the smallest with nominal coverage of the errors", {
  set.seed(5)
  # errors of one sign: the factor pulls the bounds towards the point forecast
  E <- matrix(exp(stats::rnorm(500)), 5)
  S <- t(apply(E, 1, sort))
  cover <- function(pi, lo, hi) mean(pi * lo <= E & E <= pi * hi)
  for (lv in c(0.8, 0.95)) {
    tu <- mortensor:::error_tuning(E, lv)
    M <- ncol(E)
    lo <- S[, ceiling((1 - lv) / 2 * M)]
    hi <- S[, ceiling((1 + lv) / 2 * M)]
    expect_gte(cover(tu$tuning, lo, hi), lv)
    expect_lt(cover(tu$tuning * (1 - 1e-9), lo, hi), lv)
    expect_lt(tu$tuning, 1)
    expect_equal(tu$coverage, cover(tu$tuning, lo, hi))
    # brute force over a grid of factors
    grid <- seq(0.5, 1, by = 1e-4)
    first <- grid[which(vapply(grid, cover, 0, lo = lo, hi = hi) >= lv)[1]]
    expect_equal(tu$tuning, first, tolerance = 1e-4)
  }
  # errors of both signs and fewer than 1 / (1 - level) origins: factor 1
  expect_equal(mortensor:::error_tuning(matrix(stats::rnorm(5 * 12), 5), 0.95)$tuning, 1)
})

test_that("recalibrated(method = \"errors\") refuses options of the spread method and coarse data", {
  expect_error(recalibrated(rwd(), method = "errors", bias = TRUE), class = "mt_error_argument")
  expect_error(recalibrated(rwd(), method = "errors", scale = "per_horizon"), class = "mt_error_argument")
  expect_error(recalibrated(rwd(), method = "errors", ages = 60), class = "mt_error_argument")
  expect_error(recalibrated(rwd(), method = "errors", by = "sex"), class = "mt_error_argument")
  x <- mt_example(deaths = "single", years = 1995:2016)
  expect_error(recalibrated(rwd(), years = 2, method = "errors"), class = "mt_error_argument")
  short <- mt_subset(x, periods = c(2012, 2016))
  expect_error(suppressWarnings(mt_fit(recalibrated(rwd(), years = 4, n_origins = 3, n_draws = 10,
                                                    method = "errors"), short)), class = "mt_error_calibration")
  # deaths in age groups coarser than the grid of the model
  xg <- mt_example(years = 2007:2016)
  expect_error(mt_fit(recalibrated(rwd(), years = 4, n_origins = 3, n_draws = 10, method = "errors"), xg,
                      prep = mt_prep(step_ungroup())), class = "mt_error_calibration")
})
