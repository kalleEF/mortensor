# p_splines() (Currie, Durban and Eilers 2004; MortalitySmooth) and
# cp_splines() (Camarda 2019)

# England and Wales men (HMD via StMoMo), ages 30-100, 1961-2011
ew_matrices <- function() {
  ew <- StMoMo::EWMaleData
  a <- as.character(30:100)
  y <- as.character(1961:2011)
  list(D = unname(ew$Dxt[a, y]), E = unname(ew$Ext[a, y]), age = mortensor:::mt_axis(30:100, 31:101, kind = "age"),
       years = 1961:2011)
}

# French women (HMD via demography), ages 0-99 and 100+, 1950-2001; deaths
# from the rates and exposures, rounded
fr_matrices <- function() {
  fr <- demography::fr.mort
  y <- as.character(1950:2001)
  E <- unname(fr$pop$female[1:101, y])
  list(D = round(unname(fr$rate$female[1:101, y]) * E), E = E,
       age = mortensor:::mt_axis(0:100, c(1:100, Inf), kind = "age"), years = 1950:2001)
}

ps_spec_params <- function(lambda, infant = FALSE) {
  p <- mortensor:::ps_settings(lambda, "bic", 2, 3, NULL, c(1e-4, 1e6))
  p$infant <- infant
  p
}

test_that("p_splines and cp_splines follow the model contract", {
  expect_mt_model(p_splines())
  expect_mt_model(p_splines(overdispersion = TRUE, criterion = "aic"), h = 2, n_draws = 10)
  expect_mt_model(cp_splines(n_boot = 4), h = 2, n_draws = 8)
  expect_s3_class(p_splines(), "mt_spec_p_splines")
  expect_true(p_splines()$smooth)
})

test_that("the P-spline fit reproduces MortalitySmooth 2.3.4 (stored values)", {
  skip_if_not_installed("StMoMo")
  m <- ew_matrices()
  f <- mortensor:::ps_fit_stratum(ps_spec_params(c(10, 1000)), m$D, m$E, m$age, m$years, margin = 0.01,
                                  overdispersion = FALSE)
  # values of MortalitySmooth's Mort2Dsmooth() for ages 30-100, years
  # 1961-2011, the log exposures as offset and the smoothing parameters
  # 10 and 1000 given (method 3)
  cells <- cbind(c(1, 11, 36, 56, 71), c(1, 20, 35, 45, 51))
  expect_equal(log(f$fitted[cells]), c(-6.8062299106, -6.2756732459, -3.8280188207, -2.0795798804, -0.7753790238),
               tolerance = 1e-8)
  expect_equal(f$ed, 101.0785092458, tolerance = 1e-8)
  expect_equal(f$dev, 11375.1379134501, tolerance = 1e-8)
  # predict(..., newdata = list(y = 1961:2031), se.fit = TRUE): forecasts and
  # standard errors 1, 5, 10, 15 and 20 years ahead
  fx <- mortensor:::ps2_extend_fit(f, 20, f$pord)
  fc <- cbind(c(1, 11, 36, 56, 71), 51 + c(1, 5, 10, 15, 20))
  expect_equal(fx$eta[fc], c(-7.1937404623, -6.5407476198, -4.6902026791, -2.7872635416, -1.0657242404),
               tolerance = 1e-8)
  Ainv <- chol2inv(fx$R)
  X <- t(vapply(seq_len(nrow(fc)), function(i) kronecker(fx$my$B[fc[i, 2], ], f$ma$B[fc[i, 1], ]),
                numeric(nrow(Ainv))))
  se <- sqrt(rowSums((X %*% Ainv) * X))
  expect_equal(se, c(0.0221246058, 0.0226716138, 0.0354103704, 0.0569766776, 0.1186954572), tolerance = 1e-7)
})

test_that("p_splines agrees with MortalitySmooth where it is installed", {
  skip_if_not_installed("StMoMo")
  skip_if_not_installed("MortalitySmooth")
  m <- ew_matrices()
  ms <- tryCatch(MortalitySmooth::Mort2Dsmooth(x = 30:100, y = m$years, Z = m$D, offset = log(m$E)),
                 error = function(e) NULL)
  skip_if(is.null(ms), "MortalitySmooth::Mort2Dsmooth() fails in this R version")
  f <- mortensor:::ps_fit_stratum(ps_spec_params(ms$lambdas), m$D, m$E, m$age, m$years, margin = 0.01,
                                  overdispersion = FALSE)
  expect_equal(log(f$fitted), unname(ms$logmortality), tolerance = 1e-7, ignore_attr = TRUE)
  expect_equal(f$ed, ms$df, tolerance = 1e-7)
  # selection: the criterion is at least as low as MortalitySmooth's grid search
  sel <- mortensor:::ps_fit_stratum(ps_spec_params(NULL), m$D, m$E, m$age, m$years, margin = 0.01,
                                    overdispersion = FALSE)
  expect_lt(sel$ic, ms$bic + 0.5)
  pr <- stats::predict(ms, newdata = list(x = 30:100, y = 1961:2021), se.fit = TRUE)
  fx <- mortensor:::ps2_extend_fit(f, 10, f$pord)
  expect_equal(fx$eta, unname(pr$fit), tolerance = 1e-7, ignore_attr = TRUE)
  Ainv <- chol2inv(fx$R)
  X <- kronecker(fx$my$B[52:61, ], f$ma$B)
  expect_equal(sqrt(rowSums((X %*% Ainv) * X)), as.vector(pr$se.fit[, 52:61]), tolerance = 1e-6)
})

test_that("cp_splines reproduces Camarda's (2019) code (stored values)", {
  skip_if_not_installed("demography")
  m <- fr_matrices()
  # PSinfant(Y, E, lambdas = c(1e-2, 1e3), WEI = 1) of the supplementary code
  f <- mortensor:::ps_fit_stratum(ps_spec_params(c(1e-2, 1e3), infant = TRUE), m$D, m$E, m$age, m$years,
                                  margin = 0, overdispersion = FALSE)
  expect_true(f$ma$infant)
  expect_equal(f$ndx, c(20, 10))
  expect_equal(log(f$fitted[cbind(c(1, 2, 31, 66, 101), c(1, 20, 30, 40, 52))]),
               c(-3.0661217920, -6.6182399143, -7.2737151131, -4.7422038749, -0.9060616794), tolerance = 1e-8)
  # deltasFUN(): bounds of the rates of ageing (95%) and of change (50%)
  b <- mortensor:::cps_bounds(f$alpha, f$ma, f$my, 0.95, 0.5)
  i <- c(1, 2, 31, 66, 101)
  expect_equal(b$ageing_upper[i], c(-2.2119772224, -0.4980721053, 0.0778945151, 0.1020875279, 0.1929936001),
               tolerance = 1e-8)
  expect_equal(b$ageing_lower[i], c(-2.5226828287, -1.0858055197, 0.0452780223, 0.0881035209, 0.0747932330),
               tolerance = 1e-8)
  expect_equal(b$change_upper[i], c(-0.0427033639, -0.0436896627, -0.0132574392, -0.0182827527, -0.0068991630),
               tolerance = 1e-8)
  expect_equal(b$change_lower[i], c(-0.0599639925, -0.0678984643, -0.0410216696, -0.0248072851, -0.0127689509),
               tolerance = 1e-8)
  # CPSfunction() with kappa = 1e4, forecast to 2021, run to convergence
  # (relative change 1e-12) with the basis of the fit extended over the
  # forecast years as by predict.Mort2Dsmooth() (the paper extends the basis
  # of the fit; the code builds a new basis over all years, which changes
  # these log rates by at most 0.001)
  f$kappa <- c(1e4, 1e4)
  fx <- mortensor:::cps_forecast_fit(f, 20, b)
  sel <- cbind(c(1, 2, 11, 31, 51, 66, 81, 91, 101), c(52, 53, 55, 60, 62, 64, 66, 70, 72))
  ref <- c(-5.6101082556, -7.9923261311, -9.3367183270, -7.8282506408, -6.1871773846, -5.1752271284,
           -3.4561095425, -2.1406961548, -1.1640295869)
  expect_equal(fx$eta[sel], ref, tolerance = 1e-7)
  own_basis <- c(-5.6101457018, -7.9927242448, -9.3364308590, -7.8286965886, -6.1874356471, -5.1753192885,
                 -3.4567206166, -2.1413374669, -1.1646388163)
  expect_lt(max(abs(fx$eta[sel] - own_basis)), 1e-3)
})

test_that("bootstrap deaths invert the deviance residual (Camarda 2019, eq. 22)", {
  y <- c(5, 50, 500, 5000, 1, 20, 3)
  r <- c(-1.5, 2.2, -0.3, 1.1, 0.7, -2.5, 0)
  mu <- mortensor:::boot_deaths(y, r)
  dev_res <- sign(y - mu) * sqrt(2 * (y * log(y / mu) - (y - mu)))
  expect_equal(dev_res, r, tolerance = 1e-8)
  # uniroot() as in the supplementary code
  expect_equal(mu[1:6], c(9.1424554, 36.0133889, 506.7382374, 4922.6210635, 0.4528867, 33.3537718),
               tolerance = 1e-7)
  # no deaths: r^2 / 2 for negative residuals, 0 for positive ones
  expect_equal(mortensor:::boot_deaths(c(0, 0), c(-2, 1)), c(2, 0))
})

test_that("missing cells get weight 0 and other models still refuse them", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  d <- x$obs$deaths
  drop <- d$stratum == 1 & d$age_lo %in% c(50, 51) & d$period_lo %in% c(2000, 2001)
  x2 <- x
  x2$obs$deaths <- d[!drop]
  expect_error(mt_fit(lee_carter(), x2), class = "mt_error_model_data")
  fit <- mt_fit(p_splines(), x2)
  expect_true(all(fit$diagnostics$status == "ok"))
  f <- fit$fits[[1]]
  expect_equal(sum(f$data$wt == 0), 4)
  expect_true(all(is.finite(f$fitted)))
  # the surface interpolates the missing cells: close to the fit with all cells
  full <- mt_fit(p_splines(lambda = f$lambda), x)$fits[[1]]
  expect_lt(max(abs(log(f$fitted[21:22, 6:7] / full$fitted[21:22, 6:7]))), 0.05)
  # a cell without exposure is left out as well
  wt <- mortensor:::ps_data(matrix(c(1, NA, 3, 4), 2), matrix(c(10, 10, 0, 10), 2))$wt
  expect_equal(wt, matrix(c(1, 0, 0, 1), 2))
  fc <- mt_forecast(fit, h = 3, n_draws = 20, seed = 1)
  expect_true(all(is.finite(fc$rates)))
})

test_that("p_splines handles age groups and an open interval", {
  x <- mt_regroup(mt_example(deaths = "single", years = 1995:2012), age = c(30, seq(35, 95, 5)))
  fit <- mt_fit(p_splines(), x)
  f <- fit$fits[[1]]
  # cell ages: midpoints, 97.5 for 95+
  expect_equal(f$xa[c(1, length(f$xa))], c(32.5, 97.5))
  # knots about five years of age apart: 66 years spanned
  expect_equal(f$ndx[1], floor(66 / 5))
  fc <- mt_forecast(fit, h = 5, n_draws = 30, seed = 2)
  expect_equal(dim(fc$rates), c(14, 5, 2, 30))
  expect_true(all(is.finite(fc$rates)))
})

test_that("p_splines draws follow the normal approximation and widen with the horizon", {
  x <- mt_example(deaths = "single", years = 1990:2012)
  fit <- mt_fit(p_splines(), x)
  f <- fit$fits[[1]]
  fx <- mortensor:::ps2_extend_fit(f, 15, f$pord)
  set.seed(1)
  eta <- log(simulate_model(fit$spec, f, 15, 4000))
  Ainv <- chol2inv(fx$R)
  se_at <- function(a, t) {
    v <- as.vector(kronecker(fx$my$B[length(f$xy) + t, ], f$ma$B[a, ]))
    sqrt(sum(v * as.vector(Ainv %*% v)))
  }
  for (cell in list(c(1, 1), c(30, 8), c(60, 15))) {
    expect_equal(stats::sd(eta[cell[1], cell[2], ]), se_at(cell[1], cell[2]), tolerance = 0.06)
    expect_equal(mean(eta[cell[1], cell[2], ]), fx$eta[cell[1], length(f$xy) + cell[2]], tolerance = 0.01)
  }
  s <- apply(eta, 1:2, stats::sd)
  expect_true(all(s[, 15] > s[, 1]))
  # with pord = 2 the forecast continues about linearly
  m <- fx$eta[, length(f$xy) + 1:15]
  expect_lt(max(abs(diff(t(m), differences = 2))), 0.01)
})

test_that("overdispersion widens the draws", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  f0 <- mt_fit(p_splines(lambda = c(10, 1000)), x)
  f1 <- mt_fit(p_splines(lambda = c(10, 1000), overdispersion = TRUE), x)
  psi2 <- f1$fits[[1]]$psi2
  expect_equal(psi2, f1$fits[[1]]$dispersion, tolerance = 1e-2)
  expect_gt(psi2, 0.5)
  p <- mt_parameters(f1)
  expect_true(all(c("lambda_age", "lambda_year", "ed", "bic", "dispersion") %in% p$parameter))
})

test_that("cp_splines constrains the forecast and bootstraps the fit", {
  x <- mt_example(deaths = "single", years = 1990:2012)
  set.seed(5)
  fit <- mt_fit(cp_splines(n_boot = 6), x)
  f <- fit$fits[[1]]
  expect_false(f$ma$infant)
  expect_equal(dim(f$boot$Y), c(71, 23, 6))
  expect_length(f$boot$bounds, 6)
  # the forecast keeps the rates of ageing and of change about within their
  # bounds (the weight of the asymmetric penalties is finite)
  outside <- function(fx) {
    A <- matrix(fx$alpha, ncol(f$ma$B))
    ch <- f$ma$B %*% A %*% t(fx$my$C[fx$forecast_cols, ])
    ag <- f$ma$C %*% A %*% t(fx$my$B[fx$forecast_cols, ])
    c(change = max(ch - f$bounds$change_upper, f$bounds$change_lower - ch, 0),
      ageing = max(ag - f$bounds$ageing_upper, f$bounds$ageing_lower - ag, 0))
  }
  cp <- outside(mortensor:::cps_forecast_fit(f, 20, f$bounds))
  expect_lt(cp[["change"]], 0.005)
  expect_lt(cp[["ageing"]], 0.01)
  # draws: one forecast per bootstrap sample, used in turn
  fc <- mt_forecast(fit, h = 5, n_draws = 12, seed = 1)
  expect_equal(fc$rates[, , 1, 1], fc$rates[, , 1, 7])
  expect_false(isTRUE(all.equal(fc$rates[, , 1, 1], fc$rates[, , 1, 2])))
  # reproducible with the seed of the fit
  set.seed(5)
  fit2 <- mt_fit(cp_splines(n_boot = 6), x)
  expect_identical(fit2$fits[[1]]$boot$Y, f$boot$Y)
  # without the bootstrap every draw is the forecast of the data
  f0 <- mt_fit(cp_splines(n_boot = 0), x)
  fc0 <- mt_forecast(f0, h = 3, n_draws = 4, seed = 1)
  expect_equal(fc0$rates[, , , 1], fc0$rates[, , , 4])
  p <- mt_parameters(fit)
  expect_true(all(c("ageing_lower", "ageing_upper", "change_lower", "change_upper") %in% p$parameter))
  expect_equal(sum(p$parameter == "change_upper"), 71 * 2)
})

test_that("cp_splines gives age 0 coefficients of its own", {
  skip_if_not_installed("demography")
  fr <- suppressWarnings(as_mt_data(mortensor:::pkg_data("fr.mort", "demography"), series = "female",
                                    years = 1980:2001))
  fr <- mt_regroup(fr, age = c(0, 1, seq(5, 90, 5)))
  fit <- suppressWarnings(mt_fit(cp_splines(n_boot = 2), fr))
  f <- fit$fits[[1]]
  expect_true(f$ma$infant)
  expect_equal(f$ma$B[1, ], c(1, rep(0, ncol(f$ma$B) - 1)))
  expect_equal(f$ma$B[-1, 1], rep(0, nrow(f$ma$B) - 1))
  # the rate of ageing at age 0 is the difference to the next group
  A <- matrix(f$alpha, ncol(f$ma$B))
  eta <- f$ma$B %*% A %*% t(f$my$B)
  expect_equal(as.vector(f$ma$C[1, ] %*% A %*% t(f$my$B)), (eta[2, ] - eta[1, ]) / (f$xa[2] - f$xa[1]))
  fc <- mt_forecast(fit, h = 5, n_draws = 4, seed = 1)
  expect_true(all(is.finite(fc$rates)))
})

test_that("arguments are checked", {
  expect_error(p_splines(lambda = -1), class = "mt_error_argument")
  expect_error(p_splines(ndx = 0), class = "mt_error_argument")
  expect_error(p_splines(pord = c(1, 2, 3)), class = "mt_error_argument")
  expect_error(p_splines(lambda_range = c(10, 1)), class = "mt_error_argument")
  expect_error(cp_splines(level_change = 1.5), class = "mt_error_argument")
  expect_error(cp_splines(kappa = -1), class = "mt_error_argument")
  expect_error(cp_splines(n_boot = 2.5), class = "mt_error_argument")
  expect_match(mortensor:::model_assumptions(p_splines()), "Currie, Durban and Eilers")
  expect_match(mortensor:::model_assumptions(cp_splines()), "Camarda 2019")
})
