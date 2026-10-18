# England & Wales men (HMD via StMoMo) in the 18 age groups of nordpred
ew_nordpred <- function() {
  x <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
  mt_regroup(x, age = c(seq(0, 80, 5), 85))
}

# Reference projections: the nordpred package 1.1 (Moller and Weedon-Fekjaer,
# GPL-2), run on the same data summed to the 5-year periods 1962-1966, ...,
# 2007-2011 (the year 1961 is left out), with nordpred(cases, pyr,
# startestage, startuseage, ...) and its defaults otherwise; rates (cases /
# person-years) of the five periods from 2012-2016. The package was run
# outside mortensor; only its results are stored here.
nordpred_refs <- list(
  power5 = list(spec = list(est_age = 20, use_age = 25), periods = 4, recent = TRUE,
                rows = c("25-29", "50-54", "70-74", "85-100"),
                rates = c(0.0005445302386, 0.0004552115784, 0.0004024828974, 0.0003780180962, 0.0003547576954,
                          0.0038321842080, 0.0038589651100, 0.0039776418660, 0.0040447382340, 0.0040569304810,
                          0.0230578400400, 0.0209277757700, 0.0208264255400, 0.0211386178600, 0.0219131426000,
                          0.1401941702000, 0.1211083416000, 0.1090513773000, 0.1037796141000, 0.1007035372000)),
  log = list(spec = list(link = "log", est_age = 20, use_age = 25), periods = 4, recent = TRUE,
             rows = c("25-29", "50-54", "85-100"),
             rates = c(0.0006206403878, 0.0005714105804, 0.0005407797744, 0.0005260857298, 0.0005117909511,
                       0.0038867259540, 0.0039256557330, 0.0039875046800, 0.0038901771440, 0.0036709215650,
                       0.1368959937000, 0.1146962918000, 0.0999230541500, 0.0926445328000, 0.0874278747100)),
  no_cut = list(spec = list(est_age = 20, use_age = 30, n_periods = 4, recent = FALSE, cut_trend = rep(0, 5)),
                periods = 4, recent = FALSE, rows = c("25-29", "50-54", "85-100"),
                rates = c(rep(0.0007322452974, 5),
                          0.0039233230410, 0.0038876073290, 0.0037900139730, 0.0035004192570, 0.0030440968430,
                          0.1418093837000, 0.1215588771000, 0.1063737676000, 0.0962799243000, 0.0887656833100)),
  three = list(spec = list(est_age = 30, use_age = 40, n_periods = 3), periods = 3, recent = FALSE,
               rows = c("30-34", "40-44", "55-59", "85-100"),
               rates = c(rep(0.0009417492323, 5),
                         0.0014902313710, 0.0012911139080, 0.0011706393780, 0.0011138731600, 0.0010593307640,
                         0.0059878624380, 0.0058761310370, 0.0062772578880, 0.0060584723560, 0.0058458301820,
                         0.1419386657000, 0.1228065212000, 0.1107985618000, 0.1055447970000, 0.1025694087000))
)

test_that("nordpred_apc() reproduces the projections of the nordpred package", {
  skip_if_not_installed("StMoMo")
  x <- ew_nordpred()
  for (nm in names(nordpred_refs)) {
    r <- nordpred_refs[[nm]]
    fit <- mt_fit(do.call(nordpred_apc, r$spec), x)
    f <- fit$fits[[1]]
    expect_identical(f$n_periods, as.integer(r$periods), info = nm)
    expect_identical(f$recent, r$recent, info = nm)
    proj <- mortensor:::nordpred_projection(f, 5)[r$rows, ]
    ref <- matrix(r$rates, length(r$rows), 5, byrow = TRUE)
    expect_lt(max(abs(proj / ref - 1)), 1e-7)
  }
})

test_that("nordpred_apc() gives annual draws around the projections", {
  x <- mt_regroup(mt_example(deaths = "single", years = 1993:2019), age = seq(30, 100, 5))
  res <- expect_mt_model(nordpred_apc(n_periods = 4), x, h = 12, n_draws = 50)
  f <- res$fit$fits[[1]]
  # four 5-year periods 2000-2019; earlier years have no fitted rates
  expect_identical(f$period_labels, c("2000-2004", "2005-2009", "2010-2014", "2015-2019"))
  expect_true(all(is.na(stats::fitted(res$fit)[, as.character(1993:1999), ])))
  expect_false(anyNA(stats::fitted(res$fit)[, as.character(2000:2019), ]))
  fc <- mt_forecast(res$fit, h = 13, n_draws = 2000, seed = 1)
  med <- apply(fc$rates, c(1, 2, 3), stats::median)
  # medians at the midpoints of the future periods (years 3 and 8) are the projections
  for (s in 1:2) {
    proj <- mortensor:::nordpred_projection(res$fit$fits[[s]], 2)
    expect_lt(max(abs(med[, c("2022", "2027"), s] / proj - 1)), 0.01)
  }
  # predictive draws add Poisson deaths
  fp <- mt_forecast(res$fit, h = 3, n_draws = 50, seed = 1, predictive = TRUE)
  expect_true(fp$predictive)
})

test_that("nordpred_apc() interpolates between period midpoints and projects young ages flat", {
  W <- mortensor:::nordpred_interp(10, 3)
  expect_equal(rowSums(W), rep(1, 10))
  expect_equal(W[3, ], c(0, 1, 0, 0))
  expect_equal(W[8, ], c(0, 0, 1, 0))
  expect_equal(W[1, ], c(0.4, 0.6, 0, 0))
  # a linear predictor along the midpoints gives the period value as mean of its five years
  eta <- c(-2, 3, 8, 13)
  expect_equal(mean((W %*% eta)[1:5]), 3)
  x <- mt_regroup(mt_example(deaths = "single", years = 1995:2019), age = seq(30, 100, 5))
  fit <- mt_fit(nordpred_apc(use_age = 40), x)
  fc <- mt_forecast(fit, h = 6, n_draws = 10, seed = 1)
  D <- mt_array(x, "deaths")
  E <- mt_array(x, "exposure")
  last2 <- list(as.character(2010:2014), as.character(2015:2019))
  young <- (rowSums(D[1:2, last2[[1]], 1]) / rowSums(E[1:2, last2[[1]], 1]) +
              rowSums(D[1:2, last2[[2]], 1]) / rowSums(E[1:2, last2[[2]], 1])) / 2
  expect_equal(fc$rates[1:2, , 1, ], array(young, c(2, 6, 10)), ignore_attr = TRUE)
  expect_gt(stats::sd(fc$rates[3, 6, 1, ]), 0)
})

test_that("nordpred_apc() checks its inputs and reports its parameters", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  expect_error(mt_fit(nordpred_apc(), x), class = "mt_error_model_fit")
  expect_error(nordpred_apc(n_periods = 2), class = "mt_error_argument")
  expect_error(nordpred_apc(cut_trend = c(0, 1.5)), class = "mt_error_argument")
  expect_error(nordpred_apc(est_age = 40, use_age = 30), class = "mt_error_argument")
  expect_error(nordpred_apc(recent = NA), class = "mt_error_argument")
  x5 <- mt_regroup(x, age = seq(30, 100, 5))
  expect_error(mt_fit(nordpred_apc(est_age = 33), x5), class = "mt_error_model_fit")
  fit <- mt_fit(nordpred_apc(link = "log"), x5)
  expect_identical(fit$fits[[1]]$n_periods, 3L)
  p <- mt_parameters(fit)
  expect_true(all(c("age", "drift", "period", "cohort", "drift_projected", "gof_p_value") %in% p$parameter))
  expect_equal(p[stratum == "female" & parameter == "period", value][c(1, 3)], c(0, 0))
  expect_match(summary(fit)$assumptions, "Nordpred")
  expect_match(nordpred_apc()$label, "Nordpred")
  expect_match(nordpred_apc(link = "log")$label, "log link")
})
