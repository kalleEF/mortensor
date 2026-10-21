# reconciled(method = "bottom_up") fits the base model to the bottom strata
# only and sums their draws

h3b <- function() mt_hierarchy(sex = "both", region = "all", quintile = "all")

test_that("bottom-up fits the base model to the bottom strata only and keeps the draws of the bottom strata", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  h <- h3b()
  xa <- mt_aggregate_strata(x, h)
  hs <- mortensor:::hierarchy_structure(xa$strata, h)
  b <- hs$bottom
  fit <- mt_fit(reconciled(lee_carter(), h, method = "bottom_up"), xa)
  expect_true(fit$reconcile$bottom_only)
  base <- mortensor:::unwrap_fit(fit, "mt_fit_reconciled")
  expect_identical(base$strata$label, xa$strata$label[b])
  expect_identical(base$diagnostics$stratum, xa$strata$label[b])
  expect_identical(fit$strata$label, xa$strata$label)
  expect_null(fit$fits)
  # the draws of the former approach: the base model fitted to all strata, its bottom strata summed
  fc <- mt_forecast(fit, h = 3, n_draws = 40, seed = 3)
  all_fc <- mt_forecast(mt_fit(lee_carter(), xa), h = 3, n_draws = 40, seed = 3)
  rb <- mt_reconcile(all_fc, h, method = "bottom_up")$rates
  expect_identical(fc$rates[, , b, ], rb[, , b, ])
  # the aggregates are sums by matrix products of other shapes, whose
  # rounding depends on the BLAS: equal up to rounding
  expect_equal(fc$rates, rb, tolerance = 1e-14)
  # rates times exposures divided by exposures: equal up to rounding
  expect_equal(fc$rates[, , b, ], all_fc$rates[, , b, ], tolerance = 1e-14)
  # the same from data with the bottom strata only
  expect_identical(mt_forecast(mt_fit(reconciled(lee_carter(), h, method = "bottom_up"), x), h = 3, n_draws = 40,
                               seed = 3)$rates, fc$rates)
  # the aggregates are the sums of the bottom strata, on the deaths scale
  E <- fit$exposure_last[, 1, ]
  D <- fc$rates * as.vector(fit$exposure_last[, rep(1, 3), ])
  for (i in which(!hs$is_bottom)[c(1, 10, 34)]) {
    parts <- apply(D[, , b, ], c(1, 2, 4), function(v) sum(v * hs$S[i, ]))
    expect_equal(D[, , i, ], parts, tolerance = 1e-12)
  }
  expect_true(is.na(fc$reconciliation$max_relative_gap_before))
  # the training data of the fit and of the draws cover all nodes
  ref <- mt_fit(lee_carter(), xa)
  expect_equal(fit$observed, ref$observed, tolerance = 1e-12)
  expect_equal(fit$train_D, ref$train_D, tolerance = 1e-12)
  expect_equal(fit$train_E, ref$train_E, tolerance = 1e-12)
  expect_equal(fc$observed, all_fc$observed, tolerance = 1e-12)
  fd <- fitted(fit)
  expect_equal(fd[, , b], fitted(ref)[, , b], tolerance = 1e-12)
  tot <- which(xa$strata$label == "both.all.all")
  expect_equal(fd[, , tot], apply(fd[, , b] * fit$train_E[, , b], 1:2, sum) / fit$train_E[, , tot], tolerance = 1e-12)
  expect_identical(fc$fitted, fd)
  expect_setequal(unique(mt_parameters(fit)$stratum), xa$strata$label[b])
  expect_s3_class(summary(fit), "summary.mt_fit")
  expect_match(mortensor:::model_assumptions(fit$spec), "bottom strata only")
  # predictive draws: Poisson deaths of the bottom strata given the same paths, summed to the aggregates
  pr <- mt_forecast(fit, h = 3, n_draws = 400, seed = 3, predictive = TRUE)
  u <- mt_forecast(fit, h = 3, n_draws = 400, seed = 3)
  Eb <- array(E[, b], c(dim(E)[1], 1, length(b)))[, rep(1, 3), , drop = FALSE]
  deaths <- pr$rates[, , b, ] * as.vector(Eb)
  expect_lt(max(abs(deaths - round(deaths))), 1e-6)
  mu <- u$rates[, , b, ] * as.vector(Eb)
  z <- (deaths - mu) / sqrt(mu)
  expect_lt(abs(mean(z)), 0.02)
  expect_lt(abs(stats::var(as.vector(z)) - 1), 0.05)
  expect_lt(coherence_gap(pr, h), 1e-8)
  # projection methods fit the base model to all strata
  fw <- mt_fit(reconciled(lee_carter(), h, method = "wls"), xa)
  expect_false(fw$reconcile$bottom_only)
  expect_identical(mortensor:::unwrap_fit(fw, "mt_fit_reconciled")$strata$label, xa$strata$label)
})

test_that("models whose fit of the bottom strata uses the aggregates are fitted to all strata", {
  expect_true(mortensor:::fits_bottom_alone(lee_carter()))
  expect_true(mortensor:::fits_bottom_alone(li_lee(common_by = "sex")))
  expect_true(mortensor:::fits_bottom_alone(cae()))
  expect_true(mortensor:::fits_bottom_alone(ensemble(lee_carter(), rwd())))
  expect_true(mortensor:::fits_bottom_alone(stratified("sex", female = li_lee(), .default = rwd())))
  expect_true(mortensor:::fits_bottom_alone(anchored(coupled(lee_carter(), cae()))))
  expect_true(mortensor:::fits_bottom_alone(tuned(function(w) rwd(drift_window = w), data.frame(w = c(5, 10)))))
  expect_false(mortensor:::fits_bottom_alone(tucker3()))
  expect_false(mortensor:::fits_bottom_alone(coda_multiway()))
  expect_false(mortensor:::fits_bottom_alone(saint(reference = list(region = "all"))))
  expect_false(mortensor:::fits_bottom_alone(coda_causes()))
  expect_false(mortensor:::fits_bottom_alone(ensemble(lee_carter(), tucker3())))
  expect_false(mortensor:::fits_bottom_alone(recalibrated(stratified("sex", female = tucker3(), .default = rwd()))))
  x <- mt_example_regions(years = 2004:2014, regions = 2)
  h <- mt_hierarchy(region = "all", quintile = "all")
  xa <- mt_aggregate_strata(x, h)
  fit <- mt_fit(reconciled(tucker3(), h, method = "bottom_up"), xa)
  expect_false(fit$reconcile$bottom_only)
  expect_identical(mortensor:::unwrap_fit(fit, "mt_fit_reconciled")$strata$label, xa$strata$label)
  fc <- mt_forecast(fit, h = 2, n_draws = 20, seed = 1)
  expect_identical(fc$rates, mt_reconcile(mt_forecast(mt_fit(tucker3(), xa), h = 2, n_draws = 20, seed = 1), h,
                                          method = "bottom_up")$rates)
})

test_that("joint models that pool the bottom strata keep their fits of the bottom strata", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  h <- h3b()
  xa <- mt_aggregate_strata(x, h)
  b <- mortensor:::hierarchy_structure(xa$strata, h)$bottom
  fa <- mt_fit(li_lee(common_by = "sex"), xa)
  fr <- mt_fit(reconciled(li_lee(common_by = "sex"), h, method = "bottom_up"), xa)
  fb <- mortensor:::unwrap_fit(fr, "mt_fit_reconciled")
  expect_identical(fb$strata$label, xa$strata$label[b])
  expect_equal(fitted(fb), fitted(fa)[, , b], tolerance = 1e-12)
  expect_equal(fb$fits$commons[1:2], fa$fits$commons[match(c("female", "male"), unique(xa$strata$sex))],
               tolerance = 1e-12)
  expect_equal(fb$fits$strata, fa$fits$strata[b], tolerance = 1e-12)
  expect_lt(coherence_gap(mt_forecast(fr, h = 3, n_draws = 30, seed = 1, predictive = TRUE), h), 1e-8)
  # cae(): the age effects, drifts and innovation covariance of the bottom strata
  x6 <- mt_subset(xa, ages = c(60, Inf))
  ca <- mt_fit(cae(), x6)
  cr <- mortensor:::unwrap_fit(mt_fit(reconciled(cae(), h, method = "bottom_up"), x6), "mt_fit_reconciled")
  nb <- length(b)
  expect_equal(cr$fits$beta, ca$fits$beta, tolerance = 1e-10)
  expect_equal(cr$fits$drift, ca$fits$drift[seq_len(nb)], tolerance = 1e-10)
  expect_equal(cr$fits$sigma, ca$fits$sigma[seq_len(nb), seq_len(nb)], tolerance = 1e-10)
  # coda(coherent = TRUE) and coherent_fdm(): the fitted rates of the bottom strata
  specs <- list(coda(coherent = TRUE, common_by = "sex"))
  if (requireNamespace("demography", quietly = TRUE)) specs <- c(specs, list(coherent_fdm(common_by = "sex")))
  for (sp in specs) {
    f_all <- suppressWarnings(mt_fit(sp, xa))
    f_bu <- suppressWarnings(mt_fit(reconciled(sp, h, method = "bottom_up"), xa))
    expect_true(f_bu$reconcile$bottom_only)
    expect_equal(fitted(mortensor:::unwrap_fit(f_bu, "mt_fit_reconciled")), fitted(f_all)[, , b], tolerance = 1e-10)
  }
})

test_that("wrappers inside bottom-up reconciliation are fitted to the bottom strata", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  h <- h3b()
  xa <- mt_aggregate_strata(x, h)
  base_of <- function(fit) mortensor:::unwrap_fit(fit, "mt_fit_reconciled")
  # stacking weights from inner backtests of the bottom strata
  ens <- ensemble(lc = lee_carter(), rw = rwd(), weights = "crps_stack", inner_origins = 2, inner_h = 2,
                  inner_n_draws = 30)
  fe <- mt_fit(reconciled(ens, h, method = "bottom_up"), xa)
  we <- mt_fit(ens, x)$weights
  expect_equal(base_of(fe)$weights, we)
  pw <- mt_parameters(fe)[parameter == "weight"]
  expect_equal(pw$value, unname(we))
  pe <- mt_forecast(fe, h = 2, n_draws = 30, seed = 1, predictive = TRUE)
  expect_lt(coherence_gap(pe, h), 1e-8)
  # a tuned window, chosen on the bottom strata
  tu <- tuned(function(w) rwd(drift_window = w), data.frame(w = c(5, 10)), inner_origins = 1, inner_h = 3,
              n_draws = 20)
  ft <- suppressWarnings(mt_fit(reconciled(tu, h, method = "bottom_up"), xa))
  expect_equal(base_of(ft)$tuning, suppressWarnings(mt_fit(tu, x))$tuning)
  expect_lt(coherence_gap(mt_forecast(ft, h = 2, n_draws = 20, seed = 1), h), 1e-8)
  # spread factors of recalibration from inner forecasts of the bottom strata
  rc <- recalibrated(lee_carter(), years = 4, n_origins = 2, n_draws = 20)
  fc <- suppressWarnings(mt_fit(reconciled(rc, h, method = "bottom_up"), xa))
  expect_equal(base_of(fc)$calibration, suppressWarnings(mt_fit(rc, x))$calibration)
  expect_lt(coherence_gap(mt_forecast(fc, h = 2, n_draws = 20, seed = 1, predictive = TRUE), h), 1e-8)
  # a level of a stratified model that is an aggregate of the hierarchy has no strata to fit
  st <- stratified("sex", both = rwd(), female = rwd(), .default = lee_carter())
  fs <- mt_fit(reconciled(st, h, method = "bottom_up"), xa)
  expect_setequal(names(base_of(fs)$groups), c("female", ".default"))
  ps <- mt_forecast(fs, h = 2, n_draws = 20, seed = 1, predictive = TRUE)
  us <- mt_forecast(fs, h = 2, n_draws = 20, seed = 1)
  fem <- which(xa$strata$sex == "female" & xa$strata$region != "all" & xa$strata$quintile != "all")
  mal <- which(xa$strata$sex == "male" & xa$strata$region != "all" & xa$strata$quintile != "all")
  expect_equal(ps$rates[, , fem, ], us$rates[, , fem, ], tolerance = 1e-12)
  expect_gt(mean(abs(ps$rates[, , mal, ] / us$rates[, , mal, ] - 1)), 0.01)
  expect_error(mt_fit(st, x), class = "mt_error_argument")
  # backtests score all strata of the data
  bt <- mt_backtest(reconciled(ens, h, method = "bottom_up"), xa, origins = 2011, h = 3, n_draws = 30)
  expect_setequal(unique(bt$label), xa$strata$label)
  expect_true(all(is.finite(bt$crps)))
})

test_that("wrappers around bottom-up reconciliation keep the draws coherent", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  h <- h3b()
  rb <- reconciled(stratified("sex", female = lee_carter(), .default = rwd()), h, method = "bottom_up")
  fr <- suppressWarnings(mt_fit(recalibrated(rb, years = 4, n_origins = 2, n_draws = 20, by = "sex"), x))
  for (pred in c(FALSE, TRUE)) {
    expect_lt(coherence_gap(mt_forecast(fr, h = 2, n_draws = 20, seed = 1, predictive = pred), h), 1e-8)
  }
  expect_true(any(mt_parameters(fr)$parameter == "scale_h1"))
  fa <- mt_fit(anchored(reconciled(lee_carter(), h, method = "bottom_up")), x)
  expect_equal(dim(fa$anchor$offset), c(length(fa$age), nrow(fa$strata)))
  pa <- mt_forecast(fa, h = 2, n_draws = 20, seed = 1, predictive = TRUE)
  expect_true(all(is.finite(pa$rates)))
  expect_lt(coherence_gap(pa, h), 1e-8)
  expect_equal(dim(fitted(fa)), dim(fa$train_D))
})
