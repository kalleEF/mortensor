# Regression tests for fixes in the model wrappers and the hierarchy code

h3 <- function() mt_hierarchy(sex = "both", region = "all", quintile = "all")
bottom_of <- function(st) which(st$sex != "both" & st$region != "all" & st$quintile != "all")

test_that("reconciled() reconciles the draws of wrapped base models", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  xa <- mt_aggregate_strata(x, h3())
  tot <- which(xa$strata$label == "both.all.all")
  own <- mt_forecast(mt_fit(anchored(lee_carter(), decay = Inf), xa), h = 2, n_draws = 20, seed = 1)
  fit <- mt_fit(reconciled(anchored(lee_carter(), decay = Inf), h3(), method = "proportional", fixed = "total",
                           dependence = "none"), xa)
  expect_null(mortensor:::unwrap_fit(fit, "mt_fit_reconciled")$reconcile)
  fc <- mt_forecast(fit, h = 2, n_draws = 20, seed = 1)
  # the anchored draws are not summed before they are reconciled, and the fixed total is the anchored one
  expect_gt(fc$reconciliation$max_relative_gap_before, 1e-3)
  expect_equal(fc$rates[, , tot, ], own$rates[, , tot, ], tolerance = 1e-8)
  expect_lt(coherence_gap(fc, h3()), 1e-8)
  rc <- reconciled(recalibrated(lee_carter(), years = 4, n_origins = 2, n_draws = 20), h3(), dependence = "none")
  fr <- mt_forecast(mt_fit(rc, xa), h = 2, n_draws = 20, seed = 1)
  expect_gt(fr$reconciliation$max_relative_gap_before, 1e-3)
  # outer wrappers still see the hierarchy of an inner reconciled(); nested wrappers keep their own fields
  fa <- mt_fit(anchored(reconciled(lee_carter(), h3())), xa)
  expect_false(is.null(mortensor:::unwrap_fit(fa, "mt_fit_anchored")$reconcile))
  xs <- mt_example(deaths = "single", years = 1995:2012)
  fn <- mt_fit(anchored(anchored(lee_carter(), years = 1, shrink = FALSE), years = 3), xs)
  inner <- mortensor:::unwrap_fit(fn, "mt_fit_anchored")
  expect_s3_class(inner, "mt_fit_anchored")
  expect_false(isTRUE(all.equal(inner$anchor$offset, fn$anchor$offset)))
})

test_that("predictive noise of reconciled models is independent across cells", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  xa <- mt_aggregate_strata(x, h3())
  fit <- mt_fit(reconciled(lee_carter(), h3(), method = "bottom_up", dependence = "empirical"), xa)
  p <- mt_forecast(fit, h = 2, n_draws = 300, seed = 1, predictive = TRUE)
  u <- mt_forecast(fit, h = 2, n_draws = 300, seed = 1)
  # the same seed gives the same paths, so the difference is the observation noise
  E <- fit$exposure_last[, 1, 1]
  noise <- (p$rates[, , 1, ] - u$rates[, , 1, ]) * E
  expect_lt(abs(stats::cor(noise[3, 1, ], noise[6, 1, ])), 0.3)
  expect_lt(abs(stats::cor(noise[3, 1, ], noise[3, 2, ])), 0.3)
})

test_that("reconciled() adds the observation noise of each component or group", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  xa <- mt_aggregate_strata(x, h3())
  both <- function(spec) {
    f <- mt_fit(spec, xa)
    list(u = mt_forecast(f, h = 2, n_draws = 30, seed = 1),
         p = mt_forecast(f, h = 2, n_draws = 30, seed = 1, predictive = TRUE))
  }
  expect_false(ensemble(rwd(), rwd())$smooth)
  expect_true(ensemble(lee_carter(), lee_carter(method = "svd"))$smooth)
  expect_false(tuned(function(w) rwd(drift_window = w), data.frame(w = c(5, 10)))$smooth)
  # random walks forecast observed rates: no noise is added, also in an ensemble
  r <- suppressWarnings(both(reconciled(ensemble(rwd(), rwd(drift_window = 10)), h3(), method = "bottom_up")))
  expect_equal(r$p$rates, r$u$rates, tolerance = 1e-8)
  # groups of a stratified model: Lee-Carter strata get noise, random-walk strata none
  r <- suppressWarnings(both(reconciled(stratified("sex", female = lee_carter(), male = rwd(),
                                                   .default = lee_carter()), h3(), method = "bottom_up")))
  b <- bottom_of(xa$strata)
  fem <- intersect(b, which(xa$strata$sex == "female"))
  mal <- intersect(b, which(xa$strata$sex == "male"))
  expect_equal(r$p$rates[, , mal, ], r$u$rates[, , mal, ], tolerance = 1e-8)
  expect_gt(mean(abs(r$p$rates[, , fem, ] / r$u$rates[, , fem, ] - 1)), 0.01)
  expect_lt(coherence_gap(r$p, h3()), 1e-8)
  # the spec chosen by tuned() decides, also below other wrappers
  xs <- mt_example(deaths = "single", years = 1995:2012)
  hs <- mt_hierarchy(sex = "both")
  tu <- tuned(function(w) rwd(drift_window = w), data.frame(w = c(5, 10)), inner_origins = 1, inner_h = 3,
              n_draws = 20)
  ft <- suppressWarnings(mt_fit(anchored(reconciled(tu, hs, dependence = "none")), mt_aggregate_strata(xs, hs)))
  expect_equal(mt_forecast(ft, h = 2, n_draws = 20, seed = 1, predictive = TRUE)$rates,
               mt_forecast(ft, h = 2, n_draws = 20, seed = 1)$rates, tolerance = 1e-8)
})

test_that("reconciled forecasts take exposures of the bottom strata only", {
  x <- mt_example_regions(years = 2005:2016, regions = 2)
  fit <- mt_fit(reconciled(lee_carter(), h3()), x)
  Eb <- mt_array(mt_window(x, from = 2016), "exposure")
  expect_equal(dim(Eb)[3], 20L)
  f0 <- mt_forecast(fit, h = 2, n_draws = 20, seed = 1, predictive = TRUE)
  f1 <- mt_forecast(fit, h = 2, n_draws = 20, seed = 1, predictive = TRUE, exposure = Eb)
  expect_equal(f1$rates, f0$rates)
  fc <- mt_forecast(mt_fit(lee_carter(), mt_aggregate_strata(x, h3())), h = 2, n_draws = 20, seed = 1)
  expect_equal(mt_reconcile(fc, h3(), exposure = Eb, seed = 1)$rates, mt_reconcile(fc, h3(), seed = 1)$rates)
  expect_equal(mt_coherence(fc, h3(), exposure = Eb), mt_coherence(fc, h3()))
  # wrappers around reconciled() fitted to bottom strata
  rf <- mt_fit(recalibrated(reconciled(lee_carter(), h3()), years = 4, n_origins = 2, n_draws = 20), x)
  expect_equal(rf$calibration_origins, c(2012, 2014))
  expect_lt(coherence_gap(mt_forecast(rf, h = 2, n_draws = 20, seed = 1, exposure = Eb), h3()), 1e-8)
  # components of an ensemble must forecast the same strata
  expect_error(mt_fit(ensemble(rc = reconciled(lee_carter(), h3()), lc = lee_carter()), x),
               class = "mt_error_model_fit")
})

test_that("the one-variable form of mt_reconcile() adds deaths of sexes and rates of causes", {
  x <- mt_example(deaths = "single", years = 2000:2012)
  hs <- mt_hierarchy(sex = "both")
  fc <- mt_forecast(mt_fit(lee_carter(), mt_aggregate_strata(x, hs)), h = 2, n_draws = 20, seed = 1)
  expect_equal(mt_reconcile(fc, var = "sex", total = "both", seed = 1)$rates,
               mt_reconcile(fc, hs, seed = 1)$rates)
  xc <- mt_example_causes(years = 2005:2015)
  fcc <- mt_forecast(mt_fit(lee_carter(), xc), h = 2, n_draws = 20, seed = 1)
  expect_equal(mt_reconcile(fcc, var = "cause", total = "all", seed = 1)$rates,
               mt_reconcile(fcc, mt_hierarchy(cause = "all", additive = "rates"), seed = 1)$rates)
})

test_that("stratified() keeps a level named like its default group", {
  xc <- mt_example_causes(years = 2005:2015)
  fit <- mt_fit(stratified("cause", all = lee_carter(), other = rwd(), .default = lee_carter()), xc)
  g <- match("other", names(fit$groups))
  expect_setequal(mt_strata(xc)$label[fit$groups[[g]]], c("female.other", "male.other"))
  expect_identical(fit$group_fits[[g]]$spec$label, rwd()$label)
  fc <- mt_forecast(fit, h = 2, n_draws = 10, seed = 1)
  expect_false(anyNA(fc$rates))
})

test_that("reconciled ensembles keep the component of each reordered path", {
  x <- mt_example_regions(years = 2005:2016, regions = 2)
  xa <- mt_aggregate_strata(x, h3())
  ens <- ensemble(lc = lee_carter(), rw = rwd())
  pe <- mt_forecast(mt_fit(ens, xa), h = 2, n_draws = 40, seed = 1)
  pr <- suppressWarnings(mt_forecast(mt_fit(reconciled(ens, h3(), method = "bottom_up", dependence = "empirical"), xa),
                                     h = 2, n_draws = 40, seed = 1))
  expect_equal(dim(pr$component), c(nrow(xa$strata), 40L))
  for (s in utils::head(bottom_of(xa$strata), 3)) {
    A <- matrix(pe$rates[, , s, ], ncol = 40)
    src <- vapply(seq_len(40), function(j) which.min(colSums(abs(A - as.vector(pr$rates[, , s, j])))), 1L)
    expect_equal(pr$component[s, ], pe$component[src])
  }
  expect_equal(dim(mt_reconcile(pe, h3(), seed = 1, dependence = "empirical")$component), c(nrow(xa$strata), 40L))
  # without reordering (the default), each draw keeps one component for all strata
  expect_equal(mt_reconcile(pe, h3(), seed = 1)$component, pe$component)
  fa <- suppressWarnings(mt_forecast(mt_fit(anchored(reconciled(ens, h3(), dependence = "empirical")), xa), h = 2,
                                     n_draws = 20, seed = 1, predictive = TRUE))
  expect_true(all(is.finite(fa$rates)))
  expect_lt(coherence_gap(fa, h3()), 1e-8)
})

test_that("ensembles do not take the wrapper fields of their first component", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  fit <- suppressWarnings(mt_fit(ensemble(rc = recalibrated(lee_carter(), years = 4, n_origins = 2, n_draws = 20),
                                          an = anchored(lee_carter()), rw = rwd()), x))
  expect_null(fit$calibration)
  expect_null(fit$anchor)
  expect_null(fit$inner_specs)
  expect_false(any(grepl("calibration", utils::capture.output(print(summary(fit))))))
  # components reconciled over the same hierarchy: the ensemble keeps it, and wrappers keep the draws coherent
  xc <- mt_example_causes(years = 2008:2019)
  h <- mt_hierarchy(cause = "all", additive = "rates")
  rc2 <- ensemble(a = reconciled(lee_carter(), h, method = "bottom_up"), b = reconciled(rwd(), h, method = "bottom_up"))
  fa <- mt_forecast(mt_fit(anchored(rc2), xc), h = 2, n_draws = 20, seed = 1)
  st <- mt_strata(xc)
  for (sx in c("female", "male")) {
    i <- st$stratum[st$sex == sx & st$cause == "all"]
    j <- st$stratum[st$sex == sx & st$cause != "all"]
    expect_lt(max(abs(apply(fa$rates[, , j, , drop = FALSE], c(1, 2, 4), sum) / fa$rates[, , i, ] - 1)), 1e-8)
  }
  expect_null(mt_fit(ensemble(a = reconciled(lee_carter(), h, method = "bottom_up"), b = lee_carter()), xc)$reconcile)
})

test_that("stacking and CRPS weights compare the components on common cells", {
  cells <- function(origin, crps, err) {
    data.table::data.table(origin = origin, stratum = 1L, age_lo = 1:5, age_hi = 2:6, period = origin + 1,
                           target = "rate", observed = (1:5) / 10, crps = crps, median = (1:5) / 10 + err)
  }
  # component a also forecast an origin at which b failed, with large errors there
  a <- rbind(cells(2010, 10, 5), cells(2011, 0.1, 0))
  b <- cells(2011, 0.2, 0.1)
  w <- mortensor:::ensemble_inner_weights(list(weights = "crps"), list(a, b, NULL))
  expect_equal(w, c(10, 5, 0))
  ws <- mortensor:::ensemble_inner_weights(list(weights = "stack"), list(a, b, NULL))
  expect_gt(ws[1], 0.95)
  expect_equal(ws[3], 0)
  # the draws of the inner backtests are an argument
  expect_identical(ensemble(lee_carter(), rwd())$params$inner_n_draws, 100)
  expect_error(ensemble(lee_carter(), rwd(), inner_n_draws = 2), class = "mt_error_argument")
  seen <- NULL
  local_mocked_bindings(mt_backtest = function(spec, data, origins, h, prep, n_draws, ...) {
    seen <<- c(seen, n_draws)
    NULL
  })
  x <- mt_example(deaths = "single", years = 2000:2012)
  mortensor:::ensemble_inner_backtests(ensemble(lee_carter(), rwd(), inner_n_draws = 37)$params, x, NULL)
  expect_equal(seen, c(37, 37))
})

test_that("stratified() and tuned(by = ) accept ensembles as group models", {
  x <- mt_example(deaths = "single", years = 2000:2012)
  fit <- mt_fit(stratified("sex", female = ensemble(lee_carter(), rwd()), male = lee_carter()), x)
  expect_identical(unique(fit$diagnostics$stratum), c("female", "male"))
  expect_true("component" %in% names(fit$diagnostics))
  fc <- mt_forecast(fit, h = 2, n_draws = 10, seed = 1, predictive = TRUE)
  expect_true(all(is.finite(fc$rates)))
})

test_that("coherence is checked with one forecast year or one draw", {
  x <- mt_example_regions(years = 2005:2016, regions = 2)
  xa <- mt_aggregate_strata(x, h3())
  n_agg <- nrow(xa$strata) - length(bottom_of(xa$strata))
  f <- mt_fit(lee_carter(), xa)
  co <- mt_coherence(mt_forecast(f, h = 1, n_draws = 10, seed = 1), h3())
  expect_equal(nrow(co), n_agg)
  co1 <- mt_coherence(mt_forecast(f, h = 2, n_draws = 1, seed = 1), h3())
  expect_equal(nrow(co1), 2 * n_agg)
  expect_true(all(is.finite(co1$q0.5)))
  skip_if_not_installed("ggplot2")
  expect_s3_class(mt_plot_coherence(mt_forecast(f, h = 1, n_draws = 10, seed = 1), h3()), "ggplot")
})

test_that("recalibrated() works without a seed", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  fit <- mt_fit(recalibrated(rwd(), years = 4, n_origins = 2, n_draws = 20, seed = NULL), x)
  expect_equal(nrow(fit$calibration), 1)
})

test_that("mt_parameters() finds the model below bootstrapped(tuned())", {
  x <- mt_example(deaths = "single", years = 2000:2012)
  tu <- tuned(function(method) lee_carter(method = method), data.frame(method = "poisson"), inner_origins = 1,
              inner_h = 2, n_draws = 10)
  fit <- suppressWarnings(mt_fit(bootstrapped(tu, n_boot = 1), x))
  expect_s3_class(fit, "mt_fit_bootstrap")
  expect_true(all(c("a", "b", "k") %in% mt_parameters(fit)$parameter))
})

test_that("group members of a hierarchy are checked against the data", {
  x <- mt_example_regions(years = 2010:2015, regions = 3)
  typo <- mt_hierarchy(region = list(total = "all", groups = list(North = c("R1", "R2"), South = c("R3", "R4x"))))
  expect_warning(xa <- mt_aggregate_strata(x, typo), class = "mt_warning_hierarchy")
  expect_true("female.South.Q1" %in% mt_strata(xa)$label)
  empty <- mt_hierarchy(region = list(total = "all", groups = list(North = c("R1", "R2"), East = "X1")))
  expect_error(mt_aggregate_strata(x, empty), class = "mt_error_hierarchy")
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 1, n_draws = 5, seed = 1)
  expect_error(mt_aggregate_draws(fc, hierarchy = empty), class = "mt_error_hierarchy")
  expect_error(mt_hierarchy(region = list(total = "all", groups = list("R1"))), class = "mt_error_hierarchy")
  expect_error(mt_hierarchy(region = list(total = "all", groups = list(N = c("R1", "all")))),
               class = "mt_error_hierarchy")
})
