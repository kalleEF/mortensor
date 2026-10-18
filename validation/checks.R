# Validation checks. Each check returns a data frame with columns `metric`
# and `value`; thresholds live in thresholds.csv (check, metric). Checks in
# suite "quick" run on every push, "full" adds the slower ones.

checks <- list()

checks$ungroup_ew_1d <- list(
  suite = "quick",
  description = paste("England & Wales men (HMD via StMoMo) 2000-2011, deaths grouped to 5-year groups",
                      "and 90-100, ungrouped as rates (1D PCLM); median absolute log error against",
                      "observed single-year rates (Poisson noise floor about 0.012 and 0.021)."),
  run = function() {
    r <- mortensor:::benchmark_ungroup_ew(years = 2000:2011, dims = 1)
    data.frame(metric = c("male_30_89", "male_90_100", "seconds"), value = c(r$male_30_89, r$male_90_100, r$seconds))
  })

checks$ungroup_sim <- list(
  suite = "quick",
  description = paste("Simulated data with known truth (mt_example, 2 strata, 30 years): deaths in 5-year groups",
                      "and 90+, exposures single years and 100+; median absolute log error of ungrouped rates."),
  run = function() {
    x <- mt_example()
    u <- mt_ungroup(x)
    le <- abs(log((mt_array(u, "deaths") / mt_array(u, "exposure")) / attr(x, "truth")))
    ages <- mt_grid(x)$age$lower
    data.frame(metric = c("mae_log_30_89", "mae_log_90_plus"),
               value = c(stats::median(le[ages < 90, , ]), stats::median(le[ages >= 90, , ])))
  })

coverage_check <- function(spec, n_rep) {
  r <- mortensor:::benchmark_coverage(spec, n_rep = n_rep)
  s <- attr(r, "summary")
  data.frame(metric = c("cov95_mean", "cov95_min", "cov80_mean", "cov50_mean", "crps"),
             value = c(s$cov95_mean, s$cov95_min, s$cov80_mean, s$cov50_mean, s$crps))
}

checks$coverage_lee_carter <- list(
  suite = "quick",
  description = paste("Predictive interval calibration of lee_carter() on simulated Lee-Carter-type data",
                      "(10 strata, 5-year age groups, 30 training years, horizons 1-10, 6 replicates)."),
  run = function() coverage_check(lee_carter(), 6))

checks$coverage_tucker3 <- list(
  suite = "quick",
  description = "As coverage_lee_carter for tucker3() (log rates, ranks 3/3/2, random walk factors).",
  run = function() coverage_check(tucker3(), 6))

checks$coverage_rwd <- list(
  suite = "quick",
  description = paste("As coverage_lee_carter for rwd(). The random walk treats the Poisson noise of the",
                      "observed rates as persistent, so it over-covers; the threshold only guards against",
                      "under-coverage."),
  run = function() coverage_check(rwd(), 6))

checks$coverage_local_level <- list(
  suite = "quick",
  description = paste("As coverage_lee_carter for arima_series(order = c(0, 1, 1), drift = TRUE), the reduced",
                      "form of a random walk with drift observed with noise (Harvey 1989). It separates the",
                      "Poisson noise of the observed rates from the trend."),
  run = function() coverage_check(arima_series(order = c(0, 1, 1), drift = TRUE), 6))

checks$coverage_ensemble <- list(
  suite = "full",
  description = "As coverage_lee_carter for ensemble(lee_carter(), tucker3(), rwd()).",
  run = function() coverage_check(ensemble(lee_carter(), tucker3(), rwd()), 10))

checks$backtest_ew <- list(
  suite = "quick",
  description = paste("England & Wales men (HMD via StMoMo), ages 30-100, 1961-2011: rolling-origin backtest",
                      "(origins 1996-2001, h = 10) of lee_carter() and rwd(); CRPS of log rates, coverage",
                      "and skill of Lee-Carter relative to the random walk (CRPS ratio)."),
  run = function() {
    ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
    ew <- mt_set_grid(ew, age = single_ages(30, 100, open = FALSE))
    ew <- mt_regroup(ew, age = single_ages(30, 100, open = FALSE))
    o <- 1996:2001
    lc <- mt_backtest(lee_carter(), ew, origins = o, h = 10, n_draws = 200)
    rw <- mt_backtest(rwd(), ew, origins = o, h = 10, n_draws = 200)
    data.frame(metric = c("lc_crps", "lc_cov95", "lc_cov80", "rwd_crps", "rwd_cov95", "crps_ratio_lc_rwd"),
               value = c(mean(lc$crps), mean(lc$cov95), mean(lc$cov80), mean(rw$crps), mean(rw$cov95),
                         mean(lc$crps) / mean(rw$crps)))
  })

checks$calibration_ew <- list(
  suite = "quick",
  description = paste("recalibrated() on England & Wales men (as backtest_ew: origins 1996-2001, h = 10): Lee-Carter",
                      "with its spread recalibrated on the last 10 training years at each origin, against",
                      "Lee-Carter alone; CRPS, coverage and the ratio of the CRPS."),
  run = function() {
    ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
    ew <- mt_set_grid(ew, age = single_ages(30, 100, open = FALSE))
    ew <- mt_regroup(ew, age = single_ages(30, 100, open = FALSE))
    o <- 1996:2001
    lc <- mt_backtest(lee_carter(), ew, origins = o, h = 10, n_draws = 200)
    cl <- mt_backtest(recalibrated(lee_carter()), ew, origins = o, h = 10, n_draws = 200)
    data.frame(metric = c("lc_crps", "lc_cov95", "cal_crps", "cal_cov80", "cal_cov95", "crps_ratio_cal_lc"),
               value = c(mean(lc$crps), mean(lc$cov95), mean(cl$crps), mean(cl$cov80), mean(cl$cov95),
                         mean(cl$crps) / mean(lc$crps)))
  })

checks$tuning_ew <- list(
  suite = "full",
  description = paste("Tuning and recalibration on England & Wales men (as backtest_ew: origins 1996-2001, h = 10):",
                      "hyndman_ullah() with 6 components against the number of components and the knots of the",
                      "smoother chosen at each origin by a hold-out of the last 10 training years (RMSE, tuned()),",
                      "and the tuned model calibrated (recalibrated(tuned())); CRPS, RMSE of the median and 95%",
                      "coverage of the log rates."),
  run = function() {
    ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
    ew <- mt_set_grid(ew, age = single_ages(30, 100, open = FALSE))
    ew <- mt_regroup(ew, age = single_ages(30, 100, open = FALSE))
    o <- 1996:2001
    g <- expand.grid(k = c(10, 30), order = c(1, 2, 3, 6))[, c("order", "k")]
    tu <- tuned(function(order, k) hyndman_ullah(order = order, smooth_args = list(k = k)), g, inner_origins = 1,
                inner_h = 10, metric = "rmse", label = "tuned (hold-out)")
    b0 <- mt_backtest(hyndman_ullah(), ew, origins = o, h = 10, n_draws = 200)
    b1 <- mt_backtest(tu, ew, origins = o, h = 10, n_draws = 200)
    b2 <- mt_backtest(recalibrated(tu), ew, origins = o, h = 10, n_draws = 200)
    rmse <- function(b) sqrt(mean(b$ae^2))
    data.frame(metric = c("fdm_crps", "fdm_rmse", "fdm_cov95", "tuned_crps", "tuned_rmse", "tuned_cov95",
                          "tuned_cal_crps", "tuned_cal_cov95", "rmse_ratio_tuned_fdm"),
               value = c(mean(b0$crps), rmse(b0), mean(b0$cov95), mean(b1$crps), rmse(b1), mean(b1$cov95),
                         mean(b2$crps), mean(b2$cov95), rmse(b1) / rmse(b0)))
  })

checks$smoothing_jumpoff <- list(
  suite = "quick",
  description = paste("Smoothing over time and jump-off control. England & Wales men (as backtest_ew, origins",
                      "1996-2001, h = 10): Lee-Carter fitted to rates smoothed by the age-period GAM, with and",
                      "without recalibrated() (smoothing alone narrows the intervals, because the period index is",
                      "smooth); Lee-Carter anchored to the observed rates of the last year (Lee and Miller 2001,",
                      "the default of anchored()). ONS quintiles (origins 2008-2014, h = 5): Li-Lee with and",
                      "without anchoring."),
  run = function() {
    ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
    ew <- mt_set_grid(ew, age = single_ages(30, 100, open = FALSE))
    ew <- mt_regroup(ew, age = single_ages(30, 100, open = FALSE))
    o <- 1996:2001
    gam2d <- mt_prep(step_smooth(dims = 2, method = "gam"))
    lc <- mt_backtest(lee_carter(), ew, origins = o, h = 10, n_draws = 200)
    sm <- mt_backtest(lee_carter(), ew, origins = o, h = 10, n_draws = 200, prep = gam2d)
    smc <- mt_backtest(recalibrated(lee_carter()), ew, origins = o, h = 10, n_draws = 200, prep = gam2d)
    an <- mt_backtest(anchored(lee_carter()), ew, origins = o, h = 10, n_draws = 200)
    x <- ons_quintiles()
    ll <- mt_backtest(li_lee(common_by = "sex"), x, origins = 2008:2014, h = 5, n_draws = 200)
    lla <- mt_backtest(anchored(li_lee(common_by = "sex")), x, origins = 2008:2014, h = 5, n_draws = 200)
    h12 <- function(b) mean(b$crps[b$horizon <= 2])
    data.frame(metric = c("ew_lc_crps", "ew_lc_cov95", "ew_smooth_cov95", "ew_smooth_cal_crps", "ew_smooth_cal_cov95",
                          "ew_lc_crps_h12", "ew_anchored_crps_h12", "ons_lilee_crps", "ons_lilee_anchored_crps",
                          "ons_anchored_ratio"),
               value = c(mean(lc$crps), mean(lc$cov95), mean(sm$cov95), mean(smc$crps), mean(smc$cov95), h12(lc),
                         h12(an), mean(ll$crps), mean(lla$crps), mean(lla$crps) / mean(ll$crps)))
  })

checks$apc_bayes_france <- list(
  suite = "full",
  description = paste("apc_bayes() (bamp, fast MCMC) on French women (HMD via demography), 1970-2006, deaths in",
                      "5-year groups from age 40: convergence, and a backtest from origin 1996 (h = 10)."),
  run = function() {
    fr <- suppressWarnings(as_mt_data(mortensor:::pkg_data("fr.mort", "demography"), series = "female",
                                      years = 1970:2006))
    fr <- mt_regroup(fr, age = 40:100)
    x <- fr
    x$obs$deaths <- mt_regroup(fr, age = seq(40, 100, 5))$obs$deaths
    prep <- mt_prep(step_ungroup())
    fit <- mt_fit(apc_bayes(), x, prep = prep)
    bt <- mt_backtest(apc_bayes(), x, origins = 1996, h = 10, prep = prep, n_draws = 300)
    data.frame(metric = c("converged", "crps", "cov95"),
               value = c(mean(fit$diagnostics$status == "ok"), mean(bt$crps), mean(bt$cov95)))
  })

checks$coverage_li_lee <- list(
  suite = "full",
  description = "As coverage_lee_carter for li_lee() (independent trends per stratum violate its assumption).",
  run = function() coverage_check(li_lee(), 10))

checks$coverage_fdm <- list(
  suite = "full",
  description = "As coverage_lee_carter for hyndman_ullah() (functional demographic model, demography).",
  run = function() coverage_check(hyndman_ullah(order = 3), 10))

checks$coherent_fdm_france <- list(
  suite = "full",
  description = paste("coherent_fdm() against hyndman_ullah() per sex on French women and men (HMD via demography),",
                      "ages 30-99, 1950-2006: backtest from 1986, 1991 and 1996 (h = 10), and the female/male rate",
                      "ratio at age 70 after 50 years relative to the last observed ratio (median over draws)."),
  run = function() {
    fr <- suppressWarnings(as_mt_data(mortensor:::pkg_data("fr.mort", "demography"), series = c("female", "male"),
                                      years = 1950:2006))
    fr <- mt_regroup(fr, age = 30:100)
    fr <- mt_subset(fr, ages = c(30, 100))
    org <- c(1986, 1991, 1996)
    cf <- mt_backtest(coherent_fdm(order1 = 3, order2 = 2), fr, origins = org, h = 10, n_draws = 200)
    sep <- mt_backtest(hyndman_ullah(order = 3), fr, origins = org, h = 10, n_draws = 200)
    fc <- mt_forecast(mt_fit(coherent_fdm(order1 = 3, order2 = 2), fr), h = 50, n_draws = 200, seed = 1)
    a70 <- which(fc$age$lower == 70)
    ratio_end <- stats::median(fc$rates[a70, 50, "female", ] / fc$rates[a70, 50, "male", ])
    ratio_obs <- mean(fc$observed[a70, 53:57, "female"] / fc$observed[a70, 53:57, "male"])
    data.frame(metric = c("crps", "cov95", "crps_ratio_separate", "ratio_drift"),
               value = c(mean(cf$crps), mean(cf$cov95), mean(cf$crps) / mean(sep$crps),
                         abs(log(ratio_end / ratio_obs))))
  })

checks$stanmomo_ew <- list(
  suite = "full",
  description = paste("stanmomo() (StanMoMo, 4 chains of 2000 iterations) against lee_carter() on England & Wales men",
                      "(HMD via StMoMo), single ages 55-89, 1971-2011: backtest from 1991, 1996 and 2001 (h = 10)",
                      "of the Bayesian Lee-Carter and Renshaw-Haberman models; CRPS, 95% coverage, CRPS relative to",
                      "lee_carter(), and the largest R-hat of the fits."),
  run = function() {
    ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
    ew <- mt_regroup(mt_set_grid(ew, age = single_ages(55, 90, open = FALSE)), age = single_ages(55, 90, open = FALSE))
    ew <- mt_subset(ew, periods = c(1971, 2011))
    org <- c(1991, 1996, 2001)
    set.seed(1)
    bt <- mt_compare(list(lc = lee_carter(), stan_lc = stanmomo("lc", h_max = 10, cores = 4),
                          stan_rh = stanmomo("rh", h_max = 10, cores = 4)), ew, origins = org, h = 10, n_draws = 300)
    s <- mt_score(bt)
    rhat <- vapply(c("lc", "rh"), function(m) {
      f <- suppressWarnings(mt_fit(stanmomo(m, h_max = 10, cores = 4), ew))
      f$fits[[1]]$rhat
    }, 0)
    crps <- function(m) s[model == m, crps]
    data.frame(metric = c("stan_lc_crps", "stan_lc_cov95", "stan_lc_crps_ratio", "stan_rh_crps", "stan_rh_cov95",
                          "rhat_lc", "rhat_rh"),
               value = c(crps("stan_lc"), s[model == "stan_lc", cov95], crps("stan_lc") / crps("lc"),
                         crps("stan_rh"), s[model == "stan_rh", cov95], rhat[["lc"]], rhat[["rh"]]))
  })

checks$causes_reconciled <- list(
  suite = "quick",
  description = paste("Lee-Carter per cause and for all causes, reconciled with WLS (reconciled(), hierarchy over",
                      "cause), with step_smooth() on mt_example_causes(): largest relative gap between all-cause",
                      "draws and the sum of cause draws, and CRPS in a backtest from 2012 (h = 5)."),
  run = function() {
    x <- mt_example_causes(years = 1995:2017)
    prep <- mt_prep(step_smooth())
    spec <- reconciled(lee_carter(), mt_hierarchy(cause = list(total = "all", additive = "rates")))
    fit <- mt_fit(spec, mt_subset(x, periods = c(NA, 2012)), prep = prep)
    fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
    st <- fc$strata
    gap <- max(vapply(c("female", "male"), function(sx) {
      i <- st$stratum[st$sex == sx & st$cause == "all"]
      j <- st$stratum[st$sex == sx & st$cause != "all"]
      max(abs(apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum) / fc$rates[, , i, ] - 1))
    }, 0))
    bt <- mt_backtest(spec, x, origins = 2012, h = 5, n_draws = 100, prep = prep)
    data.frame(metric = c("max_gap", "crps", "cov95"), value = c(gap, mean(bt$crps), mean(bt$cov95)))
  })

checks$regions_coherence <- list(
  suite = "quick",
  description = paste("mt_example_regions() (2 sexes x 4 regions x 5 quintiles plus all aggregates, 90 strata):",
                      "backtest from 2005 and 2008 (h = 10) of li_lee(common_by = 'sex') aggregated bottom-up and",
                      "of Lee-Carter reconciled with WLS; CRPS and 95% coverage for bottom and aggregate strata,",
                      "and the largest coherence gap of the reconciled forecast."),
  run = function() {
    x <- mt_example_regions(years = 1990:2019, regions = 4)
    h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
    xa <- mt_aggregate_strata(x, h)
    hs <- mortensor:::hierarchy_structure(mt_strata(xa), h)
    lev <- ifelse(hs$is_bottom, "bottom", "aggregate")
    score <- function(sp) {
      b <- mt_backtest(sp, xa, origins = c(2005, 2008), h = 10, n_draws = 200)
      b[, lev := lev[stratum]]
      b[, list(crps = mean(crps), cov95 = mean(cov95)), by = lev]
    }
    a <- score(reconciled(li_lee(common_by = "sex"), h, method = "bottom_up"))
    w <- score(reconciled(lee_carter(), h, method = "wls"))
    fc <- mt_forecast(mt_fit(reconciled(lee_carter(), h, method = "wls"), xa), h = 3, n_draws = 50, seed = 1)
    E <- mortensor:::expand_exposure(fc$exposure_last, dim(fc$rates))
    D <- fc$rates * as.vector(E)
    gap <- max(vapply(seq_len(dim(D)[4]), function(j) {
      full <- matrix(D[, 3, hs$bottom, j], dim(D)[1]) %*% t(hs$S)
      max(abs(full / matrix(D[, 3, , j], dim(D)[1]) - 1))
    }, 0))
    data.frame(metric = c("lilee_crps_bottom", "lilee_crps_aggregate", "lilee_cov95_bottom", "lilee_cov95_aggregate",
                          "lc_wls_cov95_bottom", "lc_wls_cov95_aggregate", "max_gap"),
               value = c(a[lev == "bottom"]$crps, a[lev == "aggregate"]$crps, a[lev == "bottom"]$cov95,
                         a[lev == "aggregate"]$cov95, w[lev == "bottom"]$cov95, w[lev == "aggregate"]$cov95, gap))
  })

# ONS data by sex and IMD quintile (all causes, or cause groups at ages 50+)
ons_quintiles <- function(causes = FALSE) {
  cols <- names(ons_data)[-(1:7)]
  x <- mt_data(deaths = ons_data, population = ons_data,
               by = c("sex", "imd_decile", if (causes) "cause"), age = "age_group",
               value = list(deaths = if (causes) c(all = "all_causes", cols) else "all_causes",
                            population = "population"))
  x <- mt_as_exposure(x, type = "midyear")
  x <- mt_group_levels(x, "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)), name = "quintile")
  if (!causes) return(x)
  x <- mt_subset(x, ages = c(50, Inf))
  mt_group_levels(x, "cause", ons_cause_groups, other = "rest", total = "all", additive = "rates")
}

checks$decomposition <- list(
  suite = "quick",
  description = paste("mt_decompose() on ONS England (sex x IMD quintile): the gap in life expectancy at birth",
                      "Q5 - Q1 in 2019 by age (Arriaga, line integral, stepwise; largest residual and largest",
                      "difference of an age contribution from Arriaga), by age and cause at ages 50+ (residual and",
                      "agreement with the all-cause decomposition), the change 2001-2019 with quintile",
                      "composition (radix: exact; pooled: line integral) and the change 2019-2035 of Li-Lee draws."),
  run = function() {
    x <- ons_quintiles()
    d <- lapply(c(arriaga = "arriaga", horiuchi = "horiuchi", stepwise = "stepwise"), function(m) {
      mt_decompose(x, var = "quintile", ref = "Q1", periods = 2019, method = m)
    })
    resid <- function(dd) max(summary(dd, by = NULL)$max_residual)
    by_age <- function(dd) summary(dd, by = "age")$mean
    xc <- ons_quintiles(causes = TRUE)
    h <- mt_hierarchy(cause = "all", additive = "rates")
    dc <- mt_decompose(xc, var = "quintile", ref = "Q1", periods = 2019, cause = "cause", hierarchy = h)
    da <- mt_decompose(mt_subset(xc, strata = cause == "all"), var = "quintile", ref = "Q1", periods = 2019)
    dr <- mt_decompose(x, periods = c(2001, 2019), groups = "quintile", composition = "radix")
    dp <- mt_decompose(x, periods = c(2001, 2019), groups = "quintile", composition = "pooled")
    fc <- mt_forecast(mt_fit(li_lee(common_by = "sex"), x), to = 2035, n_draws = 200, seed = 1)
    df <- mt_decompose(fc, periods = c(2019, 2035))
    data.frame(metric = c("arriaga_residual", "stepwise_residual", "horiuchi_residual", "horiuchi_vs_arriaga",
                          "stepwise_vs_arriaga", "causes_residual", "causes_vs_allcause", "radix_residual",
                          "pooled_residual", "forecast_residual"),
               value = c(resid(d$arriaga), resid(d$stepwise), resid(d$horiuchi),
                         max(abs(by_age(d$horiuchi) - by_age(d$arriaga))),
                         max(abs(by_age(d$stepwise) - by_age(d$arriaga))), resid(dc),
                         max(abs(by_age(dc) - by_age(da))), resid(dr), resid(dp), resid(df)))
  })

checks$ons_deprivation <- list(
  suite = "quick",
  description = paste("ONS England 2001-2019 (ons_data), sex x IMD quintile, 5-year groups: backtest from 2008 to",
                      "2014 (h = 5) of li_lee(common_by = 'sex') and of the mixture of Lee-Carter and rwd(); CRPS",
                      "and 95% coverage. Causes at ages 50+ (ons_cause_groups): Lee-Carter per cause summed",
                      "bottom-up, all-cause CRPS from origins 2008, 2011, 2014 and the largest coherence gap."),
  run = function() {
    x <- ons_quintiles()
    ll <- mt_backtest(li_lee(common_by = "sex"), x, origins = 2008:2014, h = 5, n_draws = 200)
    en <- mt_backtest(ensemble(lc = lee_carter(), rw = rwd()), x, origins = 2008:2014, h = 5, n_draws = 200)
    xc <- ons_quintiles(causes = TRUE)
    h <- mt_hierarchy(cause = "all", additive = "rates")
    sp <- reconciled(lee_carter(), h, method = "bottom_up")
    bc <- mt_backtest(sp, xc, origins = c(2008, 2011, 2014), h = 5, n_draws = 100)
    fc <- mt_forecast(mt_fit(sp, xc), h = 5, n_draws = 50, seed = 1)
    st <- mt_strata(xc)
    gap <- max(vapply(unique(paste(st$sex, st$quintile)), function(g) {
      i <- which(paste(st$sex, st$quintile) == g & st$cause == "all")
      j <- which(paste(st$sex, st$quintile) == g & st$cause != "all")
      max(abs(apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum) / fc$rates[, , i, ] - 1))
    }, 0))
    data.frame(metric = c("lilee_crps", "lilee_cov95", "ensemble_crps", "ensemble_cov95", "causes_bu_crps_all",
                          "causes_max_gap"),
               value = c(mean(ll$crps), mean(ll$cov95), mean(en$crps), mean(en$cov95),
                         mean(bc[cause == "all"]$crps), gap))
  })

checks$de_regions <- list(
  suite = "quick",
  description = paste("German Laender 1998-2019 (de_deaths, de_population), all causes: suppressed cells filled to",
                      "the national totals; Laender by sex with East, West and Germany; backtest from 2008, 2011,",
                      "2014 (h = 5) of the ensemble of li_lee(common_by = 'sex') and ARIMA(0,1,1) with drift,",
                      "bottom-up; CRPS and 95% coverage of bottom and aggregate strata; share of 2035 draws with",
                      "rates above 1."),
  run = function() {
    d <- as.data.frame(de_deaths)
    d <- d[d$cause == "all" & d$year <= 2019, setdiff(names(d), "cause")]
    x <- suppressWarnings(mt_data(deaths = d, population = de_population[de_population$year <= 2019, ],
                                  by = c("region", "sex")))
    x <- mt_regroup(mt_as_exposure(x, type = "dec31"), age = "deaths")
    x <- mt_fill_suppressed(x, mt_hierarchy(region = "Deutschland"))
    gap <- attr(x, "fill")$max_gap
    xa <- mt_subset(x, region != "Deutschland")
    laender <- unique(mt_strata(xa)$region)
    east <- c("Berlin", "Brandenburg", "Mecklenburg-Vorpommern", "Sachsen", "Sachsen-Anhalt",
              grep("^Th", laender, value = TRUE))
    h <- mt_hierarchy(region = list(total = "Deutschland", groups = list(West = setdiff(laender, east), East = east)),
                      sex = "both")
    xh <- mt_aggregate_strata(xa, h)
    sp <- reconciled(ensemble(ll = li_lee(common_by = "sex"), ll0 = arima_series(order = c(0, 1, 1), drift = TRUE)), h,
                     method = "bottom_up")
    bt <- mt_backtest(sp, xh, origins = c(2008, 2011, 2014), h = 5, n_draws = 200)
    bt <- merge(bt, mt_strata(xh, h)[, list(label, level)], by = "label")
    fc <- suppressWarnings(mt_forecast(mt_fit(sp, xh), to = 2035, n_draws = 200, seed = 1))
    closed <- is.finite(fc$age$upper)
    implausible <- mean(apply(fc$rates[closed, , , , drop = FALSE] > 1, 4, any))
    data.frame(metric = c("fill_max_gap", "crps_bottom", "crps_aggregate", "cov95_bottom", "cov95_aggregate",
                          "implausible_draws"),
               value = c(gap, mean(bt[level == "bottom"]$crps), mean(bt[level != "bottom"]$crps),
                         mean(bt[level == "bottom"]$cov95), mean(bt[level != "bottom"]$cov95), implausible))
  })

checks$ungroup_ew_2d <- list(
  suite = "full",
  description = "As ungroup_ew_1d for 1980-2011 with the 2D PCLM (age-period surface).",
  run = function() {
    r <- mortensor:::benchmark_ungroup_ew(years = 1980:2011, dims = 2)
    data.frame(metric = c("male_30_89", "male_90_100", "seconds"), value = c(r$male_30_89, r$male_90_100, r$seconds))
  })

checks$unit_tests <- list(
  suite = "full",
  description = "testthat suite including slow tests (MORTENSOR_SLOW_TESTS=true).",
  run = function() {
    old <- Sys.setenv(MORTENSOR_SLOW_TESTS = "true")
    on.exit(Sys.unsetenv("MORTENSOR_SLOW_TESTS"))
    tdir <- file.path(root, "tests", "testthat")
    if (!dir.exists(tdir)) stop("tests/testthat not found")
    r <- as.data.frame(testthat::test_dir(tdir, package = "mortensor", load_package = "installed",
                                          reporter = "silent", stop_on_failure = FALSE))
    data.frame(metric = c("failed", "errors", "passed"),
               value = c(sum(r$failed), sum(r$error), sum(r$passed)))
  })

checks <- if (exists("checks")) checks else list()

checks$coda_ew <- list(
  suite = "quick",
  description = paste("England & Wales men (HMD via StMoMo), ages 30-100, 1961-2011: rolling-origin backtest",
                      "(origins 1996-2001, h = 10) of coda() (CoDa Lee-Carter, observed jump-off) against",
                      "lee_carter(); CRPS of log rates, coverage and the ratio of the CRPS."),
  run = function() {
    ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
    ew <- mt_set_grid(ew, age = single_ages(30, 100, open = FALSE))
    ew <- mt_regroup(ew, age = single_ages(30, 100, open = FALSE))
    o <- 1996:2001
    cd <- mt_backtest(coda(), ew, origins = o, h = 10, n_draws = 200)
    lc <- mt_backtest(lee_carter(), ew, origins = o, h = 10, n_draws = 200)
    data.frame(metric = c("coda_crps", "coda_cov95", "coda_cov80", "lc_crps", "crps_ratio_coda_lc"),
               value = c(mean(cd$crps), mean(cd$cov95), mean(cd$cov80), mean(lc$crps), mean(cd$crps) / mean(lc$crps)))
  })

checks$coda_coherent_fr <- list(
  suite = "full",
  description = paste("France (HMD via demography::fr.mort), women and men, ages 0-99 and 100+, 1950-2006:",
                      "backtest (origins 1991-1996, h = 10) of coda(coherent = TRUE) against coda() and li_lee();",
                      "CRPS and coverage of log rates, CRPS of the open interval (100+), and the gap in life",
                      "expectancy at birth between women and men in 2006 and after 50 years of forecasts (medians",
                      "of the draws of the coherent and the independent model; information, no threshold)."),
  run = function() {
    fr <- as_mt_data(mortensor:::pkg_data("fr.mort", "demography"))
    fr <- mortensor:::mt_window(fr, from = 1950, to = 2006)
    fr <- mt_regroup(fr, age = single_ages(0, 100, open = TRUE))
    o <- 1991:1996
    cc <- mt_backtest(coda(coherent = TRUE), fr, origins = o, h = 10, n_draws = 200)
    cd <- mt_backtest(coda(), fr, origins = o, h = 10, n_draws = 200)
    ll <- mt_backtest(li_lee(), fr, origins = o, h = 10, n_draws = 200)
    e0_gap <- function(e, year) e[sex == "female" & period == year]$q0.5 - e[sex == "male" & period == year]$q0.5
    gap_50 <- function(spec) e0_gap(life_expectancy(mt_forecast(mt_fit(spec, fr), h = 50, n_draws = 200, seed = 1)),
                                    2056)
    data.frame(metric = c("coherent_crps", "coherent_cov95", "coherent_cov80", "open_crps", "coda_crps",
                          "lilee_crps", "crps_ratio_coherent_lilee", "e0_gap_2006", "e0_gap_2056_coherent",
                          "e0_gap_2056_independent"),
               value = c(mean(cc$crps), mean(cc$cov95), mean(cc$cov80), mean(cc[age_lo == 100]$crps), mean(cd$crps),
                         mean(ll$crps), mean(cc$crps) / mean(ll$crps), e0_gap(life_expectancy(fr), 2006),
                         gap_50(coda(coherent = TRUE)), gap_50(coda())))
  })

checks$coda_causes_coherence <- list(
  suite = "quick",
  description = paste("coda_causes() (multiple-decrement CoDa, experimental) on mt_example_causes(), 1995-2017:",
                      "largest relative gap between all-cause draws and the sum of the cause draws (expected and",
                      "predictive draws, origin 2012, h = 5), and a backtest (origins 2008-2012, h = 5) against",
                      "Lee-Carter per cause reconciled with WLS; CRPS and coverage of log rates."),
  run = function() {
    x <- mt_example_causes(years = 1995:2017)
    st <- mt_strata(x)
    fit <- mt_fit(coda_causes(), mortensor:::mt_window(x, to = 2012))
    gap <- function(fc) {
      max(vapply(c("female", "male"), function(sx) {
        i <- st$stratum[st$sex == sx & st$cause == "all"]
        j <- st$stratum[st$sex == sx & st$cause != "all"]
        tot <- fc$rates[, , i, ]
        sums <- apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum)
        max(abs(sums - tot) / pmax(tot, 1e-12))
      }, 0))
    }
    g1 <- gap(mt_forecast(fit, h = 5, n_draws = 100, seed = 1))
    g2 <- gap(mt_forecast(fit, h = 5, n_draws = 100, seed = 1, predictive = TRUE))
    o <- 2008:2012
    cc <- mt_backtest(coda_causes(), x, origins = o, h = 5, n_draws = 100)
    cs <- mt_backtest(reconciled(lee_carter(), mt_hierarchy(cause = list(total = "all", additive = "rates"))), x, origins = o, h = 5, n_draws = 100)
    data.frame(metric = c("max_gap", "max_gap_predictive", "cc_crps", "cc_cov95", "cs_crps", "crps_ratio_cc_cs"),
               value = c(g1, g2, mean(cc$crps), mean(cc$cov95), mean(cs$crps), mean(cc$crps) / mean(cs$crps)))
  })

if (!exists("checks")) checks <- list()

checks$rotation_france <- list(
  suite = "full",
  description = paste("lee_carter(rotate = TRUE) on French women (HMD via demography), ages 0-99 and 100+, 1950-2006:",
                      "forecast to 2100, where rotated and unrotated paths have the same life expectancy at birth",
                      "and the rotated ones a slower decline at ages 0-14 and a faster one at ages 15-64 (average",
                      "annual rate of decline of the death rate of each age band, 2006 exposures, median over",
                      "draws); backtests from 1986, 1991 and 1996 (h = 10) and from 1976, 1981 and 1986 (h = 20),",
                      "CRPS of log rates relative to the unrotated model."),
  run = function() {
    fr <- suppressWarnings(as_mt_data(mortensor:::pkg_data("fr.mort", "demography"), series = "female",
                                      years = 1950:2006))
    fr <- mt_regroup(fr, age = 0:100)
    f0 <- mt_fit(lee_carter(), fr)
    f1 <- mt_fit(lee_carter(rotate = TRUE), fr)
    fc0 <- mt_forecast(f0, to = 2100, n_draws = 200, seed = 1)
    t0 <- proc.time()[["elapsed"]]
    fc1 <- mt_forecast(f1, to = 2100, n_draws = 200, seed = 1)
    secs <- proc.time()[["elapsed"]] - t0
    e0 <- life_expectancy(fc0, paths = TRUE)
    e1 <- life_expectancy(fc1, paths = TRUE)
    lo <- fc0$age$lower
    band <- cut(lo, c(0, 15, 65, 85, Inf), right = FALSE, labels = c("0-14", "15-64", "65-84", "85+"))
    ex <- mt_array(fr, "exposure")[, length(f0$periods), 1]
    group_rates <- function(m) tapply(ex * m, band, sum) / tapply(ex, band, sum)
    jump <- group_rates(f0$fits[[1]]$fitted[, length(f0$periods)])
    nh <- length(fc0$periods)
    decline <- function(fc) {
      r <- vapply(seq_len(dim(fc$rates)[4]), function(j) -log(group_rates(fc$rates[, nh, 1, j]) / jump) / nh,
                  numeric(4))
      100 * apply(r, 1, stats::median)
    }
    d0 <- decline(fc0)
    d1 <- decline(fc1)
    bt <- function(spec, origins, h) mt_backtest(spec, fr, origins = origins, h = h, n_draws = 200)
    b10 <- bt(lee_carter(), c(1986, 1991, 1996), 10)
    r10 <- bt(lee_carter(rotate = TRUE), c(1986, 1991, 1996), 10)
    b20 <- bt(lee_carter(), c(1976, 1981, 1986), 20)
    r20 <- bt(lee_carter(rotate = TRUE), c(1976, 1981, 1986), 20)
    data.frame(metric = c("e0_max_diff", "decline_unrot_0_14", "decline_rot_0_14", "decline_diff_0_14",
                          "decline_diff_15_64", "decline_diff_65_84", "decline_diff_85_plus", "weight_2100",
                          "crps_unrot_h10", "crps_rot_h10", "crps_ratio_h10", "crps_unrot_h20", "crps_rot_h20",
                          "crps_ratio_h20", "rot_cov95_h20", "seconds_forecast"),
               value = c(max(abs(e1 - e0)), d0[1], d1[1], d1[1] - d0[1], d1[2] - d0[2], d1[3] - d0[3],
                         d1[4] - d0[4], stats::median(mortensor:::rotation_weights(e1[nh, 1, ])),
                         mean(b10$crps), mean(r10$crps), mean(r10$crps) / mean(b10$crps), mean(b20$crps),
                         mean(r20$crps), mean(r20$crps) / mean(b20$crps), mean(r20$cov95), secs))
  })

if (!exists("checks")) checks <- list()

# England & Wales men (HMD via StMoMo) on single years of age
ew_single_ages <- function(from, to) {
  ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
  a <- single_ages(from, to, open = FALSE)
  mt_regroup(mt_set_grid(ew, age = a), age = a)
}

checks$dm_ew <- list(
  suite = "quick",
  description = paste("Diebold-Mariano test (mt_dm_test()) on the backtest of backtest_ew (England & Wales men, ages",
                      "30-100, origins 1996-2001, h = 10): Lee-Carter against rwd(), CRPS of log rates. Agreement",
                      "of varestimator = 'acf' with forecast::dm.test() on the losses per origin (largest absolute",
                      "difference of the statistic and the p-value over the tested horizons, with the lag used;",
                      "horizons beyond 5 and all horizons together are not tested with 6 consecutive origins);",
                      "p-values of the default fixed-smoothing (EWC) test at horizons 1 and 5, and over all",
                      "horizons for 7 origins 5 years apart (1971-2001), as info."),
  run = function() {
    ew <- ew_single_ages(30, 100)
    models <- list(rwd = rwd(), lee_carter = lee_carter())
    bt <- suppressWarnings(mt_compare(models, ew, origins = 1996:2001, h = 10, n_draws = 200))
    r <- suppressWarnings(mt_dm_test(bt, ref = "rwd", varestimator = "acf"))
    rw <- suppressWarnings(mt_dm_test(bt, ref = "rwd"))
    b <- as.data.frame(bt)
    loss <- function(mod, h = NULL) {
      s <- b[b$model == mod, ]
      if (!is.null(h)) return(tapply(s$crps[s$horizon == h], s$origin[s$horizon == h], mean))
      rowMeans(tapply(s$crps, list(s$origin, s$horizon), mean))
    }
    hs <- as.character(r$horizon)
    # rows without a test (too few origins for the lags of overlapping years) are not compared
    ref <- t(vapply(seq_along(hs), function(i) {
      if (is.na(r$lag[i])) return(c(NA_real_, NA_real_))
      h <- if (hs[i] == "all") NULL else as.numeric(hs[i])
      d <- suppressWarnings(forecast::dm.test(loss("lee_carter", h), loss("rwd", h), h = r$lag[i] + 1, power = 1))
      c(unname(d$statistic), unname(d$p.value))
    }, numeric(2)))
    bt5 <- suppressWarnings(mt_compare(models, ew, origins = seq(1971, 2001, by = 5), h = 10, n_draws = 200))
    r5 <- suppressWarnings(mt_dm_test(bt5, ref = "rwd"))
    data.frame(metric = c("max_abs_diff_stat", "max_abs_diff_p", "mean_diff_all", "p_h1", "p_h5", "p_all_spaced"),
               value = c(max(abs(r$statistic - ref[, 1]), na.rm = TRUE),
                         max(abs(r$p_value - ref[, 2]), na.rm = TRUE), r$mean_diff[hs == "all"],
                         rw$p_value[rw$horizon == "1"], rw$p_value[rw$horizon == "5"],
                         r5$p_value[r5$horizon == "all"]))
  })

checks$dm_size <- list(
  suite = "quick",
  description = paste("Size of mt_dm_test() at nominal 5% under the null: loss differences of forecasts h years",
                      "ahead from n consecutive origins simulated as a moving average of order h - 1 with mean",
                      "zero (2000 replicates for each n in 7, 16 and h in 1, 3, 5). Rejection rates of the default",
                      "fixed-smoothing (EWC) test, of the Daniell-kernel test of Coroneo and Iacone (2020;",
                      "varestimator = 'daniell') at h = 1 and over all h (info), and of the HLN-corrected test",
                      "(varestimator = 'acf', info)."),
  run = function() {
    set.seed(5)
    grid <- expand.grid(n = c(7, 16), h = c(1, 3, 5))
    size <- function(n, h, ve) {
      mean(vapply(seq_len(2000), function(i) {
        e <- stats::rnorm(n + h - 1)
        d <- stats::filter(e, rep(1, h), sides = 1)[h:(n + h - 1)]
        p <- mortensor:::dm_stat(d, h, varestimator = ve)$p_value
        !is.na(p) && p < 0.05
      }, TRUE))
    }
    ewc <- mapply(size, grid$n, grid$h, MoreArgs = list(ve = "ewc"))
    acf <- mapply(size, grid$n, grid$h, MoreArgs = list(ve = "acf"))
    dan <- mapply(size, grid$n, grid$h, MoreArgs = list(ve = "daniell"))
    data.frame(metric = c("ewc_size_max", "ewc_size_min", "acf_size_max", "daniell_size_h1_max", "daniell_size_max"),
               value = c(max(ewc), min(ewc), max(acf), max(dan[grid$h == 1]), max(dan)))
  })

checks$stmomo_cbd_ew <- list(
  suite = "full",
  description = paste("CBD family of stmomo() on England & Wales men (HMD via StMoMo), ages 55-89, 1961-2011:",
                      "rolling-origin backtest (origins 1996-2001, h = 10) of CBD, M6, M7 and M8 (xc = 89); CRPS",
                      "and 95% coverage of log rates, CRPS ratio M7 / CBD and the Diebold-Mariano p-value at horizon 1",
                      "(info). Fitted rates of M7 and M8 on 1961-2011 against StMoMo's own example",
                      "(central2initial(), clip = 3)."),
  run = function() {
    ew <- ew_single_ages(55, 89)
    models <- list(cbd = stmomo("cbd"), m6 = stmomo("m6"), m7 = stmomo("m7"), m8 = stmomo("m8"))
    bt <- mt_compare(models, ew, origins = 1996:2001, h = 10, n_draws = 200)
    s <- mt_score(bt)
    r <- suppressWarnings(mt_dm_test(bt, ref = "cbd"))
    crps <- function(m) s$crps[s$model == m]
    # fitted rates against StMoMo called as in its documentation
    sd <- StMoMo::central2initial(mortensor:::pkg_data("EWMaleData", "StMoMo"))
    wxt <- StMoMo::genWeightMat(55:89, sd$years, clip = 3)
    rel <- function(m, mod) {
      direct <- suppressWarnings(StMoMo::fit(mod, data = sd, ages.fit = 55:89, wxt = wxt, verbose = FALSE))
      q <- stats::fitted(direct, type = "rates")
      f <- stats::fitted(mt_fit(stmomo(m), ew))[, , 1]
      ok <- !is.na(q)
      max(abs(f[ok] / -log(1 - q[ok]) - 1))
    }
    data.frame(metric = c("cbd_crps", "m6_crps", "m7_crps", "m8_crps", "cbd_cov95", "m7_cov95", "m8_cov95",
                          "crps_ratio_m7_cbd", "p_m7_cbd_h1", "m7_max_rel_diff_stmomo", "m8_max_rel_diff_stmomo"),
               value = c(crps("cbd"), crps("m6"), crps("m7"), crps("m8"), s$cov95[s$model == "cbd"],
                         s$cov95[s$model == "m7"], s$cov95[s$model == "m8"], crps("m7") / crps("cbd"),
                         r$p_value[r$model == "m7" & r$horizon == "1"], rel("m7", StMoMo::m7()),
                         rel("m8", StMoMo::m8(xc = 89))))
  })

# France (HMD via demography::fr.mort) life tables: single ages 0-109 and
# 110+, or deaths and exposures summed to the age groups starting at `breaks`
fr_life_tables <- function() {
  fr <- mortensor:::pkg_data("fr.mort", "demography")
  grouped <- function(s, yr, breaks) {
    m <- fr$rate[[s]][, as.character(yr)]
    p <- fr$pop[[s]][, as.character(yr)]
    g <- findInterval(0:110, breaks)
    D <- tapply(ifelse(is.na(m), 0, m * p), g, sum)
    E <- tapply(p, g, sum)
    list(m = as.vector(D / E), age = age_groups(breaks))
  }
  list(female_2000_single = list(m = fr$rate$female[, "2000"], age = single_ages(0, 110)),
       male_1987_single = list(m = fr$rate$male[, "1987"], age = single_ages(0, 110)),
       female_1950_5y = grouped("female", 1950, seq(0, 100, 5)))
}

checks$lifetable_thatcher <- list(
  suite = "quick",
  description = paste("Kannisto and Gompertz closures of life_table() fitted by maximum likelihood as in Thatcher,",
                      "Kannisto and Vaupel (1998): from the probabilities of dying at ages 80-98 of the cohorts born",
                      "in 1871-1880 (their Table 6.1, four decimals), for which the life-table survivors are the",
                      "numbers reaching each age, the fitted a and b against their Appendix D (Tables 5 and 1);",
                      "largest relative differences, and those of the former least-squares fit (info)."),
  run = function() {
    q <- list(male = c(0.1224, 0.1326, 0.1452, 0.1565, 0.1702, 0.1822, 0.1965, 0.2116, 0.2276, 0.2424, 0.2580,
                       0.2749, 0.2907, 0.3073, 0.3253, 0.3434, 0.3571, 0.3691, 0.3817),
              female = c(0.0987, 0.1074, 0.1180, 0.1291, 0.1409, 0.1520, 0.1653, 0.1787, 0.1926, 0.2070, 0.2225,
                         0.2354, 0.2524, 0.2687, 0.2867, 0.3050, 0.3228, 0.3358, 0.3504))
    pub <- list(kannisto = list(male = c(3.39306e-05, 1.04287e-01), female = c(2.31232e-05, 1.05869e-01)),
                gompertz = list(male = c(2.17126e-04, 7.97629e-02), female = c(1.22095e-04, 8.41400e-02)))
    out <- list()
    for (law in names(pub)) {
      da <- db <- la <- lb <- 0
      for (sx in names(q)) {
        m <- c(q[[sx]] / (1 - q[[sx]] / 2), 0.6)
        tl <- attr(life_table(m, single_ages(80, 99), open = law, k = 19), "tail")
        da <- max(da, abs(tl$a / pub[[law]][[sx]][1] - 1))
        db <- max(db, abs(tl$b / pub[[law]][[sx]][2] - 1))
        # the former fit: least squares on the log (logit) rates at mid-year
        y <- if (law == "kannisto") stats::qlogis(m[1:19]) else log(m[1:19])
        cf <- stats::coef(stats::lm(y ~ I(80:98 + 0.5)))
        la <- max(la, abs(exp(cf[[1]]) / pub[[law]][[sx]][1] - 1))
        lb <- max(lb, abs(cf[[2]] / pub[[law]][[sx]][2] - 1))
      }
      out[[law]] <- data.frame(metric = paste0(law, c("_max_rel_diff_a", "_max_rel_diff_b", "_ls_rel_diff_a",
                                                      "_ls_rel_diff_b")), value = c(da, db, la, lb))
    }
    do.call(rbind, out)
  })

checks$reconcile_nonneg <- list(
  suite = "quick",
  description = paste("Non-negative reconciliation (Wickramasuriya, Turlach and Hyndman 2020): on 400 random",
                      "base forecasts of a hierarchy of 8 bottom series (total, 3 groups) with aggregates too low,",
                      "diagonal or full weights and 0, 1 or 2 fixed nodes, the largest difference (relative to",
                      "the largest base forecast) between the block principal pivoting solution and",
                      "quadprog::solve.QP() of their Eq. 3 (fixed nodes as equality constraints), the largest",
                      "relative change of a fixed node, and the mean number of exchanges (info)."),
  run = function() {
    if (!requireNamespace("quadprog", quietly = TRUE)) return(data.frame(metric = "skipped", value = 1))
    groups <- list(1:3, 4:5, 6:8)
    sagg <- rbind(rep(1, 8), t(vapply(groups, function(g) replace(numeric(8), g, 1), numeric(8))))
    S <- rbind(sagg, diag(8))
    Ut <- cbind(diag(4), -sagg)
    b <- 5:12
    set.seed(11)
    diff_qp <- diff_fixed <- 0
    its <- integer()
    for (i in 1:400) {
      y <- as.vector(S %*% (stats::rgamma(8, 2) * 10)) * exp(stats::rnorm(12, 0, 0.6))
      y[1:4] <- y[1:4] * stats::runif(4, 0.2, 1)
      full <- i %% 2 == 0
      Lam <- if (full) crossprod(matrix(stats::rnorm(144), 12)) / 12 + diag(stats::runif(12, 0.5, 2)) else
        diag(stats::runif(12, 0.5, 3))
      keep <- list(integer(), 1L, c(2L, 9L))[[i %% 3 + 1]]
      W <- mortensor:::keep_weights(if (full) Lam else diag(Lam), keep)
      x0 <- mortensor:::project_rows(matrix(y, 1), W, Ut)
      if (all(x0[, b] >= 0)) next
      nn <- mortensor:::nnls_reconcile(matrix(y, 1), W, Ut, b, x0)
      its <- c(its, nn$iterations)
      li <- solve(Lam)
      A <- diag(8)
      b0 <- numeric(8)
      if (length(keep)) {
        A <- cbind(t(S[keep, , drop = FALSE]), A)
        b0 <- c(y[keep], b0)
      }
      ref <- quadprog::solve.QP(t(S) %*% li %*% S, as.vector(t(S) %*% li %*% y), A, b0, meq = length(keep))$solution
      diff_qp <- max(diff_qp, max(abs(nn$X[1, b] - pmax(ref, 0))) / max(abs(y)))
      if (length(keep)) diff_fixed <- max(diff_fixed, max(abs((S %*% nn$X[1, b])[keep] / y[keep] - 1)))
    }
    data.frame(metric = c("max_rel_diff_quadprog", "max_rel_change_fixed", "cases", "mean_iterations"),
               value = c(diff_qp, diff_fixed, length(its), mean(its)))
  })

checks$lifespan_reference <- list(
  suite = "quick",
  description = paste("lifespan_variation() on France life tables (fr.mort: women 2000 and men 1987 in single ages to",
                      "110+, women 1950 in 5-year groups to 100+) at ages 0 and 65: largest relative difference",
                      "from LifeIneq 00.04.09 (e-dagger, entropy, Gini, standard deviation; reference values as in",
                      "tests/testthat/test-lifespan.R), largest difference of the interquartile range (single ages),",
                      "and the difference of the standard deviation, Gini coefficient and interquartile range from",
                      "1e6 ages at death sampled from the life table of women 2000 (uniform within intervals,",
                      "exponential in the open interval)."),
  run = function() {
    ref <- list(
      female_2000_single = rbind(edagger = c(9.663024709, 7.008585465), sd = c(14.18425557, 8.211834891),
                                 entropy = c(0.116674102, 0.3298611941), gini = c(0.08453618794, 0.0537459495),
                                 iqr = c(13.41011018, 11.22595719)),
      male_1987_single = rbind(edagger = c(12.4459491, 7.202768889), sd = c(17.07315081, 7.945863166),
                               entropy = c(0.1727633064, 0.480640979), gini = c(0.1237094046, 0.0568842088),
                               iqr = c(19.32465572, 12.05458686)),
      female_1950_5y = rbind(edagger = c(14.64681576, 6.98491648), sd = c(22.05606492, 7.485704412),
                             entropy = c(0.2120541073, 0.4770597137), gini = c(0.1577954255, 0.05295861841),
                             iqr = c(19.20655623, 11.21954679)))
    lts <- fr_life_tables()
    rel <- 0
    iqr <- 0
    for (nm in names(ref)) {
      lt <- lts[[nm]]
      for (k in 1:2) {
        for (ms in c("edagger", "sd", "entropy", "gini")) {
          v <- lifespan_variation(lt$m, ms, lt$age, at = c(0, 65)[k])
          rel <- max(rel, abs(v / ref[[nm]][ms, k] - 1))
        }
        if (grepl("single", nm)) {
          iqr <- max(iqr, abs(lifespan_variation(lt$m, "iqr", lt$age, at = c(0, 65)[k]) - ref[[nm]]["iqr", k]))
        }
      }
    }
    lt <- lts$female_2000_single
    r <- mortensor:::lt_columns(matrix(lt$m, ncol = 1), lt$age)
    na <- length(lt$m)
    n <- 1e6
    set.seed(1)
    j <- sample.int(na, n, replace = TRUE, prob = r$dx[, 1])
    y <- lt$age$lower[j] + ifelse(j < na, stats::runif(n), stats::rexp(n, r$dx[na, 1] / r$Lx[na, 1]))
    ys <- sort(y)
    gini <- 2 * sum(ys * seq_len(n)) / (n * sum(ys)) - (n + 1) / n
    data.frame(metric = c("max_rel_diff_lifeineq", "max_iqr_diff_lifeineq", "mc_sd_diff", "mc_gini_diff",
                          "mc_iqr_diff"),
               value = c(rel, iqr, abs(lifespan_variation(lt$m, "sd", lt$age) - stats::sd(y)),
                         abs(lifespan_variation(lt$m, "gini", lt$age) - gini),
                         abs(lifespan_variation(lt$m, "iqr", lt$age) - diff(stats::quantile(y, c(0.25, 0.75))))))
  })

# One replicate of the excess simulation: mt_example() with deaths at ages
# 70+ in 2016 and 2017 raised by 10%, Lee-Carter baseline 1990-2014, excess
# by sex and age group (30-69, 70+) and in total for 2015-2019, with the
# injected deaths as truth
excess_replicate <- function(r, n_draws = 500) {
  x <- mt_example(years = 1990:2019, deaths = "single", seed = r)
  d <- x$obs$deaths
  hit <- d$period_lo %in% c(2016, 2017) & d$age_lo >= 70
  add <- round(1.1 * d$value[hit]) - d$value[hit]
  inj <- data.table::data.table(sex = x$strata$sex[d$stratum[hit]], period = d$period_lo[hit], add = add)
  x$obs$deaths$value[hit] <- d$value[hit] + add
  ex <- mt_excess(lee_carter(), x, baseline = c(1990, 2014), ages = c(30, 70), n_draws = n_draws, seed = r)
  tot <- mt_excess(lee_carter(), x, baseline = c(1990, 2014), by = character(), n_draws = n_draws, seed = r)
  tr <- inj[, list(true = sum(add)), by = c("sex", "period")]
  ex <- merge(ex, tr, by = c("sex", "period"), all.x = TRUE)
  ex$true[ex$age == "30-69" | is.na(ex$true)] <- 0
  tt <- inj[, list(true = sum(add)), by = "period"]
  tot <- merge(tot, tt, by = "period", all.x = TRUE)
  tot$true[is.na(tot$true)] <- 0
  rbind(data.frame(level = "sex_age", true = ex$true, med = ex$excess_q0.5, lo = ex$excess_q0.025,
                   hi = ex$excess_q0.975),
        data.frame(level = "total", true = tot$true, med = tot$excess_q0.5, lo = tot$excess_q0.025,
                   hi = tot$excess_q0.975))
}

checks$excess_simulation <- list(
  suite = "quick",
  description = paste("mt_excess() on simulated data with known truth (mt_example, single ages, 1990-2019, 20",
                      "replicates): deaths at ages 70+ in 2016 and 2017 raised by 10%, Lee-Carter baseline",
                      "1990-2014, excess 2015-2019 by sex and age group (30-69, 70+) and in total; 95% coverage of",
                      "zero where there is no shock and of the injected deaths where there is one, and the mean",
                      "relative error of the median excess."),
  run = function() {
    t0 <- proc.time()[["elapsed"]]
    res <- do.call(rbind, lapply(1:20, excess_replicate))
    secs <- proc.time()[["elapsed"]] - t0
    covered <- res$lo <= res$true & res$true <= res$hi
    null <- res$true == 0
    shock <- res$true > 0
    data.frame(metric = c("cov95_null", "cov95_shock", "cov95_total_null", "rel_bias_shock", "seconds"),
               value = c(mean(covered[null & res$level == "sex_age"]), mean(covered[shock]),
                         mean(covered[null & res$level == "total"]),
                         abs(mean((res$med[shock] - res$true[shock]) / res$true[shock])), secs))
  })

checks$excess_france <- list(
  suite = "quick",
  description = paste("mt_excess() on France (HMD via demography::fr.mort), ages 0-99 and 100+: excess deaths in",
                      "2003 (August heatwave) and 2004 against a Lee-Carter baseline 1990-2002 (both sexes, median",
                      "and 95% interval), the loss of life expectancy at birth of women in 2003, and the P-scores of",
                      "1918 against a 1900-1913 baseline (ages 0-94 and 95+)."),
  run = function() {
    fr <- as_mt_data(mortensor:::pkg_data("fr.mort", "demography"))
    x <- mt_regroup(mt_subset(fr, periods = c(1990, 2004)), age = 0:100)
    tot <- mt_excess(lee_carter(), x, baseline = c(1990, 2002), by = character(), n_draws = 1000, seed = 1)
    le <- mt_excess(lee_carter(), x, baseline = c(1990, 2002), periods = 2003, quantity = "life_expectancy",
                    n_draws = 1000, seed = 1)
    y <- mt_regroup(mt_subset(fr, periods = c(1900, 1918)), age = 0:95)
    e18 <- mt_excess(lee_carter(), y, baseline = c(1900, 1913), periods = 1918, n_draws = 1000, seed = 1)
    data.frame(metric = c("excess_2003", "excess_2003_lower", "excess_2003_upper", "excess_2004",
                          "le_excess_2003_female", "p_score_1918_female", "p_score_1918_male"),
               value = c(tot$excess_q0.5[tot$period == 2003], tot$excess_q0.025[tot$period == 2003],
                         tot$excess_q0.975[tot$period == 2003], tot$excess_q0.5[tot$period == 2004],
                         le$excess_q0.5[le$sex == "female"], e18$p_score_q0.5[e18$sex == "female"],
                         e18$p_score_q0.5[e18$sex == "male"]))
  })

# Score-optimal reconciliation (reconciled(method = "energy_score"), R/hierarchy-scoreopt.R). The
# thresholds are regression guards around the values at 2026-10: with the defaults the method is
# less accurate than WLS on both data sets.
checks$scoreopt_probreco <- list(
  suite = "quick",
  description = paste("Score-optimal reconciliation against ProbReco::scoreopt() (Panagiotelis et al. 2023) on a",
                      "toy hierarchy of 7 nodes: 60 Adam steps from OLS with the same base draws; largest",
                      "difference of d and G."),
  run = function() {
    if (!requireNamespace("ProbReco", quietly = TRUE)) return(data.frame(metric = "skipped", value = 1))
    S <- rbind(c(1, 1, 1, 1), c(1, 1, 0, 0), c(0, 0, 1, 1), diag(4))
    n <- nrow(S)
    m <- ncol(S)
    nt <- 6
    q <- 20
    set.seed(2)
    Y <- vapply(seq_len(nt), function(t) as.vector(S %*% (c(5, 3, 2, 4) + stats::rnorm(4))), numeric(n))
    mu <- as.vector(S %*% c(5, 3, 2, 4)) * c(1.2, 0.9, 1.1, 1, 1, 1, 1)
    prob <- lapply(seq_len(nt), function(t) function() matrix(stats::rnorm(n * q, mu, 0.7), n))
    set.seed(3)
    ref <- suppressWarnings(ProbReco::scoreopt(lapply(seq_len(nt), function(t) Y[, t]), prob, S,
                                               control = list(maxIter = 60, eta = 0.01)))
    # ProbReco calls every function three times to check its inputs
    set.seed(3)
    for (k in 1:3) for (f in prob) f()
    draw <- function() {
      x1 <- lapply(prob, function(f) f())
      x2 <- lapply(prob, function(f) f())
      list(x1 = do.call(cbind, x1), x2 = do.call(cbind, x2), per = rep(seq_len(nt), each = q))
    }
    ctl <- list(max_iter = 60, eta = 0.01, beta1 = 0.9, beta2 = 0.999, tol = 1e-4, epsilon = 1e-8)
    G0 <- solve(crossprod(S), t(S))
    mine <- mortensor:::scoreopt_sgd(Y, draw, S, c(rep(0, m), as.vector(G0)), ctl)
    data.frame(metric = c("max_diff_G", "max_diff_d", "max_change_from_ols"),
               value = c(max(abs(mine$G - ref$G)), max(abs(mine$d - ref$d)), max(abs(ref$G - G0))))
  })

# Backtest of score-optimal reconciliation against WLS: CRPS ratios and 95%
# coverage by level of the hierarchy, the coherence gap of a forecast and the
# runtime of one fit
scoreopt_backtest <- function(xa, h, origins, spec_es) {
  t0 <- proc.time()[["elapsed"]]
  fit <- suppressWarnings(mt_fit(spec_es, xa))
  secs <- proc.time()[["elapsed"]] - t0
  fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1, predictive = TRUE)
  hs <- mortensor:::hierarchy_structure(fc$strata, h)
  E <- mortensor:::expand_exposure(fc$exposure_last, dim(fc$rates))
  D <- fc$rates * as.vector(E)
  gap <- max(vapply(seq_len(dim(D)[4]), function(j) {
    max(vapply(seq_len(dim(D)[2]), function(t) {
      full <- matrix(D[, t, hs$bottom, j], dim(D)[1]) %*% t(hs$S)
      max(abs(full / matrix(D[, t, , j], dim(D)[1]) - 1))
    }, 0))
  }, 0))
  bt <- suppressWarnings(mt_compare(list(wls = reconciled(lee_carter(), h, method = "wls"), es = spec_es), xa,
                                    origins = origins, h = 5, n_draws = 200))
  s <- mt_score_table(bt, metrics = c("crps", "cov95"), hierarchy = h, digits = 6)
  v <- function(m, l, col) s[model == m & level == l][[col]]
  lev <- c("bottom", "aggregate", "total")
  data.frame(metric = c(paste0("crps_ratio_", lev), paste0("cov95_", lev), "max_gap", "seconds_fit"),
             value = c(vapply(lev, function(l) v("es", l, "crps") / v("wls", l, "crps"), 0),
                       vapply(lev, function(l) v("es", l, "cov95"), 0), gap, secs))
}

checks$scoreopt_regions <- list(
  suite = "full",
  description = paste("Score-optimal reconciliation (reconciled(lee_carter(), method = 'energy_score'), defaults:",
                      "10 training years, 5 inner origins, horizons 1-5) on mt_example_regions() with 2 sexes x 2",
                      "regions x 5 quintiles and all aggregates (54 strata), 1995-2019: backtest from 2009, 2011",
                      "and 2013 (h = 5); CRPS relative to WLS and 95% coverage by level, largest coherence gap of",
                      "predictive draws, seconds for one fit."),
  run = function() {
    x <- mt_example_regions(years = 1995:2019, regions = 2)
    h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
    xa <- mt_aggregate_strata(x, h)
    scoreopt_backtest(xa, h, c(2009, 2011, 2013), reconciled(lee_carter(), h, method = "energy_score"))
  })

checks$scoreopt_ons <- list(
  suite = "full",
  description = paste("As scoreopt_regions on ONS England (ons_data), sex x IMD quintile with England (18 strata),",
                      "training on 6 years (5 inner origins, horizons 1-3): backtest from 2011-2014 (h = 5)."),
  run = function() {
    x <- ons_quintiles()
    h <- mt_hierarchy(sex = "both", quintile = list(total = "England"))
    xa <- mt_aggregate_strata(x, h)
    scoreopt_backtest(xa, h, 2011:2014, reconciled(lee_carter(), h, method = "energy_score", years = 6,
                                                   n_origins = 5, horizon = 3))
  })
