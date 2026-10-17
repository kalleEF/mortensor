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
  description = paste("As coverage_lee_carter for rwd(). The random walk treats the iid noise of the",
                      "generator as persistent, so it over-covers; the threshold only guards against",
                      "under-coverage."),
  run = function() coverage_check(rwd(), 6))

checks$coverage_rwd_noise <- list(
  suite = "quick",
  description = paste("As coverage_lee_carter for rwd(observation_noise = TRUE), the random walk plus noise.",
                      "It separates the iid noise of the generator from the trend, so its intervals should be",
                      "close to nominal at every level (see dev/calibration.md)."),
  run = function() coverage_check(rwd(observation_noise = TRUE), 6))

checks$coverage_tucker3_coda <- list(
  suite = "full",
  description = "As coverage_lee_carter for tucker3(transform = \"coda\").",
  run = function() coverage_check(tucker3(transform = "coda"), 10))

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
                      "smooth); Lee-Carter anchored to the observed deaths of the last 3 years. ONS quintiles",
                      "(origins 2008-2014, h = 5): Li-Lee with and without anchoring."),
  run = function() {
    ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
    ew <- mt_set_grid(ew, age = single_ages(30, 100, open = FALSE))
    ew <- mt_regroup(ew, age = single_ages(30, 100, open = FALSE))
    o <- 1996:2001
    gam2d <- mt_prep(step_smooth(dims = 2, method = "gam"))
    lc <- mt_backtest(lee_carter(), ew, origins = o, h = 10, n_draws = 200)
    sm <- mt_backtest(lee_carter(), ew, origins = o, h = 10, n_draws = 200, prep = gam2d)
    smc <- mt_backtest(recalibrated(lee_carter()), ew, origins = o, h = 10, n_draws = 200, prep = gam2d)
    an <- mt_backtest(anchored(lee_carter(), years = 3, decay = 10), ew, origins = o, h = 10, n_draws = 200)
    x <- ons_quintiles()
    ll <- mt_backtest(li_lee(common_by = "sex"), x, origins = 2008:2014, h = 5, n_draws = 200)
    lla <- mt_backtest(anchored(li_lee(common_by = "sex"), years = 3, decay = 10), x, origins = 2008:2014, h = 5,
                       n_draws = 200)
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

checks$cause_shares_coherence <- list(
  suite = "quick",
  description = paste("cause_shares(lee_carter()) with step_smooth() on mt_example_causes(): largest relative gap",
                      "between all-cause",
                      "draws and the sum of cause draws, and CRPS by cause in a backtest from 2012 (h = 5)."),
  run = function() {
    x <- mt_example_causes(years = 1995:2017)
    prep <- mt_prep(step_smooth())
    fit <- mt_fit(cause_shares(total = lee_carter()), mt_subset(x, periods = c(NA, 2012)), prep = prep)
    fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
    st <- fc$strata
    gap <- max(vapply(c("female", "male"), function(sx) {
      i <- st$stratum[st$sex == sx & st$cause == "all"]
      j <- st$stratum[st$sex == sx & st$cause != "all"]
      max(abs(apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum) / fc$rates[, , i, ] - 1))
    }, 0))
    bt <- mt_backtest(cause_shares(total = lee_carter()), x, origins = 2012, h = 5, n_draws = 100, prep = prep)
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
                      "2014 (h = 5) of the ensemble of li_lee(common_by = 'sex') and rwd(observation_noise = TRUE),",
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
    sp <- reconciled(ensemble(ll = li_lee(common_by = "sex"), rw = rwd(observation_noise = TRUE)), h,
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
