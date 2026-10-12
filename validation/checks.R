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

checks$coverage_tucker3_coda <- list(
  suite = "full",
  description = "As coverage_lee_carter for tucker3(transform = \"coda\").",
  run = function() coverage_check(tucker3(transform = "coda"), 10))

checks$coverage_ensemble <- list(
  suite = "full",
  description = "As coverage_lee_carter for mt_ensemble(lee_carter(), tucker3(), rwd()).",
  run = function() coverage_check(mt_ensemble(lee_carter(), tucker3(), rwd()), 10))

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
