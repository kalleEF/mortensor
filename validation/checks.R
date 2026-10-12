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
