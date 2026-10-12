# Capture point forecasts of every mortensor 0.1.1 model with fixed hyperparameters.
# Used as the regression baseline for the 0.2 rewrite (implementation plan WP P0.5).
# Run from the package root with mortensor 0.1.1 installed:
#   Rscript dev/baseline/capture_baseline.R
suppressPackageStartupMessages(library(mortensor))
set.seed(20261012)
data(example_mortality)
sm <- suppressMessages(smooth_mortality(example_mortality, method = "gam", covid_years = 2020L))
calib <- list(R_A = 3, R_B = 3, R_C = 2, ts_method = "ets")
models <- c("Tucker3-CoDa", "Tucker3-Direct", "Lee-Carter", "FDM", "StMoMo-APC",
            "RWD", "ARIMA-perSeries", "ETS-Damped")
h <- 10L
out <- list()
for (m in models) {
  res <- tryCatch(
    suppressWarnings(suppressMessages(utils::capture.output(
      fc <- dispatch_forecast(m, sm$mx, sm$pop, sm$years, sm$ages, sm$groups, h,
                              calib, calib, width = 5L, calib_all = NULL,
                              return_fitted = FALSE)))),
    error = function(e) e)
  out[[m]] <- if (inherits(res, "error")) list(error = conditionMessage(res)) else
    list(mx_forecast = fc$mx_forecast)
  cat(m, if (inherits(res, "error")) paste("ERROR:", conditionMessage(res)) else "ok", "\n")
}
baseline <- list(
  created = format(Sys.time(), tz = "UTC"),
  mortensor_version = as.character(utils::packageVersion("mortensor")),
  r_version = R.version.string,
  settings = list(smoothing = "gam", covid_years = 2020L, h = h,
                  tucker3 = calib, hyperparameters = "0.1.1 defaults in dispatch_forecast()"),
  ages = sm$ages, years = sm$years, groups = sm$groups,
  smoothed_input = sm$mx,
  forecasts = out)
saveRDS(baseline, "dev/baseline/baseline_0.1.rds", version = 3)
cat("saved dev/baseline/baseline_0.1.rds\n")
