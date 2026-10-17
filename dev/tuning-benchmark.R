# Benchmark of tuning and calibration (2026-10-17): what choosing settings by
# forecast accuracy (tuned(), the hold-out workflow of demography::fdm()
# users) and recalibrating the spread (recalibrated()) change in outer
# rolling-origin backtests. Every tuning and calibration uses only the
# training data of its outer origin.
#
#   Rscript dev/tuning-benchmark.R            (about 40 minutes)
#
# Writes dev/tuning-benchmark.rds; the results are summarised in
# dev/tuning.md.

suppressMessages({
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(quiet = TRUE)
  } else {
    library(mortensor)
  }
  library(data.table)
})

log_step <- function(...) message(format(Sys.time(), "%H:%M:%S "), ...)

ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
ew <- mt_set_grid(ew, age = single_ages(30, 100, open = FALSE))
ew <- mt_regroup(ew, age = single_ages(30, 100, open = FALSE))
fr <- suppressWarnings(as_mt_data(mortensor:::pkg_data("fr.mort", "demography"), series = "male"))
fr <- mt_subset(fr, ages = c(0, 96), periods = c(1950, 2006))
fr <- mt_set_grid(fr, age = single_ages(0, 95, open = FALSE))

# fdm: number of components and knots of the smoother over age, simplest first
fdm_fun <- function(order, k) hyndman_ullah(order = order, smooth_args = list(k = k))
fdm_grid <- expand.grid(k = c(10, 30), order = c(1, 2, 3, 6))[, c("order", "k")]
# Lee-Carter: fitting window (the last `window` training years)
win_fun <- function(window) list(spec = lee_carter(), prep = mt_prep(step_window(window)))
win_grid <- data.frame(window = c(15, 20, 30, 60))

holdout <- function(f, g, metric = "rmse") {
  tuned(f, g, inner_origins = 1, inner_h = 10, metric = metric, n_draws = 100, label = "tuned (hold-out)")
}
specs <- list(
  "Lee-Carter" = lee_carter(),
  "Lee-Carter, calibrated" = recalibrated(lee_carter()),
  "Lee-Carter, window tuned (hold-out 10, RMSE)" = holdout(win_fun, win_grid),
  "fdm (6 components, smoothed)" = hyndman_ullah(),
  "fdm (6 components, unsmoothed)" = hyndman_ullah(smooth = FALSE),
  "fdm tuned (hold-out 10, RMSE)" = holdout(fdm_fun, fdm_grid),
  "fdm tuned (3 origins, CRPS, one-SE)" = tuned(fdm_fun, fdm_grid, inner_origins = 3, inner_h = 5,
                                                select = "one_se", n_draws = 100, label = "tuned (rolling)"),
  "fdm tuned (hold-out 10, RMSE), calibrated" = recalibrated(holdout(fdm_fun, fdm_grid))
)

data_sets <- list(
  ew = list(label = "England & Wales men (HMD via StMoMo), ages 30-100, 1961-2011", x = ew, origins = 1996:2001),
  fr = list(label = "France men (HMD via demography), ages 0-95, 1950-2006", x = fr, origins = c(1986, 1991, 1996))
)

res <- list()
for (dn in names(data_sets)) {
  ds <- data_sets[[dn]]
  for (sn in names(specs)) {
    log_step(dn, ": ", sn)
    t0 <- proc.time()[["elapsed"]]
    bt <- tryCatch(mt_backtest(specs[[sn]], ds$x, origins = ds$origins, h = 10, n_draws = 200,
                               quantities = "life_expectancy"), error = function(e) e)
    if (inherits(bt, "error")) {
      log_step("  failed: ", conditionMessage(bt))
      next
    }
    bt[, `:=`(spec = sn, data = dn, seconds = proc.time()[["elapsed"]] - t0)]
    res[[length(res) + 1L]] <- bt
  }
}
bt <- rbindlist(res, fill = TRUE)

# settings chosen at each outer origin by the hold-out tuning
choices <- rbindlist(lapply(names(data_sets), function(dn) {
  ds <- data_sets[[dn]]
  rbindlist(lapply(ds$origins, function(o) {
    tr <- mt_window(ds$x, to = o)
    a <- mt_tune(fdm_fun, fdm_grid, tr, inner_origins = 1, inner_h = 10, metric = "rmse")$chosen
    b <- mt_tune(win_fun, win_grid, tr, inner_origins = 1, inner_h = 10, metric = "rmse")$chosen
    data.table(data = dn, origin = o, fdm_order = a$order, fdm_k = a$k, lc_window = b$window)
  }))
}))

score <- function(b, tgt) {
  b[target == tgt, list(crps = mean(crps), rmse = sqrt(mean(ae^2)), cov80 = mean(cov80), cov95 = mean(cov95),
                        n = .N), by = list(data, spec)]
}
summary_rates <- score(bt, "rate")
summary_e <- score(bt, "life_expectancy")
by_h <- bt[target == "rate", list(crps = mean(crps), cov95 = mean(cov95)), by = list(data, spec, horizon)]
saveRDS(list(rates = summary_rates, life_expectancy = summary_e, by_horizon = by_h, choices = choices,
             seconds = bt[, list(seconds = max(seconds)), by = list(data, spec)], date = Sys.Date()),
        "dev/tuning-benchmark.rds")
print(summary_rates)
print(summary_e)
print(choices)
