# Interval calibration on the simulation benchmark and HMD data (see dev/calibration.md).
# Run from the package root: Rscript dev/calibration-benchmarks.R
suppressMessages(devtools::load_all(quiet = TRUE))
library(data.table)
# Fits `spec` on 30 years of simulate_lc_dgp() data and scores the draws of
# the next `h` years, over `n_rep` replicates
benchmark_coverage <- function(spec, n_rep = 10, h = 10, n_draws = 300, n_train = 30) {
  bt <- rbindlist(lapply(seq_len(n_rep), function(r) {
    x <- simulate_lc_dgp(n_t = n_train + h, seed = 1000 + r)
    mt_backtest(spec, x, origins = 2000 + n_train, h = h, n_draws = n_draws, seed = r)
  }))
  by_h <- mt_score(bt, by = "horizon")
  attr(by_h, "summary") <- list(cov95_mean = mean(by_h$cov95), cov80_mean = mean(by_h$cov80),
                                cov50_mean = mean(by_h$cov50), crps = mean(bt$crps))
  by_h
}
specs <- list(rwd = rwd(), rwdn = rwd(observation_noise = TRUE), rwdn_g = rwd(observation_noise = TRUE, innovations = "gaussian"), lc = lee_carter(), ensn = ensemble(lc = lee_carter(), rw = rwd(observation_noise = TRUE)),
              ens = ensemble(lc = lee_carter(), rw = rwd()))
cat("simulation benchmark (LC DGP, 10 reps, h = 10)\n")
sim <- rbindlist(lapply(names(specs), function(n) {
  r <- attr(benchmark_coverage(specs[[n]], n_rep = 10, n_draws = 300), "summary")
  data.table(model = n, crps = r$crps, cov50 = r$cov50_mean, cov80 = r$cov80_mean, cov95 = r$cov95_mean)
}))
print(sim, digits = 3)
ew <- as_mt_data(pkg_data("EWMaleData", "StMoMo"))
ew <- mt_subset(ew, periods = c(1971, 2011))
ew <- mt_set_grid(ew, age = single_ages(0, 100, open = FALSE))
cat("England & Wales men (HMD via StMoMo), ages 0-100, origins 1996-2001, h = 10\n")
hmd <- rbindlist(lapply(names(specs), function(n) {
  b <- mt_backtest(specs[[n]], ew, origins = 1996:2001, h = 10, n_draws = 300)
  b[, .(model = n, crps = mean(crps), cov50 = mean(cov50), cov80 = mean(cov80), cov95 = mean(cov95))]
}))
print(hmd, digits = 3)
fr <- as_mt_data(pkg_data("fr.mort", "demography"), years = 1950:2006)
fr <- mt_set_grid(fr, age = single_ages(0, 100, open = FALSE))
cat("France (HMD via demography), ages 0-100, origins 1991-1996, h = 10\n")
frr <- rbindlist(lapply(names(specs), function(n) {
  b <- mt_backtest(specs[[n]], fr, origins = 1991:1996, h = 10, n_draws = 300)
  b[, .(model = n, crps = mean(crps), cov50 = mean(cov50), cov80 = mean(cov80), cov95 = mean(cov95))]
}))
print(frr, digits = 3)

