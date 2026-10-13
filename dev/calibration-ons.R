# Interval calibration of random walks, Lee-Carter and ensembles (see dev/calibration.md).
# Run from the package root: Rscript dev/calibration-ons.R
suppressMessages(devtools::load_all(quiet = TRUE))
library(data.table)
x <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile"), age = "age_group",
             value = c(deaths = "all_causes", population = "population"))
x <- mt_group_levels(mt_as_exposure(x, type = "midyear"), "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)),
                     name = "quintile")
specs <- list(rwd = rwd(), rwdn = rwd(observation_noise = TRUE), ensn = ensemble(lc = lee_carter(), rw = rwd(observation_noise = TRUE)),
              lc = lee_carter(),
              ens = ensemble(lc = lee_carter(), rw = rwd()),
              ens_crps = ensemble(lc = lee_carter(), rw = rwd(observation_noise = TRUE), weights = "crps"))
res <- rbindlist(lapply(names(specs), function(n) {
  t0 <- Sys.time()
  b <- mt_backtest(specs[[n]], x, origins = 2006:2014, h = 5, n_draws = 300)
  s <- b[, .(crps = mean(crps), cov50 = mean(cov50), cov80 = mean(cov80), cov95 = mean(cov95), pit_sd = sd(pit))]
  s[, `:=`(model = n, secs = round(as.numeric(Sys.time() - t0, units = "secs")))]
  hz <- b[, .(cov95 = mean(cov95)), by = horizon]
  s[, cov95_h1 := hz[horizon == 1, cov95]][, cov95_h5 := hz[horizon == 5, cov95]]
  s
}))
print(res, digits = 3)

