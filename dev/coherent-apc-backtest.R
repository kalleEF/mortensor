# Backtest of apc_bayes(coherent = TRUE) against li_lee() and separate
# Lee-Carter models on England, women and men (ons_data pooled over the
# deprivation deciles, ages 30-89 in 5-year groups, 2001-2019): gaps between
# the sexes and joint scores (coherence plan, A1). Results:
# dev/coherent-apc-check.md.

pkgload::load_all(".", quiet = TRUE)
library(data.table)
set.seed(1)
ages <- c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84",
          "85-89")
d <- as.data.table(ons_data)[age_group %in% ages,
                             list(all_causes = sum(all_causes), population = sum(population)),
                             by = c("year", "sex", "age_group")]
x <- mt_data(deaths = d, population = d, by = "sex", age = "age_group",
             value = c(deaths = "all_causes", population = "population"))
x <- mt_as_exposure(x, type = "midyear")
specs <- list(lee_carter = lee_carter(), li_lee = li_lee(),
              coherent_apc = apc_bayes(coherent = TRUE, seed = 1, mcmc = mcmc_control("fast")))
t0 <- proc.time()[["elapsed"]]
bt <- mt_compare(specs, x, origins = 2009:2014, h = 5, n_draws = 500, quantities = c("life_expectancy", "asr"),
                 contrast = c(var = "sex", ref = "male"), multivariate = TRUE, cores = 1)
secs <- proc.time()[["elapsed"]] - t0
res <- list(
  seconds = secs,
  rates = mt_score_table(bt, metrics = c("crps", "cov95")),
  joint = mt_score_table(bt, metrics = c("energy", "variogram"), targets = c("life_expectancy_joint", "asr_joint")),
  gap = mt_score_table(bt, metrics = c("crps", "cov80", "cov95"), targets = c("life_expectancy_gap", "asr_gap")))
print(res)
saveRDS(res, file.path(tempdir(), "coherent-apc-backtest.rds"))
