# Comparison of apc_causes() with bamp_multicause() of the bamp development
# branch "coherent_forecasting" (commit 4363a78), its reference
# implementation. Results: dev/apc-causes-check.md.
#
# Install the branch into a separate library (see dev/check-coherent-apc.R)
# and set BAMP_CF_LIB. Data: ons_data, England, women, ages 50-89 in 5-year
# groups, 2001-2019: deaths from cancer (all cancers), circulatory causes
# (ischaemic heart disease, stroke) and the rest of all causes.
#
# Both fit stick-breaking binomial-logit APC models of the cause shares with
# Wishart-coupled period and cohort effects (RW1) and the same priors. The
# branch fits the all-cause total with bamp() and the shares with Polya-Gamma
# Gibbs steps; mortensor fits both with its own sampler. Shares should agree
# up to Monte Carlo error.

lib <- Sys.getenv("BAMP_CF_LIB")
stopifnot(nzchar(lib))
library(bamp, lib.loc = lib)
stopifnot(packageVersion("bamp", lib.loc = lib) == "2.2.0", exists("bamp_multicause"))
pkgload::load_all(".", quiet = TRUE)
library(data.table)

ages <- c("50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89")
d <- as.data.table(ons_data)[sex == "female" & age_group %in% ages]
cancer_cols <- c("breast_cancer", "colorectal_cancer", "lung_cancer", "prostate_cancer", "all_other_cancers")
d[, cancer := Reduce(`+`, .SD), .SDcols = cancer_cols]
d[, circulatory := chronic_ischaemic_heart_disease + stroke_incl_transient_ischaemic_attack]
d <- d[, list(all = sum(all_causes), cancer = sum(cancer), circulatory = sum(circulatory),
              pop = sum(population)), by = c("year", "age_group")]
d[, other := all - cancer - circulatory]
mat <- function(v) as.matrix(dcast(d, year ~ age_group, value.var = v)[, ages, with = FALSE])
cases <- list(cancer = mat("cancer"), circulatory = mat("circulatory"), other = mat("other"))
pop <- mat("pop")
I <- length(ages)
J <- nrow(pop)
h <- 10

t0 <- proc.time()[["elapsed"]]
fb <- bamp_multicause(cases, pop, age = "rw1", period = "rw1", cohort = "rw1", periods_per_agegroup = 5,
                      mcmc = list(iterations = 8000, burn_in = 2000, thin = 4), seed = 1)
tb <- proc.time()[["elapsed"]] - t0
pb <- predict_multicause(fb, periods = h)

# mortensor: the same data as an mt_data object of one sex with four strata
long <- melt(d, id.vars = c("year", "age_group", "pop"), measure.vars = c("all", "cancer", "circulatory", "other"),
             variable.name = "cause", value.name = "deaths")
x <- mt_data(deaths = long, population = long, by = "cause", age = "age_group",
             value = c(deaths = "deaths", population = "pop"))
x <- mt_as_exposure(x, type = "midyear")
t0 <- proc.time()[["elapsed"]]
fm <- mt_fit(apc_causes(age = "rw1", residual = FALSE, seed = 1,
                        mcmc = mcmc_control(iterations = 5000, burn_in = 1000, thin = 5)), x, cores = 4)
tm <- proc.time()[["elapsed"]] - t0
fcm <- mt_forecast(fm, h = h, n_draws = 2000, seed = 2)

lab <- fm$strata$label
share_m <- function(cause, years) {
  r <- fcm$rates
  apply(r[, years, match(cause, lab), , drop = FALSE] / r[, years, match("all", lab), , drop = FALSE], 1, median)
}
share_b <- function(cause, k) {
  z <- pb[[cause]]$samples$rate[J + k, , ] / pb$total$samples$rate[J + k, , ]
  apply(z, 1, median)
}
fitted_share_m <- function(cause) fitted(fm)[, J, match(cause, lab)] / fitted(fm)[, J, match("all", lab)]
fitted_share_b <- function(cause) {
  z <- pb[[cause]]$samples$rate[J, , ] / pb$total$samples$rate[J, , ]
  apply(z, 1, mean)
}
res <- list(
  seconds = c(branch = tb, mortensor = tm),
  share_last_year = rbind(branch = sapply(names(cases), fitted_share_b),
                          mortensor = sapply(names(cases), function(c) {
                            # fitted() is on the grid (5-year groups here)
                            fitted_share_m(c)
                          })),
  share_h = lapply(names(cases), function(c) rbind(branch = share_b(c, h), mortensor = share_m(c, h))),
  cor_period_branch = pb$cor_omega,
  cor_period_mortensor = mt_parameters(fm)[parameter == "period_correlation"],
  rhat = fm$diagnostics)
names(res$share_h) <- names(cases)
print(res)
saveRDS(res, file.path(tempdir(), "apc-causes-check.rds"))
