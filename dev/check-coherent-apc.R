# Comparison of the native coherent APC sampler of mortensor
# (apc_bayes(coherent = TRUE)) with bamp_coherent() of the bamp development
# branch "coherent_forecasting" (commit 4363a78), its reference
# implementation. Results: dev/coherent-apc-check.md.
#
# The branch is not on CRAN: install it into a separate library, e.g.
#   git clone -b coherent_forecasting https://github.com/ChristK/bamp bamp-cf
#   R CMD INSTALL -l <lib> bamp-cf
# and set BAMP_CF_LIB=<lib>. Data: ons_data (England, all deprivation deciles
# pooled), women and men, ages 30-89 in 5-year groups, 2001-2019.
#
# The two parameterisations differ in one respect: for two strata the
# branch writes the deviations as +delta and -delta with sum_j delta_j = 0
# (the levels a_s carry the mean gap), mortensor lets them sum to zero over
# the strata without a time constraint (Cairns et al. 2011, Section 4.5).
# The likelihood, the other priors and the forecast processes are the same,
# so fitted rates, gaps and forecasts should agree up to Monte Carlo error
# and the slightly different prior on the mean gap. The branch's innovation
# precision of delta is twice mortensor's lambda_d (d_1 = delta = -d_2), so
# the branch's Gamma(1, 0.05) prior of that precision corresponds to
# Gamma(1, 0.1) here; mortensor's default, Gamma(1, 5e-4), is fitted too.

lib <- Sys.getenv("BAMP_CF_LIB")
stopifnot(nzchar(lib))
library(bamp, lib.loc = lib)
stopifnot(packageVersion("bamp", lib.loc = lib) == "2.2.0", exists("bamp_coherent"))
pkgload::load_all(".", quiet = TRUE)
library(data.table)

d <- as.data.table(ons_data)
ages <- c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84",
          "85-89")
d <- d[age_group %in% ages, list(deaths = sum(all_causes), pop = sum(population)), by = c("sex", "year", "age_group")]
mat <- function(sx, v) {
  w <- dcast(d[sex == sx], year ~ age_group, value.var = v)
  as.matrix(w[, ages, with = FALSE])
}
cases <- list(female = mat("female", "deaths"), male = mat("male", "deaths"))
pop <- list(female = mat("female", "pop"), male = mat("male", "pop"))
I <- length(ages)
J <- nrow(cases$female)
h <- 10

# branch: S = 2 path, AR(1) deviation with sampled rho, RW1 age/period/cohort
t0 <- proc.time()[["elapsed"]]
fb <- bamp_coherent(cases, pop, age = "rw1", period = "rw1", cohort = "rw1", periods_per_agegroup = 5,
                    deviation = "ar1", mcmc = list(iterations = 20000, burn_in = 5000, thin = 10), seed = 1)
tb <- proc.time()[["elapsed"]] - t0
pb <- predict_coherent(fb, periods = h)

# mortensor: arrays [age group, period, stratum]
D <- array(c(t(cases$female), t(cases$male)), c(I, J, 2))
E <- array(c(t(pop$female), t(pop$male)), c(I, J, 2))
t0 <- proc.time()[["elapsed"]]
mc <- list(options = list(number_of_iterations = 5000, burn_in = 1000, step = 5))
hy <- coh_hyper()
hy$dev <- c(1, 0.1)
fm <- with_seed(1, coh_fit(D, E, 5, c(1L, 1L, 1L), FALSE, NULL, mc, cores = 4L, max_refits = 1L, hyper = hy))
tm <- proc.time()[["elapsed"]] - t0
sim <- with_seed(2, coh_simulate(fm, h, 2000))
# mortensor's default prior of the deviation precision
fd <- with_seed(1, coh_fit(D, E, 5, c(1L, 1L, 1L), FALSE, NULL, mc, cores = 4L, max_refits = 1L))
sim_d <- with_seed(2, coh_simulate(fd, h, 2000))

# fitted rates: [age, period, stratum]
fit_b <- array(0, c(I, J, 2))
for (s in 1:2) fit_b[, , s] <- t(pb[[c("female", "male")[s]]]$rate[2, seq_len(J), ])
rel_fit <- abs(fm$fitted / fit_b - 1)
# observed in-sample gap (mean over ages of the logit difference) per year
gap_in <- function(P) colMeans(qlogis(P[, , 1]) - qlogis(P[, , 2]))
# forecast: medians and 90% intervals of the rates and of the gap at h
qb <- function(s) pb[[s]]$samples$rate[J + h, , ]                       # [age, draw]
rate_b <- list(female = qb("female"), male = qb("male"))
rate_m <- list(female = sim[, h, 1, ], male = sim[, h, 2, ])
gap_draws <- function(r) colMeans(qlogis(r$female) - qlogis(r$male))
summ <- function(v) stats::quantile(v, c(0.05, 0.5, 0.95))
res <- list(
  seconds = c(branch = tb, mortensor = tm),
  fitted_rel_diff = summary(as.vector(rel_fit)),
  gap_in_sample = rbind(branch = gap_in(fit_b), mortensor = gap_in(fm$fitted)),
  rho = rbind(branch = summ(fb$samples$rho), mortensor = summ(fm$h[, "rho"])),
  gap_h = rbind(branch = summ(gap_draws(rate_b)), mortensor = summ(gap_draws(rate_m)),
                mortensor_default_prior = summ(gap_draws(list(female = sim_d[, h, 1, ], male = sim_d[, h, 2, ])))),
  default_prior = rbind(rho = summ(fd$h[, "rho"]), lambda_d = summ(fd$h[, "lambda_d"])),
  rate_h_median_rel_diff = rbind(
    female = range(apply(rate_m$female, 1, median) / apply(rate_b$female, 1, median) - 1),
    male = range(apply(rate_m$male, 1, median) / apply(rate_b$male, 1, median) - 1)),
  rate_h_width_ratio = rbind(
    female = range(apply(rate_m$female, 1, function(v) diff(quantile(v, c(0.05, 0.95)))) /
                     apply(rate_b$female, 1, function(v) diff(quantile(v, c(0.05, 0.95))))),
    male = range(apply(rate_m$male, 1, function(v) diff(quantile(v, c(0.05, 0.95)))) /
                   apply(rate_b$male, 1, function(v) diff(quantile(v, c(0.05, 0.95)))))),
  # posterior medians and 90% intervals of the hyperparameters; the
  # branch's precision of delta is compared as lambda_d / 2
  hyper = rbind(
    branch_kappa = summ(fb$samples$kappa_theta), mortensor_kappa = summ(fm$h[, "kappa"]),
    branch_lambda = summ(fb$samples$lambda_phi), mortensor_lambda = summ(fm$h[, "lambda"]),
    branch_nu = summ(fb$samples$nu_psi), mortensor_nu = summ(fm$h[, "nu"]),
    branch_lambda_d_half = summ(fb$samples$lambda_d / 2), mortensor_lambda_d = summ(fm$h[, "lambda_d"])),
  # posterior sd of the last period effect and of the last cohort effect
  sd_last = rbind(branch = c(phi = sd(fb$samples$phi[, J]), psi = sd(fb$samples$psi[, ncol(fb$samples$psi)])),
                  mortensor = c(phi = sd(fm$x[, fm$setup$idx$ph[J]]), psi = sd(fm$x[, fm$setup$idx$ps[fm$setup$K]]))),
  mortensor_rhat = fm$psrf, mortensor_accept = fm$accept)
print(res)
saveRDS(res, file.path(tempdir(), "coherent-apc-check.rds"))
