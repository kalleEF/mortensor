# Audit (coherence plan, item B3): do the stratum deviations of li_lee()
# move together, and does drawing their innovations jointly change the
# forecasts? Results: dev/audit-li-lee-innovations.md.
#
# li_lee() draws the innovations of the AR(1) stratum indexes independently
# (Li and Lee 2005 fit one AR(1) per population). The alternative tested
# here resamples the residual vectors of all strata of one year together,
# which keeps their correlation without estimating a covariance matrix
# (ons_data has 20 strata and 18 residual years). It is not part of the
# package: the script replaces the forecast method of li_lee() for
# specifications marked `audit_joint`.
#
# Data: ons_data (England, sex x deprivation decile, 5-year age groups,
# 2001-2019, common factor by sex) and HMD France (demography::fr.mort,
# women and men, ages 0-99, 1950-2006).

pkgload::load_all(".", quiet = TRUE)
library(data.table)
set.seed(1)

# forecast method of li_lee() with joint resampling of the residuals of the
# stratum indexes (centred, scaled to the variance of the fit) for
# specifications with params$audit_joint = TRUE; the package's method
# otherwise
sim_package <- simulate_model.mt_spec_multipop
sim_audit <- function(spec, fit, h, n_draws, ...) {
  if (!isTRUE(spec$params$audit_joint)) return(sim_package(spec, fit, h, n_draws, ...))
  na <- length(fit$commons[[1]]$a)
  ns <- length(fit$strata)
  kps <- lapply(fit$commons, function(cm) {
    rw_paths(cm$k[length(cm$k)], cm$drift, matrix(cm$dk - cm$drift, nrow = 1), h, n_draws, "gaussian",
             sigma = cm$sigma_k, exp_out = FALSE)
  })
  pool <- lapply(fit$strata, function(st) lapply(seq_along(st$ar), function(f) {
    k <- st$k[, f]
    e <- k[-1] - st$ar[[f]]$phi * k[-length(k)]
    e <- e - mean(e)
    e * st$ar[[f]]$sigma / sqrt(mean(e^2))
  }))
  yr <- matrix(sample.int(length(pool[[1]][[1]]), h * n_draws, replace = TRUE), h, n_draws)
  out <- array(NA_real_, c(na, h, ns, n_draws))
  for (j in seq_len(n_draws)) {
    for (s in seq_len(ns)) {
      cm <- fit$commons[[fit$group[s]]]
      st <- fit$strata[[s]]
      kk <- vapply(seq_along(st$ar), function(f) {
        a <- st$ar[[f]]
        eps <- pool[[s]][[f]][yr[, j]]
        v <- numeric(h)
        prev <- a$last
        for (t in seq_len(h)) {
          prev <- a$phi * prev + eps[t]
          v[t] <- prev
        }
        v
      }, numeric(h))
      lp <- outer(cm$b, kps[[fit$group[s]]][1, , j]) + st$a + st$b %*% t(matrix(kk, h))
      if (!is.null(fit$jump)) lp <- lp + fit$jump[, s]
      out[, , s, j] <- exp(lp)
    }
  }
  out
}
ns_env <- asNamespace("mortensor")
unlockBinding("simulate_model.mt_spec_multipop", ns_env)
assign("simulate_model.mt_spec_multipop", sim_audit, envir = ns_env)
registerS3method("simulate_model", "mt_spec_multipop", sim_audit, envir = ns_env)
spec_ll <- function(innovations, ...) {
  sp <- li_lee(...)
  if (innovations == "joint") {
    sp$params$audit_joint <- TRUE
    sp$label <- paste(sp$label, "(joint innovations)")
  }
  sp
}

ons <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile"), age = "age_group",
               value = c(deaths = "all_causes", population = "population"))
ons <- mt_as_exposure(ons, type = "midyear")
fr <- suppressWarnings(as_mt_data(demography::fr.mort, years = 1950:2006))
fr <- mt_regroup(fr, age = 0:100)

# correlation of the AR(1) residuals of the stratum indexes, signs aligned so
# that a positive index raises the rates on average over age
resid_cor <- function(fit) {
  E <- vapply(fit$fits$strata, function(st) {
    k <- st$k[, 1]
    sign(sum(st$b[, 1])) * (k[-1] - st$ar[[1]]$phi * k[-length(k)])
  }, numeric(nrow(fit$fits$strata[[1]]$k) - 1L))
  R <- stats::cor(E)
  off <- R[upper.tri(R)]
  n <- nrow(E)
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  list(n_years = n, n_strata = ncol(E), mean_r = mean(off), mean_abs_r = mean(abs(off)),
       share_beyond_2se = mean(abs(off) > 2 / sqrt(n)), first_eigen_share = ev[1] / sum(ev), R = R)
}

# spread of forecasts at horizon h: life expectancy at the first age of each
# stratum, the largest gap, and the total deaths of all strata (last
# exposures)
spread <- function(fit, h, contrast) {
  out <- list()
  for (inn in c("independent", "joint")) {
    fit$spec$params$audit_joint <- inn == "joint"
    fc <- mt_forecast(fit, h = h, n_draws = 2000, seed = 1)
    le <- life_expectancy(fc, probs = c(0.05, 0.95))[period == max(period)]
    cf <- mt_contrast(fc, var = contrast[["var"]], ref = contrast[["ref"]], probs = c(0.05, 0.95))
    cf <- cf[period == max(period)]
    E <- fit$exposure_last
    tot <- apply(fc$rates[, h, , , drop = FALSE] * as.vector(E[, 1, ]), 4, sum)
    out[[inn]] <- c(le_width_mean = mean(le$q0.95 - le$q0.05),
                    gap_width_mean = mean(cf$q0.95 - cf$q0.05),
                    total_deaths_width_rel = unname(diff(stats::quantile(tot, c(0.05, 0.95))) / stats::median(tot)))
  }
  do.call(rbind, out)
}

backtest <- function(data, spec_fun, origins, h, contrast) {
  bt <- mt_compare(list(independent = spec_fun("independent"), joint = spec_fun("joint")), data, origins = origins,
                   h = h, n_draws = 500, quantities = c("life_expectancy", "asr"), contrast = contrast,
                   multivariate = TRUE)
  list(joint = mt_score_table(bt, metrics = c("energy", "variogram"),
                              targets = c("life_expectancy_joint", "asr_joint")),
       gap = mt_score_table(bt, metrics = c("crps", "cov80", "cov95"), targets = c("life_expectancy_gap", "asr_gap")),
       total_rate = mt_score_table(bt, metrics = c("crps", "cov95"), targets = "rate"))
}

res <- list()
# England by sex and deprivation decile
f_ons <- mt_fit(li_lee(common_by = "sex"), ons)
res$ons_cor <- resid_cor(f_ons)
res$ons_spread <- spread(f_ons, 10, c(var = "imd_decile", ref = "10"))
res$ons_backtest <- backtest(ons, function(i) spec_ll(i, common_by = "sex"), origins = 2009:2014, h = 5,
                             contrast = c(var = "imd_decile", ref = "10"))
# France, women and men
f_fr <- mt_fit(li_lee(), fr)
res$fr_cor <- resid_cor(f_fr)
res$fr_spread <- spread(f_fr, 20, c(var = "sex", ref = "male"))
res$fr_backtest <- backtest(fr, function(i) spec_ll(i), origins = 1987:1996, h = 10,
                            contrast = c(var = "sex", ref = "male"))

str(res[c("ons_cor", "fr_cor")], max.level = 1)
print(lapply(res[c("ons_cor", "fr_cor")], function(r) r[setdiff(names(r), "R")]))
print(res[c("ons_spread", "fr_spread")])
print(res[c("ons_backtest", "fr_backtest")])
saveRDS(res, file.path(tempdir(), "audit-li-lee-innovations.rds"))
