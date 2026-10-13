# Comparison of mortensor with the original packages (demography, bamp,
# StMoMo): fidelity of the wrappers and forecast accuracy of the native
# workflows against mortensor pipelines. Results are stored in
# vignettes/comparison-results.rds and shown in vignettes/comparison.Rmd.
#
# Run from the package root (about one hour, most of it bamp chains):
#   Rscript dev/compare-packages.R
# Set COMPARE_QUICK=true for a short run (fewer origins and draws).

suppressMessages(devtools::load_all(quiet = TRUE))
suppressMessages({
  library(data.table)
  library(demography)
  library(bamp)
  library(StMoMo)
})
quick <- identical(Sys.getenv("COMPARE_QUICK"), "true")
n_draws <- if (quick) 100 else 300
log_step <- function(...) message(format(Sys.time(), "%H:%M:%S "), sprintf(...))
ns <- asNamespace("mortensor")

# ---- Native workflows as mortensor specifications ---------------------------
# Each adapter calls the original package as its documentation recommends and
# returns that package's own predictive distribution, so mt_backtest() scores
# all forecasts in the same way. smooth = FALSE: no noise is added.

native_fdm <- function(order = 6) {
  new_mt_spec("native_fdm", list(order = order), min_periods = 8, label = "demography fdm (native)", smooth = FALSE)
}
registerS3method("fit_model", "mt_spec_native_fdm", function(spec, D, E, ages, periods, ...) {
  rates <- pmax(D, 0.5) / E
  dd <- demography::demogdata(rates, pop = E, ages = ages$lower, years = periods, type = "mortality",
                              label = "native", name = "s")
  # smooth.demogdata() is the recommended first step for single years of age;
  # it re-grids grouped ages, so grouped data are modelled unsmoothed
  single <- all(utils::head(ages$upper - ages$lower, -1) == 1)
  if (single) dd <- suppressWarnings(demography::smooth.demogdata(dd))
  list(fit = demography::fdm(dd, series = "s", order = min(spec$params$order, nrow(D) - 1)))
}, envir = ns)
registerS3method("simulate_model", "mt_spec_native_fdm", function(spec, fit, h, n_draws, ...) {
  # forecast.fdm: point forecast and total variance (model, error and
  # observational) of the log rates; draws from the normal it implies
  fc <- forecast::forecast(fit$fit, h = h)
  mu <- log(fc$rate$s)
  sd <- sqrt(fc$var$total)
  array(exp(rep(mu, n_draws) + rep(sd, n_draws) * stats::rnorm(length(mu) * n_draws)), c(dim(mu), n_draws))
}, envir = ns)

native_lca <- function() {
  new_mt_spec("native_lca", list(), min_periods = 5, label = "demography lca (native)", smooth = FALSE)
}
registerS3method("fit_model", "mt_spec_native_lca", function(spec, D, E, ages, periods, ...) {
  dd <- demography::demogdata(pmax(D, 0.5) / E, pop = E, ages = ages$lower, years = periods, type = "mortality",
                              label = "native", name = "s")
  list(fit = demography::lca(dd, series = "s", adjust = "dt", interpolate = FALSE))
}, envir = ns)
registerS3method("simulate_model", "mt_spec_native_lca", function(spec, fit, h, n_draws, ...) {
  # forecast.lca: random walk with drift for k with its 95% bounds; rates
  # move with b times the deviation of k (one index for all ages)
  fc <- forecast::forecast(fit$fit, h = h, level = 95)
  sd_k <- (as.numeric(fc$kt.f$upper) - as.numeric(fc$kt.f$mean)) / stats::qnorm(0.975)
  mu <- log(fc$rate$s)
  out <- array(NA_real_, c(dim(mu), n_draws))
  for (j in seq_len(n_draws)) out[, , j] <- exp(mu + outer(fit$fit$bx, sd_k * stats::rnorm(1)))
  out
}, envir = ns)

native_bamp <- function(mcmc = mcmc_control("fast")) {
  new_mt_spec("native_bamp", list(mcmc = mcmc), min_periods = 5, label = "bamp (native)", smooth = FALSE)
}
registerS3method("fit_model", "mt_spec_native_bamp", function(spec, D, E, ages, periods, ...) {
  # bamp on the age groups: random walks of first order, cohorts from
  # 5-year age groups and annual periods
  w <- unique(utils::head(ages$upper - ages$lower, -1))
  utils::capture.output(m <- bamp::bamp(round(t(D)), round(t(E)), age = "rw1", period = "rw1", cohort = "rw1",
                                        periods_per_agegroup = w, mcmc.options = spec$params$mcmc$options,
                                        parallel = TRUE, verbose = FALSE))
  list(model = m)
}, envir = ns)
registerS3method("simulate_model", "mt_spec_native_bamp", function(spec, fit, h, n_draws, ...) {
  pr <- bamp::predict_apc(fit$model, periods = h)$samples$pr   # [period, age, sample]
  np <- dim(pr)[1]
  fut <- pr[(np - h + 1):np, , , drop = FALSE]
  pick <- sample.int(dim(fut)[3], n_draws, replace = n_draws > dim(fut)[3])
  aperm(fut[, , pick, drop = FALSE], c(2, 1, 3))
}, envir = ns)

# ---- Data --------------------------------------------------------------------

ons <- mt_regroup(mt_read_ons(ons_data), age = c(0, seq(5, 90, 5)))   # equal 5-year groups for bamp
fr <- mt_set_grid(as_mt_data(pkg_data("fr.mort", "demography"), years = 1950:2006),
                  age = single_ages(0, 95, open = FALSE))
ew <- mt_set_grid(mt_window(as_mt_data(pkg_data("EWMaleData", "StMoMo")), 1961, 2011),
                  age = single_ages(55, 89, open = FALSE))
data("apc", package = "bamp", envir = environment())
bamp_ex <- local({
  ages <- seq(40, 75, 5)
  long <- data.table(year = rep(2001:2010, times = 8), age_lo = rep(ages, each = 10), age_hi = rep(ages + 5, each = 10),
                     deaths = as.vector(cases), exposure = as.vector(population))
  x <- mt_data(deaths = long, exposure = long, age = c("age_lo", "age_hi"))
  mt_set_grid(x, age = mt_axis(ages, ages + 5, kind = "age"))
})

# ---- A. Fidelity ---------------------------------------------------------------

log_step("fidelity: Lee-Carter")
fid_lc <- local({
  D <- mt_array(ew, "deaths")[, , 1]
  E <- mt_array(ew, "exposure")[, , 1]
  ages <- 55:89
  yrs <- 1961:2011
  ref <- StMoMo::fit(StMoMo::lc(), data = pkg_data("EWMaleData", "StMoMo"), ages.fit = ages, years.fit = yrs,
                     verbose = FALSE)
  ref_fc <- forecast::forecast(ref, h = 10)
  mine <- ns$lc_poisson(D, E)
  drift <- mean(diff(mine$k))
  mine_fc <- mine$a + outer(mine$b, mine$k[length(mine$k)] + drift * 1:10)
  st_fitted <- ref$ax + outer(ref$bx[, 1], ref$kt[1, ])
  frD <- mt_array(fr, "deaths")[, , 2]
  frE <- mt_array(fr, "exposure")[, , 2]
  dd <- demography::demogdata(pmax(frD, 0.5) / frE, pop = frE, ages = 0:95, years = 1950:2006, type = "mortality",
                              label = "FR", name = "male")
  l <- demography::lca(dd, adjust = "dt", interpolate = FALSE)
  lf <- forecast::forecast(l, h = 10)
  s <- ns$lc_svd(frD, frE)
  s_fc <- s$a + outer(s$b, s$k[length(s$k)] + mean(diff(s$k)) * 1:10)
  data.table(comparison = c("lee_carter() vs StMoMo lc(): fitted log rates",
                            "lee_carter() vs StMoMo lc(): central forecast, log rates, h = 1-10",
                            "lee_carter(method = \"svd\") vs demography lca(): fitted log rates",
                            "lee_carter(method = \"svd\") vs demography lca(): central forecast, log rates, h = 1-10"),
             data = c("England & Wales men 55-89, 1961-2011", "England & Wales men 55-89, 1961-2011",
                      "France men 0-95, 1950-2006", "France men 0-95, 1950-2006"),
             max_abs_difference = c(max(abs(mine$a + outer(mine$b, mine$k) - st_fitted)),
                                    max(abs(mine_fc - log(ref_fc$rates))),
                                    max(abs(s$a + outer(s$b, s$k) - (l$ax + outer(l$bx, as.numeric(l$kt))))),
                                    max(abs(s_fc - log(lf$rate$male)))))
})

log_step("fidelity: fdm")
fid_fdm <- local({
  x1 <- mt_window(fr, 1950, 2006)
  D <- mt_array(x1, "deaths")[, , 2]
  E <- mt_array(x1, "exposure")[, , 2]
  dd <- demography::demogdata(pmax(D, 0.5) / E, pop = E, ages = 0:95, years = 1950:2006, type = "mortality",
                              label = "FR", name = "s")
  nat <- forecast::forecast(demography::fdm(dd, series = "s", order = 6), h = 20, method = "arima", level = 80)
  fit <- mt_fit(hyndman_ullah(order = 6), mt_window(fr, 1950, 2006))
  fc <- mt_forecast(fit, h = 20, n_draws = 2000, seed = 1)
  r <- log(fc$rates[, , 2, ])
  q <- apply(r, 1:2, stats::quantile, c(0.1, 0.5, 0.9))
  hz <- c(1, 10, 20)
  data.table(horizon = hz,
             median_abs_log_diff_point = sapply(hz, function(k) stats::median(abs(q[2, , k] - log(nat$rate$s[, k])))),
             max_abs_log_diff_point = sapply(hz, function(k) max(abs(q[2, , k] - log(nat$rate$s[, k])))),
             width80_ratio_median = sapply(hz, function(k) stats::median((q[3, , k] - q[1, , k]) /
                                                                            (log(nat$rate$upper[, k]) -
                                                                               log(nat$rate$lower[, k])))))
})

log_step("fidelity: bamp")
fid_bamp <- local({
  opts <- mcmc_control("fast")
  run_native <- function(seed) {
    set.seed(seed)
    utils::capture.output(m <- bamp::bamp(cases, population, age = "rw1", period = "rw1", cohort = "rw1",
                                          periods_per_agegroup = 5, mcmc.options = opts$options, parallel = TRUE,
                                          verbose = FALSE))
    pr <- bamp::predict_apc(m, periods = 3)$samples$pr
    np <- dim(pr)[1]
    log(pr[(np - 2):np, , , drop = FALSE])   # [period, age, sample]
  }
  n1 <- run_native(1)
  n2 <- run_native(2)
  fit <- mt_fit(apc_bayes(age = "rw1", period = "rw1", cohort = "rw1", mcmc = opts), bamp_ex)
  fc <- mt_forecast(fit, h = 3, n_draws = 3200, seed = 3)
  mine <- aperm(log(fc$rates[, , 1, ]), c(2, 1, 3))
  summ <- function(a) list(med = apply(a, 1:2, stats::median), sd = apply(a, 1:2, stats::sd),
                           lo = apply(a, 1:2, stats::quantile, 0.05), hi = apply(a, 1:2, stats::quantile, 0.95))
  s1 <- summ(n1)
  s2 <- summ(n2)
  sm <- summ(mine)
  data.table(comparison = c("apc_bayes() vs bamp (seed 1)", "bamp (seed 2) vs bamp (seed 1)"),
             max_abs_diff_median_in_sd = c(max(abs(sm$med - s1$med) / s1$sd), max(abs(s2$med - s1$med) / s1$sd)),
             max_abs_diff_90_bounds_in_sd = c(max(abs(c(sm$lo - s1$lo, sm$hi - s1$hi)) / c(s1$sd, s1$sd)),
                                              max(abs(c(s2$lo - s1$lo, s2$hi - s1$hi)) / c(s1$sd, s1$sd))))
})

# ---- B. Forecast accuracy --------------------------------------------------------

run_bt <- function(name, data, specs, origins, h) {
  rbindlist(lapply(specs, function(sp) {
    log_step("%s: %s", name, sp$label)
    t0 <- proc.time()[["elapsed"]]
    b <- tryCatch(mt_backtest(sp, data, origins = origins, h = h, n_draws = n_draws, seed = 1),
                  error = function(e) {
                    message("  failed: ", conditionMessage(e))
                    NULL
                  })
    if (is.null(b)) return(NULL)
    b[, `:=`(data = name, seconds = proc.time()[["elapsed"]] - t0)]
    b
  }), fill = TRUE)
}

acc <- list()
acc$ons <- run_bt("ONS England, sex x IMD quintile, 5-year groups", ons,
                  list(native_fdm(), native_lca(), native_bamp(), hyndman_ullah(), lee_carter(),
                       li_lee(common_by = "sex"), apc_bayes(age = "rw1", period = "rw1", cohort = "rw1"),
                       mt_ensemble(lc = lee_carter(), rw = rwd())),
                  origins = if (quick) 2014 else c(2009, 2011, 2014), h = 5)
acc$france <- run_bt("France (HMD via demography), single ages 0-95", fr,
                     list(native_fdm(), native_lca(), hyndman_ullah(), lee_carter(), rwd(),
                          mt_ensemble(lc = lee_carter(), rw = rwd())),
                     origins = if (quick) 1996 else c(1986, 1991, 1996), h = 10)
acc$ew <- run_bt("England & Wales men (HMD via StMoMo), single ages 55-89", ew,
                 list(stmomo("lc"), stmomo("cbd"), stmomo("apc"), native_lca(), lee_carter(), rwd(),
                      mt_ensemble(lc = lee_carter(), rw = rwd())),
                 origins = if (quick) 2001 else c(1991, 1996, 2001), h = 10)
acc$bamp_example <- run_bt("bamp example data (8 age groups, 10 years)", bamp_ex,
                           list(native_bamp(), apc_bayes(age = "rw1", period = "rw1", cohort = "rw1"), lee_carter(),
                                rwd()),
                           origins = if (quick) 2008 else 2006:2008, h = 2)

score <- function(b) {
  s <- b[, list(crps = mean(crps), wis = mean(wis), mae = mean(ae), cov50 = mean(cov50), cov80 = mean(cov80),
                cov95 = mean(cov95), cells = .N, seconds = max(seconds)), by = list(data, model)]
  setorderv(s, c("data", "crps"))
  s
}
results <- list(
  fidelity = list(lee_carter = fid_lc, fdm = fid_fdm, bamp = fid_bamp),
  accuracy = rbindlist(lapply(acc, score)),
  by_horizon = rbindlist(lapply(acc, function(b) b[, list(crps = mean(crps), cov95 = mean(cov95)),
                                                   by = list(data, model, horizon)])),
  settings = list(n_draws = n_draws, quick = quick, date = Sys.Date(),
                  versions = vapply(c("mortensor", "demography", "bamp", "StMoMo", "forecast"),
                                    function(p) as.character(utils::packageVersion(p)), "")))
saveRDS(results, "vignettes/comparison-results.rds")
log_step("done")
print(results$fidelity)
print(results$accuracy, digits = 3)
