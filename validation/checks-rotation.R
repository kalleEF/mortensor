# Proposed validation check for lee_carter(rotate = TRUE), the rotated
# Lee-Carter model of Li, Lee and Gerland (2013). To be merged into
# checks.R (same format: suite, description, run returning metric and value).
#
# Proposed lines for thresholds.csv (check,metric,comparator,threshold,note);
# values measured 2026-10-18 with 200 draws:
#
# rotation_france,e0_max_diff,<=,1e-6,"rotated paths keep the e0 of the unrotated paths; 1.0e-10 at 2026-10"
# rotation_france,decline_diff_0_14,<=,-0.8,"rotated minus unrotated decline (pp a year) 2006-2100; -1.25 at 2026-10"
# rotation_france,decline_diff_15_64,>=,0.15,"0.29 at 2026-10"
# rotation_france,decline_diff_65_84,>=,-0.15,"-0.03 at 2026-10 (b at 65 is above the mean over 15-64 in France)"
# rotation_france,crps_ratio_h10,<=,1.01,"rotated / unrotated CRPS of log rates, origins 1986-1996; 0.990 at 2026-10"
# rotation_france,crps_ratio_h20,<=,1,"origins 1976-1986, h = 20; 0.976 at 2026-10"
# rotation_france,rot_cov95_h20,>=,0.9,"0.947 at 2026-10"
# rotation_france,seconds_forecast,<=,40,"runtime budget; rotated forecast to 2100, 101 ages x 94 years x 200 draws"

if (!exists("checks")) checks <- list()

checks$rotation_france <- list(
  suite = "full",
  description = paste("lee_carter(rotate = TRUE) on French women (HMD via demography), ages 0-99 and 100+, 1950-2006:",
                      "forecast to 2100, where rotated and unrotated paths have the same life expectancy at birth",
                      "and the rotated ones a slower decline at ages 0-14 and a faster one at ages 15-64 (average",
                      "annual rate of decline of the death rate of each age band, 2006 exposures, median over",
                      "draws); backtests from 1986, 1991 and 1996 (h = 10) and from 1976, 1981 and 1986 (h = 20),",
                      "CRPS of log rates relative to the unrotated model."),
  run = function() {
    fr <- suppressWarnings(as_mt_data(mortensor:::pkg_data("fr.mort", "demography"), series = "female",
                                      years = 1950:2006))
    fr <- mt_regroup(fr, age = 0:100)
    f0 <- mt_fit(lee_carter(), fr)
    f1 <- mt_fit(lee_carter(rotate = TRUE), fr)
    fc0 <- mt_forecast(f0, to = 2100, n_draws = 200, seed = 1)
    t0 <- proc.time()[["elapsed"]]
    fc1 <- mt_forecast(f1, to = 2100, n_draws = 200, seed = 1)
    secs <- proc.time()[["elapsed"]] - t0
    e0 <- life_expectancy(fc0, paths = TRUE)
    e1 <- life_expectancy(fc1, paths = TRUE)
    lo <- fc0$age$lower
    band <- cut(lo, c(0, 15, 65, 85, Inf), right = FALSE, labels = c("0-14", "15-64", "65-84", "85+"))
    ex <- mt_array(fr, "exposure")[, length(f0$periods), 1]
    group_rates <- function(m) tapply(ex * m, band, sum) / tapply(ex, band, sum)
    jump <- group_rates(f0$fits[[1]]$fitted[, length(f0$periods)])
    nh <- length(fc0$periods)
    decline <- function(fc) {
      r <- vapply(seq_len(dim(fc$rates)[4]), function(j) -log(group_rates(fc$rates[, nh, 1, j]) / jump) / nh,
                  numeric(4))
      100 * apply(r, 1, stats::median)
    }
    d0 <- decline(fc0)
    d1 <- decline(fc1)
    bt <- function(spec, origins, h) mt_backtest(spec, fr, origins = origins, h = h, n_draws = 200)
    b10 <- bt(lee_carter(), c(1986, 1991, 1996), 10)
    r10 <- bt(lee_carter(rotate = TRUE), c(1986, 1991, 1996), 10)
    b20 <- bt(lee_carter(), c(1976, 1981, 1986), 20)
    r20 <- bt(lee_carter(rotate = TRUE), c(1976, 1981, 1986), 20)
    data.frame(metric = c("e0_max_diff", "decline_unrot_0_14", "decline_rot_0_14", "decline_diff_0_14",
                          "decline_diff_15_64", "decline_diff_65_84", "decline_diff_85_plus", "weight_2100",
                          "crps_unrot_h10", "crps_rot_h10", "crps_ratio_h10", "crps_unrot_h20", "crps_rot_h20",
                          "crps_ratio_h20", "rot_cov95_h20", "seconds_forecast"),
               value = c(max(abs(e1 - e0)), d0[1], d1[1], d1[1] - d0[1], d1[2] - d0[2], d1[3] - d0[3],
                         d1[4] - d0[4], stats::median(mortensor:::rotation_weights(e1[nh, 1, ])),
                         mean(b10$crps), mean(r10$crps), mean(r10$crps) / mean(b10$crps), mean(b20$crps),
                         mean(r20$crps), mean(r20$crps) / mean(b20$crps), mean(r20$cov95), secs))
  })
