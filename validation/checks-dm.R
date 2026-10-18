# Proposed validation checks for mt_dm_test() and the M7/M8 models of
# stmomo(), in the format of checks.R (each check returns a data frame with
# columns `metric` and `value`). To adopt them, move them into checks.R and
# the thresholds below into thresholds.csv. The file can also be sourced on
# its own after pkgload::load_all() or library(mortensor).
#
# Proposed lines for thresholds.csv (check,metric,comparator,threshold,note):
#
# dm_ew,max_abs_diff_stat,<=,1e-8,"agreement with forecast::dm.test() on the losses per origin (3e-14 at 2026-10)"
# dm_ew,max_abs_diff_p,<=,1e-8,"agreement with forecast::dm.test() (7e-16 at 2026-10)"
# stmomo_cbd_ew,m7_crps,<=,0.062,"regression guard; 0.0558 at 2026-10 (cbd 0.0671, m6 0.0384, m8 0.0728)"
# stmomo_cbd_ew,m7_cov95,>=,0.9,"0.956 at 2026-10 (cbd 0.790)"
# stmomo_cbd_ew,m8_crps,<=,0.08,"regression guard; 0.0728 at 2026-10 with xc = 89 (xc = 100: 0.0586)"
# stmomo_cbd_ew,crps_ratio_m7_cbd,<=,0.95,"0.832 at 2026-10: the cohort effect helps on E&W men"
# stmomo_cbd_ew,m7_max_rel_diff_stmomo,<=,1e-6,"fitted rates as StMoMo's own example (0 at 2026-10)"
# stmomo_cbd_ew,m8_max_rel_diff_stmomo,<=,1e-6,"fitted rates as StMoMo's own example (0 at 2026-10)"
#
# The other metrics are info only. Values at 2026-10: dm_ew stat_all 5.61, p_h1 0.020,
# p_h5 0.137, p_all 0.0025 (lag 0 fallback: 6 origins, h = 10), p_all_spaced 0.419;
# stmomo_cbd_ew p_m7_cbd_all 0.0061 (lag 0 fallback). stmomo_cbd_ew takes about 16 s
# and dm_ew about 8 s.

if (!exists("checks")) checks <- list()

# England & Wales men (HMD via StMoMo) on single years of age
ew_single_ages <- function(from, to) {
  ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
  a <- single_ages(from, to, open = FALSE)
  mt_regroup(mt_set_grid(ew, age = a), age = a)
}

checks$dm_ew <- list(
  suite = "quick",
  description = paste("Diebold-Mariano test (mt_dm_test()) on the backtest of backtest_ew (England & Wales men, ages",
                      "30-100, origins 1996-2001, h = 10): Lee-Carter against rwd(), CRPS of log rates. Agreement",
                      "with forecast::dm.test() on the losses per origin (largest absolute difference of the",
                      "statistic and the p-value over horizons 1-10 and all, with the lag used); p-values at",
                      "horizons 1 and 5 and over all horizons, and over all horizons for 7 origins 5 years apart",
                      "(1971-2001) as info."),
  run = function() {
    ew <- ew_single_ages(30, 100)
    models <- list(rwd = rwd(), lee_carter = lee_carter())
    bt <- suppressWarnings(mt_compare(models, ew, origins = 1996:2001, h = 10, n_draws = 200))
    r <- suppressWarnings(mt_dm_test(bt, ref = "rwd"))
    b <- as.data.frame(bt)
    loss <- function(mod, h = NULL) {
      s <- b[b$model == mod, ]
      if (!is.null(h)) return(tapply(s$crps[s$horizon == h], s$origin[s$horizon == h], mean))
      rowMeans(tapply(s$crps, list(s$origin, s$horizon), mean))
    }
    hs <- as.character(r$horizon)
    ref <- t(vapply(seq_along(hs), function(i) {
      h <- if (hs[i] == "all") NULL else as.numeric(hs[i])
      d <- suppressWarnings(forecast::dm.test(loss("lee_carter", h), loss("rwd", h), h = r$lag[i] + 1, power = 1))
      c(unname(d$statistic), unname(d$p.value))
    }, numeric(2)))
    bt5 <- suppressWarnings(mt_compare(models, ew, origins = seq(1971, 2001, by = 5), h = 10, n_draws = 200))
    r5 <- suppressWarnings(mt_dm_test(bt5, ref = "rwd"))
    data.frame(metric = c("max_abs_diff_stat", "max_abs_diff_p", "mean_diff_all", "stat_all", "p_h1", "p_h5",
                          "p_all", "p_all_spaced"),
               value = c(max(abs(r$statistic - ref[, 1])), max(abs(r$p_value - ref[, 2])),
                         r$mean_diff[hs == "all"], r$statistic[hs == "all"], r$p_value[hs == "1"],
                         r$p_value[hs == "5"], r$p_value[hs == "all"], r5$p_value[r5$horizon == "all"]))
  })

checks$stmomo_cbd_ew <- list(
  suite = "full",
  description = paste("CBD family of stmomo() on England & Wales men (HMD via StMoMo), ages 55-89, 1961-2011:",
                      "rolling-origin backtest (origins 1996-2001, h = 10) of CBD, M6, M7 and M8 (xc = 89); CRPS",
                      "and 95% coverage of log rates, CRPS ratio M7 / CBD and the Diebold-Mariano p-value over all",
                      "horizons (info). Fitted rates of M7 and M8 on 1961-2011 against StMoMo's own example",
                      "(central2initial(), clip = 3)."),
  run = function() {
    ew <- ew_single_ages(55, 89)
    models <- list(cbd = stmomo("cbd"), m6 = stmomo("m6"), m7 = stmomo("m7"), m8 = stmomo("m8"))
    bt <- mt_compare(models, ew, origins = 1996:2001, h = 10, n_draws = 200)
    s <- mt_score(bt)
    r <- suppressWarnings(mt_dm_test(bt, ref = "cbd"))
    crps <- function(m) s$crps[s$model == m]
    # fitted rates against StMoMo called as in its documentation
    sd <- StMoMo::central2initial(mortensor:::pkg_data("EWMaleData", "StMoMo"))
    wxt <- StMoMo::genWeightMat(55:89, sd$years, clip = 3)
    rel <- function(m, mod) {
      direct <- suppressWarnings(StMoMo::fit(mod, data = sd, ages.fit = 55:89, wxt = wxt, verbose = FALSE))
      q <- stats::fitted(direct, type = "rates")
      f <- stats::fitted(mt_fit(stmomo(m), ew))[, , 1]
      ok <- !is.na(q)
      max(abs(f[ok] / -log(1 - q[ok]) - 1))
    }
    data.frame(metric = c("cbd_crps", "m6_crps", "m7_crps", "m8_crps", "cbd_cov95", "m7_cov95", "m8_cov95",
                          "crps_ratio_m7_cbd", "p_m7_cbd_all", "m7_max_rel_diff_stmomo", "m8_max_rel_diff_stmomo"),
               value = c(crps("cbd"), crps("m6"), crps("m7"), crps("m8"), s$cov95[s$model == "cbd"],
                         s$cov95[s$model == "m7"], s$cov95[s$model == "m8"], crps("m7") / crps("cbd"),
                         r$p_value[r$model == "m7" & r$horizon == "all"], rel("m7", StMoMo::m7()),
                         rel("m8", StMoMo::m8(xc = 89))))
  })
