# Validation checks of the compositional data (CoDa) models: coda(),
# coda(coherent = TRUE) and coda_causes(). To be merged into checks.R. Run
# times at 2026-10: coda_ew 9 s, coda_coherent_fr 24 s (suite full),
# coda_causes_coherence 7 s. Proposed lines for thresholds.csv
# (check,metric,comparator,threshold,note):
#
# coda_ew,coda_crps,<=,0.06,"0.0548 at 2026-10 (lee_carter() 0.0687, rwd() 0.0586)"
# coda_ew,coda_cov95,>=,0.9,"0.991 at 2026-10 (over-covers: clr noise of ages with few deaths enters k)"
# coda_ew,crps_ratio_coda_lc,<=,0.88,"0.797 at 2026-10"
# coda_coherent_fr,coherent_crps,<=,0.084,"0.0767 at 2026-10 (coda() 0.0786, li_lee() 0.0802)"
# coda_coherent_fr,coherent_cov95,>=,0.9,"0.980 at 2026-10"
# coda_coherent_fr,open_crps,<=,0.07,"0.0604 at 2026-10 (100+; a random walk on the open rate gave 0.117)"
# coda_coherent_fr,crps_ratio_coherent_lilee,<=,1.02,"0.957 at 2026-10"
# (coda_coherent_fr also reports e0_gap_* without thresholds: female minus male life expectancy at birth,
#  observed 6.94 in 2006, after 50 years 4.55 coherent and 5.42 independent at 2026-10)
# coda_causes_coherence,max_gap,<=,1e-12,"causes add up to all causes in expected draws"
# coda_causes_coherence,max_gap_predictive,<=,1e-12,"and in predictive draws (multinomial split)"
# coda_causes_coherence,cc_crps,<=,0.09,"0.0814 at 2026-10 (cause_shares(lee_carter()) 0.0629; experimental)"
# coda_causes_coherence,cc_cov95,>=,0.8,"0.862 at 2026-10"

checks <- if (exists("checks")) checks else list()

checks$coda_ew <- list(
  suite = "quick",
  description = paste("England & Wales men (HMD via StMoMo), ages 30-100, 1961-2011: rolling-origin backtest",
                      "(origins 1996-2001, h = 10) of coda() (CoDa Lee-Carter, observed jump-off) against",
                      "lee_carter(); CRPS of log rates, coverage and the ratio of the CRPS."),
  run = function() {
    ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
    ew <- mt_set_grid(ew, age = single_ages(30, 100, open = FALSE))
    ew <- mt_regroup(ew, age = single_ages(30, 100, open = FALSE))
    o <- 1996:2001
    cd <- mt_backtest(coda(), ew, origins = o, h = 10, n_draws = 200)
    lc <- mt_backtest(lee_carter(), ew, origins = o, h = 10, n_draws = 200)
    data.frame(metric = c("coda_crps", "coda_cov95", "coda_cov80", "lc_crps", "crps_ratio_coda_lc"),
               value = c(mean(cd$crps), mean(cd$cov95), mean(cd$cov80), mean(lc$crps), mean(cd$crps) / mean(lc$crps)))
  })

checks$coda_coherent_fr <- list(
  suite = "full",
  description = paste("France (HMD via demography::fr.mort), women and men, ages 0-99 and 100+, 1950-2006:",
                      "backtest (origins 1991-1996, h = 10) of coda(coherent = TRUE) against coda() and li_lee();",
                      "CRPS and coverage of log rates, CRPS of the open interval (100+), and the gap in life",
                      "expectancy at birth between women and men in 2006 and after 50 years of forecasts (medians",
                      "of the draws of the coherent and the independent model; information, no threshold)."),
  run = function() {
    fr <- as_mt_data(mortensor:::pkg_data("fr.mort", "demography"))
    fr <- mt_window(fr, from = 1950, to = 2006)
    fr <- mt_regroup(fr, age = single_ages(0, 100, open = TRUE))
    o <- 1991:1996
    cc <- mt_backtest(coda(coherent = TRUE), fr, origins = o, h = 10, n_draws = 200)
    cd <- mt_backtest(coda(), fr, origins = o, h = 10, n_draws = 200)
    ll <- mt_backtest(li_lee(), fr, origins = o, h = 10, n_draws = 200)
    e0_gap <- function(e, year) e[sex == "female" & period == year]$q0.5 - e[sex == "male" & period == year]$q0.5
    gap_50 <- function(spec) e0_gap(life_expectancy(mt_forecast(mt_fit(spec, fr), h = 50, n_draws = 200, seed = 1)),
                                    2056)
    data.frame(metric = c("coherent_crps", "coherent_cov95", "coherent_cov80", "open_crps", "coda_crps",
                          "lilee_crps", "crps_ratio_coherent_lilee", "e0_gap_2006", "e0_gap_2056_coherent",
                          "e0_gap_2056_independent"),
               value = c(mean(cc$crps), mean(cc$cov95), mean(cc$cov80), mean(cc[age_lo == 100]$crps), mean(cd$crps),
                         mean(ll$crps), mean(cc$crps) / mean(ll$crps), e0_gap(life_expectancy(fr), 2006),
                         gap_50(coda(coherent = TRUE)), gap_50(coda())))
  })

checks$coda_causes_coherence <- list(
  suite = "quick",
  description = paste("coda_causes() (multiple-decrement CoDa, experimental) on mt_example_causes(), 1995-2017:",
                      "largest relative gap between all-cause draws and the sum of the cause draws (expected and",
                      "predictive draws, origin 2012, h = 5), and a backtest (origins 2008-2012, h = 5) against",
                      "cause_shares(lee_carter()); CRPS and coverage of log rates."),
  run = function() {
    x <- mt_example_causes(years = 1995:2017)
    st <- mt_strata(x)
    fit <- mt_fit(coda_causes(), mt_window(x, to = 2012))
    gap <- function(fc) {
      max(vapply(c("female", "male"), function(sx) {
        i <- st$stratum[st$sex == sx & st$cause == "all"]
        j <- st$stratum[st$sex == sx & st$cause != "all"]
        tot <- fc$rates[, , i, ]
        sums <- apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum)
        max(abs(sums - tot) / pmax(tot, 1e-12))
      }, 0))
    }
    g1 <- gap(mt_forecast(fit, h = 5, n_draws = 100, seed = 1))
    g2 <- gap(mt_forecast(fit, h = 5, n_draws = 100, seed = 1, predictive = TRUE))
    o <- 2008:2012
    cc <- mt_backtest(coda_causes(), x, origins = o, h = 5, n_draws = 100)
    cs <- mt_backtest(cause_shares(total = lee_carter()), x, origins = o, h = 5, n_draws = 100)
    data.frame(metric = c("max_gap", "max_gap_predictive", "cc_crps", "cc_cov95", "cs_crps", "crps_ratio_cc_cs"),
               value = c(g1, g2, mean(cc$crps), mean(cc$cov95), mean(cs$crps), mean(cc$crps) / mean(cs$crps)))
  })
