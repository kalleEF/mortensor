# Proposed validation checks for score-optimal reconciliation
# (reconciled(method = "energy_score"), R/hierarchy-scoreopt.R). They use the
# format of checks.R and are meant to be moved there (this file is not
# sourced by run_validation.R). Proposed lines for thresholds.csv:
#
# check,metric,comparator,threshold,note
# scoreopt_probreco,max_diff_G,<=,1e-6,"Adam path of ProbReco::scoreopt() reproduced (same draws; 2e-14 at 2026-10)"
# scoreopt_regions,max_gap,<=,1e-6,"reconciled forecasts are coherent"
# scoreopt_regions,crps_ratio_bottom,<=,TBD,""
# scoreopt_regions,crps_ratio_aggregate,<=,TBD,""
# scoreopt_regions,crps_ratio_total,<=,TBD,""
# scoreopt_regions,cov95_bottom,>=,TBD,""
# scoreopt_regions,seconds_fit,<=,TBD,"runtime budget of one fit"

checks$scoreopt_probreco <- list(
  suite = "quick",
  description = paste("Score-optimal reconciliation against ProbReco::scoreopt() (Panagiotelis et al. 2023) on a",
                      "toy hierarchy of 7 nodes: 60 Adam steps from OLS with the same base draws; largest",
                      "difference of d and G."),
  run = function() {
    if (!requireNamespace("ProbReco", quietly = TRUE)) return(data.frame(metric = "skipped", value = 1))
    S <- rbind(c(1, 1, 1, 1), c(1, 1, 0, 0), c(0, 0, 1, 1), diag(4))
    n <- nrow(S)
    m <- ncol(S)
    nt <- 6
    q <- 20
    set.seed(2)
    Y <- vapply(seq_len(nt), function(t) as.vector(S %*% (c(5, 3, 2, 4) + stats::rnorm(4))), numeric(n))
    mu <- as.vector(S %*% c(5, 3, 2, 4)) * c(1.2, 0.9, 1.1, 1, 1, 1, 1)
    prob <- lapply(seq_len(nt), function(t) function() matrix(stats::rnorm(n * q, mu, 0.7), n))
    set.seed(3)
    ref <- suppressWarnings(ProbReco::scoreopt(lapply(seq_len(nt), function(t) Y[, t]), prob, S,
                                               control = list(maxIter = 60, eta = 0.01)))
    # ProbReco calls every function three times to check its inputs
    set.seed(3)
    for (k in 1:3) for (f in prob) f()
    draw <- function() {
      x1 <- lapply(prob, function(f) f())
      x2 <- lapply(prob, function(f) f())
      list(x1 = do.call(cbind, x1), x2 = do.call(cbind, x2), per = rep(seq_len(nt), each = q))
    }
    ctl <- list(max_iter = 60, eta = 0.01, beta1 = 0.9, beta2 = 0.999, tol = 1e-4, epsilon = 1e-8)
    G0 <- solve(crossprod(S), t(S))
    mine <- mortensor:::scoreopt_sgd(Y, draw, S, c(rep(0, m), as.vector(G0)), ctl)
    data.frame(metric = c("max_diff_G", "max_diff_d", "max_change_from_ols"),
               value = c(max(abs(mine$G - ref$G)), max(abs(mine$d - ref$d)), max(abs(ref$G - G0))))
  })

# Backtest of score-optimal reconciliation against WLS: CRPS ratios and 95%
# coverage by level of the hierarchy, the coherence gap of a forecast and the
# runtime of one fit
scoreopt_backtest <- function(xa, h, origins, spec_es) {
  t0 <- proc.time()[["elapsed"]]
  fit <- suppressWarnings(mt_fit(spec_es, xa))
  secs <- proc.time()[["elapsed"]] - t0
  fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1, predictive = TRUE)
  hs <- mortensor:::hierarchy_structure(fc$strata, h)
  E <- mortensor:::expand_exposure(fc$exposure_last, dim(fc$rates))
  D <- fc$rates * as.vector(E)
  gap <- max(vapply(seq_len(dim(D)[4]), function(j) {
    max(vapply(seq_len(dim(D)[2]), function(t) {
      full <- matrix(D[, t, hs$bottom, j], dim(D)[1]) %*% t(hs$S)
      max(abs(full / matrix(D[, t, , j], dim(D)[1]) - 1))
    }, 0))
  }, 0))
  bt <- suppressWarnings(mt_compare(list(wls = reconciled(lee_carter(), h, method = "wls"), es = spec_es), xa,
                                    origins = origins, h = 5, n_draws = 200))
  s <- mt_score_table(bt, metrics = c("crps", "cov95"), hierarchy = h, digits = 6)
  v <- function(m, l, col) s[model == m & level == l][[col]]
  lev <- c("bottom", "aggregate", "total")
  data.frame(metric = c(paste0("crps_ratio_", lev), paste0("cov95_", lev), "max_gap", "seconds_fit"),
             value = c(vapply(lev, function(l) v("es", l, "crps") / v("wls", l, "crps"), 0),
                       vapply(lev, function(l) v("es", l, "cov95"), 0), gap, secs))
}

checks$scoreopt_regions <- list(
  suite = "full",
  description = paste("Score-optimal reconciliation (reconciled(lee_carter(), method = 'energy_score'), defaults:",
                      "10 training years, 5 inner origins, horizons 1-5) on mt_example_regions() with 2 sexes x 2",
                      "regions x 5 quintiles and all aggregates (54 strata), 1995-2019: backtest from 2009, 2011",
                      "and 2013 (h = 5); CRPS relative to WLS and 95% coverage by level, largest coherence gap of",
                      "predictive draws, seconds for one fit."),
  run = function() {
    x <- mt_example_regions(years = 1995:2019, regions = 2)
    h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
    xa <- mt_aggregate_strata(x, h)
    scoreopt_backtest(xa, h, c(2009, 2011, 2013), reconciled(lee_carter(), h, method = "energy_score"))
  })

checks$scoreopt_ons <- list(
  suite = "full",
  description = paste("As scoreopt_regions on ONS England (ons_data), sex x IMD quintile with England (18 strata),",
                      "training on 6 years (5 inner origins, horizons 1-3): backtest from 2011-2014 (h = 5)."),
  run = function() {
    x <- ons_quintiles()
    h <- mt_hierarchy(sex = "both", quintile = list(total = "England"))
    xa <- mt_aggregate_strata(x, h)
    scoreopt_backtest(xa, h, 2011:2014, reconciled(lee_carter(), h, method = "energy_score", years = 6,
                                                   n_origins = 5, horizon = 3))
  })
