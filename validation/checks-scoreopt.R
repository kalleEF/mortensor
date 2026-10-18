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

checks$scoreopt_regions <- list(
  suite = "full",
  description = "TBD",
  run = function() data.frame(metric = "todo", value = NA_real_))
