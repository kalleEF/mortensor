# Benchmarks shared by the slow tests and the validation pipeline (validation/)

#' Ungrouping benchmark on England & Wales (HMD via StMoMo)
#'
#' Groups observed single-year deaths of men aged 30-100 to 30-34, ..., 85-89
#' and 90-100, ungroups them with single-year exposures as offset and compares
#' the result with the observed single-year rates. Errors are median absolute
#' log errors; the Poisson noise of the observed rates is a floor (about 0.012
#' at 30-89 and 0.021 at 90-100 for 1980-2011).
#'
#' @param years Years used.
#' @param dims 1 or 2.
#' @param ... Passed to [mt_ungroup()].
#' @return A list with `male_30_89`, `male_90_100`, `seconds`.
#' @keywords internal
benchmark_ungroup_ew <- function(years = 1980:2011, dims = 1, ...) {
  need_pkg("StMoMo", "for the England & Wales benchmark")
  ew <- pkg_data("EWMaleData", "StMoMo")
  ages <- 30:100
  D <- ew$Dxt[as.character(ages), as.character(years), drop = FALSE]
  E <- ew$Ext[as.character(ages), as.character(years), drop = FALSE]
  lo <- c(seq(30, 85, 5), 90)
  gid <- findInterval(ages, lo)
  Dg <- rowsum(D, gid)
  dg <- data.table(year = rep(years, each = length(lo)), age_lo = rep(lo, length(years)),
                   age_hi = rep(c(lo[-1], 101), length(years)), deaths = as.vector(Dg))
  ex <- data.table(year = rep(years, each = length(ages)), age_lo = rep(ages, length(years)),
                   age_hi = rep(ages + 1, length(years)), exposure = as.vector(E))
  x <- mt_data(deaths = dg, exposure = ex, age = c("age_lo", "age_hi"))
  t0 <- proc.time()[["elapsed"]]
  u <- mt_ungroup(x, dims = dims, ...)
  secs <- proc.time()[["elapsed"]] - t0
  est <- mt_array(u, "deaths")[, , 1] / mt_array(u, "exposure")[, , 1]
  le <- abs(log(est) - log(D / E))
  list(male_30_89 = stats::median(le[ages < 90, ]), male_90_100 = stats::median(le[ages >= 90, ]),
       seconds = secs, estimate = est)
}

#' Simulated Lee-Carter data for calibration checks
#'
#' The Lee-Carter model with Poisson deaths (Brouhns, Denuit and Vermunt
#' 2002): `log m = a_x + b_x k_t + shift_g`, deaths Poisson with mean
#' `m * exposure`, with `k_t` a random walk with drift (-0.015, sd 0.02) per
#' stratum, 5-year age groups 30-34 ... 90-94 and 10 strata (2 sexes x 5
#' deprivation quintiles). `sd_obs > 0` adds iid normal noise to the log
#' rates before the Poisson draw (the generator of the 0.1 review used 0.03),
#' for studies of misspecified models.
#'
#' @param n_t Number of years.
#' @param seed Random seed.
#' @param pop Population at age 30 per stratum.
#' @param sd_k,drift Standard deviation and drift of the random walk `k_t`.
#' @param sd_obs Standard deviation of iid noise on the log rates (0: none,
#'   the Lee-Carter model with Poisson deaths).
#' @return An [mt_data] object.
#' @keywords internal
simulate_lc_dgp <- function(n_t = 40, seed = 1, pop = 2e5, sd_k = 0.02, drift = -0.015, sd_obs = 0) {
  with_seed(seed, {
    ages <- seq(30, 90, 5)
    na <- length(ages)
    groups <- expand.grid(sex = c("male", "female"), quintile = paste0("Q", 1:5), stringsAsFactors = FALSE)
    a_x <- -10.5 + 0.085 * ages
    b_x <- seq(1.3, 0.7, length.out = na)
    rows <- vector("list", nrow(groups))
    for (g in seq_len(nrow(groups))) {
      shift <- ifelse(groups$sex[g] == "male", 0.3, 0) + (g - 1) %/% 2 * 0.08
      k <- cumsum(drift + stats::rnorm(n_t, 0, sd_k))
      m <- exp(a_x + outer(b_x, k) + shift + matrix(stats::rnorm(na * n_t, 0, sd_obs), na))
      e <- matrix(round(pop * 5 * exp(-0.03 * (ages - 30))), na, n_t)
      rows[[g]] <- data.table(sex = groups$sex[g], quintile = groups$quintile[g],
                              year = rep(2000 + seq_len(n_t), each = na), age_lo = rep(ages, n_t),
                              age_hi = rep(ages + 5, n_t), deaths = stats::rpois(na * n_t, m * e),
                              exposure = as.vector(e))
    }
    long <- rbindlist(rows)
    x <- mt_data(deaths = long, exposure = long, by = c("sex", "quintile"), age = c("age_lo", "age_hi"))
    mt_set_grid(x, age = age_groups(ages, open = FALSE, last_upper = 95))
  })
}

#' Interval calibration benchmark on simulated data
#'
#' Fits `spec` on 30 years of [simulate_lc_dgp()] data and scores predictive
#' draws for the next `h` years, over `n_rep` replicates.
#'
#' @param spec Model specification.
#' @param n_rep Replicates.
#' @param h Horizon.
#' @param n_draws Draws per forecast.
#' @param n_train Training years.
#' @param prep Optional [mt_prep()] pipeline.
#' @param sd_obs Noise of the generator (see [simulate_lc_dgp()]).
#' @return A `data.table` with coverage by horizon and a summary list in
#'   attribute `"summary"` (`cov95_mean`, `cov95_min`, `cov80_mean`,
#'   `cov50_mean`).
#' @keywords internal
benchmark_coverage <- function(spec, n_rep = 10, h = 10, n_draws = 300, n_train = 30, prep = NULL, sd_obs = 0) {
  bt <- rbindlist(lapply(seq_len(n_rep), function(r) {
    x <- simulate_lc_dgp(n_t = n_train + h, seed = 1000 + r, sd_obs = sd_obs)
    b <- mt_backtest(spec, x, origins = 2000 + n_train, h = h, n_draws = n_draws, prep = prep, seed = r)
    b[, rep := r]
  }))
  by_h <- mt_score(bt, by = "horizon")
  attr(by_h, "summary") <- list(cov95_mean = mean(by_h$cov95), cov95_min = min(by_h$cov95),
                                cov80_mean = mean(by_h$cov80), cov50_mean = mean(by_h$cov50),
                                crps = mean(bt$crps))
  by_h
}
