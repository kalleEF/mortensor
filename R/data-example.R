#' Simulated mortality data with known truth
#'
#' Generates a two-stratum dataset (women and men) from a Lee-Carter-type
#' model with Gompertz-Makeham age pattern, a random walk with drift for the
#' period index, and a population with cohort-size waves. Deaths are Poisson
#' draws. By default deaths are reported in 5-year age groups with an open
#' 90+ interval and exposures in single years of age with an open 100+
#' interval, as in many official statistics. A model fits these data after
#' the deaths are ungrouped ([mt_ungroup()] or `step_ungroup()`) or all
#' tables are aggregated to common age groups ([mt_regroup()]); the
#' examples of the other help pages pass `deaths = "single"` (deaths in
#' single years of age) so that every model fits them directly.
#'
#' The true rates on the default grid are attached as attribute `"truth"`
#' (array `[age, period, stratum]`), which makes the dataset suitable for
#' testing ungrouping and forecast calibration.
#'
#' @param years Calendar years.
#' @param deaths `"grouped"` (5-year groups and 90+) or `"single"` (single
#'   ages and 100+).
#' @param size Population size at age 30 per year and stratum.
#' @param seed Random seed.
#' @return An [mt_data] object with attribute `"truth"`.
#' @examples
#' x <- mt_example()
#' x
#' dim(attr(x, "truth"))
#' @export
mt_example <- function(years = 1990:2019, deaths = c("grouped", "single"), size = 1e5, seed = 1) {
  deaths <- match.arg(deaths)
  sim <- with_seed(seed, simulate_lc_population(years, size))
  ages <- sim$ages
  open_e <- 100
  rows <- list()
  for (s in seq_along(sim$sex)) {
    for (j in seq_along(years)) {
      e <- sim$E[, j, s]
      d <- sim$D[, j, s]
      ae <- pmin(ages, open_e)
      rows[[length(rows) + 1L]] <- data.table(
        sex = sim$sex[s], year = years[j], age_lo = ae, d = d, e = e)
    }
  }
  long <- rbindlist(rows)
  long[, age_hi := ifelse(age_lo >= open_e, Inf, age_lo + 1)]
  E <- long[, list(exposure = sum(e)), by = c("sex", "year", "age_lo", "age_hi")]
  if (deaths == "grouped") {
    long[, `:=`(g_lo = ifelse(age_lo >= 90, 90, 30 + 5 * floor((age_lo - 30) / 5)))]
    long[, g_hi := ifelse(g_lo >= 90, Inf, g_lo + 5)]
    D <- long[, list(deaths = sum(d)), by = c("sex", "year", "g_lo", "g_hi")]
    setnames(D, c("g_lo", "g_hi"), c("age_lo", "age_hi"))
  } else {
    D <- long[, list(deaths = sum(d)), by = c("sex", "year", "age_lo", "age_hi")]
  }
  x <- mt_data(deaths = D, exposure = E, by = "sex", age = c("age_lo", "age_hi"), omega = 110)
  g <- mt_grid(x)
  truth_dx <- sim$m * sim$E
  in_open <- ages >= open_e
  tr <- array(NA_real_, c(length(g$age), length(years), length(sim$sex)),
              dimnames = list(age = g$age$labels, period = g$period$labels, stratum = mt_strata(x)$label))
  for (s in seq_along(sim$sex)) {
    tr[!grepl("\\+", g$age$labels), , s] <- sim$m[!in_open, , s]
    tr[length(g$age), , s] <- colSums(truth_dx[in_open, , s, drop = FALSE]) / colSums(sim$E[in_open, , s, drop = FALSE])
  }
  attr(x, "truth") <- tr
  attr(x, "truth_single") <- list(ages = ages, m = sim$m, E = sim$E)
  x
}

simulate_lc_population <- function(years, size) {
  ages <- 30:109
  sex <- c("female", "male")
  nt <- length(years)
  bx <- 0.02 * (1 - (ages - 30) / 100)
  m <- E <- D <- array(NA_real_, c(length(ages), nt, 2))
  for (s in 1:2) {
    lvl <- if (s == 1) 0.6 else 1
    ax <- log(lvl * (2e-5 * exp(0.1 * ages)) + 3e-4)
    kt <- cumsum(c(0, stats::rnorm(nt - 1, -1, 0.6)))
    kt <- kt - mean(kt)
    m[, , s] <- exp(ax + outer(bx, kt))
    surv <- exp(-cumsum(c(0, head(exp(ax + bx * kt[1]), -1))))
    for (j in seq_len(nt)) {
      cohort <- years[j] - ages
      wave <- 1 + 0.25 * sin(2 * pi * (cohort - 1900) / 25) - 0.3 * (cohort %in% c(1918, 1919, 1945, 1946))
      E[, j, s] <- size * surv * wave * (1 + 0.004 * (years[j] - years[1]))
    }
    D[, , s] <- stats::rpois(length(ages) * nt, m[, , s] * E[, , s])
  }
  list(ages = ages, sex = sex, m = m, E = E, D = D)
}

with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (has_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv())
  })
  seed_default(seed)
  expr
}

# set.seed() with R's default generators, so that a seed gives the same
# numbers also where another generator is active (the L'Ecuyer-CMRG streams of
# par_lapply())
seed_default <- function(seed) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
}

#' Simulated deaths by sex, region and deprivation quintile
#'
#' A grouped structure as in small-area mortality data: `sex` x `region` x
#' `quintile` bottom strata. Log rates follow a Lee-Carter-type model with a
#' period index per sex, region-specific deviations (AR(1)), a deprivation
#' gradient that widens slowly over time, Gompertz-Makeham age patterns and
#' regional population sizes that differ by a factor of up to four. Deaths
#' are Poisson. True rates of the bottom strata are attached as attribute
#' `"truth"`.
#'
#' @param years Calendar years.
#' @param regions Number of regions.
#' @param ages `"grouped"` (5-year groups 40-44 ... 85-89 and 90+, also used
#'   as grid) or `"single"` (single years 40-99 and 100+).
#' @param size Population at age 40 per bottom stratum in the largest region.
#' @param seed Random seed.
#' @return An [mt_data] object with strata `sex` x `region` x `quintile`
#'   (bottom level only; add aggregates with [mt_aggregate_strata()]).
#' @examples
#' x <- mt_example_regions(years = 2010:2019, regions = 3)
#' x
#' @export
mt_example_regions <- function(years = 1995:2019, regions = 4, ages = c("grouped", "single"), size = 2e4,
                               seed = 1) {
  ages <- match.arg(ages)
  with_seed(seed, {
    a <- 40:109
    nt <- length(years)
    reg <- paste0("R", seq_len(regions))
    reg_size <- seq(1, 0.25, length.out = regions)
    reg_level <- stats::rnorm(regions, 0, 0.08)
    bx <- 0.02 * (1 - (a - 40) / 90)
    rows <- list()
    for (sx in c("female", "male")) {
      lvl <- if (sx == "female") 0.6 else 1
      kt <- cumsum(c(0, stats::rnorm(nt - 1, -1, 0.5)))
      kt <- kt - mean(kt)
      for (r in seq_len(regions)) {
        dev <- numeric(nt)
        for (t in seq_len(nt)[-1]) dev[t] <- 0.7 * dev[t - 1] + stats::rnorm(1, 0, 0.02)
        for (q in 1:5) {
          grad <- 0.12 * (q - 3) * (1 + 0.01 * (seq_len(nt) - 1))
          ax <- log(lvl * 2e-5 * exp(0.1 * a) + 3e-4)
          logm <- ax + outer(bx, kt) + reg_level[r] + matrix(grad + dev, length(a), nt, byrow = TRUE)
          m <- exp(logm)
          e <- size * reg_size[r] * exp(-cumsum(c(0, utils::head(exp(ax), -1)))) %o% (1 + 0.003 * seq_len(nt))
          d <- matrix(stats::rpois(length(m), m * e), length(a))
          rows[[length(rows) + 1L]] <- data.table(sex = sx, region = reg[r], quintile = paste0("Q", q),
                                                  year = rep(years, each = length(a)), age = rep(a, nt),
                                                  deaths = as.vector(d), exposure = as.vector(e), m = as.vector(m))
        }
      }
    }
    long <- rbindlist(rows)
    if (ages == "grouped") {
      long[, `:=`(age_lo = ifelse(age >= 90, 90, 40 + 5 * floor((age - 40) / 5)))]
      long[, age_hi := ifelse(age_lo >= 90, Inf, age_lo + 5)]
    } else {
      long[, `:=`(age_lo = pmin(age, 100))]
      long[, age_hi := ifelse(age_lo >= 100, Inf, age_lo + 1)]
    }
    agg <- long[, list(deaths = sum(deaths), exposure = sum(exposure), m = sum(m * exposure) / sum(exposure)),
                by = c("sex", "region", "quintile", "year", "age_lo", "age_hi")]
    x <- mt_data(deaths = agg, exposure = agg, by = c("sex", "region", "quintile"), age = c("age_lo", "age_hi"))
    if (ages == "grouped") x <- mt_set_grid(x, age = age_groups(c(seq(40, 85, 5), 90)))
    g <- mt_grid(x)
    st <- mt_strata(x)
    tr <- array(NA_real_, c(length(g$age), length(g$period), nrow(st)),
                dimnames = list(age = g$age$labels, period = g$period$labels, stratum = st$label))
    agg[, label := paste(sex, region, quintile, sep = ".")]
    idx <- cbind(match(agg$age_lo, g$age$lower), match(agg$year, g$period$lower), match(agg$label, st$label))
    tr[idx] <- agg$m
    attr(x, "truth") <- tr
    x
  })
}

#' Simulated Lee-Carter data with known truth
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

#' Simulated deaths by cause
#'
#' Extends [mt_example()] with deaths by cause (`"cvd"`, `"cancer"`,
#' `"other"`) and the all-cause total (`"all"`), single years of age. Cause
#' shares change smoothly with age and time; deaths are split multinomially,
#' so causes add up to all causes exactly.
#'
#' @inheritParams mt_example
#' @return An [mt_data] object with strata `sex` x `cause`.
#' @examples
#' mt_example_causes(years = 2010:2019)
#' @export
mt_example_causes <- function(years = 1990:2019, size = 1e5, seed = 1) {
  x <- mt_example(years = years, deaths = "single", size = size, seed = seed)
  D <- mt_obs(x, "deaths")
  E <- mt_obs(x, "exposure")
  with_seed(seed + 1, {
    t0 <- (D$period_lo - min(years)) / 10
    a <- (pmin(D$age_lo, 100) - 30) / 10
    eta_cvd <- -1.5 + 0.35 * a - 0.25 * t0
    eta_can <- 0.3 - 0.15 * a + 0.05 * t0
    p <- cbind(exp(eta_cvd), exp(eta_can), 1)
    p <- p / rowSums(p)
    split <- t(vapply(seq_len(nrow(D)), function(i) as.numeric(stats::rmultinom(1, D$value[i], p[i, ])), numeric(3)))
  })
  causes <- c("cvd", "cancer", "other")
  base <- D[, list(sex, period_lo, age_lo, age_hi)]
  long <- rbind(cbind(base, cause = "all", deaths = D$value),
                rbindlist(lapply(1:3, function(k) cbind(base, cause = causes[k], deaths = split[, k]))))
  ex <- E[, list(sex, period_lo, age_lo, age_hi, exposure = value)]
  ex <- rbindlist(lapply(c("all", causes), function(cc) cbind(ex, cause = cc)))
  mt_data(deaths = long, exposure = ex, by = c("sex", "cause"), age = c("age_lo", "age_hi"), period = "period_lo")
}
