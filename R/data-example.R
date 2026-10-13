#' Simulated mortality data with known truth
#'
#' Generates a two-stratum dataset (women and men) from a Lee-Carter-type
#' model with Gompertz-Makeham age pattern, a random walk with drift for the
#' period index, and a population with cohort-size waves. Deaths are Poisson
#' draws. By default deaths are reported in 5-year age groups with an open
#' 90+ interval and exposures in single years of age with an open 100+
#' interval, as in many official statistics.
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
  set.seed(seed)
  expr
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
        for (t in 2:nt) dev[t] <- 0.7 * dev[t - 1] + stats::rnorm(1, 0, 0.02)
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
