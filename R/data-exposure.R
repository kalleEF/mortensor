#' Exposures from population counts
#'
#' Approximates person-years at risk from population stocks and stores them
#' as the exposure table:
#'
#' * `"midyear"`: the mid-year population of a calendar year is used as the
#'   exposure of that year;
#' * `"jan1"`: the exposure of year t is the mean of the January 1
#'   populations of years t and t + 1 in the same age interval. The last year
#'   has no successor and is dropped.
#' * `"dec31"`: populations on 31 December of year t (as published by the
#'   German Federal Statistical Office); the exposure of year t is the mean
#'   of the populations on 31 December of years t - 1 and t. The first year
#'   has no predecessor and is dropped.
#'
#' Both are standard approximations; the HMD uses more refined cohort-based
#' formulas. The approximation used is recorded in the provenance column.
#'
#' @param data An [mt_data] object with a population table.
#' @param type `"midyear"`, `"jan1"` or `"dec31"`.
#' @return The `mt_data` object with an exposure table.
#' @examples
#' p <- data.frame(year = 2001:2003, age = 70, population = c(1000, 990, 985))
#' x <- mt_as_exposure(mt_add_population(mt_data(), p), type = "jan1")
#' mt_obs(x, "exposure")
#' @export
mt_as_exposure <- function(data, type = c("midyear", "jan1", "dec31")) {
  check_class(data, "mt_data", "data")
  type <- match.arg(type)
  P <- data$obs$population
  if (type == "dec31" && !is.null(P)) {
    # the population on 31 December of year t is the population on 1 January of t + 1
    P <- copy(P)[, `:=`(period_lo = period_lo + 1, period_hi = period_hi + 1)]
  }
  if (is.null(P)) mt_abort("No population table.", class = "missing_measure")
  if (any(P$period_hi - P$period_lo != 1)) {
    mt_abort("Population must be given for single calendar years.", class = "exposure")
  }
  E <- if (type == "midyear") {
    copy(P)
  } else {
    nxt <- copy(P)[, period_lo := period_lo - 1][, period_hi := period_hi - 1]
    j <- merge(P, nxt, by = c("stratum", "age_lo", "age_hi", "period_lo", "period_hi"), suffixes = c("", ".next"))
    j[, value := (value + value.next) / 2][, c("value.next", "provenance.next") := NULL]
    j
  }
  E[, provenance := paste0("population_", type)]
  E[, value := as.numeric(value)]
  setorderv(E, c("stratum", "period_lo", "age_lo"))
  data$obs$exposure <- E[]
  data$grid <- NULL
  data
}
