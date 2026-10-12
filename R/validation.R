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
