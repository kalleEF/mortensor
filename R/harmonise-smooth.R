#' Smooth death rates on the grid
#'
#' Poisson P-spline smoothing of death rates with the exposures as offset
#' (Currie, Durban and Eilers 2004 for the 2D case): the same engine as
#' [mt_ungroup()] with every grid cell as its own observation. Deaths are
#' replaced by exposure times the smoothed rate; exposures are unchanged.
#' Deaths and exposures must be available on the grid (use [mt_ungroup()]
#' first for grouped data).
#'
#' @param data An [mt_data] object.
#' @param dims 1 (over age, separately per period) or 2 (age-period
#'   surface).
#' @inheritParams mt_ungroup
#' @return An `mt_data` object with smoothed deaths; diagnostics as attribute
#'   `"smooth"`.
#' @references
#' Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
#' mortality rates. *Statistical Modelling* 4(4): 279-298.
#' \doi{10.1191/1471082X04st080oa}
#' @examples
#' x <- mt_example(years = 2001:2003, deaths = "single")
#' s <- mt_smooth(x)
#' head(attr(s, "smooth"))
#' @export
mt_smooth <- function(data, dims = 1, lambda = NULL, criterion = c("bic", "aic"), kr = 2, deg = 3,
                      pord = 2, lambda_range = c(1e-2, 1e8)) {
  check_class(data, "mt_data", "data")
  criterion <- match.arg(criterion)
  check_number(dims, "dims", lower = 1, upper = 2, integerish = TRUE)
  g <- mt_grid(data)
  for (m in c("deaths", "exposure")) {
    if (is.null(data$obs[[m]])) mt_abort(sprintf("Smoothing needs a %s table.", m), class = "missing_measure")
    if (needs_ungroup(data$obs[[m]], g)) {
      mt_abort(sprintf("The %s table is coarser than the grid; use mt_ungroup() first.", m), class = "resolution")
    }
  }
  ctl <- list(lambda = lambda, criterion = criterion, kr = kr, deg = deg, pord = pord,
              lambda_range = lambda_range, preserve_totals = FALSE, keep_observed = FALSE,
              provenance = "smoothed", method_label = "pspline")
  res <- ungroup_rates(data$obs$deaths, data$obs$exposure, g, data$omega, ctl, dims)
  data$obs$deaths <- res$obs
  data$grid <- g
  attr(data, "smooth") <- res$diag
  data
}
