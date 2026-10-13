#' Convert data to demography and StMoMo objects
#'
#' Exports observed deaths and exposures for cross-checks with the original
#' packages: `as_demogdata()` returns a `demogdata` object of the
#' demography package with one series per stratum (rates and exposures as
#' population), and `as_stmomodata()` a `StMoMoData` object of the StMoMo
#' package for one stratum (central exposures). The data must be on the
#' grid (deaths and exposures in the same age groups and single years).
#' [as_mt_data()] converts in the other direction.
#'
#' @param data An [mt_data] object.
#' @param stratum Label of the stratum for `as_stmomodata()` (default: the
#'   first).
#' @param name Name of the population.
#' @return A `demogdata` or `StMoMoData` object.
#' @examples
#' if (requireNamespace("demography", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 2001:2010)
#'   dd <- as_demogdata(x)
#'   names(dd$rate)
#' }
#' @export
as_demogdata <- function(data, name = "mortensor") {
  check_class(data, "mt_data", "data")
  need_pkg("demography", "to create demogdata objects")
  md <- model_inputs(data, list(min_periods = 1, label = "Conversion"))
  rates <- lapply(seq_len(dim(md$D)[3]), function(s) {
    m <- md$D[, , s] / md$E[, , s]
    matrix(m, dim(md$D)[1], dimnames = list(md$age$lower, md$periods))
  })
  pops <- lapply(seq_len(dim(md$E)[3]), function(s) {
    matrix(md$E[, , s], dim(md$E)[1], dimnames = list(md$age$lower, md$periods))
  })
  labs <- data$strata$label
  out <- demography::demogdata(rates[[1]], pops[[1]], ages = md$age$lower, years = md$periods,
                               type = "mortality", label = name, name = labs[1])
  for (s in seq_along(labs)[-1]) {
    out$rate[[labs[s]]] <- rates[[s]]
    out$pop[[labs[s]]] <- pops[[s]]
  }
  out
}

#' @rdname as_demogdata
#' @export
as_stmomodata <- function(data, stratum = NULL, name = "mortensor") {
  check_class(data, "mt_data", "data")
  need_pkg("StMoMo", "to create StMoMoData objects")
  md <- model_inputs(data, list(min_periods = 1, label = "Conversion"))
  labs <- data$strata$label
  s <- if (is.null(stratum)) 1L else match(stratum, labs)
  if (is.na(s)) mt_abort(sprintf("Unknown stratum; available: %s.", paste(labs, collapse = ", ")), class = "argument")
  D <- matrix(md$D[, , s], dim(md$D)[1], dimnames = list(md$age$lower, md$periods))
  E <- matrix(md$E[, , s], dim(md$E)[1], dimnames = list(md$age$lower, md$periods))
  structure(list(Dxt = D, Ext = E, ages = md$age$lower, years = md$periods, type = "central",
                 series = labs[s], label = name), class = "StMoMoData")
}
