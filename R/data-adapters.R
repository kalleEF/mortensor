#' Convert other mortality data formats
#'
#' Methods for long data frames, `StMoMoData` objects from \pkg{StMoMo} (for
#' example `StMoMo::EWMaleData`, Human Mortality Database data for England
#' and Wales) and `demogdata` objects from \pkg{demography} (for example
#' `demography::fr.mort`, HMD data for France). For `demogdata`, deaths are
#' reconstructed as rate times exposure and rounded.
#'
#' @param x Object to convert.
#' @param ... Passed to methods; for data frames, arguments of [mt_data()].
#' @return An [mt_data] object.
#' @examples
#' if (requireNamespace("StMoMo", quietly = TRUE)) {
#'   ew <- as_mt_data(StMoMo::EWMaleData)
#'   ew
#' }
#' @export
as_mt_data <- function(x, ...) UseMethod("as_mt_data")

#' @export
as_mt_data.mt_data <- function(x, ...) x

#' @rdname as_mt_data
#' @param deaths,exposure Column names of deaths and exposures in a data
#'   frame.
#' @export
as_mt_data.data.frame <- function(x, deaths = "deaths", exposure = "exposure", ...) {
  mt_data(deaths = x, exposure = if (exposure %in% names(x)) x else NULL,
          value = c(deaths = deaths, exposure = exposure), ...)
}

#' @rdname as_mt_data
#' @export
as_mt_data.StMoMoData <- function(x, ...) {
  if (!identical(x$type, "central")) {
    mt_abort("Only central exposures are supported; convert with StMoMo::initial2central().", class = "adapter")
  }
  long <- matrices_to_long(x$Dxt, x$Ext, x$ages, x$years)
  long[, sex := x$series]
  mt_data(deaths = long, exposure = long, by = "sex", ...)
}

#' @rdname as_mt_data
#' @param series Series of a `demogdata` object to keep (default all except
#'   `"total"`).
#' @param years Optional subset of years.
#' @export
as_mt_data.demogdata <- function(x, series = NULL, years = NULL, ...) {
  series <- series %||% setdiff(names(x$rate), "total")
  yrs <- years %||% x$year
  out <- rbindlist(lapply(series, function(s) {
    r <- x$rate[[s]][, as.character(yrs), drop = FALSE]
    p <- x$pop[[s]][, as.character(yrs), drop = FALSE]
    r[is.na(r) & p == 0] <- 0
    long <- matrices_to_long(round(r * p), p, x$age, yrs)
    long[, sex := s]
  }))
  open <- any(grepl("\\+", rownames(x$rate[[series[1]]])))
  mt_data(deaths = out, exposure = out, by = "sex", age_open = if (open) max(x$age) else NULL, ...)
}

matrices_to_long <- function(D, E, ages, years) {
  data.table(age = rep(as.numeric(ages), length(years)), year = rep(as.numeric(years), each = length(ages)),
             deaths = as.vector(D), exposure = as.vector(E))
}
