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

#' Human Mortality Database tables
#'
#' `mt_from_hmd()` converts HMD tables in the layout returned by
#' `HMDHFDplus::readHMDweb()` (columns `Year`, `Age`, `Female`, `Male`,
#' `Total`, and `OpenInterval`) into an [mt_data] object with strata `sex`.
#' Without column `OpenInterval`, an age written with `+` (`"110+"`) is the
#' open interval.
#' `mt_read_hmd()` downloads `Deaths_1x1` and `Exposures_1x1` for one country
#' with \pkg{HMDHFDplus}; an HMD account is needed (<https://www.mortality.org>).
#'
#' @param deaths,exposures Data frames from `HMDHFDplus::readHMDweb()`.
#' @param sexes Columns to use.
#' @param years Optional subset of years.
#' @param country HMD country code, for example `"SWE"`.
#' @param username,password HMD credentials (default: environment variables
#'   `HMD_USER` and `HMD_PASSWORD`).
#' @return An [mt_data] object.
#' @examples
#' d <- data.frame(Year = 2000, Age = 0:2, Female = c(10, 1, 1), Male = c(12, 2, 1),
#'                 Total = c(22, 3, 2), OpenInterval = c(FALSE, FALSE, TRUE))
#' e <- transform(d, Female = 1000, Male = 1000, Total = 2000)
#' mt_from_hmd(d, e)
#' @export
mt_from_hmd <- function(deaths, exposures, sexes = c("Female", "Male"), years = NULL) {
  conv <- function(x, nm) {
    x <- as.data.table(x)
    if (!is.null(years)) x <- x[x$Year %in% years]
    open <- if ("OpenInterval" %in% names(x)) x$OpenInterval else grepl("+", as.character(x$Age), fixed = TRUE)
    age <- suppressWarnings(as.numeric(gsub("\\+", "", as.character(x$Age))))
    rbindlist(lapply(sexes, function(s) {
      data.table(sex = tolower(s), year = x$Year, age_lo = age, age_hi = ifelse(open, Inf, age + 1),
                 value = as.numeric(x[[s]]))
    }))[, (nm) := value][, value := NULL][]
  }
  d <- conv(deaths, "deaths")
  e <- conv(exposures, "exposure")
  mt_data(deaths = d, exposure = e, by = "sex", age = c("age_lo", "age_hi"))
}

#' @rdname mt_from_hmd
#' @export
mt_read_hmd <- function(country, username = Sys.getenv("HMD_USER"), password = Sys.getenv("HMD_PASSWORD"),
                        years = NULL) {
  need_pkg("HMDHFDplus", "to download Human Mortality Database data")
  if (!nzchar(username) || !nzchar(password)) {
    mt_abort("HMD credentials are needed (arguments or HMD_USER and HMD_PASSWORD).", class = "adapter")
  }
  d <- HMDHFDplus::readHMDweb(country, "Deaths_1x1", username, password)
  e <- HMDHFDplus::readHMDweb(country, "Exposures_1x1", username, password)
  mt_from_hmd(d, e, years = years)
}
