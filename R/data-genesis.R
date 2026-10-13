#' Read deaths and populations from GENESIS (Destatis) flat files
#'
#' Reads tables of the German Federal Statistical Office in the "flat" CSV
#' format of GENESIS-Online and the regional database (columns `time`,
#' `1_variable_code`, `1_variable_attribute_label`, ..., `value`): deaths
#' (for example table 23211-0002) and population on 31 December (for example
#' table 12411-0013).
#'
#' Variables are recognised by their codes: `GES` (sex), codes starting with
#' `ALT` (age, German labels such as `"unter 1 Jahr"`, `"15 bis unter 20
#' Jahre"`, `"85 Jahre und mehr"`, `"68-Jaehrige"` with umlaut), `DLAND`,
#' `KREISE` or `REGBEZ` (region), codes starting with `TODUR` (cause of
#' death), and `DINSG` (Germany as a whole, dropped). Other variables are
#' kept as stratum variables named by their code. Values `"-"` (nothing) are
#' zero; `"."`, `"..."`, `"x"` and `"/"` are missing.
#'
#' Exposures for year t are the mean of the populations on 31 December of
#' years t - 1 and t. Stratum variables that only the population table has
#' (for example regions, when the deaths are national) are summed out; those
#' that only the deaths table has (causes) share the exposure.
#'
#' @param deaths Path of a flat CSV file with deaths, or a data frame read
#'   from one.
#' @param population Path or data frame with populations on 31 December.
#' @param causes `FALSE` (all causes only) or `TRUE` (all causes as level
#'   `"all"`, the listed causes, and `"rest"`, all causes minus the listed
#'   ones) when the deaths table has a cause variable.
#' @param years Optional subset of years.
#' @return An [mt_data] object; the grid has the age groups of the deaths.
#'   With `causes = TRUE`, the levels of `cause` are the GENESIS codes (for
#'   example `"TDU-021"`), and attribute `"cause_labels"` holds their labels.
#' @examples
#' d <- data.frame(time = 2020, `1_variable_code` = "GES", `1_variable_attribute_label` = "weiblich",
#'                 `2_variable_code` = "ALT027", `2_variable_attribute_label` = "85 Jahre und mehr",
#'                 value = "120", check.names = FALSE)
#' p <- data.frame(time = c("2019-12-31", "2020-12-31"), `1_variable_code` = "GES",
#'                 `1_variable_attribute_label` = "weiblich", `2_variable_code` = "ALT103",
#'                 `2_variable_attribute_label` = "85 Jahre und mehr", value = c(900, 950),
#'                 check.names = FALSE)
#' x <- mt_read_genesis(d, p)
#' mt_obs(x, "exposure")
#' @export
mt_read_genesis <- function(deaths, population, causes = FALSE, years = NULL) {
  check_flag(causes, "causes")
  d <- genesis_long(genesis_table(deaths))
  cause_labels <- attr(d, "cause_labels")
  p <- genesis_long(genesis_table(population))
  if (!is.numeric(d$time)) d[, time := as.integer(substr(time, 1, 4))]
  # population on 31 December of year y is the population at the start of y + 1
  p[, time := as.integer(substr(as.character(time), 1, 4)) + 1L]
  d_vars <- setdiff(names(d), c("time", "age", "value"))
  p_vars <- setdiff(names(p), c("time", "age", "value"))
  by <- intersect(d_vars, p_vars)
  extra <- setdiff(p_vars, d_vars)
  if (length(extra)) p <- p[, list(value = sum(value)), by = c("time", "age", by)]
  if (!is.null(years)) d <- d[time %in% years]
  yrs <- unique(d$time)
  yrs <- yrs[yrs %in% p$time & (yrs + 1L) %in% p$time]
  if (!length(yrs)) {
    mt_abort("No year has deaths and populations at the start and end of the year.", class = "genesis")
  }
  d <- d[time %in% yrs]
  p <- p[time %in% c(yrs, yrs + 1L)]
  if ("cause" %in% d_vars) {
    if (!any(d$cause == "all")) mt_abort("The cause variable has no total (\"Insgesamt\").", class = "genesis")
    if (!causes) d <- d[cause == "all"][, cause := NULL]
  }
  if (anyNA(d$value)) {
    mt_warn(sprintf("%d death cell(s) are missing in the source and are dropped.", sum(is.na(d$value))),
            class = "genesis")
    d <- d[!is.na(value)]
  }
  if ("cause" %in% names(d)) {
    keys <- c("time", "age", setdiff(d_vars, "cause"))
    rest <- d[, list(value = sum(value[cause == "all"]) - sum(value[cause != "all"])), by = keys]
    if (any(rest$value < -0.5)) mt_warn("Listed causes exceed all causes in some cells.", class = "genesis")
    rest[, `:=`(cause = "rest", value = pmax(value, 0))]
    d <- rbind(d, rest, use.names = TRUE)
  }
  d_vars <- setdiff(names(d), c("time", "age", "value"))
  if ("cause" %in% d_vars) {
    lev <- unique(d$cause)
    p <- rbindlist(lapply(lev, function(l) copy(p)[, cause := l]))
  }
  setnames(d, "value", "deaths")
  setnames(p, "value", "population")
  x <- mt_data(deaths = d, population = p, by = if (length(d_vars)) d_vars else NULL, age = "age", period = "time")
  x <- mt_as_exposure(x, type = "jan1")
  x$obs$population <- NULL
  ag <- unique(x$obs$deaths[, c("age_lo", "age_hi")])[order(age_lo)]
  # exposures (often single years) are summed to the age groups of the deaths
  x <- mt_regroup(x, age = ag$age_lo)
  if (causes && !is.null(cause_labels)) attr(x, "cause_labels") <- cause_labels
  x
}

genesis_table <- function(x) {
  if (is.data.frame(x)) return(as.data.table(x))
  if (!is.character(x) || length(x) != 1 || !file.exists(x)) {
    mt_abort("Expected a data frame or the path of an existing CSV file.", class = "genesis")
  }
  fread(x, encoding = "UTF-8", colClasses = list(character = "value"), showProgress = FALSE)
}

genesis_long <- function(dt) {
  if (!all(c("time", "value") %in% names(dt))) {
    mt_abort("Not a GENESIS flat table (no `time`/`value`).", class = "genesis")
  }
  k <- sort(as.integer(sub("_variable_code$", "", grep("^[0-9]+_variable_code$", names(dt), value = TRUE))))
  out <- data.table(time = dt$time)
  v <- trimws(as.character(dt$value))
  num <- suppressWarnings(as.numeric(gsub(",", ".", v)))
  num[v == "-"] <- 0
  num[v %in% c(".", "...", "x", "/")] <- NA
  for (i in k) {
    code <- unique(dt[[sprintf("%d_variable_code", i)]])
    lab <- dt[[sprintf("%d_variable_attribute_label", i)]]
    if (length(code) != 1) mt_abort("A variable column holds several variables.", class = "genesis")
    if (code == "DINSG") next
    if (code == "GES") {
      l <- tolower(enc2utf8(lab))
      out[, sex := ifelse(grepl("^m", l, useBytes = TRUE), "male", ifelse(grepl("^w", l, useBytes = TRUE), "female",
                                                                         "total"))]
    } else if (startsWith(code, "ALT")) {
      out[, age := enc2utf8(lab)]
    } else if (code %in% c("DLAND", "KREISE", "REGBEZ")) {
      out[, region := enc2utf8(lab)]
    } else if (startsWith(code, "TODUR")) {
      cc <- dt[[sprintf("%d_variable_attribute_code", i)]]
      out[, cause := ifelse(lab == "Insgesamt" | !nzchar(cc), "all", cc)]
      labs <- unique(data.frame(cause = out$cause, label = enc2utf8(lab)))
      setattr(out, "cause_labels", labs[labs$cause != "all", ])
    } else {
      out[, (tolower(code)) := enc2utf8(lab)]
    }
  }
  if (!"age" %in% names(out)) mt_abort("The table has no age variable.", class = "genesis")
  out[, value := num]
  out[]
}
