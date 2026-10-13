#' Read deaths and populations from GENESIS (Destatis) flat files
#'
#' Reads tables of the German Federal Statistical Office in the "flat" CSV
#' format of GENESIS-Online and the regional database (columns `time`,
#' `1_variable_code`, `1_variable_attribute_label`, ..., `value`): deaths
#' (for example table 23211-0002) and populations (for example table
#' 12411-0013, on 31 December).
#'
#' Variables are recognised by their codes: `GES` (sex), codes starting with
#' `ALT` (age, German labels such as `"unter 1 Jahr"`, `"15 bis unter 20
#' Jahre"`, `"85 Jahre und mehr"`, `"68-Jaehrige"` with umlaut), `DLAND`,
#' `KREISE` or `REGBEZ` (region), codes starting with `TODUR` (cause of
#' death; the total "Insgesamt" becomes level `"all"`, the other levels are
#' the GENESIS codes, their labels are in attribute `"cause_labels"`), and
#' `DINSG` (Germany as a whole, dropped). Other variables are kept as
#' stratum variables named by their code. Values `"-"` (nothing) are zero;
#' `"."`, `"..."`, `"x"` and `"/"` are missing. Stratum variables that only
#' the population table has (for example regions, when the deaths are
#' national) are summed out of the populations.
#'
#' The reader makes no analysis decisions: populations stay populations, and
#' the periods are those of the source (the year of the reference date).
#' The next steps are usually [mt_as_exposure()] with `type = "dec31"`,
#' [mt_regroup()] with `age = "deaths"` (single-year populations to the age
#' groups of the deaths), [mt_fill_suppressed()] and [mt_group_levels()]
#' (cause groups, with `total = "all"` for a residual group).
#'
#' @param deaths Path of a flat CSV file with deaths, or a data frame read
#'   from one.
#' @param population Path or data frame with populations, or `NULL`.
#' @param years Optional subset of years of the deaths (populations of the
#'   year before are kept for exposures from 31 December populations).
#' @return An [mt_data] object with deaths and populations.
#' @examples
#' d <- data.frame(time = 2020, `1_variable_code` = "GES", `1_variable_attribute_label` = "weiblich",
#'                 `2_variable_code` = "ALT027", `2_variable_attribute_label` = "85 Jahre und mehr",
#'                 value = "120", check.names = FALSE)
#' p <- data.frame(time = c("2019-12-31", "2020-12-31"), `1_variable_code` = "GES",
#'                 `1_variable_attribute_label` = "weiblich", `2_variable_code` = "ALT103",
#'                 `2_variable_attribute_label` = "85 Jahre und mehr", value = c(900, 950),
#'                 check.names = FALSE)
#' x <- mt_read_genesis(d, p)
#' x <- mt_as_exposure(x, type = "dec31")
#' mt_obs(x, "exposure")
#' @export
mt_read_genesis <- function(deaths, population = NULL, years = NULL) {
  d <- genesis_long(genesis_table(deaths))
  cause_labels <- attr(d, "cause_labels")
  if (!is.numeric(d$time)) d[, time := as.integer(substr(time, 1, 4))]
  if (!is.null(years)) d <- d[time %in% years]
  d_vars <- setdiff(names(d), c("time", "age", "value"))
  setnames(d, "value", "deaths")
  p <- NULL
  if (!is.null(population)) {
    p <- genesis_long(genesis_table(population))
    p[, time := as.integer(substr(as.character(time), 1, 4))]
    p_vars <- setdiff(names(p), c("time", "age", "value"))
    extra <- setdiff(p_vars, d_vars)
    if (length(extra)) p <- p[, list(value = sum(value)), by = c("time", "age", intersect(d_vars, p_vars))]
    if (!is.null(years)) p <- p[time %in% c(years, min(years) - 1L)]
    setnames(p, "value", "population")
  }
  x <- mt_data(deaths = d, population = p, by = if (length(d_vars)) d_vars else NULL, age = "age",
               period = "time")
  if (!is.null(cause_labels)) attr(x, "cause_labels") <- cause_labels
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
