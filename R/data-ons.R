#' Read ONS deaths and populations by deprivation
#'
#' Reads death registrations and population estimates in the layout of the
#' Office for National Statistics (ONS) ad hoc tables by sex, age group and
#' Index of Multiple Deprivation (IMD) decile, the reader of mortensor 0.1.
#' Deaths are on sheet `"1"` and populations on sheet `"2"`; the header row
#' is the first row that contains `"Sex"`. Columns are recognised by name:
#' year, sex, age, the IMD decile (`"imd"`, `"decile"`, `"deprivation"`) and,
#' optionally, a region column. All other columns of the deaths sheet are
#' causes of death; one of them is all causes.
#'
#' Deciles are combined to quintiles by default (deciles 1 and 2 to quintile
#' 1, and so on; in the English IMD decile 1 is the most deprived). A
#' deprivation column whose name contains "quintile" is read as quintiles. The
#' mid-year population is used as exposure. Population years that are
#' missing are filled with the last earlier year (with a warning), as in
#' mortensor 0.1.
#'
#' @param deaths Path to the Excel file with deaths, or a data frame in the
#'   same layout.
#' @param population Path(s) to Excel files with populations (stacked), or a
#'   data frame.
#' @param causes `NULL` (all causes only) or a named list: each element lists
#'   the cause columns that are summed into one level of the stratum variable
#'   `cause`; the value `"rest"` gives all causes minus the other levels.
#'   Level `"all"` (all causes) is added. Column names are matched after
#'   conversion to lower case with underscores (`"Lung cancer"` becomes
#'   `"lung_cancer"`).
#' @param deprivation `"quintile"` or `"decile"`.
#' @param region Name of a region column, or `NULL`.
#' @param years Optional subset of years.
#' @param fill_population Fill missing population years with the last
#'   earlier year.
#' @param sheets Sheet names of deaths and populations.
#' @return An [mt_data] object with strata `sex`, the deprivation variable
#'   and, if requested, `region` and `cause`.
#' @examples
#' d <- data.frame(Year = 2020, Sex = rep(c("Males", "Females"), each = 4),
#'                 Age = rep(c("85-89", "90+"), 4), IMD = rep(rep(1:2, each = 2), 2),
#'                 `All causes` = 10:17, `Lung cancer` = 1, check.names = FALSE)
#' p <- data.frame(Year = 2020, Sex = d$Sex, Age = d$Age, IMD = d$IMD, Population = 1000)
#' x <- mt_read_ons(d, p, causes = list(lung = "lung_cancer", other = "rest"))
#' mt_strata(x)
#' @export
mt_read_ons <- function(deaths, population, causes = NULL, deprivation = c("quintile", "decile"), region = NULL,
                        years = NULL, fill_population = TRUE, sheets = c(deaths = "1", population = "2")) {
  deprivation <- match.arg(deprivation)
  dt <- ons_keys(ons_table(deaths, sheets[["deaths"]]), region, "deaths")
  pops <- if (is.data.frame(population)) list(population) else as.list(population)
  pop <- rbindlist(lapply(pops, function(p) ons_keys(ons_table(p, sheets[["population"]]), region, "population")))
  keys <- c("year", "sex", "age", "dep", if (!is.null(region)) "region")
  pop <- unique(pop, by = keys)
  if (!is.null(years)) {
    dt <- dt[year %in% years]
    pop <- pop[year %in% years]
  }
  vals <- setdiff(names(pop), keys)
  if (!length(vals)) mt_abort("The population table has no population column.", class = "ons")
  setnames(pop, vals[length(vals)], "population")
  pop[, population := as.numeric(population)]
  causes_dt <- ons_causes(dt, keys, causes)
  if (attr(dt, "dep_type") == "quintile") {
    if (deprivation == "decile") mt_abort("The data have quintiles, not deciles.", class = "ons")
  } else if (deprivation == "quintile") {
    causes_dt[, dep := ceiling(dep / 2)]
    pop[, dep := ceiling(dep / 2)]
  }
  by_d <- c(keys, "cause")
  d_agg <- causes_dt[, list(deaths = sum(deaths)), by = by_d]
  p_agg <- pop[, list(population = sum(population, na.rm = TRUE)), by = keys]
  out <- merge(d_agg, p_agg, by = keys, all.x = TRUE)
  miss <- sort(unique(out[is.na(population), year]))
  if (length(miss)) {
    if (!fill_population) {
      mt_abort(sprintf("Population missing for %s.", paste(miss, collapse = ", ")), class = "ons")
    }
    for (yr in miss) {
      src <- suppressWarnings(max(p_agg$year[p_agg$year <= yr]))
      if (!is.finite(src)) mt_abort(sprintf("No population at or before %d.", yr), class = "ons")
      fill <- p_agg[year == src][, year := yr]
      out[fill, on = c(keys), population := ifelse(is.na(population), i.population, population)]
    }
    mt_warn(sprintf("Population for %s filled with the last earlier year.", paste(miss, collapse = ", ")),
            class = "ons_fill")
  }
  out[, sex := ifelse(grepl("^f", tolower(sex)), "female", ifelse(grepl("^m", tolower(sex)), "male", tolower(sex)))]
  dep_name <- deprivation
  out[, (dep_name) := paste0(if (deprivation == "quintile") "Q" else "D", dep)]
  by <- c("sex", if (!is.null(region)) "region", dep_name, if (!is.null(causes)) "cause")
  out[, exposure := population]
  x <- mt_data(deaths = out, exposure = out, by = by, age = "age", period = "year")
  ag <- unique(x$obs$deaths[, c("age_lo", "age_hi")])[order(age_lo)]
  mt_set_grid(x, age = mt_axis(ag$age_lo, ag$age_hi, kind = "age"))
}

ons_table <- function(x, sheet) {
  if (is.data.frame(x)) return(as.data.table(x))
  if (!is.character(x) || length(x) != 1 || !file.exists(x)) {
    mt_abort("Expected a data frame or the path of an existing Excel file.", class = "ons")
  }
  need_pkg("readxl", "to read Excel files")
  peek <- suppressMessages(as.data.frame(readxl::read_excel(x, sheet = sheet, col_names = FALSE, n_max = 30)))
  hit <- which(apply(peek, 1, function(r) any(tolower(trimws(as.character(r))) == "sex", na.rm = TRUE)))
  if (!length(hit)) {
    mt_abort(sprintf("No header row with \"Sex\" on sheet %s of %s.", sheet, basename(x)), class = "ons")
  }
  as.data.table(suppressMessages(readxl::read_excel(x, sheet = sheet, skip = hit[1] - 1)))
}

ons_clean_names <- function(x) gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(trimws(x))))

ons_keys <- function(dt, region, what) {
  nm <- ons_clean_names(names(dt))
  find <- function(pattern, label) {
    i <- grep(pattern, nm)
    if (!length(i)) mt_abort(sprintf("No %s column in the %s table.", label, what), class = "ons")
    i[1]
  }
  idx <- c(year = find("^year", "year"), sex = find("^sex", "sex"), age = find("^age", "age"),
           dep = find("imd|decile|deprivation|quintile", "deprivation"))
  if (!is.null(region)) {
    i <- which(names(dt) == region | nm == ons_clean_names(region))
    if (!length(i)) mt_abort(sprintf("Region column `%s` not found in the %s table.", region, what), class = "ons")
    idx <- c(idx, region = i[1])
  }
  dep_type <- if (grepl("quintile", nm[idx[["dep"]]])) "quintile" else "decile"
  setnames(dt, names(dt)[idx], names(idx))
  rest <- setdiff(seq_along(dt), idx)
  setnames(dt, rest, nm[rest])
  dt <- dt[!is.na(year) & !is.na(sex)]
  dt[, year := as.integer(year)]
  dt[, dep := as.integer(dep)]
  dt[, age := as.character(age)]
  dt[, sex := as.character(sex)]
  if (!is.null(region)) dt[, region := as.character(region)]
  setattr(dt, "dep_type", dep_type)
  dt
}

ons_causes <- function(dt, keys, causes) {
  cols <- setdiff(names(dt), keys)
  all_col <- grep("^all_?causes?$", cols, value = TRUE)
  if (!length(all_col)) mt_abort("No all-cause column (\"All causes\") in the deaths table.", class = "ons")
  num <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    v[is.na(v)] <- 0
    v
  }
  all <- num(dt[[all_col[1]]])
  base <- dt[, keys, with = FALSE]
  if (is.null(causes)) return(base[, `:=`(cause = "all", deaths = all)])
  if (is.null(names(causes)) || any(!nzchar(names(causes)))) causes <- stats::setNames(causes, unlist(causes))
  if ("all" %in% names(causes)) mt_abort("Level \"all\" is added automatically.", class = "ons")
  sums <- list()
  for (lev in names(causes)) {
    cc <- ons_clean_names(causes[[lev]])
    if (identical(cc, "rest")) next
    miss <- setdiff(cc, cols)
    if (length(miss)) mt_abort(sprintf("Cause column(s) not found: %s.", paste(miss, collapse = ", ")), class = "ons")
    sums[[lev]] <- Reduce(`+`, lapply(cc, function(k) num(dt[[k]])))
  }
  rest <- names(causes)[vapply(causes, function(cc) identical(ons_clean_names(cc), "rest"), TRUE)]
  if (length(rest) > 1) mt_abort("Only one level can be \"rest\".", class = "ons")
  if (length(rest)) {
    r <- all - Reduce(`+`, sums, 0)
    if (any(r < -0.5)) mt_warn("Listed causes exceed all causes in some cells; the rest is set to zero there.",
                               class = "ons_rest")
    sums[[rest]] <- pmax(r, 0)
  }
  rbindlist(c(list(copy(base)[, `:=`(cause = "all", deaths = all)]),
              lapply(names(sums), function(lev) copy(base)[, `:=`(cause = lev, deaths = sums[[lev]])])))
}
