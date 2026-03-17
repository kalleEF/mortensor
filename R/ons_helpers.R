#' Read mortality data from ONS-format Excel files
#'
#' Reads UK Office for National Statistics (ONS) death registration and
#' population data from Excel files and returns a [mort_data] object.
#' Handles IMD decile-to-quintile mapping, cause column aggregation,
#' population forward-filling, and array construction.
#'
#' @note This function is specific to the UK ONS data format. For generic
#'   mortality data, use [read_mortality()] instead.
#'
#' @param deaths_file Path to the Excel file with death registration data
#'   (Sheet 1).
#' @param pop_files Character vector of paths to population Excel files
#'   (Sheet 2, one file per year). Multiple files are stacked.
#' @param config A [mort_config] object specifying cause columns, age groups,
#'   year range, etc.
#' @return A [mort_data] object.
#' @export
read_ons_excel <- function(deaths_file, pop_files, config) {
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("Install the 'readxl' package to read Excel files.")

  deaths <- load_deaths_ons(deaths_file)
  pop_list <- lapply(pop_files, load_population_ons)
  pop <- data.table::rbindlist(pop_list)
  pop <- unique(pop, by = c("year", "sex", "age_group", "imd_decile"))

  dt <- prepare_ons_analysis_dt(deaths, pop, config)

  cause_short <- config$cause_short
  mx_col <- paste0("mx_", cause_short)
  ages   <- config$age_starts
  years  <- sort(unique(dt$year))
  sexes  <- c("Male", "Female")
  quints <- seq_len(5L)
  groups <- data.table::CJ(sex = sexes, quintile = quints)

  mx_arr  <- build_array(dt, mx_col, ages, years, groups)
  pop_arr <- build_array(dt, "population", ages, years, groups, floor = NULL)

  mort_data(
    mx     = mx_arr,
    pop    = pop_arr,
    ages   = ages,
    years  = years,
    groups = as.data.frame(groups),
    cause  = config$cause_label,
    age_labels = config$age_groups
  )
}

#' @rdname read_ons_excel
#' @usage # Deprecated alias
#' @export
read_mortality_excel <- function(deaths_file, pop_files, config) {
  .Deprecated("read_ons_excel")
  read_ons_excel(deaths_file, pop_files, config)
}


# ---- ONS-specific internal loading helpers ----

#' Find the header row in an ONS Excel sheet
#' @note ONS-format specific.
#' @keywords internal
find_header_skip <- function(path, sheet, search = "Sex") {
  peek <- suppressWarnings(
    as.data.frame(readxl::read_excel(path, sheet = sheet,
                                     col_names = FALSE, n_max = 20)))
  for (i in seq_len(nrow(peek))) {
    vals <- as.character(unlist(peek[i, ]))
    if (search %in% vals) return(i - 1L)
  }
  stop("Could not find header row containing '", search, "' in ",
       basename(path), " sheet ", sheet)
}

#' Load deaths from an ONS-format Excel file
#' @note ONS-format specific.
#' @keywords internal
load_deaths_ons <- function(path) {
  skip <- find_header_skip(path, sheet = "1")
  dt <- suppressWarnings(
    data.table::as.data.table(readxl::read_excel(path, sheet = "1", skip = skip)))

  std_names <- c(
    "year", "sex", "age_group", "imd_decile",
    "all_causes", "breast_cancer", "colorectal", "lung_cancer",
    "prostate", "other_cancer", "dementia", "ihd", "stroke",
    "t1_diabetes", "t2_diabetes", "diabetes_other", "copd",
    "afib", "heart_failure", "epilepsy", "ibs",
    "alc_liver", "nonalc_liver", "ckd", "anxiety_dep",
    "ra_other", "asthma", "hypertension", "hearing",
    "psychosis", "alc_specific", "constipation", "chronic_pain",
    "covid19", "suicide"
  )

  if (ncol(dt) == length(std_names)) {
    data.table::setnames(dt, std_names)
  } else {
    data.table::setnames(dt, 1:4, c("year", "sex", "age_group", "imd_decile"))
    orig <- names(dt)
    map <- list(
      ihd = "Ischaemic Heart", stroke = "Stroke",
      all_causes = "^All causes$", dementia = "Dementia",
      copd = "COPD", lung_cancer = "Lung cancer",
      breast_cancer = "Breast cancer", colorectal = "Colorectal",
      prostate = "Prostate", other_cancer = "All other cancers",
      afib = "Atrial fibrillation", heart_failure = "Heart failure",
      hypertension = "Hypertension",
      t1_diabetes = "Type 1 Diabetes", t2_diabetes = "Type 2 Diabetes",
      diabetes_other = "Diabetes other",
      alc_liver = "Alcohol-related Liver", nonalc_liver = "Non-alcohol",
      alc_specific = "Alcohol specific", suicide = "Suicide"
    )
    for (nm in names(map)) {
      idx <- grep(map[[nm]], orig, ignore.case = TRUE)[1]
      if (!is.na(idx)) data.table::setnames(dt, idx, nm)
    }
  }

  dt[, year := as.integer(year)]
  dt[, imd_decile := as.integer(imd_decile)]
  dt[, sex := data.table::fifelse(tolower(sex) == "males", "Male", "Female")]

  cause_cols <- setdiff(names(dt), c("year", "sex", "age_group", "imd_decile"))
  for (cc in cause_cols) dt[, (cc) := as.numeric(get(cc))]
  dt
}

#' Load population from an ONS-format Excel file
#' @note ONS-format specific.
#' @keywords internal
load_population_ons <- function(path) {
  skip <- find_header_skip(path, sheet = "2")
  dt <- suppressWarnings(
    data.table::as.data.table(readxl::read_excel(path, sheet = "2", skip = skip)))
  data.table::setnames(dt, c("year", "sex", "age_group", "imd_decile", "population"))
  dt[, year := as.integer(year)]
  dt[, imd_decile := as.integer(imd_decile)]
  dt[, population := as.numeric(population)]
  dt[, sex := data.table::fifelse(tolower(sex) == "males", "Male", "Female")]
  dt
}

#' Prepare an ONS analysis data.table
#'
#' Aggregates IMD deciles to quintiles, computes cause-specific and
#' all-cause death counts, merges with population, forward-fills missing
#' population years, and computes mortality rates.
#'
#' @note ONS-format specific.
#' @keywords internal
prepare_ons_analysis_dt <- function(deaths, pop, config) {
  deaths <- deaths[age_group %in% config$age_groups &
                     year %in% config$year_range]
  pop <- pop[age_group %in% config$age_groups &
               year %in% config$year_range]

  deaths[, quintile := ceiling(imd_decile / 2)]
  pop[,    quintile := ceiling(imd_decile / 2)]

  cause_cols <- config$cause_cols
  all_disease_cols <- default_disease_cols()

  if (length(cause_cols) == 1L && cause_cols == "RESIDUAL") {
    disease_cols <- intersect(all_disease_cols, names(deaths))
    deaths[, deaths_disease_sum := Reduce(`+`, lapply(.SD, function(x)
      data.table::fifelse(is.na(x), 0, x))), .SDcols = disease_cols]
    deaths[, deaths_cause := pmax(
      data.table::fifelse(is.na(all_causes), 0, all_causes) -
        deaths_disease_sum, 0)]
    deaths[, deaths_disease_sum := NULL]
  } else {
    missing_cols <- setdiff(cause_cols, names(deaths))
    if (length(missing_cols) > 0)
      stop("Missing cause columns in data: ",
           paste(missing_cols, collapse = ", "))
    deaths[, deaths_cause := Reduce(`+`, lapply(.SD, function(x)
      data.table::fifelse(is.na(x), 0, x))), .SDcols = cause_cols]
  }
  deaths[, deaths_all := data.table::fifelse(is.na(all_causes), 0, all_causes)]

  d_agg <- deaths[, .(deaths_cause = sum(deaths_cause, na.rm = TRUE),
                       deaths_all   = sum(deaths_all, na.rm = TRUE)),
                   by = .(year, sex, quintile, age_group)]

  p_agg <- pop[, .(population = sum(population, na.rm = TRUE)),
                by = .(year, sex, quintile, age_group)]

  dt <- merge(d_agg, p_agg, by = c("year", "sex", "quintile", "age_group"),
              all.x = TRUE)

  if (anyNA(dt$population)) {
    for (fill_yr in sort(unique(dt[is.na(population), year]))) {
      src_yr <- max(p_agg$year[p_agg$year <= fill_yr])
      src <- p_agg[year == src_yr,
                    .(sex, quintile, age_group, pop_fill = population)]
      dt[year == fill_yr & is.na(population),
         population := src[.SD, on = .(sex, quintile, age_group), pop_fill]]
    }
  }

  dt[, population := pmax(population, 1)]
  cause_short <- config$cause_short
  dt[, (paste0("mx_", cause_short)) := deaths_cause / population]
  dt[, mx_all := deaths_all / population]

  dt[, age_group := factor(age_group, levels = config$age_groups)]
  data.table::setorder(dt, year, sex, quintile, age_group)
  dt
}

#' Load and prepare ONS data for a single cause
#'
#' Internal helper used by the pipeline functions when a `mort_config` with
#' ONS file paths is provided. Reads deaths + population Excel files, then
#' calls [prepare_ons_analysis_dt()].
#'
#' @note ONS-format specific.
#' @param config A [mort_config] object with `deaths_file` and `pop_files`
#'   fields.
#' @return A `data.table` ready for array construction.
#' @keywords internal
load_and_prepare_ons <- function(config) {
  cat("=== Loading ONS data ===\n")

  deaths <- load_deaths_ons(config$deaths_file)
  cat("  Deaths:", nrow(deaths), "rows,",
      min(deaths$year), "-", max(deaths$year), "\n")

  pop_list <- lapply(config$pop_files, load_population_ons)
  pop <- data.table::rbindlist(pop_list)
  pop <- unique(pop, by = c("year", "sex", "age_group", "imd_decile"))
  cat("  Population:", nrow(pop), "rows,",
      min(pop$year), "-", max(pop$year), "\n")

  dt <- prepare_ons_analysis_dt(deaths, pop, config)

  cat("  Cause:", config$cause_label,
      "(", paste(config$cause_cols, collapse = " + "), ")\n")
  cat("  Analysis data:", nrow(dt), "rows\n")
  cat("  Years:", min(dt$year), "-", max(dt$year), "\n")
  dt
}
