# Builds data/ons_data.rda and data/ons_cause_groups.rda from the ONS ad hoc
# release 3064 "Death registrations and populations by Index of Multiple
# Deprivation, England" (released 13 October 2025), first workbook.
#
# Source: Office for National Statistics licensed under the Open Government
# Licence v3.0 (Crown copyright 2025).
# https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/adhocs/3064deathregistrationsandpopulationsbyindexofmultipledeprivationengland
#
# Usage: Rscript data-raw/ons_data.R [path to deathsandpopsimdfinal.xlsx]
# Without a path the workbook is downloaded to tempdir().

library(data.table)
args <- commandArgs(trailingOnly = TRUE)
url <- paste0("https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/",
              "adhocs/3064deathregistrationsandpopulationsbyindexofmultipledeprivationengland/",
              "deathsandpopsimdfinal.xlsx")
path <- if (length(args)) args[1] else {
  f <- file.path(tempdir(), "deathsandpopsimdfinal.xlsx")
  utils::download.file(url, f, mode = "wb", quiet = TRUE)
  f
}
stopifnot(tools::md5sum(path)[[1]] == "2dee9698887b508b9a539ad3f92ee470")

clean <- function(x) gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", tolower(trimws(x))))
read_sheet <- function(sheet) {
  x <- as.data.table(suppressMessages(readxl::read_excel(path, sheet = sheet, skip = 5)))
  setnames(x, clean(names(x)))
  x
}
deaths <- read_sheet("1")
pop <- read_sheet("2")
setnames(pop, "all_causes", "population")   # the value column of sheet 2 is headed "All causes"
years <- 2001:2019                          # before the COVID-19 pandemic
deaths <- deaths[year %in% years]
pop <- pop[year %in% years]
stopifnot(all(deaths$deaths_involving_covid_19 == 0))
deaths[, deaths_involving_covid_19 := NULL]

imd <- data.table(year = years, imd_version = cut(years, c(2000, 2005, 2008, 2012, 2017, 2024),
                                                 labels = c(2004, 2007, 2010, 2015, 2019)))
imd[, imd_version := as.integer(as.character(imd_version))]

keys <- c("year", "sex", "age_group", "imd_decile")
ons_data <- merge(pop, deaths, by = keys, sort = FALSE)
ons_data <- merge(imd, ons_data, by = "year")
ages <- c("<1", "01-04", sprintf("%02d-%02d", seq(5, 85, 5), seq(9, 89, 5)), "90+")
ons_data[, age_group := factor(age_group, levels = ages)]
setorderv(ons_data, c("year", "sex", "age_group", "imd_decile"))
ons_data[, age_group := as.character(age_group)]
num <- setdiff(names(ons_data), c("sex", "age_group"))
ons_data[, (num) := lapply(.SD, as.integer), .SDcols = num]
setcolorder(ons_data, c(keys, "imd_version", "population", "all_causes"))
stopifnot(nrow(ons_data) == length(years) * 2 * 20 * 10, !anyNA(ons_data))
ons_data <- as.data.frame(ons_data)

# Cause groups used in the vignettes (deaths at ages 50 and over): the causes
# of the table are mutually exclusive (ICD-10 lookup of the workbook). Common
# cancers and causes stay separate; small causes are combined within their
# ICD-10 chapter; everything else, including suicide, is in "rest".
ons_cause_groups <- list(
  lung_cancer = "lung_cancer",
  colorectal_cancer = "colorectal_cancer",
  breast_prostate = c("breast_cancer", "prostate_cancer"),
  other_cancers = "all_other_cancers",
  ihd = "chronic_ischaemic_heart_disease",
  stroke = "stroke_incl_transient_ischaemic_attack",
  other_circulatory = c("heart_failure", "hypertension", "atrial_fibrillation"),
  dementia = "dementia",
  respiratory = c("copd", "asthma"),
  diabetes = c("type_1_diabetes", "type_2_diabetes", "diabetes_other"),
  liver_alcohol = c("alcohol_related_liver_disease", "non_alcohol_related_liver_disease", "alcohol_specific"),
  rest = "rest")

save(ons_data, file = "data/ons_data.rda", compress = "xz")
save(ons_cause_groups, file = "data/ons_cause_groups.rda", compress = "xz")
