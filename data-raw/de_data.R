# Builds data/de_deaths.rda, data/de_population.rda and data/de_cause_groups.rda.
#
# Deaths: Gesundheitsberichterstattung des Bundes (GBE-Bund), table "Sterbefaelle,
# Sterbeziffern (je 100.000 Einwohner, altersstandardisiert) (ab 1998)", deaths by
# year, Land, age group, sex and ICD-10 chapter, all nationalities, exported as CSV
# from https://www.gbe-bund.de on 13 October 2026. Source: Statistisches Bundesamt,
# Todesursachenstatistik.
# Population: Statistisches Bundesamt (Destatis), GENESIS-Online table 12411-0013,
# population on 31 December by Land, sex and single year of age (flat CSV).
# Licence: Datenlizenz Deutschland - Namensnennung - Version 2.0 (dl-de/by-2-0);
# the data are reshaped and suppressed cells kept as missing.
#
# Usage: Rscript data-raw/de_data.R <deaths CSV from GBE> <population CSV from GENESIS>

suppressMessages(devtools::load_all(quiet = TRUE))
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)

raw <- mt_read_gbe(args[1], row_names = c("region", "age_group", "sex"), col_names = c("year", "cause_label"))
raw[, cause := ifelse(startsWith(cause_label, "Alle"), "all", sub(" .*", "", cause_label))]
raw[, sex := ifelse(grepl("^M", sex), "male", "female")]
# "X" (not applicable) occurs only for pregnancy and childbirth (O00-O99) in men: zero by definition
stopifnot(all(raw[symbol == "X", unique(cause)] == "O00-O99"), all(raw[symbol == "X", unique(sex)] == "male"))
raw[symbol == "X", value := 0]
# 2025: cause chapters do not add up to all causes (provisional); kept out
de_deaths <- raw[year <= 2024, list(year, region, sex, age = age_group, cause, deaths = as.integer(round(value)))]
ages <- unique(raw$age_group)
de_deaths[, age := factor(age, levels = ages)]
setorderv(de_deaths, c("year", "region", "sex", "age", "cause"))
de_deaths[, age := as.character(age)]
de_cause_labels <- unique(raw[, list(cause, cause_label)])

pop <- mortensor:::genesis_long(mortensor:::genesis_table(args[2]))
pop[, year := as.integer(substr(as.character(time), 1, 4))]
de_population <- pop[year >= 1997 & year <= 2024, list(year, region, sex, age, population = as.integer(value))]
setorderv(de_population, c("year", "region", "sex"))
de_population[, age_lo := parse_age_labels(age)$age_lo]
setorderv(de_population, c("year", "region", "sex", "age_lo"))
de_population[, age_lo := NULL]

# consistency checks
stopifnot(nrow(de_population) == 28 * 16 * 2 * 91)

de_deaths <- as.data.frame(de_deaths)
de_population <- as.data.frame(de_population)
de_cause_groups <- list(
  neoplasms = "C00-D48",
  circulatory = "I00-I99",
  respiratory = "J00-J99",
  digestive = "K00-K93",
  mental_nervous = c("F00-F99", "G00-H95"),
  endocrine = "E00-E90",
  external = "S00-T98")
attr(de_deaths, "cause_labels") <- as.data.frame(de_cause_labels)
save(de_deaths, file = "data/de_deaths.rda", compress = "xz")
save(de_population, file = "data/de_population.rda", compress = "xz")
save(de_cause_groups, file = "data/de_cause_groups.rda", compress = "xz")
