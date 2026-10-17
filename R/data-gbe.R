#' Read tables exported from GBE-Bund
#'
#' Parses the CSV export of tables of the Federal Health Monitoring
#' (Gesundheitsberichterstattung des Bundes, <https://www.gbe-bund.de>): a
#' few title lines, two header lines (for example years and causes of
#' death) and data rows whose first columns hold the row variables (for
#' example region, age and sex), each written only where it changes. The
#' file is Latin-1 encoded with `;` as separator.
#'
#' Symbols of official statistics are kept in column `symbol`: `"-"`
#' (nothing, value 0), `"."` (unknown or confidential, value `NA`), `"X"`
#' (not applicable, value `NA`), `"..."` (not yet available, value `NA`).
#'
#' @param file Path of the exported CSV file.
#' @param row_names Names of the row variables, in their order.
#' @param col_names Names of the two column variables (outer, inner).
#' @return A `data.table` with the row and column variables, `value` and
#'   `symbol`.
#' @export
mt_read_gbe <- function(file, row_names = c("region", "age", "sex"), col_names = c("year", "cause")) {
  if (!is.character(file) || length(file) != 1 || !file.exists(file)) {
    mt_abort("`file` must be the path of an existing file.", class = "gbe")
  }
  txt <- readChar(file, file.size(file), useBytes = TRUE)
  txt <- iconv(txt, "latin1", "UTF-8")
  lines <- strsplit(gsub("\r\n?", "\n", txt), "\n", fixed = TRUE)[[1]]
  split <- function(l) {
    v <- strsplit(paste0(l, ";<end>"), ";", fixed = TRUE)[[1]]
    gsub('^"|"$', "", v[-length(v)])
  }
  nr <- length(row_names)
  lead <- paste0("^", strrep(";", nr), "\"")
  hdr <- grep(lead, lines)
  if (length(hdr) < 2) mt_abort("No two header lines found (check `row_names`).", class = "gbe")
  outer <- split(lines[hdr[1]])
  inner <- split(lines[hdr[2]])
  nc <- length(inner)
  length(outer) <- nc
  for (i in seq_len(nc)[-1]) if (is.na(outer[i]) || !nzchar(outer[i])) outer[i] <- outer[i - 1]
  first <- hdr[2] + 1
  stop_at <- grep("^Fu.{1,2}noten|^Quelle|^Die Tabelle wurde", lines, useBytes = TRUE)
  stop_at <- if (length(stop_at)) min(stop_at[stop_at > first]) - 1 else length(lines)
  rows <- lines[first:stop_at]
  rows <- rows[nzchar(gsub(";", "", rows))]
  M <- do.call(rbind, lapply(rows, function(l) {
    v <- split(l)
    length(v) <- nc
    v
  }))
  carry <- function(v) {
    v <- trimws(v)
    for (i in seq_along(v)[-1]) if (is.na(v[i]) || !nzchar(v[i])) v[i] <- v[i - 1]
    v
  }
  keys <- lapply(seq_len(nr), function(j) carry(M[, j]))
  vals <- trimws(M[, (nr + 1):nc, drop = FALSE])
  ncell <- ncol(vals)
  out <- data.table(idx = rep(seq_len(nrow(M)), ncell))
  for (j in seq_len(nr)) out[, (row_names[j]) := rep(keys[[j]], ncell)]
  out[, (col_names[1]) := rep(trimws(outer[(nr + 1):nc]), each = nrow(M))]
  out[, (col_names[2]) := rep(trimws(inner[(nr + 1):nc]), each = nrow(M))]
  raw <- as.vector(vals)
  num <- suppressWarnings(as.numeric(gsub(",", ".", gsub("\\.(?=[0-9]{3}(,|$))", "", raw, perl = TRUE))))
  out[, value := num]
  out[, symbol := ifelse(is.na(num), raw, "")]
  out[symbol == "-", `:=`(value = 0)]
  out[, idx := NULL]
  out <- out[nzchar(out[[col_names[2]]])]   # trailing separators give empty columns
  if (identical(col_names[1], "year") && all(grepl("^[0-9]{4}$", out$year))) out[, year := as.integer(year)]
  out[]
}

# nolint start: line_length_linter.
#' Deaths and populations by Land, Germany
#'
#' Deaths by year (1998-2024), Land, sex, age group and ICD-10 chapter, and
#' populations on 31 December (1997-2024) by Land, sex and single year of
#' age. The vignette "Differences between the German Laender" reads them
#' with [mt_data()] and prepares them with [mt_as_exposure()],
#' [mt_regroup()], [mt_subset()] and [mt_fill_suppressed()].
#'
#' `de_deaths` holds the national total (`region = "Deutschland"`) and the 16
#' Laender, 18 age groups (column `age`: `"Unter 1 Jahr"`, `"1 Jahr bis unter 15 Jahre"`,
#' 5-year groups, `"90 Jahre und aelter"` with umlaut), all causes (`cause =
#' "all"`) and the 18 ICD-10 chapters of the source table (codes such as
#' `"C00-D48"`; the labels are in attribute `"cause_labels"`). Deaths that
#' the source suppresses for confidentiality (`"."`) are `NA`: 13% of the
#' cells, mostly rare chapters in small Laender; 33 all-cause cells
#' (Bremen, Saarland, Hamburg, Sachsen-Anhalt, ages 1-34). Deaths from
#' pregnancy and childbirth in men (`"X"`, not applicable) are 0. The causes
#' add up to all causes and the Laender to Germany wherever all cells are
#' published; 2025 is left out because its chapters did not add up. Injuries
#' and poisonings (S00-T98) are the chapter of external causes.
#'
#' `de_cause_groups` lists cause groups for cause-specific analyses at ages
#' 50 and over: neoplasms, diseases of the circulatory, respiratory and
#' digestive systems, mental and nervous system disorders combined (dementia
#' is coded in both chapters), endocrine diseases and external causes; the
#' remaining chapters form `"rest"` with [mt_group_levels()].
#'
#' @format `de_deaths`: a data frame with columns `year`, `region`, `sex`
#'   (`"male"`, `"female"`), `age`, `cause`, `deaths`.
#'   `de_population`: a data frame with columns `year` (population on 31
#'   December), `region`, `sex`, `age` (`"unter 1 Jahr"`, single years, `"90
#'   Jahre und mehr"`), `population`. `de_cause_groups`: a named list of
#'   chapter codes.
#' @source Deaths: Statistisches Bundesamt (Destatis), Todesursachenstatistik, table "Sterbefaelle, Sterbeziffern (ab 1998)" of the Gesundheitsberichterstattung des Bundes, <https://www.gbe-bund.de>, retrieved 13 October 2026. Population: Statistisches Bundesamt (Destatis), Genesis-Online, table 12411-0013, retrieved 13 October 2026. Datenlizenz Deutschland - Namensnennung - Version 2.0 (<https://www.govdata.de/dl-de/by-2-0>); reshaped, suppressed cells kept as missing (own processing, `data-raw/de_data.R`).
#' @examples
#' head(de_deaths)
#' x <- mt_data(deaths = subset(de_deaths, cause == "all" & year >= 2018),
#'              population = subset(de_population, year >= 2017),
#'              by = c("region", "sex"))
#' x <- mt_as_exposure(x, type = "dec31")
#' x
#' @name de_deaths
"de_deaths"

#' @rdname de_deaths
"de_population"

#' @rdname de_deaths
"de_cause_groups"
# nolint end
