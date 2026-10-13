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

#' Read the German regional data
#'
#' Builds an [mt_data] object from deaths by Land, sex, age group and cause
#' (in the layout of [de_deaths]) and populations on 31 December (in the
#' layout of [de_population]). Exposures are the mean of the populations on
#' 31 December of the previous and the current year, summed to the age
#' groups of the deaths.
#'
#' Suppressed death counts (`NA`) are estimated by [mt_fill_suppressed()]
#' from the published totals: deaths of each Land over all causes and deaths
#' of Germany by cause. Cause chapters are then summed into the groups given
#' in `causes`; the rest of all causes forms group `"rest"`.
#'
#' @param deaths A data frame in the layout of [de_deaths].
#' @param population A data frame in the layout of [de_population].
#' @param causes `NULL` (all causes only) or a named list of ICD-10 chapter
#'   codes per cause group (for example [de_cause_groups]).
#' @param years Optional subset of years.
#' @param min_age Optional lower age bound.
#' @param regions Regions to keep (default all Laender; `"Deutschland"` is
#'   the national total and is used for the estimation of suppressed cells).
#' @param keep_total Keep the national total as strata (region
#'   `"Deutschland"`), for example to check the sums with [mt_validate()].
#' @return An [mt_data] object with strata `region`, `sex` and, with
#'   `causes`, `cause`; attribute `"fill"` summarises the estimation of
#'   suppressed cells.
#' @examples
#' x <- mt_read_destatis(de_deaths, de_population, years = 2015:2019)
#' x
#' @export
mt_read_destatis <- function(deaths, population, causes = NULL, years = NULL, min_age = NULL, regions = NULL,
                             keep_total = FALSE) {
  d <- as.data.table(deaths)
  p <- as.data.table(population)
  need_d <- c("year", "region", "sex", "age_group", "cause", "deaths")
  need_p <- c("year", "region", "sex", "age", "population")
  if (!all(need_d %in% names(d)) || !all(need_p %in% names(p))) {
    mt_abort("`deaths` and `population` must have the columns of de_deaths and de_population.", class = "destatis")
  }
  if (!is.null(years)) d <- d[year %in% years]
  ag <- parse_age_labels(unique(d$age_group))
  d <- merge(d, ag, by.x = "age_group", by.y = "label")
  if (!is.null(min_age)) d <- d[age_lo >= min_age]
  if (!nrow(d)) mt_abort("No deaths left after the selection.", class = "destatis")
  p <- p[year %in% c(min(d$year) - 1, unique(d$year))]
  pa <- parse_age_labels(unique(p$age))
  p <- merge(p, pa, by.x = "age", by.y = "label")
  # population ages summed to the age groups of the deaths
  grp_lo <- sort(unique(d$age_lo))
  gi <- findInterval(p$age_lo, grp_lo)
  p[, age_lo := ifelse(gi > 0, grp_lo[pmax(gi, 1)], NA_real_)]
  p <- p[!is.na(age_lo)]
  p <- p[, list(population = sum(population)), by = c("year", "region", "sex", "age_lo")]
  bounds <- unique(d[, c("age_lo", "age_hi")])
  p <- merge(p, bounds, by = "age_lo")
  tot <- "Deutschland"
  if (!tot %in% p$region) p <- rbind(p, p[, list(region = tot, population = sum(population)),
                                          by = c("year", "sex", "age_lo", "age_hi")], use.names = TRUE)
  # all causes only: the suppressed all-cause cells are estimated from the national total
  if (is.null(causes)) d <- d[cause == "all"]
  cause_levels <- unique(d$cause)
  dd <- d[, c("year", "region", "sex", "cause", "age_lo", "age_hi", "deaths")]
  pp <- p[rep(seq_len(nrow(p)), each = length(cause_levels))][, cause := rep(cause_levels, nrow(p))]
  x <- suppressWarnings(mt_data(deaths = dd, population = pp, by = c("region", "sex", "cause"),
                                age = c("age_lo", "age_hi"), period = "year", value = c(deaths = "deaths",
                                                                                        population = "population")))
  x <- mt_as_exposure(x, type = "dec31")
  x$obs$population <- NULL
  x <- mt_fill_suppressed(x, if (is.null(causes)) mt_hierarchy(region = tot) else
    mt_hierarchy(region = tot, cause = "all"))
  fill <- attr(x, "fill")
  if (!is.null(causes)) {
    groups <- c(causes, list(rest = setdiff(setdiff(cause_levels, "all"), unlist(causes))))
    x <- regroup_strata(x, "cause", groups, keep = "all")
  } else {
    x <- regroup_strata(x, "cause", list(), keep = "all")
  }
  keep_regions <- regions %||% setdiff(unique(x$strata$region), tot)
  if (keep_total) keep_regions <- union(keep_regions, tot)
  x <- strata_subset(x, x$strata$region %in% keep_regions)
  if (is.null(causes)) x <- drop_stratum_var(x, "cause")
  attr(x, "fill") <- fill
  mt_set_grid(x, age = mt_axis(bounds[order(age_lo)]$age_lo, bounds[order(age_lo)]$age_hi, kind = "age"))
}

# Sum levels of one stratum variable into groups (deaths and exposures are
# summed for deaths; for a cause variable the exposure is shared, so it is
# taken from the first member). Levels in `keep` stay as they are.
regroup_strata <- function(x, var, groups, keep = character()) {
  st <- x$strata
  others <- setdiff(x$by, var)
  key <- function(s) {
    if (length(others)) do.call(paste, c(lapply(others, function(o) s[[o]]), sep = "\r")) else rep("", nrow(s))
  }
  obs_out <- lapply(x$obs, function(o) NULL)
  combos <- unique(st[, others, with = FALSE])
  levels_out <- c(keep, names(groups))
  rows <- list()
  k <- 0L
  for (lev in levels_out) {
    mem_levels <- if (lev %in% keep) lev else groups[[lev]]
    for (ci in seq_len(max(nrow(combos), 1))) {
      sel <- st[[var]] %in% mem_levels
      if (length(others)) sel <- sel & key(st) == key(combos[ci])
      members <- st$stratum[sel]
      if (!length(members)) next
      k <- k + 1L
      row <- if (length(others)) copy(combos[ci]) else data.table()
      row[, (var) := lev]
      rows[[k]] <- row
      for (m in names(x$obs)) {
        o <- x$obs[[m]][stratum %in% members]
        if (m == "deaths") {
          a <- o[, list(value = sum(value), n_members = .N), by = c("age_lo", "age_hi", "period_lo", "period_hi")]
          a <- a[a$n_members == length(members)]
          a[, n_members := NULL]
        } else {
          a <- o[stratum == members[1], c("age_lo", "age_hi", "period_lo", "period_hi", "value"), with = FALSE]
        }
        a[, `:=`(stratum = k, provenance = if (length(members) > 1) "aggregated" else "observed")]
        obs_out[[m]] <- rbind(obs_out[[m]], a)
      }
    }
  }
  new <- rbindlist(rows, use.names = TRUE)
  setcolorder(new, x$by)
  new[, label := do.call(paste, c(lapply(x$by, function(b) new[[b]]), sep = "."))]
  new[, stratum := seq_len(.N)]
  setcolorder(new, c("stratum", x$by, "label"))
  x$strata <- new
  for (m in names(obs_out)) {
    o <- obs_out[[m]]
    setcolorder(o, c("stratum", "age_lo", "age_hi", "period_lo", "period_hi", "value", "provenance"))
    setorderv(o, c("stratum", "period_lo", "age_lo"))
    x$obs[[m]] <- o
  }
  x
}

strata_subset <- function(x, keep) {
  st <- x$strata[keep]
  map <- stats::setNames(seq_len(nrow(st)), st$stratum)
  for (m in names(x$obs)) {
    o <- x$obs[[m]][stratum %in% st$stratum]
    o[, stratum := as.integer(map[as.character(stratum)])]
    x$obs[[m]] <- o
  }
  st[, stratum := seq_len(.N)]
  x$strata <- st
  x
}

drop_stratum_var <- function(x, var) {
  x$by <- setdiff(x$by, var)
  x$strata[, (var) := NULL]
  x$strata[, label := do.call(paste, c(lapply(x$by, function(b) x$strata[[b]]), sep = "."))]
  x
}

# nolint start: line_length_linter.
#' Deaths and populations by Land, Germany
#'
#' Deaths by year (1998-2024), Land, sex, age group and ICD-10 chapter, and
#' populations on 31 December (1997-2024) by Land, sex and single year of
#' age, for the regional examples. Read them with [mt_read_destatis()].
#'
#' `de_deaths` holds the national total (`region = "Deutschland"`) and the 16
#' Laender, 18 age groups (`"Unter 1 Jahr"`, `"1 Jahr bis unter 15 Jahre"`,
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
#' `de_cause_groups` lists the cause groups of the regional vignette for ages
#' 50 and over: neoplasms, diseases of the circulatory, respiratory and
#' digestive systems, mental and nervous system disorders combined (dementia
#' is coded in both chapters), endocrine diseases and external causes; the
#' remaining chapters form `"rest"` in [mt_read_destatis()].
#'
#' @format `de_deaths`: a data frame with columns `year`, `region`, `sex`
#'   (`"male"`, `"female"`), `age_group`, `cause`, `deaths`.
#'   `de_population`: a data frame with columns `year` (population on 31
#'   December), `region`, `sex`, `age` (`"unter 1 Jahr"`, single years, `"90
#'   Jahre und mehr"`), `population`. `de_cause_groups`: a named list of
#'   chapter codes.
#' @source Deaths: Statistisches Bundesamt (Destatis), Todesursachenstatistik, table "Sterbefaelle, Sterbeziffern (ab 1998)" of the Gesundheitsberichterstattung des Bundes, <https://www.gbe-bund.de>, retrieved 13 October 2026. Population: Statistisches Bundesamt (Destatis), Genesis-Online, table 12411-0013, retrieved 13 October 2026. Datenlizenz Deutschland - Namensnennung - Version 2.0 (<https://www.govdata.de/dl-de/by-2-0>); reshaped, suppressed cells kept as missing (own processing, `data-raw/de_data.R`).
#' @examples
#' head(de_deaths)
#' x <- mt_read_destatis(de_deaths, de_population, years = 2017:2019)
#' x
#' @name de_deaths
"de_deaths"

#' @rdname de_deaths
"de_population"

#' @rdname de_deaths
"de_cause_groups"
# nolint end
