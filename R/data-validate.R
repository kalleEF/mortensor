#' Check mortality data for common problems
#'
#' Runs plausibility checks and returns one row per finding:
#'
#' * deaths with zero exposure in the same cell (error);
#' * central death rates above 1 below age 100 (warning);
#' * zero exposures (info);
#' * grid cells without deaths or exposures (info);
#' * hierarchy sums, when `hierarchy` is given ([mt_hierarchy()], or the
#'   older `list(var, total)`): every aggregate stratum in the data is
#'   compared with the sum of its bottom strata per cell, for deaths and, if
#'   the hierarchy is additive on deaths, exposures; ratios near 2 are
#'   reported as "total equals twice the sum of parts", a pattern seen in
#'   real data.
#'
#' @param data An [mt_data] object.
#' @param hierarchy Optional [mt_hierarchy()], or a list with `var` (stratum
#'   column) and `total` (value of the total).
#' @param tol Relative tolerance for hierarchy sums.
#' @return A `data.table` of class `mt_validation` with columns `check`,
#'   `severity`, `measure`, `stratum`, `n`, `detail`.
#' @examples
#' mt_validate(mt_example())
#' @export
mt_validate <- function(data, hierarchy = NULL, tol = 0.005) {
  check_class(data, "mt_data", "data")
  found <- new.env(parent = emptyenv())
  found$rows <- list()
  add <- function(check, severity, measure, stratum, n, detail) {
    found$rows[[length(found$rows) + 1L]] <- data.table(check = check, severity = severity, measure = measure,
                                                        stratum = stratum, n = n, detail = detail)
  }
  D <- data$obs$deaths
  E <- data$obs$exposure
  lab <- function(s) data$strata$label[match(s, data$strata$stratum)]
  if (!is.null(E)) {
    z <- E[value == 0, .N, by = "stratum"]
    for (i in seq_len(nrow(z))) {
      add("zero exposure", "info", "exposure", lab(z$stratum[i]), z$N[i], "cells with zero exposure")
    }
  }
  if (!is.null(D) && !is.null(E)) {
    j <- merge(D, E, by = c("stratum", "age_lo", "age_hi", "period_lo", "period_hi"), suffixes = c(".d", ".e"))
    z <- j[value.d > 0 & value.e == 0, .N, by = "stratum"]
    for (i in seq_len(nrow(z))) add("deaths without exposure", "error", "deaths", lab(z$stratum[i]), z$N[i],
                                    "deaths > 0 where exposure is 0")
    z <- j[value.e > 0 & value.d / value.e > 1 & age_lo < 100,
           list(N = .N, ex = paste(head(unique(age_lo), 3), collapse = ",")), by = "stratum"]
    for (i in seq_len(nrow(z))) add("rate above 1", "warning", "deaths", lab(z$stratum[i]), z$N[i],
                                    sprintf("central death rate > 1 below age 100 (ages %s)", z$ex[i]))
  }
  g <- tryCatch(mt_grid(data), error = function(e) NULL)
  if (!is.null(g)) {
    for (m in intersect(c("deaths", "exposure"), names(data$obs))) {
      idx <- obs_grid_index(data$obs[[m]], g, m)
      ok <- !is.na(idx$a1)
      area <- sum((idx$a2[ok] - idx$a1[ok] + 1) * (idx$p2[ok] - idx$p1[ok] + 1))
      total <- length(g$age) * length(g$period) * nrow(data$strata)
      if (area < total) add("grid coverage", "info", m, "all", total - area, "grid cells not covered by observations")
    }
  }
  if (!is.null(hierarchy)) found$rows <- c(found$rows, list(check_hierarchy(data, hierarchy, tol)))
  res <- rbindlist(found$rows)
  if (!nrow(res)) {
    res <- data.table(check = character(), severity = character(), measure = character(),
                      stratum = character(), n = integer(), detail = character())
  }
  setattr(res, "class", c("mt_validation", class(res)))
  res
}

check_hierarchy <- function(data, hierarchy, tol) {
  h <- as_hierarchy(hierarchy)
  st <- data$strata
  hs <- tryCatch(hierarchy_structure(st, h), mt_error_hierarchy = function(e) {
    mt_abort(conditionMessage(e), class = "hierarchy")
  })
  agg <- which(!hs$is_bottom)
  if (!length(agg)) {
    return(data.table(check = "hierarchy", severity = "info", measure = "deaths", stratum = "all", n = 0L,
                      detail = "no aggregate strata in the data (see mt_aggregate_strata())"))
  }
  measures <- intersect(if (h$additive == "deaths") c("deaths", "exposure") else "deaths", names(data$obs))
  bottom_ids <- st$stratum[hs$bottom]
  rows <- list()
  for (m in measures) {
    o <- data$obs[[m]]
    keys <- c("age_lo", "age_hi", "period_lo", "period_hi")
    for (i in agg) {
      mem <- bottom_ids[hs$S[i, ] > 0]
      parts <- o[o$stratum %in% mem, list(parts = sum(value), k = .N), by = keys]
      parts <- parts[k == length(mem)]
      tot <- o[o$stratum == st$stratum[i], c(keys, "value"), with = FALSE]
      j <- merge(tot, parts, by = keys)
      j <- j[parts > 0]
      if (!nrow(j)) next
      ratio <- j$value / j$parts
      bad <- abs(ratio - 1) > tol
      med <- stats::median(ratio)
      detail <- sprintf("total / sum(parts): median %.3f, range %.3f-%.3f; %d of %d cells outside tolerance",
                        med, min(ratio), max(ratio), sum(bad), length(ratio))
      if (abs(med - 2) < 0.02) detail <- paste0(detail, "; total equals twice the sum of parts")
      rows[[length(rows) + 1L]] <- data.table(check = "hierarchy", severity = if (any(bad)) "error" else "info",
                                              measure = m, stratum = st$label[i], n = sum(bad), detail = detail)
    }
  }
  if (!length(rows)) {
    return(data.table(check = "hierarchy", severity = "info", measure = "deaths", stratum = "all", n = 0L,
                      detail = "no cells with the total and all parts observed"))
  }
  out <- rbindlist(rows)
  ok <- out[severity == "info"]
  bad <- out[severity != "info"]
  if (nrow(ok) > 3) {
    ok <- data.table(check = "hierarchy", severity = "info", measure = paste(unique(ok$measure), collapse = ","),
                     stratum = sprintf("%d aggregate strata", length(unique(ok$stratum))), n = 0L,
                     detail = "sums of parts agree with the totals within tolerance")
  }
  rbind(bad, ok)
}

#' @export
print.mt_validation <- function(x, ...) {
  if (!nrow(x)) {
    cat("<mt_validation> no findings\n")
    return(invisible(x))
  }
  cat(sprintf("<mt_validation> %d finding(s): %d error, %d warning, %d info\n", nrow(x),
              sum(x$severity == "error"), sum(x$severity == "warning"), sum(x$severity == "info")))
  for (i in seq_len(nrow(x))) {
    cat(sprintf("  [%s] %s (%s, %s, n = %s): %s\n", x$severity[i], x$check[i], x$measure[i],
                x$stratum[i], x$n[i], x$detail[i]))
  }
  invisible(x)
}
