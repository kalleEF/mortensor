#' Check mortality data for common problems
#'
#' Runs plausibility checks and returns one row per finding:
#'
#' * deaths with zero exposure in the same cell (error);
#' * central death rates above 1 below age 100 (warning);
#' * zero exposures (info);
#' * grid cells without deaths or exposures (info);
#' * hierarchy sums, when `hierarchy` is given: the total stratum is compared
#'   with the sum of its parts per cell; ratios near 2 are reported as
#'   "total equals twice the sum of parts", a pattern seen in real data.
#'
#' @param data An [mt_data] object.
#' @param hierarchy Optional list with `var` (stratum column), `total` (value
#'   of the total) and `parts` (values of the parts; default all other
#'   values).
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
  var <- hierarchy$var
  if (is.null(var) || !var %in% data$by) mt_abort("`hierarchy$var` must be a stratum column.", class = "hierarchy")
  st <- data$strata
  total <- hierarchy$total
  part_levels <- hierarchy$parts %||% setdiff(unique(st[[var]]), total)
  n_parts <- length(part_levels)
  others <- setdiff(data$by, var)
  D <- merge(data$obs$deaths, st, by = "stratum")
  keys <- c(others, "age_lo", "age_hi", "period_lo", "period_hi")
  tot <- D[get(var) == total, c(keys, "value"), with = FALSE]
  prt <- D[get(var) %in% part_levels, list(parts = sum(value), k = .N), by = keys]
  j <- merge(tot, prt, by = keys)
  j <- j[k == n_parts & parts > 0]
  if (!nrow(j)) {
    return(data.table(check = "hierarchy", severity = "info", measure = "deaths", stratum = total, n = 0L,
                      detail = "no cells with the total and all parts observed"))
  }
  j[, ratio := value / parts]
  bad <- abs(j$ratio - 1) > tol
  med <- stats::median(j$ratio)
  detail <- sprintf("total / sum(parts): median %.3f, range %.3f-%.3f; %d of %d cells outside tolerance",
                    med, min(j$ratio), max(j$ratio), sum(bad), nrow(j))
  if (abs(med - 2) < 0.02) detail <- paste0(detail, "; total equals twice the sum of parts")
  data.table(check = "hierarchy", severity = if (any(bad)) "error" else "info", measure = "deaths",
             stratum = total, n = sum(bad), detail = detail)
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
