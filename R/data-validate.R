#' Check mortality data for common problems
#'
#' Runs plausibility checks and returns one row per finding:
#'
#' * deaths with zero exposure in the same cell (error);
#' * central death rates above 1 below age 100 (warning);
#' * zero exposures (info);
#' * grid cells without deaths or exposures (info);
#' * hierarchy sums, when `hierarchy` is given ([mt_hierarchy()]): every
#'   aggregate stratum in the data is compared with the sum of its bottom
#'   strata per cell, for deaths and exposures (exposures are summed over
#'   the variables whose exposures add up and taken once over those whose
#'   levels share an exposure, such as causes of death).
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
  for (m in names(data$obs)) {
    z <- data$obs[[m]][is.na(value), .N, by = "stratum"]
    for (i in seq_len(nrow(z))) {
      add("missing values", "warning", m, lab(z$stratum[i]), z$N[i],
          "cells without a value (suppressed or unknown); see mt_fill_suppressed()")
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
  # strata without any observation of a measure (for example populations
  # that lack the national total)
  for (m in intersect(c("deaths", "exposure"), names(data$obs))) {
    none <- setdiff(data$strata$stratum, unique(data$obs[[m]]$stratum))
    if (length(none)) {
      add(sprintf("no %s", m), "warning", m, paste(utils::head(lab(none), 4), collapse = ", "), length(none),
          sprintf("%d stratum/strata without any %s; remove them with mt_subset() or add the table", length(none),
                  m))
    }
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
  # levels that are the sum of the other levels of their variable (for
  # example "Persons" next to "female" and "male") but are not declared as
  # totals of the hierarchy
  declared <- if (!is.null(hierarchy)) {
    h <- as_hierarchy(hierarchy)
    unlist(lapply(names(h$vars), function(v) paste(v, c(h$vars[[v]]$total, names(h$vars[[v]]$groups)))))
  } else {
    character()
  }
  tl <- find_total_levels(data)
  tl <- tl[!paste(var, level) %in% declared]
  for (i in seq_len(nrow(tl))) {
    add("total among levels", "warning", "deaths", sprintf("%s = %s", tl$var[i], tl$level[i]), tl$n[i],
        sprintf(paste("equals the sum of the other levels of `%s` in %.0f%% of cells: declare it as the total in",
                      "mt_hierarchy(%s = \"%s\") or remove it with mt_subset(), so that it is not counted twice"),
                tl$var[i], 100 * tl$share[i], tl$var[i], tl$level[i]))
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
  measures <- intersect(c("deaths", "exposure"), names(data$obs))
  bottom_ids <- st$stratum[hs$bottom]
  rows <- list()
  for (m in measures) {
    o <- data$obs[[m]]
    keys <- c("age_lo", "age_hi", "period_lo", "period_hi")
    Sm <- if (m == "deaths") hs$S else hs$Sx
    for (i in agg) {
      mem <- bottom_ids[Sm[i, ] > 0]
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

# Levels of a stratum variable whose deaths equal the sum of the other levels
# (within 0.5%) in at least 90% of the cells where all levels are observed
find_total_levels <- function(data, vars = data$by) {
  o <- data$obs$deaths
  empty <- data.table(var = character(), level = character(), share = numeric(), n = integer())
  if (is.null(o) || !length(vars)) return(empty)
  st <- data$strata
  o <- merge(o[, c("stratum", "age_lo", "period_lo", "value"), with = FALSE], st, by = "stratum")
  out <- list()
  for (v in vars) {
    levs <- unique(st[[v]])
    if (length(levs) < 3) next
    keys <- c(setdiff(data$by, v), "age_lo", "period_lo")
    w <- o[, c(keys, v, "value"), with = FALSE]
    tot <- w[, list(sum_all = sum(value), n_lev = .N), by = keys]
    for (l in levs) {
      x <- merge(w[w[[v]] == l, c(keys, "value"), with = FALSE], tot, by = keys)
      x <- x[n_lev == length(levs) & value > 0]
      if (nrow(x) < 5) next
      rest <- x$sum_all - x$value
      share <- mean(abs(x$value - rest) <= pmax(0.005 * x$value, 0.5))
      if (share >= 0.9) out[[length(out) + 1L]] <- data.table(var = v, level = l, share = share, n = nrow(x))
    }
  }
  if (length(out)) rbindlist(out) else empty
}
