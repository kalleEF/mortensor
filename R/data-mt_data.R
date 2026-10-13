#' Mortality data observed at any granularity
#'
#' `mt_data()` stores death counts, exposures (person-years at risk) and
#' populations for one or more strata. Each table may use its own age and
#' period resolution, for example deaths in 5-year age groups with an open
#' 90+ interval and exposures in single years of age. The container keeps
#' the observations as they were supplied; harmonisation onto a common grid
#' happens explicitly with [mt_regroup()] and [mt_ungroup()].
#'
#' @param deaths,exposure,population Data frames in long format, one row per
#'   cell, or `NULL`.
#' @param by Character vector of stratum columns (for example `"sex"`).
#' @param age Age column: one column with single ages (numeric) or labels such
#'   as `"30-34"` and `"90+"` (see [parse_age_labels()]), or two columns with
#'   the lower and upper bound (`Inf` for an open interval).
#' @param period Period column: calendar years (numeric), labels such as
#'   `"2010-2014"`, or two columns with lower and upper bounds.
#' @param value Name of the value column. Either one name used for all tables,
#'   or a named vector or list such as `c(deaths = "D", exposure = "E")`. By
#'   default the column named after the measure, otherwise `"value"`. A table
#'   in wide format, with one value column per level of a stratum variable
#'   (for example columns `Female` and `Male`, or one column per cause of
#'   death), is read by giving all its value columns, for example
#'   `value = list(deaths = c(female = "Female", male = "Male"))`: the
#'   columns become the levels (named by the names of the vector, else the
#'   column names) of the one variable in `by` that the table lacks.
#' @param age_open Lower bound of an open age interval coded as a single
#'   numeric age (for example 110 for the Human Mortality Database "110+").
#' @param omega Assumed maximum age, used wherever an open interval needs a
#'   finite upper bound.
#' @param missing Values that mark missing cells in character value columns,
#'   for example the symbols of official statistics for confidential or
#'   unknown counts (`"."`, `"x"`). Missing cells are left out; suppressed
#'   death counts can then be estimated with [mt_fill_suppressed()].
#' @param zero Values that mark a count of zero (for example `"-"`).
#' @details A table may leave out stratum columns over which its values are
#'   shared. For example, deaths by region, sex and cause of death and
#'   populations by region and sex: every cause of death then uses the
#'   population of its region and sex. Tables with all stratum columns are
#'   read first and define the strata.
#' @return An object of class `mt_data`.
#' @seealso [mt_add_deaths()], [mt_grid()], [mt_regroup()], [mt_ungroup()],
#'   [as_mt_data()].
#' @examples
#' d <- data.frame(sex = "female", year = 2001:2002, age = "90+",
#'                 deaths = c(120, 130), exposure = c(800, 820))
#' mt_data(deaths = d, exposure = d, by = "sex")
#' @export
mt_data <- function(deaths = NULL, exposure = NULL, population = NULL, by = NULL,
                    age = "age", period = "year", value = NULL, age_open = NULL,
                    omega = 110, missing = c(".", "..", "...", ":", "x", "X", "/", "NA", ""), zero = "-") {
  check_number(omega, "omega", lower = 1)
  out <- new_mt_data(omega = omega)
  inputs <- list(deaths = deaths, exposure = exposure, population = population)
  inputs <- inputs[!vapply(inputs, is.null, logical(1))]
  values <- list()
  for (m in names(inputs)) {
    w <- melt_wide(inputs[[m]], pick_value_col(value, m, inputs[[m]]), by, m)
    inputs[[m]] <- w$x
    values[[m]] <- w$value
  }
  # tables with every stratum column first: they define the strata
  complete <- vapply(inputs, function(x) all(by %in% names(x)), logical(1))
  for (m in names(inputs)[order(!complete)]) {
    out <- mt_add(out, m, inputs[[m]], by = by, age = age, period = period, value = values[[m]],
                  age_open = age_open, missing = missing, zero = zero)
  }
  out
}

# A wide table (one value column per level of a stratum variable) in long format
melt_wide <- function(x, value, by, measure) {
  if (length(value) <= 1L) return(list(x = x, value = value))
  if (!is.data.frame(x)) mt_abort(sprintf("The %s table must be a data frame.", measure), class = "obs")
  x <- as.data.table(x)
  miss <- setdiff(value, names(x))
  if (length(miss)) {
    mt_abort(sprintf("Column(s) not found in the %s table: %s.", measure, paste(miss, collapse = ", ")),
             class = "obs")
  }
  var <- setdiff(by, names(x))
  if (length(var) != 1L) {
    mt_abort(sprintf(paste("Several value columns of the %s table are read as the levels of a stratum variable,",
                           "which must be the one variable in `by` that the table lacks."), measure), class = "obs")
  }
  labs <- names(value) %||% value
  labs[!nzchar(labs)] <- value[!nzchar(labs)]
  if (!all(vapply(x[, value, with = FALSE], is.numeric, logical(1)))) {
    x[, (value) := lapply(.SD, as.character), .SDcols = value]
  }
  long <- melt(x, id.vars = setdiff(names(x), value), measure.vars = unname(value), variable.name = var,
               value.name = ".value", variable.factor = FALSE)
  set(long, j = var, value = labs[match(long[[var]], value)])
  list(x = long, value = ".value")
}

new_mt_data <- function(omega = 110) {
  structure(list(obs = list(), strata = data.table(stratum = integer(), label = character()),
                 by = character(), grid = NULL, omega = omega),
            class = "mt_data")
}

pick_value_col <- function(value, measure, x) {
  if (is.null(value)) return(if (measure %in% names(x)) measure else "value")
  if (!is.null(names(value)) && measure %in% names(value)) return(value[[measure]])
  if (is.null(names(value)) && length(value) == 1L) return(value)
  if (measure %in% names(x)) measure else "value"
}

mt_measures <- c("deaths", "exposure", "population")

#' Add observations to an `mt_data` object
#'
#' @param data An [mt_data] object.
#' @param x A data frame in long format.
#' @inheritParams mt_data
#' @param value Name of the value column.
#' @param measure One of `"deaths"`, `"exposure"`, `"population"`.
#' @return The updated `mt_data` object. An existing table for the same
#'   measure is replaced.
#' @examples
#' d <- data.frame(year = 2001, age = 30:31, deaths = c(5, 7), exposure = c(9e3, 9e3))
#' x <- mt_add_deaths(mt_data(), d)
#' x <- mt_add_exposure(x, d)
#' x
#' @export
mt_add <- function(data, measure, x, by = NULL, age = "age", period = "year",
                   value = NULL, age_open = NULL, missing = c(".", "..", "...", ":", "x", "X", "/", "NA", ""),
                   zero = "-") {
  check_class(data, "mt_data", "data")
  check_choice(measure, mt_measures, "measure")
  if (!is.data.frame(x)) mt_abort("`x` must be a data frame.", class = "obs")
  by <- by %||% data$by
  if (length(data$obs) > 0L && !setequal(by, data$by)) {
    mt_abort(sprintf("All tables need the same stratum columns (%s).",
                     paste(data$by, collapse = ", ")), class = "obs")
  }
  w <- melt_wide(x, value %||% pick_value_col(NULL, measure, x), by, measure)
  x <- as.data.table(w$x)
  value <- w$value
  absent <- setdiff(by, names(x))
  if (length(absent)) x <- share_over_strata(x, data$strata, by, absent, measure)
  value <- value %||% pick_value_col(NULL, measure, x)
  obs <- ingest_obs(x, by, age, period, value, measure, age_open, missing = missing, zero = zero)
  st <- register_strata(data$strata, obs, by)
  data$strata <- st
  data$by <- by
  if (length(by)) {
    obs <- merge(obs, st[, c("stratum", by), with = FALSE], by = by, sort = FALSE)
    obs[, (by) := NULL]
  } else {
    obs[, stratum := 1L]
  }
  setcolorder(obs, c("stratum", "age_lo", "age_hi", "period_lo", "period_hi", "value"))
  setorderv(obs, c("stratum", "period_lo", "age_lo"))
  validate_obs(obs, measure)
  obs[, provenance := "observed"]
  data$obs[[measure]] <- obs
  data$grid <- NULL
  data
}

# A table without some stratum columns is shared over their levels (for
# example populations shared by all causes of death)
share_over_strata <- function(x, strata, by, absent, measure) {
  if (!nrow(strata) || !all(absent %in% names(strata))) {
    mt_abort(sprintf(paste("The %s table has no column(s) %s. Tables without some stratum columns are shared over",
                           "their levels, which must come from another table (for example the deaths)."),
                     measure, paste(absent, collapse = ", ")), class = "obs")
  }
  present <- setdiff(by, absent)
  combos <- unique(strata[, by, with = FALSE])
  if (!length(present)) {
    idx <- CJ(i = seq_len(nrow(x)), j = seq_len(nrow(combos)))
    return(cbind(x[idx$i], combos[idx$j]))
  }
  for (b in present) set(x, j = b, value = as.character(x[[b]]))
  merge(x, combos, by = present, allow.cartesian = TRUE)
}

#' @rdname mt_add
#' @param ... Arguments passed to `mt_add()`.
#' @export
mt_add_deaths <- function(data, x, ...) mt_add(data, "deaths", x, ...)

#' @rdname mt_add
#' @export
mt_add_exposure <- function(data, x, ...) mt_add(data, "exposure", x, ...)

#' @rdname mt_add
#' @export
mt_add_population <- function(data, x, ...) mt_add(data, "population", x, ...)

ingest_obs <- function(x, by, age, period, value, measure, age_open = NULL, missing = character(),
                       zero = character()) {
  need <- unique(c(by, age, period, value))
  miss <- setdiff(need, names(x))
  if (length(miss)) {
    mt_abort(sprintf("Column(s) not found in the %s table: %s.", measure,
                     paste(miss, collapse = ", ")), class = "obs")
  }
  a <- interval_cols(x, age, "age", age_open)
  p <- interval_cols(x, period, "period")
  out <- data.table(age_lo = a$lo, age_hi = a$hi, period_lo = p$lo, period_hi = p$hi,
                    value = parse_values(x[[value]], missing, zero, measure))
  if (length(by)) {
    for (b in by) set(out, j = b, value = as.character(x[[b]]))
  }
  nas <- is.na(out$value)
  if (any(nas)) {
    mt_warn(sprintf("%d missing %s value(s) left out%s.", sum(nas), measure,
                    if (measure == "deaths") "; suppressed counts can be estimated with mt_fill_suppressed()" else ""),
            class = "missing")
    out <- out[!nas]
  }
  out
}

# Numeric values from numbers or text with markers of missing and zero cells
parse_values <- function(v, missing, zero, measure) {
  if (is.numeric(v)) return(as.numeric(v))
  v <- trimws(as.character(v))
  num <- suppressWarnings(as.numeric(v))
  num[v %in% zero] <- 0
  unknown <- unique(v[is.na(num) & !is.na(v) & !v %in% c(missing, zero)])
  if (length(unknown)) {
    mt_warn(sprintf(paste("Values of %s that are neither numbers nor declared markers are treated as missing: %s.",
                          "Declare them with `missing` or `zero`."), measure,
                    paste0("\"", utils::head(unknown, 5), "\"", collapse = ", ")), class = "markers")
  }
  num
}

interval_cols <- function(x, cols, kind, open = NULL) {
  if (length(cols) == 2L) {
    return(list(lo = as.numeric(x[[cols[1]]]), hi = as.numeric(x[[cols[2]]])))
  }
  if (length(cols) != 1L) mt_abort(sprintf("`%s` must name one or two columns.", kind), class = "obs")
  v <- x[[cols]]
  if (is.factor(v)) v <- as.character(v)
  if (is.numeric(v)) {
    lo <- as.numeric(v)
    hi <- lo + 1
    if (!is.null(open)) hi[lo >= open] <- Inf
    if (any(abs(lo - round(lo)) > 1e-9)) {
      mt_abort(sprintf("Numeric %s values must be whole numbers (single years).", kind), class = "obs")
    }
    return(list(lo = lo, hi = hi))
  }
  u <- unique(v)
  p <- parse_age_labels(u)
  i <- match(v, u)
  lo <- p$age_lo[i]
  hi <- p$age_hi[i]
  if (kind == "period" && any(is.infinite(hi))) {
    mt_abort("Period intervals must be closed.", class = "obs")
  }
  list(lo = lo, hi = hi)
}

register_strata <- function(strata, obs, by) {
  if (!length(by)) {
    return(data.table(stratum = 1L, label = "all"))
  }
  new <- unique(obs[, by, with = FALSE])
  if (nrow(strata) && all(by %in% names(strata))) {
    old <- strata[, c("stratum", by), with = FALSE]
    new <- new[!old, on = by]
  } else {
    old <- new[0L][, stratum := integer()]
    setcolorder(old, c("stratum", by))
  }
  if (nrow(new)) {
    new[, stratum := nrow(old) + seq_len(.N)]
    old <- rbind(old, new[, c("stratum", by), with = FALSE])
  }
  old[, label := do.call(paste, c(lapply(by, function(b) old[[b]]), sep = "."))]
  setcolorder(old, c("stratum", by, "label"))
  old[]
}

validate_obs <- function(obs, measure) {
  if (any(obs$value < 0, na.rm = TRUE)) mt_abort(sprintf("Negative %s values.", measure), class = "obs")
  if (any(obs$age_hi <= obs$age_lo) || any(obs$period_hi <= obs$period_lo)) {
    mt_abort("Every interval needs upper > lower.", class = "obs")
  }
  dup <- duplicated(obs, by = c("stratum", "age_lo", "period_lo"))
  if (any(dup)) {
    mt_abort(sprintf("Duplicate %s cells (same stratum, age and period), e.g. age %s in %s.",
                     measure, obs$age_lo[dup][1], obs$period_lo[dup][1]), class = "obs")
  }
  ov <- obs[order(stratum, period_lo, age_lo), list(bad = any(head(age_hi, -1) > tail(age_lo, -1) + 1e-9)),
            by = c("stratum", "period_lo")]
  if (any(ov$bad)) mt_abort(sprintf("Overlapping age intervals in the %s table.", measure), class = "obs")
  pp <- unique(obs[, c("stratum", "period_lo", "period_hi")])[order(stratum, period_lo)]
  pov <- pp[, list(bad = any(head(period_hi, -1) > tail(period_lo, -1) + 1e-9) |
                     anyDuplicated(period_lo) > 0), by = "stratum"]
  if (any(pov$bad)) mt_abort(sprintf("Overlapping period intervals in the %s table.", measure), class = "obs")
  if (measure == "deaths" && any(abs(obs$value - round(obs$value)) > 1e-6, na.rm = TRUE)) {
    mt_warn("Death counts are not whole numbers; they are kept as supplied.", class = "noninteger")
  }
  invisible(TRUE)
}

#' Strata, grid and tables of an `mt_data` object
#'
#' `mt_strata()` returns the stratum table, `mt_grid()` the modelling grid
#' (single years of age and calendar years by default, see [mt_set_grid()]),
#' `mt_obs()` one observation table.
#'
#' @param data An [mt_data] object.
#' @param measure One of `"deaths"`, `"exposure"`, `"population"`.
#' @return A `data.table` or, for `mt_grid()`, a list with `age` and `period`
#'   axes.
#' @examples
#' x <- mt_example()
#' mt_strata(x)
#' mt_grid(x)
#' head(mt_obs(x, "deaths"))
#' @export
mt_strata <- function(data) {
  check_class(data, "mt_data", "data")
  copy(data$strata)
}

#' @rdname mt_strata
#' @export
mt_obs <- function(data, measure = "deaths") {
  check_class(data, "mt_data", "data")
  obs <- data$obs[[measure]]
  if (is.null(obs)) mt_abort(sprintf("No %s table.", measure), class = "missing_measure")
  out <- merge(data$strata, obs, by = "stratum", sort = FALSE)
  setorderv(out, c("stratum", "period_lo", "age_lo"))
  out[]
}

#' @rdname mt_strata
#' @export
mt_grid <- function(data) {
  check_class(data, "mt_data", "data")
  data$grid %||% default_grid(data)
}

# The coarsest grid on which every table is exact: all interval bounds of
# all tables (age groups when every table has the same groups, single years
# when one table has single years)
default_grid <- function(data) {
  if (!length(data$obs)) mt_abort("No observations.", class = "obs")
  all <- rbindlist(lapply(data$obs, function(o) o[, c("age_lo", "age_hi", "period_lo", "period_hi")]))
  ab <- sort(unique(c(all$age_lo, all$age_hi)))
  fin <- ab[is.finite(ab)]
  age <- if (any(is.infinite(ab))) {
    mt_axis(fin, c(fin[-1], Inf), kind = "age")
  } else {
    mt_axis(fin[-length(fin)], fin[-1], kind = "age")
  }
  pb <- sort(unique(c(all$period_lo, all$period_hi)))
  list(age = age, period = mt_axis(pb[-length(pb)], pb[-1], kind = "period"))
}

# Single years of age from the lowest age to the start of the open interval
# (or the highest age), the target of ungrouping
single_year_grid <- function(data) {
  all <- rbindlist(lapply(data$obs, function(o) o[, c("age_lo", "age_hi")]))
  open_lo <- all$age_lo[is.infinite(all$age_hi)]
  a0 <- min(all$age_lo)
  if (length(open_lo)) {
    top <- max(c(open_lo, all$age_hi[is.finite(all$age_hi)]))
    if (top > a0) single_ages(a0, top, open = TRUE) else age_groups(a0)
  } else {
    single_ages(a0, max(all$age_hi) - 1, open = FALSE)
  }
}

#' Set the modelling grid
#'
#' The grid is the resolution on which models work. Every observed cell must
#' be a union of grid cells, so the grid is at least as fine as the finest
#' table. Without a grid set, the grid is the coarsest one on which every
#' table is exact: the age groups of the data when all tables have the same
#' groups, single years when one table (for example populations) has single
#' years. Coarser grids are obtained with [mt_regroup()], single years from
#' grouped deaths with [mt_ungroup()].
#'
#' @param data An [mt_data] object.
#' @param age,period [mt_axis] objects (`NULL` keeps the current axis).
#' @return The updated `mt_data` object.
#' @examples
#' x <- mt_example()
#' x <- mt_set_grid(x, age = single_ages(30, 100))
#' mt_grid(x)$age
#' @export
mt_set_grid <- function(data, age = NULL, period = NULL) {
  check_class(data, "mt_data", "data")
  g <- mt_grid(data)
  if (!is.null(age)) {
    check_class(age, "mt_axis", "age")
    g$age <- age
  }
  if (!is.null(period)) {
    check_class(period, "mt_axis", "period")
    g$period <- period
  }
  for (m in names(data$obs)) obs_grid_index(data$obs[[m]], g, m)
  data$grid <- g
  data
}

#' Map observed cells to grid indices
#'
#' @return A `data.table` with `a1`, `a2`, `p1`, `p2` (first and last grid
#'   cell covered on each axis), aligned with the rows of `obs`. Cells outside
#'   the grid get `NA`.
#' @keywords internal
obs_grid_index <- function(obs, grid, measure = "observations", strict = TRUE) {
  a1 <- match_bound(obs$age_lo, grid$age$lower)
  a2 <- match_bound(obs$age_hi, grid$age$upper)
  p1 <- match_bound(obs$period_lo, grid$period$lower)
  p2 <- match_bound(obs$period_hi, grid$period$upper)
  inside <- obs$age_lo >= grid$age$lower[1] - 1e-9 &
    obs$age_hi <= grid$age$upper[length(grid$age)] + 1e-9 &
    obs$period_lo >= grid$period$lower[1] - 1e-9 &
    obs$period_hi <= grid$period$upper[length(grid$period)] + 1e-9
  bad <- inside & (is.na(a1) | is.na(a2) | is.na(p1) | is.na(p2))
  if (strict && any(bad)) {
    i <- which(bad)[1]
    mt_abort(sprintf(paste0("The %s cell [%s, %s) x [%s, %s) is not a union of grid cells. ",
                            "Choose a finer grid with mt_set_grid() or a coarser one with mt_regroup()."),
                     measure, obs$age_lo[i], obs$age_hi[i], obs$period_lo[i], obs$period_hi[i]),
             class = "grid")
  }
  a1[!inside] <- a2[!inside] <- p1[!inside] <- p2[!inside] <- NA_integer_
  data.table(a1 = a1, a2 = a2, p1 = p1, p2 = p2)
}

match_bound <- function(x, b) {
  i <- match(round(x, 8), round(b, 8))
  i
}

#' Arrays on the modelling grid
#'
#' Returns a measure as an array `[age, period, stratum]` on the grid. Cells
#' observed at grid resolution are filled; cells that are not observed are
#' `NA`. Observed cells that are coarser than the grid raise an error unless
#' `coarse = "na"`.
#'
#' @param data An [mt_data] object.
#' @param measure Measure name.
#' @param coarse `"error"` or `"na"`.
#' @return A numeric array with dimnames.
#' @examples
#' x <- mt_example()
#' dim(mt_array(x, "exposure"))
#' @export
mt_array <- function(data, measure = "deaths", coarse = c("error", "na")) {
  check_class(data, "mt_data", "data")
  coarse <- match.arg(coarse)
  obs <- data$obs[[measure]]
  if (is.null(obs)) mt_abort(sprintf("No %s table.", measure), class = "missing_measure")
  g <- mt_grid(data)
  idx <- obs_grid_index(obs, g, measure)
  arr <- array(NA_real_, c(length(g$age), length(g$period), nrow(data$strata)),
               dimnames = list(age = g$age$labels, period = g$period$labels,
                               stratum = data$strata$label))
  fine <- !is.na(idx$a1) & idx$a1 == idx$a2 & idx$p1 == idx$p2
  is_coarse <- !is.na(idx$a1) & !fine
  if (any(is_coarse) && coarse == "error") {
    mt_abort(sprintf(paste0("%d %s cell(s) are coarser than the grid. ",
                            "Use mt_ungroup() to estimate them on the grid or mt_regroup() to coarsen the grid."),
                     sum(is_coarse), measure), class = "resolution")
  }
  arr[cbind(idx$a1[fine], idx$p1[fine], obs$stratum[fine])] <- obs$value[fine]
  arr
}

#' Restrict data to a window of periods
#'
#' Used by backtests to give models and preprocessing steps only the training
#' data. Cells that are not fully inside the window are dropped.
#'
#' @param data An [mt_data] object.
#' @param from,to First and last calendar year to keep (inclusive).
#' @return An `mt_data` object.
#' @examples
#' x <- mt_window(mt_example(), to = 2010)
#' mt_grid(x)$period
#' @export
mt_window <- function(data, from = -Inf, to = Inf) {
  check_class(data, "mt_data", "data")
  keep <- function(o) o[period_lo >= from - 1e-9 & period_hi <= to + 1 + 1e-9]
  data$obs <- lapply(data$obs, keep)
  data$obs <- data$obs[vapply(data$obs, nrow, 1L) > 0L]
  if (!is.null(data$grid)) {
    p <- data$grid$period
    k <- p$lower >= from - 1e-9 & p$upper <= to + 1 + 1e-9
    if (!any(k)) mt_abort("The window contains no grid period.", class = "grid")
    data$grid$period <- mt_axis(p$lower[k], p$upper[k], p$labels[k], kind = "period")
  }
  data
}

#' Describe the granularity of each table
#'
#' @param data An [mt_data] object.
#' @return A `data.table` with one row per measure and stratum: number of
#'   cells, age widths, start of the open interval, period widths and range.
#' @examples
#' mt_granularity(mt_example())
#' @export
mt_granularity <- function(data) {
  check_class(data, "mt_data", "data")
  rows <- lapply(names(data$obs), function(m) {
    o <- data$obs[[m]]
    o[, list(measure = m, cells = .N,
             ages = describe_widths(age_hi - age_lo, age_lo, age_hi),
             open_from = open_from_label(age_lo, age_hi),
             periods = sprintf("%s-%s", min(period_lo), max(period_hi) - 1),
             period_width = paste(unique(period_hi - period_lo), collapse = ",")),
      by = "stratum"]
  })
  out <- rbindlist(rows)
  out <- merge(data$strata[, c("stratum", "label")], out, by = "stratum")
  setcolorder(out, c("measure", "stratum", "label"))
  out[order(measure, stratum)]
}

open_from_label <- function(lo, hi) {
  if (any(is.infinite(hi))) paste(unique(lo[is.infinite(hi)]), collapse = ",") else "none"
}

describe_widths <- function(w, lo, hi) {
  w <- w[is.finite(w)]
  if (!length(w)) return(paste0("open interval ", min(lo), "+"))
  rng <- sprintf("%s-%s%s", min(lo), max(hi[is.finite(hi)]) - 1,
                 if (any(is.infinite(hi))) sprintf(", %s+", max(lo[is.infinite(hi)])) else "")
  u <- sort(unique(w))
  kind <- if (length(u) == 1L && u == 1) "single years" else if (length(u) == 1L) sprintf("%s-year groups", u) else
    sprintf("mixed widths (%s)", paste(u, collapse = ","))
  paste(kind, rng)
}

#' @export
print.mt_data <- function(x, ...) {
  ns <- nrow(x$strata)
  cat(sprintf("<mt_data> %d %s%s, omega = %s\n", ns, if (ns == 1) "stratum" else "strata",
              if (length(x$by)) sprintf(" by %s", paste(x$by, collapse = ", ")) else "", x$omega))
  if (!length(x$obs)) {
    cat("No observations.\n")
    return(invisible(x))
  }
  gr <- mt_granularity(x)
  for (m in unique(gr$measure)) {
    gm <- gr[measure == m]
    cat(sprintf("  %-10s %s; periods %s (width %s); %d cells\n", m,
                paste(unique(gm$ages), collapse = " | "),
                paste(unique(gm$periods), collapse = "|"), paste(unique(gm$period_width), collapse = "|"),
                sum(gm$cells)))
  }
  g <- mt_grid(x)
  cat(sprintf("  grid       %s\n             %s\n", format(g$age), format(g$period)))
  invisible(x)
}

#' @rdname mt_strata
#' @param x An [mt_data] object.
#' @param keep.rownames,... Ignored.
#' @method as.data.table mt_data
#' @export
as.data.table.mt_data <- function(x, keep.rownames = FALSE, ...) {
  rbindlist(lapply(names(x$obs), function(m) {
    o <- mt_obs(x, m)
    o[, measure := m]
    o
  }), fill = TRUE)
}

#' @export
as.data.frame.mt_data <- function(x, ...) as.data.frame(as.data.table.mt_data(x))
