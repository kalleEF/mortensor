#' Regroup data onto coarser age or period intervals
#'
#' Aggregates every table exactly onto the target axes. Counts (deaths and
#' exposures) are summed; populations are summed over age and averaged over
#' periods (weighted by period width). An observed cell that straddles a
#' target boundary raises an error, because exact aggregation is not possible.
#' Target cells that are only partly covered by observations are dropped with
#' a warning. The grid becomes the target axes, restricted to the range of the
#' data.
#'
#' @param data An [mt_data] object.
#' @param age Target age axis ([mt_axis]), numeric lower bounds, or the
#'   name of a table (`"deaths"`, `"exposure"`, `"population"`) whose age
#'   intervals are the target (for example deaths in 5-year groups and
#'   populations in single years: `age = "deaths"`); `NULL` keeps the age
#'   intervals.
#' @param period Target period axis, numeric lower bounds or the name of a
#'   table; `NULL` keeps the periods.
#' @return An `mt_data` object.
#' @examples
#' x <- mt_example()
#' mt_regroup(x, age = c(30, seq(40, 90, 10)))
#' mt_regroup(x, age = "deaths")
#' @export
mt_regroup <- function(data, age = NULL, period = NULL) {
  check_class(data, "mt_data", "data")
  g <- mt_grid(data)
  if (is.character(age)) age <- table_axis(data, age, "age")
  if (is.character(period)) period <- table_axis(data, period, "period")
  if (is.numeric(age)) {
    age <- if (is_open(g$age)) age_groups(age) else age_groups(age, open = FALSE, last_upper = max(g$age$upper))
  }
  if (is.numeric(period)) period <- period_groups(period, last_upper = max(g$period$upper))
  if (!is.null(age)) check_class(age, "mt_axis", "age")
  if (!is.null(period)) check_class(period, "mt_axis", "period")
  for (m in names(data$obs)) {
    data$obs[[m]] <- regroup_obs(data$obs[[m]], m, age, period, data$omega)
  }
  data$obs <- data$obs[vapply(data$obs, nrow, 1L) > 0L]
  all <- rbindlist(lapply(data$obs, function(o) o[, c("age_lo", "age_hi", "period_lo", "period_hi")]))
  new_age <- if (is.null(age)) g$age else axis_trim(age, min(all$age_lo), max(all$age_hi))
  new_period <- if (is.null(period)) g$period else axis_trim(period, min(all$period_lo), max(all$period_hi))
  data$grid <- NULL
  mt_set_grid(data, age = new_age, period = new_period)
}

# The intervals of one table as an axis
table_axis <- function(data, measure, kind) {
  if (length(measure) != 1 || !measure %in% names(data$obs)) {
    mt_abort(sprintf("`%s` must be an axis, numeric bounds or one of the tables (%s).", kind,
                     paste(names(data$obs), collapse = ", ")), class = "argument")
  }
  cols <- paste0(kind, c("_lo", "_hi"))
  iv <- unique(data$obs[[measure]][, cols, with = FALSE])
  setorderv(iv, cols)
  if (anyDuplicated(iv[[1]]) || any(iv[[2]][-nrow(iv)] != iv[[1]][-1])) {
    mt_abort(sprintf(paste("The %s intervals of the %s table overlap or have gaps (for example different groups",
                           "in different years); give the target intervals explicitly."), kind, measure),
             class = "regroup")
  }
  mt_axis(iv[[1]], iv[[2]], kind = kind)
}

axis_trim <- function(axis, lo, hi) {
  k <- axis$lower >= lo - 1e-9 & axis$upper <= hi + 1e-9
  if (!any(k)) mt_abort("No target interval is covered by the data.", class = "grid")
  k <- seq(min(which(k)), max(which(k)))
  mt_axis(axis$lower[k], axis$upper[k], axis$labels[k], kind = axis$kind)
}

regroup_obs <- function(o, measure, age, period, omega) {
  o <- copy(o)
  o[, `:=`(nalo = age_lo, nahi = age_hi, nplo = period_lo, nphi = period_hi)]
  if (!is.null(age)) {
    ia <- axis_locate(age, o$age_lo, o$age_hi)
    straddle <- is.na(ia) & o$age_lo < max(age$upper) & o$age_hi > min(age$lower)
    if (any(straddle)) {
      mt_abort(sprintf("The %s age interval [%s, %s) straddles a target boundary.", measure,
                       o$age_lo[straddle][1], o$age_hi[straddle][1]), class = "regroup")
    }
    o[, `:=`(nalo = age$lower[ia], nahi = age$upper[ia])]
  }
  if (!is.null(period)) {
    ip <- axis_locate(period, o$period_lo, o$period_hi)
    straddle <- is.na(ip) & o$period_lo < max(period$upper) & o$period_hi > min(period$lower)
    if (any(straddle)) {
      mt_abort(sprintf("The %s period [%s, %s) straddles a target boundary.", measure,
                       o$period_lo[straddle][1], o$period_hi[straddle][1]), class = "regroup")
    }
    o[, `:=`(nplo = period$lower[ip], nphi = period$upper[ip])]
  }
  o <- o[!is.na(nalo) & !is.na(nplo)]
  aw <- function(lo, hi) pmin(hi, omega) - lo
  o[, area := aw(age_lo, age_hi) * (period_hi - period_lo)]
  # stage 1: over age within original periods
  s1 <- o[, list(value = sum(value), area = sum(area), pw = period_hi[1] - period_lo[1]),
          by = c("stratum", "nalo", "nahi", "period_lo", "period_hi", "nplo", "nphi")]
  # stage 2: over periods
  s2 <- if (measure == "population") {
    s1[, list(value = sum(value * pw) / sum(pw), area = sum(area)), by = c("stratum", "nalo", "nahi", "nplo", "nphi")]
  } else {
    s1[, list(value = sum(value), area = sum(area)), by = c("stratum", "nalo", "nahi", "nplo", "nphi")]
  }
  s2[, full := aw(nalo, nahi) * (nphi - nplo)]
  partial <- s2$area < s2$full * (1 - 1e-9)
  if (any(partial)) {
    mt_warn(sprintf("%d %s cell(s) only partly covered by observations were dropped.", sum(partial), measure),
            class = "regroup_partial")
    s2 <- s2[!partial]
  }
  out <- s2[, list(stratum, age_lo = nalo, age_hi = nahi, period_lo = nplo, period_hi = nphi, value,
                   provenance = "aggregated")]
  setorderv(out, c("stratum", "period_lo", "age_lo"))
  out[]
}

#' Aggregation operator from grid cells to observed cells
#'
#' A sparse matrix with one row per observed cell of one stratum and one
#' column per grid cell (age varies fastest), with entry 1 where the grid cell
#' lies inside the observed cell. It links latent quantities on the grid to
#' what is observed, for example in the penalized composite link model.
#'
#' @param data An [mt_data] object.
#' @param measure Measure name.
#' @param stratum Stratum id.
#' @return A `dgCMatrix`. Attribute `rows` gives the observed cells.
#' @examples
#' x <- mt_example()
#' C <- mt_operator(x, "deaths", 1)
#' dim(C)
#' @export
mt_operator <- function(data, measure, stratum) {
  check_class(data, "mt_data", "data")
  g <- mt_grid(data)
  o <- data$obs[[measure]]
  rows_in_stratum <- which(o[["stratum"]] == stratum)
  o <- o[rows_in_stratum]
  idx <- obs_grid_index(o, g, measure)
  keep <- !is.na(idx$a1)
  o <- o[keep]
  idx <- idx[keep]
  build_operator(idx, length(g$age), length(g$period), rows = o)
}

build_operator <- function(idx, n_age, n_period, rows = NULL) {
  na <- idx$a2 - idx$a1 + 1L
  np <- idx$p2 - idx$p1 + 1L
  n <- na * np
  r <- rep(seq_len(nrow(idx)), n)
  a <- unlist(Map(function(a1, a2, p1, p2) rep(a1:a2, times = p2 - p1 + 1L), idx$a1, idx$a2, idx$p1, idx$p2))
  p <- unlist(Map(function(a1, a2, p1, p2) rep(p1:p2, each = a2 - a1 + 1L), idx$a1, idx$a2, idx$p1, idx$p2))
  C <- Matrix::sparseMatrix(i = r, j = a + (p - 1L) * n_age, x = 1, dims = c(nrow(idx), n_age * n_period))
  attr(C, "rows") <- rows
  C
}
