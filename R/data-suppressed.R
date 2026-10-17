#' Estimate suppressed death counts
#'
#' Official statistics often suppress small cells for confidentiality (for
#' example `"."` in tables of the German Federal Statistical Office, `"<3"`
#' or `"<5"` elsewhere) while publishing totals: deaths of a region over all
#' causes, deaths of the whole country by cause. `mt_fill_suppressed()`
#' estimates the missing death counts of the bottom strata from two sources
#' of information:
#'
#' * **published totals** (`hierarchy`): the missing cells of every
#'   aggregate with a known total must add up to that total minus its known
#'   cells. They are estimated by iterative proportional fitting (raking):
#'   starting values, the expected deaths at the pooled death rate of the
#'   same age and period, are scaled constraint by constraint until every
#'   total is matched. With two crossed variables (regions and causes) the
#'   result is the log-linear independence fit of the missing cells given
#'   both margins. Missing aggregates are then filled as sums of their parts.
#' * **bounds of the suppressed counts** (`lower`, `upper`): the range that
#'   the suppression rule implies, for example `lower = 1, upper = 2` when
#'   counts of 1 and 2 are suppressed, or `upper = 4` for `"<5"`. Estimates
#'   are kept inside the bounds during raking. Cells that no published total
#'   covers get the expected deaths at the pooled rate, limited to the
#'   bounds; this needs a finite `upper`.
#'
#' Without bounds, cells that no published total covers stay missing (with a
#' warning). The estimates are plausible values that agree with the
#' published information; they carry no information beyond it, so forecasts
#' of strata with many estimated cells are less certain than their intervals
#' suggest. Filled cells are marked with provenance `"estimated"`. When the
#' known cells of an aggregate already exceed its total (inconsistent data),
#' the missing cells of that aggregate are set to the lower bound and the
#' case is reported.
#'
#' @param data An [mt_data] object in which some death cells are missing
#'   (cells that other strata have and that no observed cell of the stratum
#'   overlaps, so a stratum reported in other age groups is not filled;
#'   suppressed values are left out when the data are read).
#' @param hierarchy An [mt_hierarchy()] that describes the published totals
#'   in the data, or `NULL` when there are none.
#' @param lower,upper Bounds of the suppressed counts.
#' @param maxit,tol Iterations and tolerance (in deaths) of the raking; cells
#'   are raked until no estimate changes by more than `tol`.
#' @return The `mt_data` object with filled deaths; a summary in attribute
#'   `"fill"` (cells filled, cells left missing, inconsistent cells,
#'   iterations, largest relative gap to a published total).
#' @examples
#' d <- data.frame(region = c("A", "B", "C", "all"), year = 2020, age = "90+",
#'                 deaths = c(10, NA, NA, 40), exposure = c(100, 200, 200, 500))
#' x <- mt_data(deaths = d, exposure = d, by = "region")
#' y <- mt_fill_suppressed(x, mt_hierarchy(region = "all"))
#' mt_obs(y)
#' attr(y, "fill")
#' # counts below 3 suppressed, no totals: expected deaths limited to 1-2
#' d2 <- data.frame(region = c("A", "B"), year = 2020, age = "90+", deaths = c(30, NA),
#'                  exposure = c(1000, 50))
#' z <- mt_fill_suppressed(mt_data(deaths = d2, exposure = d2, by = "region"), lower = 1, upper = 2)
#' mt_obs(z)
#' @export
mt_fill_suppressed <- function(data, hierarchy = NULL, lower = 0, upper = Inf, maxit = 10000, tol = 1e-6) {
  check_class(data, "mt_data", "data")
  check_number(lower, "lower", lower = 0)
  check_number(upper, "upper", lower = lower)
  if (is.null(hierarchy) && !is.finite(upper)) {
    mt_abort(paste("Give the published totals (`hierarchy`) or the range of the suppressed counts (`upper`, and",
                   "`lower`)."), class = "fill")
  }
  D <- data$obs$deaths
  if (is.null(D)) mt_abort("No deaths table.", class = "missing_measure")
  ns <- nrow(data$strata)
  if (!is.null(hierarchy)) {
    h <- as_hierarchy(hierarchy)
    hs <- hierarchy_structure(data$strata, h)
    if (all(hs$is_bottom)) {
      tot <- vapply(names(h$vars), function(v) sprintf("%s = \"%s\"", v, h$vars[[v]]$total), "")
      mt_abort(sprintf(paste("None of the totals of the hierarchy (%s) is in the data, so they cannot constrain",
                             "the suppressed cells. Read the published totals as strata, or give the range of the",
                             "suppressed counts with `upper` (and `lower`) and no hierarchy."),
                       paste(tot, collapse = ", ")), class = "fill")
    }
  } else {
    hs <- list(bottom = seq_len(ns), is_bottom = rep(TRUE, ns), S = diag(ns))
  }
  D <- D[!is.na(value)]
  cells <- unique(D[, c("age_lo", "age_hi", "period_lo", "period_hi")])
  cells[, cell := seq_len(.N)]
  Dk <- merge(D, cells, by = c("age_lo", "age_hi", "period_lo", "period_hi"))
  # a cell that some strata have and others lack is missing (suppressed)
  # there, unless an observed cell of the stratum overlaps it: then the
  # stratum reports these ages and periods at another resolution
  Y <- matrix(NA_real_, nrow(cells), ns)
  Y[cbind(Dk$cell, Dk$stratum)] <- Dk$value
  ov <- D[cells, on = c("age_lo<age_hi", "age_hi>age_lo", "period_lo<period_hi", "period_hi>period_lo"),
          list(cell, stratum), nomatch = NULL, allow.cartesian = TRUE]
  other_res <- matrix(FALSE, nrow(cells), ns)
  other_res[cbind(ov$cell, ov$stratum)] <- TRUE
  other_res <- other_res & is.na(Y)
  if (!any(is.na(Y) & !other_res)) {
    attr(data, "fill") <- data.table(filled = 0L, missing = 0L, inconsistent = 0L, iterations = 0L, max_gap = 0)
    return(data)
  }
  # starting values: expected deaths at the pooled rate of the cell (observed
  # bottom strata), the exposure being the only information on the cell
  Ex <- matrix(NA_real_, nrow(cells), ns)
  E <- data$obs$exposure
  if (!is.null(E)) {
    Ek <- merge(E, cells, by = c("age_lo", "age_hi", "period_lo", "period_hi"))
    ok <- !is.na(Ek$value) & Ek$value > 0
    Ex[cbind(Ek$cell[ok], Ek$stratum[ok])] <- Ek$value[ok]
  }
  bb <- hs$bottom
  obs_b <- !is.na(Y[, bb, drop = FALSE]) & !is.na(Ex[, bb, drop = FALSE])
  ref <- rowSums(ifelse(obs_b, Y[, bb, drop = FALSE], 0)) / pmax(rowSums(ifelse(obs_b, Ex[, bb, drop = FALSE], 0)),
                                                                   1e-12)
  start <- Ex * ref
  start[is.na(start)] <- 1
  start <- pmin(pmax(start, max(lower, 1e-6)), upper)
  b <- hs$bottom
  agg <- which(!hs$is_bottom)
  clip <- function(v) pmin(pmax(v, lower), upper)
  members <- lapply(seq_along(agg), function(i) b[hs$S[agg[i], ] > 0])
  miss_b <- is.na(Y[, b, drop = FALSE]) & !other_res[, b, drop = FALSE]
  colnames(miss_b) <- b
  # rows where the known total and the members can be compared (no member at another resolution)
  usable <- function(i) {
    which(!is.na(Y[, agg[i]]) & rowSums(miss_b[, as.character(members[[i]]), drop = FALSE]) > 0 &
            rowSums(other_res[, members[[i]], drop = FALSE]) == 0)
  }
  est <- Y
  est[, b][miss_b] <- start[, b][miss_b]
  inconsistent <- 0L
  # constraints with a known total and missing members, and the cells concerned
  cons <- lapply(seq_along(agg), function(i) {
    mk <- miss_b[, as.character(members[[i]]), drop = FALSE]
    rows <- usable(i)
    if (!length(rows)) return(NULL)
    known <- rowSums(ifelse(mk[rows, , drop = FALSE], 0, Y[rows, members[[i]], drop = FALSE]), na.rm = TRUE)
    resid <- Y[rows, agg[i]] - known
    inconsistent <<- inconsistent + sum(resid < -0.5)
    list(mem = members[[i]], rows = rows, mk = mk[rows, , drop = FALSE], target = pmax(resid, 0))
  })
  cons <- Filter(Negate(is.null), cons)
  active <- rep(FALSE, nrow(Y))
  active[unique(unlist(lapply(cons, `[[`, "rows")))] <- TRUE
  it <- 0L
  while (any(active) && it < maxit) {
    it <- it + 1L
    change <- numeric(nrow(Y))
    for (cn in cons) {
      use <- active[cn$rows]
      if (!any(use)) next
      rows <- cn$rows[use]
      mk <- cn$mk[use, , drop = FALSE]
      old <- est[rows, cn$mem, drop = FALSE]
      cur <- rowSums(old * mk)
      f <- ifelse(cur > 0, cn$target[use] / cur, 0)
      upd <- ifelse(mk, clip(old * f), old)
      est[rows, cn$mem] <- upd
      change[rows] <- pmax(change[rows], apply(abs(upd - old), 1, max))
    }
    active <- active & change > tol
  }
  # cells that no published total covers: expected deaths limited to the
  # bounds, or left missing without a finite upper bound
  covered <- matrix(FALSE, nrow(Y), length(b))
  colnames(covered) <- b
  for (cn in cons) {
    idx <- cbind(rep(cn$rows, ncol(cn$mk)), rep(match(cn$mem, b), each = length(cn$rows)))
    covered[idx[as.vector(cn$mk), , drop = FALSE]] <- TRUE
  }
  uncovered <- miss_b & !covered
  if (any(uncovered)) {
    eb <- est[, b, drop = FALSE]
    eb[uncovered] <- if (is.finite(upper)) start[, b, drop = FALSE][uncovered] else NA
    est[, b] <- eb
  }
  # missing aggregates as sums of their parts
  for (i in seq_along(agg)) {
    a <- agg[i]
    rows <- which(is.na(Y[, a]) & !other_res[, a])
    if (length(rows)) est[rows, a] <- rowSums(est[rows, members[[i]], drop = FALSE])
  }
  gap <- 0
  for (i in seq_along(agg)) {
    a <- agg[i]
    rows <- usable(i)
    if (length(rows)) {
      s <- rowSums(est[rows, members[[i]], drop = FALSE])
      gap <- max(gap, abs(s - Y[rows, a]) / pmax(Y[rows, a], 1))
    }
  }
  out <- CJ(cell = cells$cell, stratum = seq_len(ns))
  out <- out[!other_res[cbind(cell, stratum)]]
  out <- merge(out, cells, by = "cell")
  out[, value := est[cbind(cell, stratum)]]
  out <- merge(out, Dk[, c("cell", "stratum", "provenance")], by = c("cell", "stratum"), all.x = TRUE)
  filled <- is.na(out$provenance)
  out[filled, provenance := "estimated"]
  still <- sum(is.na(out$value))
  n_filled <- sum(filled & !is.na(out$value))
  out <- out[!is.na(value)][, cell := NULL]
  setcolorder(out, c("stratum", "age_lo", "age_hi", "period_lo", "period_hi", "value", "provenance"))
  setorderv(out, c("stratum", "period_lo", "age_lo"))
  data$obs$deaths <- out[]
  if (still) {
    mt_warn(sprintf(paste("%d cell(s) stay missing because no published total covers them; give the range of the",
                          "suppressed counts with `upper` to estimate them."), still), class = "fill")
  }
  if (inconsistent) {
    mt_warn(sprintf(paste("In %d cell(s) the known parts exceed the published total; their missing parts are set",
                          "to the lower bound."), inconsistent), class = "fill")
  }
  if (gap > 0.01) {
    mt_warn(sprintf(paste("The estimates miss a published total by up to %.1f%%: the bounds cannot be met together",
                          "with the totals."), 100 * gap), class = "fill")
  }
  attr(data, "fill") <- data.table(filled = n_filled, missing = still,
                                   inconsistent = inconsistent, iterations = it, max_gap = gap)
  data
}
