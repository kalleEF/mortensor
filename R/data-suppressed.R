#' Estimate suppressed death counts from published totals
#'
#' Official statistics often suppress small cells for confidentiality (for
#' example `"."` in tables of the German Federal Statistical Office) while
#' publishing the totals: deaths of a region over all causes, deaths of the
#' whole country by cause. `mt_fill_suppressed()` estimates the missing death
#' counts of the bottom strata so that every published total of the hierarchy
#' is matched, by iterative proportional fitting (raking): starting values
#' proportional to the exposures are scaled, constraint by constraint, until
#' the missing cells of every aggregate with a known total add up to that
#' total minus its known cells. With two crossed variables (regions and
#' causes) the result is the log-linear independence fit of the missing
#' cells given both margins. Missing aggregates are then filled as sums of
#' their (filled) parts. Death counts sum over every hierarchy variable,
#' also over causes, whatever the `additive` setting of the hierarchy.
#'
#' Filled cells are marked with provenance `"estimated"`. When the known
#' cells of an aggregate already exceed its total (inconsistent data), the
#' missing cells of that aggregate are set to zero and the case is reported.
#'
#' @param data An [mt_data] object in which some death cells are missing
#'   (cells that other strata have; rows with `NA` values are dropped when
#'   the data are created).
#' @param hierarchy An [mt_hierarchy()] that describes the totals in the data.
#' @param maxit,tol Iterations and tolerance (in deaths) of the raking; cells
#'   are raked until no estimate changes by more than `tol`.
#' @return The `mt_data` object with filled deaths; a summary in attribute
#'   `"fill"`.
#' @examples
#' d <- data.frame(region = c("A", "B", "C", "all"), year = 2020, age = "90+",
#'                 deaths = c(10, NA, NA, 40), exposure = c(100, 200, 200, 500))
#' x <- mt_data(deaths = d, exposure = d, by = "region")
#' y <- mt_fill_suppressed(x, mt_hierarchy(region = "all"))
#' mt_obs(y)
#' attr(y, "fill")
#' @export
mt_fill_suppressed <- function(data, hierarchy, maxit = 10000, tol = 1e-6) {
  check_class(data, "mt_data", "data")
  h <- as_hierarchy(hierarchy)
  hs <- hierarchy_structure(data$strata, h)
  D <- data$obs$deaths
  if (is.null(D)) mt_abort("No deaths table.", class = "missing_measure")
  ns <- nrow(data$strata)
  D <- D[!is.na(value)]
  cells <- unique(D[, c("age_lo", "age_hi", "period_lo", "period_hi")])
  cells[, cell := seq_len(.N)]
  Dk <- merge(D, cells, by = c("age_lo", "age_hi", "period_lo", "period_hi"))
  # a cell that some strata have and others lack is missing (suppressed) there
  Y <- matrix(NA_real_, nrow(cells), ns)
  Y[cbind(Dk$cell, Dk$stratum)] <- Dk$value
  if (!anyNA(Y)) {
    attr(data, "fill") <- data.table(filled = 0L, inconsistent = 0L, iterations = 0L, max_gap = 0)
    return(data)
  }
  start <- matrix(1, nrow(cells), ns)
  E <- data$obs$exposure
  if (!is.null(E)) {
    Ek <- merge(E, cells, by = c("age_lo", "age_hi", "period_lo", "period_hi"))
    ok <- !is.na(Ek$value) & Ek$value > 0
    start[cbind(Ek$cell[ok], Ek$stratum[ok])] <- Ek$value[ok]
  }
  b <- hs$bottom
  agg <- which(!hs$is_bottom)
  members <- lapply(seq_along(agg), function(i) b[hs$S[agg[i], ] > 0])
  miss_b <- is.na(Y[, b, drop = FALSE])
  colnames(miss_b) <- b
  est <- Y
  est[, b][miss_b] <- start[, b][miss_b]
  inconsistent <- 0L
  # constraints with a known total and missing members, and the cells concerned
  cons <- lapply(seq_along(agg), function(i) {
    mk <- miss_b[, as.character(members[[i]]), drop = FALSE]
    rows <- which(!is.na(Y[, agg[i]]) & rowSums(mk) > 0)
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
      upd <- ifelse(mk, old * f, old)
      est[rows, cn$mem] <- upd
      change[rows] <- pmax(change[rows], apply(abs(upd - old), 1, max))
    }
    active <- active & change > tol
  }
  # missing aggregates as sums of their parts
  for (i in seq_along(agg)) {
    a <- agg[i]
    rows <- which(is.na(Y[, a]))
    if (length(rows)) est[rows, a] <- rowSums(est[rows, members[[i]], drop = FALSE])
  }
  gap <- 0
  for (i in seq_along(agg)) {
    a <- agg[i]
    rows <- which(!is.na(Y[, a]) & rowSums(miss_b[, as.character(members[[i]]), drop = FALSE]) > 0)
    if (length(rows)) {
      s <- rowSums(est[rows, members[[i]], drop = FALSE])
      gap <- max(gap, abs(s - Y[rows, a]) / pmax(Y[rows, a], 1))
    }
  }
  out <- CJ(cell = cells$cell, stratum = seq_len(ns))
  out <- merge(out, cells, by = "cell")
  out[, value := est[cbind(cell, stratum)]]
  out <- merge(out, Dk[, c("cell", "stratum", "provenance")], by = c("cell", "stratum"), all.x = TRUE)
  filled <- is.na(out$provenance)
  out[filled, provenance := "estimated"]
  still <- sum(is.na(out$value))
  out <- out[!is.na(value)][, cell := NULL]
  setcolorder(out, c("stratum", "age_lo", "age_hi", "period_lo", "period_hi", "value", "provenance"))
  setorderv(out, c("stratum", "period_lo", "age_lo"))
  data$obs$deaths <- out[]
  if (still) {
    mt_warn(sprintf("%d cell(s) could not be estimated (no total covers them).", still), class = "fill")
  }
  if (inconsistent) {
    mt_warn(sprintf("In %d cell(s) the known parts exceed the published total; their missing parts are set to 0.",
                    inconsistent), class = "fill")
  }
  attr(data, "fill") <- data.table(filled = sum(filled) - still, inconsistent = inconsistent, iterations = it,
                                   max_gap = gap)
  data
}
