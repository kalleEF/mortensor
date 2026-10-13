#' Estimate grid-resolution data from grouped observations
#'
#' Fills the modelling grid (single years of age by default) from deaths and
#' exposures observed in coarser age groups, using the penalized composite
#' link model (PCLM; Eilers 2007, Rizzi et al. 2015).
#'
#' * **Exposures** observed coarser than the grid are ungrouped first, as
#'   counts (`exposure_method = "pclm"`) or by a monotone spline of the
#'   cumulative counts (`"monotone"`, Hyman filter), which preserves group
#'   totals exactly.
#' * **Deaths** are ungrouped as **rates** when exposures are available on
#'   the grid: the exposures enter as offset, so cohort-size waves stay in the
#'   exposures and the estimated rate schedule is smooth. Without exposures,
#'   deaths are ungrouped as counts per year of age (1D only; 2D smoothing of
#'   death counts cannot represent cohort waves and is refused).
#' * The open age interval is expanded into single years up to `omega`; for
#'   rates the exposure of the open interval is spread with
#'   stationary-population weights implied by the fitted rates.
#' * `dims = 1` fits each stratum and period separately; `dims = 2` fits a
#'   smooth age-period surface per stratum (two smoothing parameters).
#'
#' With `preserve_totals = TRUE`, the estimated deaths within each observed
#' group are rescaled to the observed total, so aggregation reproduces the
#' data exactly. Cells already observed at grid resolution are kept.
#'
#' In backtests use `step_ungroup()` so that ungrouping, in particular in 2D,
#' sees only the training window.
#'
#' @param data An [mt_data] object.
#' @param age Target age axis ([mt_axis]). By default the grid when it was
#'   set (with [mt_set_grid()] or `step_grid()`) and is finer than the data,
#'   otherwise single years of age from the lowest age to the start of the
#'   open interval.
#' @param dims 1 or 2.
#' @param method Method for deaths: `"pclm"` (default); `"monotone"`
#'   (cumulative deaths interpolated by a monotone Hyman spline, counts only,
#'   1D, ignores exposures); or rate-based spline methods, which
#'   need exposures on the grid: `"hyman_kannisto"` (monotone spline of log
#'   group rates with a Kannisto tail for the open interval),
#'   `"natural_gompertz"` (natural spline with a log-linear tail) and `"gam"`
#'   (2D tensor-product GAM of log group rates over age and year with a
#'   Kannisto tail, package mgcv). All three preserve the exposure-weighted
#'   group rates exactly. Median absolute log errors against single-year
#'   rates (England & Wales men 1980-2011, 5-year groups to 90-100; and the
#'   simulated [mt_example()] with deaths to 90+): PCLM 0.015/0.024 and
#'   0.012/0.006 (ages 30-89/90+), `"hyman_kannisto"` 0.015/0.029 and
#'   0.013/0.028, `"natural_gompertz"` 0.015/0.027 and 0.013/0.030, `"gam"`
#'   0.015/0.031 and 0.012/0.031. PCLM is the most accurate in the open
#'   interval; the spline methods are 5 to 10 times faster.
#' @param exposure_method `"pclm"` or `"monotone"`.
#' @param lambda Fixed smoothing parameter(s); `NULL` selects them by
#'   `criterion`.
#' @param criterion `"bic"` or `"aic"`.
#' @param kr Knot spacing in years (age and period).
#' @param deg,pord Degree of the B-splines and order of the difference
#'   penalty.
#' @param lambda_range Search range for the smoothing parameters.
#' @param preserve_totals Rescale estimates to observed group totals.
#' @return An `mt_data` object with deaths (and exposures) on the grid.
#'   Diagnostics per fit are attached as attribute `"ungroup"`: smoothing
#'   parameters, effective dimension, the largest relative difference between
#'   fitted and observed group totals before rescaling, and a flag when a
#'   smoothing parameter is at the boundary of the search range.
#' @references
#' Eilers PHC (2007). Ill-posed problems with counts, the composite link model
#' and penalized likelihood. *Statistical Modelling* 7(3): 239-254.
#'
#' Rizzi S, Gampe J, Eilers PHC (2015). Efficient estimation of smooth
#' distributions from coarsely grouped data. *American Journal of
#' Epidemiology* 182(2): 138-147. \doi{10.1093/aje/kwv020}
#' @examples
#' x <- mt_example(years = 2010:2014)
#' u <- mt_ungroup(x)
#' u
#' head(attr(u, "ungroup"))
#' @export
mt_ungroup <- function(data, age = NULL, dims = 1,
                       method = c("pclm", "monotone", "hyman_kannisto", "natural_gompertz", "gam"),
                       exposure_method = c("pclm", "monotone"),
                       lambda = NULL,
                       criterion = c("bic", "aic"), kr = 2, deg = 3, pord = 2,
                       lambda_range = c(1e-2, 1e8), preserve_totals = TRUE) {
  check_class(data, "mt_data", "data")
  method <- match.arg(method)
  exposure_method <- match.arg(exposure_method)
  criterion <- match.arg(criterion)
  check_number(dims, "dims", lower = 1, upper = 2, integerish = TRUE)
  g <- mt_grid(data)
  if (!is.null(age)) {
    check_class(age, "mt_axis", "age")
    g$age <- age
  } else if (is.null(data$grid) || !any(vapply(data$obs, needs_ungroup, logical(1), g = g))) {
    g$age <- single_year_grid(data)
  }
  ctl <- list(lambda = lambda, criterion = criterion, kr = kr, deg = deg, pord = pord,
              lambda_range = lambda_range, preserve_totals = preserve_totals)
  diag_rows <- list()
  if (!is.null(data$obs$exposure) && needs_ungroup(data$obs$exposure, g)) {
    res <- ungroup_counts_1d(data$obs$exposure, g, data$omega, ctl, exposure_method, "exposure")
    data$obs$exposure <- res$obs
    diag_rows <- c(diag_rows, list(res$diag))
  }
  if (!is.null(data$obs$deaths) && needs_ungroup(data$obs$deaths, g)) {
    E <- data$obs$exposure
    if (method == "monotone") {
      res <- ungroup_counts_1d(data$obs$deaths, g, data$omega, ctl, "monotone", "deaths")
    } else if (method %in% c("hyman_kannisto", "natural_gompertz", "gam")) {
      if (is.null(E) || needs_ungroup(E, g)) {
        mt_abort(sprintf("method = \"%s\" needs exposures on the grid.", method), class = "ungroup")
      }
      res <- list(obs = ungroup_rates_spline(data$obs$deaths, E, g, data$omega, method),
                  diag = data.table(measure = "deaths", stratum = NA_integer_, period = "all", method = method))
    } else if (!is.null(E) && !needs_ungroup(E, g)) {
      res <- ungroup_rates(data$obs$deaths, E, g, data$omega, ctl, dims)
    } else {
      if (dims == 2) {
        mt_abort(paste0("2D ungrouping of death counts without exposures is not supported: ",
                        "smoothing counts over time cannot represent cohort-size waves. ",
                        "Add exposures on the grid or use dims = 1."), class = "ungroup")
      }
      mt_warn(paste0("Deaths are ungrouped as counts without exposures; the split of the open ",
                     "interval is weakly identified."), class = "ungroup_counts")
      res <- ungroup_counts_1d(data$obs$deaths, g, data$omega, ctl, "pclm", "deaths")
    }
    data$obs$deaths <- res$obs
    diag_rows <- c(diag_rows, list(res$diag))
  }
  dg <- rbindlist(diag_rows, fill = TRUE)
  if (nrow(dg) && any(dg$at_bound)) {
    mt_warn(sprintf("The smoothing parameter is at the boundary of the search range in %d of %d fits.",
                    sum(dg$at_bound), nrow(dg)), class = "lambda_bound")
  }
  data$grid <- g
  attr(data, "ungroup") <- dg
  data
}

needs_ungroup <- function(obs, g) {
  idx <- obs_grid_index(obs, g)
  any(!is.na(idx$a1) & (idx$a1 != idx$a2 | idx$p1 != idx$p2))
}

# Counts (exposures, or deaths without exposures), 1D per stratum and period
ungroup_counts_1d <- function(obs, g, omega, ctl, method, measure) {
  idx <- obs_grid_index(obs, g, measure)
  obs <- cbind(obs, idx)[!is.na(a1)]
  if (any(obs$p1 != obs$p2)) {
    mt_abort(sprintf("The %s table has periods coarser than the grid; 1D ungrouping needs grid periods.", measure),
             class = "ungroup")
  }
  fa <- fine_axis(g$age, omega)
  x <- (fa$lower + fa$upper) / 2
  bas <- pclm_basis(x, min(fa$lower), max(fa$upper), ctl$kr, ctl$deg, ctl$pord)
  width <- fa$upper - fa$lower
  open_split <- is.infinite(obs$age_hi) & obs$a1 != obs$a2
  if (any(open_split)) {
    mt_warn(sprintf(paste0("The open %s interval from age %s is split into %d grid cells from its total alone; ",
                           "the age pattern inside it is extrapolated and can be badly wrong. Prefer data with ",
                           "a higher open age or a grid whose open interval starts at %s."),
                    measure, min(obs$age_lo[open_split]), max(obs$a2 - obs$a1 + 1), min(obs$age_lo[open_split])),
            class = "ungroup_open_counts", call = NULL)
  }
  out <- list()
  dg <- list()
  for (key in split(seq_len(nrow(obs)), list(obs$stratum, obs$p1), drop = TRUE)) {
    o <- obs[key]
    hi <- pmin(o$age_hi, omega)
    if (all(o$a1 == o$a2)) {
      out[[length(out) + 1L]] <- o[, list(stratum, a = a1, p = p1, value, provenance)]
      next
    }
    C <- link_matrix(o$age_lo, hi, fa$lower, fa$upper)
    cov <- colSums(C) > 0
    if (method == "monotone") {
      est <- monotone_split(o$age_lo, hi, o$value, fa$lower[cov], fa$upper[cov])
      fit <- list(lambda = NA_real_, ed = NA_real_, at_bound = FALSE, mu = o$value)
    } else {
      fit <- pclm_fit(o$value, C[, cov, drop = FALSE], width[cov], bas$B[cov, , drop = FALSE], list(bas$P),
                      ctl$lambda, ctl$criterion, ctl$lambda_range)
      est <- width[cov] * fit$gamma
    }
    rel <- max(abs(fit$mu / pmax(o$value, 1e-12) - 1)[o$value > 0], 0)
    if (ctl$preserve_totals) est <- rescale_to_groups(est, C[, cov, drop = FALSE], o$value)
    gi <- fa$grid[cov]
    v <- tapply(est, gi, sum)
    exact <- o[a1 == a2]
    res <- data.table(stratum = o$stratum[1], a = as.integer(names(v)), p = o$p1[1], value = as.numeric(v),
                      provenance = method)
    res[exact, on = c(a = "a1"), `:=`(value = i.value, provenance = i.provenance)]
    out[[length(out) + 1L]] <- res
    dg[[length(dg) + 1L]] <- data.table(measure = measure, stratum = o$stratum[1], period = g$period$labels[o$p1[1]],
                                        method = method, lambda_age = fit$lambda[1], lambda_period = NA_real_,
                                        ed = fit$ed,
                                        max_rel_diff = rel, at_bound = isTRUE(fit$at_bound))
  }
  list(obs = grid_table(rbindlist(out), g), diag = rbindlist(dg))
}

# Deaths as rates with grid exposures as offset
ungroup_rates <- function(D, E, g, omega, ctl, dims) {
  idx <- obs_grid_index(D, g, "deaths")
  D <- cbind(D, idx)[!is.na(a1)]
  Earr <- obs_to_matrix(E, g)
  fa <- fine_axis(g$age, omega)
  xa <- (fa$lower + fa$upper) / 2
  n_open <- sum(fa$grid == length(g$age) & fa$open)
  bas <- pclm_basis(xa, min(fa$lower), max(fa$upper), ctl$kr, ctl$deg, ctl$pord)
  out <- list()
  dg <- list()
  for (s in unique(D$stratum)) {
    Ds <- D[stratum == s]
    if (dims == 1) {
      if (any(Ds$p1 != Ds$p2)) {
        mt_abort("Deaths have periods coarser than the grid; use dims = 2.", class = "ungroup")
      }
      for (p in sort(unique(Ds$p1))) {
        o <- Ds[p1 == p]
        e_grid <- Earr[, p, s]
        r <- ungroup_rate_block(o, list(e_grid), fa, bas$B, list(bas$P), ctl, omega, n_open)
        out[[length(out) + 1L]] <- data.table(stratum = s, a = r$a, p = p, value = r$deaths,
                                              provenance = r$provenance)
        dg[[length(dg) + 1L]] <- data.table(measure = "deaths", stratum = s, period = g$period$labels[p],
                                            method = paste0(ctl$method_label %||% "pclm_rate", "_1d"),
                                            lambda_age = r$fit$lambda[1],
                                            lambda_period = NA_real_, ed = r$fit$ed,
                                            max_rel_diff = r$rel, at_bound = isTRUE(r$fit$at_bound))
      }
    } else {
      ps <- seq(min(Ds$p1), max(Ds$p2))
      xp <- (g$period$lower[ps] + g$period$upper[ps]) / 2
      bp <- pclm_basis(xp, g$period$lower[min(ps)], g$period$upper[max(ps)], ctl$kr, ctl$deg, ctl$pord)
      B2 <- kronecker(bp$B, bas$B)
      pen <- list(kronecker(diag(ncol(bp$B)), bas$P), kronecker(bp$P, diag(ncol(bas$B))))
      r <- ungroup_rate_block(Ds, lapply(ps, function(p) Earr[, p, s]), fa, B2, pen, ctl, omega, n_open, ps)
      out[[length(out) + 1L]] <- data.table(stratum = s, a = r$a, p = r$p, value = r$deaths,
                                            provenance = r$provenance)
      dg[[length(dg) + 1L]] <- data.table(measure = "deaths", stratum = s, period = "all",
                                          method = paste0(ctl$method_label %||% "pclm_rate", "_2d"),
                                          lambda_age = r$fit$lambda[1],
                                          lambda_period = r$fit$lambda[2],
                                          ed = r$fit$ed, max_rel_diff = r$rel, at_bound = isTRUE(r$fit$at_bound))
    }
  }
  list(obs = grid_table(rbindlist(out), g), diag = rbindlist(dg, fill = TRUE))
}

# One stratum: a block of periods (one for 1D). e_list: grid exposures per period.
ungroup_rate_block <- function(o, e_list, fa, B, pen, ctl, omega, n_open, ps = NULL) {
  ps <- ps %||% o$p1[1]
  nf <- length(fa$lower)
  ng <- max(fa$grid)
  # composite link: rows = observed cells; columns = fine cells x periods
  rows <- lapply(seq_len(nrow(o)), function(i) {
    cols_a <- which(fa$lower >= o$age_lo[i] - 1e-9 & fa$upper <= pmin(o$age_hi[i], omega) + 1e-9)
    cols_p <- match(seq(o$p1[i], o$p2[i]), ps)
    as.vector(outer(cols_a, (cols_p - 1L) * nf, `+`))
  })
  C <- Matrix::sparseMatrix(i = rep(seq_along(rows), lengths(rows)), j = unlist(rows), x = 1,
                            dims = c(nrow(o), nf * length(ps)))
  e_grid <- do.call(cbind, e_list)
  if (anyNA(e_grid)) mt_abort("Exposures are missing for grid cells needed by the ungrouping.", class = "ungroup")
  open_idx <- which(fa$grid == ng & fa$open)
  w <- matrix(1, nf, length(ps))
  init_rate <- sum(o$value) / sum(e_grid)
  if (length(open_idx)) {
    w[open_idx, ] <- open_weights(init_rate * exp(0.1 * seq_along(open_idx)))
  }
  fit <- NULL
  for (outer_it in seq_len(if (length(open_idx)) 3L else 1L)) {
    off <- as.vector(e_grid[fa$grid, , drop = FALSE] * w)
    cov <- Matrix::colSums(C) > 0 & off > 0
    fit <- pclm_fit(o$value, C[, cov, drop = FALSE], off[cov], B[cov, , drop = FALSE], pen, ctl$lambda,
                    ctl$criterion, ctl$lambda_range, beta = fit$beta)
    rate <- rep(NA_real_, length(off))
    rate[cov] <- fit$gamma
    if (!length(open_idx)) break
    rate_m <- matrix(exp(as.numeric(B %*% fit$beta)), nf, length(ps))
    for (j in seq_along(ps)) w[open_idx, j] <- open_weights(rate_m[open_idx, j])
  }
  off <- as.vector(e_grid[fa$grid, , drop = FALSE] * w)
  cov <- Matrix::colSums(C) > 0 & off > 0
  deaths_f <- rep(0, length(off))
  deaths_f[cov] <- off[cov] * fit$gamma
  rel <- max(abs(fit$mu / pmax(o$value, 1e-12) - 1)[o$value > 0], 0)
  if (ctl$preserve_totals) deaths_f[cov] <- rescale_to_groups(deaths_f[cov], C[, cov, drop = FALSE], o$value)
  dm <- matrix(deaths_f, nf, length(ps))
  dgrid <- rowsum(dm, fa$grid)
  covered <- rowsum(matrix(as.numeric(Matrix::colSums(C) > 0), nf, length(ps)), fa$grid) > 0
  a <- rep(seq_len(ng), length(ps))
  p <- rep(ps, each = ng)
  keep <- as.vector(covered)
  prov <- rep(ctl$provenance %||% "pclm", length(a))
  vals <- as.vector(dgrid)
  if (!isFALSE(ctl$keep_observed)) {
    exact <- o[a1 == a2 & p1 == p2]
    m <- match(paste(exact$a1, exact$p1), paste(a, p))
    vals[m] <- exact$value
    prov[m] <- "observed"
  }
  list(a = a[keep], p = p[keep], deaths = vals[keep], provenance = prov[keep], fit = fit, rel = rel)
}

rescale_to_groups <- function(est, C, y) {
  C <- as.matrix(C)
  fitted <- as.numeric(C %*% est)
  f <- ifelse(fitted > 0, y / fitted, 0)
  grp <- apply(C, 2, function(col) if (any(col > 0)) which(col > 0)[1] else NA_integer_)
  est * ifelse(is.na(grp), 1, f[grp])
}

monotone_split <- function(lo, hi, value, fine_lo, fine_hi) {
  ord <- order(lo)
  lo <- lo[ord]
  hi <- hi[ord]
  value <- value[ord]
  if (any(abs(hi[-length(hi)] - lo[-1]) > 1e-9)) {
    mt_abort("Monotone ungrouping needs contiguous groups.", class = "ungroup")
  }
  xs <- c(lo, hi[length(hi)])
  cs <- c(0, cumsum(value))
  f <- stats::splinefun(xs, cs, method = "hyman")
  pmax(diff(f(c(fine_lo, fine_hi[length(fine_hi)]))), 0)
}

obs_to_matrix <- function(obs, g) {
  idx <- obs_grid_index(obs, g)
  ok <- !is.na(idx$a1) & idx$a1 == idx$a2 & idx$p1 == idx$p2
  arr <- array(NA_real_, c(length(g$age), length(g$period), max(obs$stratum)))
  arr[cbind(idx$a1[ok], idx$p1[ok], obs$stratum[ok])] <- obs$value[ok]
  arr
}

grid_table <- function(x, g) {
  out <- x[, list(stratum, age_lo = g$age$lower[a], age_hi = g$age$upper[a], period_lo = g$period$lower[p],
                  period_hi = g$period$upper[p], value, provenance)]
  setorderv(out, c("stratum", "period_lo", "age_lo"))
  out[]
}
