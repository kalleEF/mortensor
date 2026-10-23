#' Pessimism of forecasts by cause of death
#'
#' Compares a forecast of all-cause mortality with the sum of separate
#' forecasts of the causes of death, as Wilmoth (1995) does. Wilmoth calls
#' the first the "total" forecast and the second the "all-cause" forecast;
#' here they are `"total"` and `"causes"`. The sum of cause forecasts tends
#' to be higher, because it is dominated over time by the causes that fall
#' slowest or rise. Every measure is computed path by path:
#'
#' * `by = "total"`: life expectancy at `at` of the total forecast
#'   (`measure = "total"`) and of the sum of the causes (`"causes"`), and
#'   their `"difference"`, total minus causes, the difference
#'   `e0(Tot) - e0(All)` that Wilmoth decomposes (Section 3, equation 15 and
#'   Table 5). It is positive where the forecast by cause is more
#'   pessimistic.
#' * `by = "age"`: the rates of the total forecast (`"total"`) and of the
#'   sum of the causes (`"causes"`), and their ratio, the m-ratio of
#'   Wilmoth's equation (6), `M(x, t) = sum_i m_i(x, t) / m(x, t)`
#'   (`"m_ratio"`), his measure of the pessimism of forecasts by cause. It is
#'   above 1 where the forecast by cause is more pessimistic.
#' * `by = "cause"`: the contribution of each cause to the m-ratio minus 1
#'   (equation 13), `m_i(x, t) / m(x, t) - p_i(x)`, where `p_i(x)` is the
#'   proportion of the deaths at age `x` that are due to cause `i` in the last
#'   observed year (equation 5). The contributions of an age add up to its
#'   m-ratio minus 1 (Table 3).
#'
#' @section Forecasts:
#' For [mt_draws], the forecast of all causes comes from `x` and the
#' forecasts of the causes from `causes`. Without `causes`, `x` holds both:
#' the strata of `var` at `total_level` and those of the causes, for
#' example the forecasts of [mt_fit()] with a separate model for each
#' stratum. With `causes`, for example the draws of a model of causes such
#' as [coda_causes()], `x` holds the total forecast: its strata at
#' `total_level`, or all its strata if `var` is not a stratum variable of
#' `x`. The causes are the levels of `var` other than `total_level` and the
#' groups of a hierarchy stored in the object (see [mt_aggregate_strata()]),
#' and the strata of the total and of the causes are matched on the other
#' stratum variables (for example sex and deprivation quintile). The causes
#' should cover all deaths: a warning reports observed causes that do not
#' add up to the observed total in the last training year.
#'
#' Wilmoth derives equations (6) and (13) for proportional rate of change
#' models, which project the observed rates of the last year `t0` with
#' constant annual factors. Equation (6) is defined for any forecast; for
#' other models, and for models whose forecasts start from fitted rates
#' (his footnote 2), the m-ratio is still the ratio of the forecasts and the
#' contributions still add up to the m-ratio minus 1, but the m-ratio of
#' the first forecast years also contains the jump-off of the models (see
#' [mt_jump_off()]).
#'
#' Draws of two objects come from different models and are independent: by
#' default (`pairing = "independent"`) the draws of `causes` are shuffled
#' before they are paired with those of `x`, as in [mt_decompose()]. Use
#' `pairing = "draw"` when `causes` is derived from the draws of `x`.
#'
#' @section Approximation from two periods:
#' For an [mt_data] object, `mt_pessimism()` computes Wilmoth's
#' approximation from the distributions of deaths by cause in two periods
#' (Section 3, equations 12 and 14). The rates of each cause and period are
#' its deaths divided by the exposures, both summed over the years of the
#' period (Appendix C), and the total rate is the sum of the causes (a
#' warning reports causes that do not add up to the level `total_level` of
#' the data, where it exists). Both
#' are projected from the midpoint `t2` of the second period at the
#' geometric mean rate of change between the midpoints `t1` and `t2`, the
#' "first-to-last" method of his equations (8) to (11) applied to periods:
#' `m(x, t) = m2(x) (m2(x) / m1(x))^((t - t2) / (t2 - t1))`. The m-ratio of
#' these forecasts is equation (12), `sum_i p_i (p_i / q_i)^((t - t2) / (t2 -
#' t1))`, with `q_i` and `p_i` the proportions of deaths by cause in the
#' first and the second period, and the contributions are equation (14).
#' Wilmoth shows (Appendix B) that this m-ratio exceeds 1 at every horizon
#' unless the two distributions are equal. The difference of life
#' expectancies is the quantity that his Table 5 decomposes by age and
#' cause.
#'
#' The default `periods` are the first and the second half of the observed
#' years, the "first-to-second" method that Wilmoth uses for his
#' decompositions; with an odd number of years the middle year is left out,
#' so that the halves have the same length (package-specific, reference
#' unknown). His "first-to-last" method uses the first and the last three
#' years, for example `periods = list(1951:1953, 1988:1990)`. The forecast
#' years are those after the last observed year (`h`, `to`).
#'
#' With `n_draws = 1` (default) the observed deaths are used. With more
#' draws the deaths are Poisson draws, as in [mt_observed()], and the
#' intervals show the sampling variation of the approximation
#' (package-specific, reference unknown). A cause without deaths at an age
#' in both periods contributes nothing there (the dashes of Wilmoth's Table
#' 3). A cause without deaths in the first period and with deaths in the
#' second would grow without bound: its deaths in the first period are set
#' to 0.5, as [rwd()] does with zero counts, with a warning
#' (package-specific, reference unknown).
#'
#' @param x An [mt_draws] object (forecast) or an [mt_data] object (observed
#'   deaths by cause, for the approximation from two periods).
#' @param causes Optional [mt_draws] object with the forecasts of the causes
#'   (see Forecasts).
#' @param var Stratum variable with the causes.
#' @param total_level Level of `var` for all causes.
#' @param by `"total"` (life expectancy), `"age"` (rates and m-ratios) or
#'   `"cause"` (contributions of the causes to the m-ratios).
#' @param at Age at which life expectancy is computed (start of an age
#'   interval; default the first age).
#' @param periods For data: a list of two vectors of observed years, the
#'   first and the second period (default: the two halves of the observed
#'   years).
#' @param h,to For data: the number of forecast years after the last
#'   observed year, or the last forecast year.
#' @param pairing `"independent"` or `"draw"`: how the draws of `x` and
#'   `causes` are paired.
#' @param n_draws For data: 1 (the observed deaths) or the number of Poisson
#'   draws.
#' @param probs Quantile levels.
#' @param paths Return the value of every path (columns `draw` and `value`)
#'   instead of the mean and quantiles.
#' @param seed Seed of the shuffling of `causes` and of Poisson draws.
#' @param ... Passed to the life table (`ax`, `open`, `k`, `omega`, see
#'   [life_table()]) for life expectancy.
#' @return A `data.table` with the stratum variables other than `var`,
#'   `period`, for `by = "age"` and `"cause"` also `age_lo` and `age`, for
#'   `"cause"` the column `var` with the cause, then `measure`, for `"total"`
#'   the age `at`, and the `mean` and quantiles over paths (`q0.025`, ...).
#'   With `paths = TRUE`, `draw` and `value` replace the mean and quantiles.
#' @references
#' Wilmoth JR (1995). Are mortality projections always more pessimistic
#' when disaggregated by cause of death? *Mathematical Population Studies*
#' 5(4): 293-319. \doi{10.1080/08898489509525409}
#' @seealso [mt_coherence()] for the ratio of summed strata to aggregates in
#'   any hierarchy, [mt_decompose()] for decompositions of life expectancy
#'   by age and cause, [reconciled()] and [coda_causes()] for forecasts by
#'   cause that add up.
#' @examples
#' x <- mt_example_causes(years = 2000:2019)
#' # separate Lee-Carter models for all causes and for each cause
#' fc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 50, seed = 1)
#' mt_pessimism(fc, at = 65)[period == 2029]
#' m <- mt_pessimism(fc, by = "age")
#' m[measure == "m_ratio" & period == 2029 & age %in% c("50", "70", "90")]
#' # Wilmoth's approximation from the halves 2000-2009 and 2010-2019
#' a <- mt_pessimism(x, by = "cause", to = 2039)
#' a[period == 2039 & age == "70"]
#' @export
mt_pessimism <- function(x, causes = NULL, var = "cause", total_level = "all", by = c("total", "age", "cause"),
                         at = NULL, periods = NULL, h = 10, to = NULL, pairing = c("independent", "draw"),
                         n_draws = 1, probs = c(0.025, 0.5, 0.975), paths = FALSE, seed = 1, ...) {
  by <- match.arg(by)
  pairing <- match.arg(pairing)
  if (!is.character(var) || length(var) != 1) mt_abort("`var` must be the name of one stratum variable.",
                                                       class = "argument")
  if (!is.character(total_level) || length(total_level) != 1) {
    mt_abort("`total_level` must be one level of `var`.", class = "argument")
  }
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) || any(probs < 0 | probs > 1)) {
    mt_abort("`probs` must be numbers between 0 and 1.", class = "argument")
  }
  check_flag(paths, "paths")
  src <- if (inherits(x, "mt_data")) {
    if (!is.null(causes)) mt_abort("`causes` takes forecasts; for data, `x` holds the causes.", class = "argument")
    if (!is.null(to) && !missing(h)) mt_warn("Both `h` and `to` are given; `to` is used.", class = "argument")
    pess_from_data(x, var, total_level, periods, h, to, n_draws, seed)
  } else {
    check_class(x, "mt_draws", "x")
    pess_from_draws(x, causes, var, total_level, pairing, seed, by == "cause")
  }
  pess_measures(src, by, at, probs, paths, ...)
}

pess_key <- function(st, vars) {
  if (!length(vars)) return(rep("", nrow(st)))
  do.call(paste, c(lapply(vars, function(v) as.character(st[[v]])), sep = "\r"))
}

# Positions of the strata of the total and of the causes in each group (the
# combinations of the stratum variables other than `var`); causes are the
# levels other than the total and the groups of a stored hierarchy
pess_strata <- function(obj, var, total_level, name) {
  st <- obj$strata
  if (!var %in% names(st)) {
    mt_abort(sprintf("`%s` is not a stratum variable of `%s` (%s).", var, name,
                     paste(setdiff(names(st), c("stratum", "label")), collapse = ", ")), class = "argument")
  }
  gv <- setdiff(names(st), c("stratum", "label", var))
  hv <- object_hierarchy(obj)$vars[[var]]
  lev <- as.character(st[[var]])
  cl <- setdiff(unique(lev), c(total_level, names(hv$groups)))
  if (!length(cl)) mt_abort(sprintf("`%s` has no cause levels besides \"%s\".", var, total_level), class = "argument")
  key <- pess_key(st, gv)
  gk <- unique(key)
  first <- match(gk, key)
  ci <- matrix(NA_integer_, length(gk), length(cl))
  for (g in seq_along(gk)) {
    for (k in seq_along(cl)) {
      i <- which(key == gk[g] & lev == cl[k])
      if (length(i) != 1) {
        mt_abort(sprintf("Stratum %s of `%s` lacks the cause \"%s\".", gsub("\r", ".", gk[g]), name, cl[k]),
                 class = "argument")
      }
      ci[g, k] <- i
    }
  }
  ti <- vapply(gk, function(k) {
    i <- which(key == k & lev == total_level)
    if (length(i) == 1) i else NA_integer_
  }, 1L, USE.NAMES = FALSE)
  list(gv = gv, groups = st[first, gv, with = FALSE], key = gk, total = ti, causes = ci, cause_levels = cl)
}

# Total and cause forecasts [age, period, group, (cause,) draw] of draws, and
# the shares of the causes in the last observed year [age, group, cause]
pess_from_draws <- function(x, causes, var, total_level, pairing, seed, need_shares) {
  cs <- causes %||% x
  if (!is.null(causes)) check_class(causes, "mt_draws", "causes")
  sc <- pess_strata(cs, var, total_level, if (is.null(causes)) "x" else "causes")
  if (is.null(causes)) {
    ti <- sc$total
    if (anyNA(ti)) {
      mt_abort(sprintf(paste("`x` has no stratum with %s = \"%s\" for %s. Give the forecast of all causes in `x`",
                             "and those of the causes in `causes`, or set `total_level`."), var, total_level,
                       gsub("\r", ".", sc$key[which(is.na(ti))[1]])), class = "argument")
    }
  } else {
    if (!identical(x$age$labels, cs$age$labels)) {
      mt_abort("`x` and `causes` must have the same age grid.", class = "argument")
    }
    st <- x$strata
    keep <- if (var %in% names(st)) which(st[[var]] == total_level) else seq_len(nrow(st))
    if (!setequal(setdiff(names(st), c("stratum", "label", var)), sc$gv)) {
      mt_abort(sprintf("`x` and `causes` need the same stratum variables besides `%s`.", var), class = "argument")
    }
    kx <- pess_key(st, sc$gv)
    ti <- vapply(sc$key, function(k) {
      i <- keep[kx[keep] == k]
      if (length(i) == 1) i else NA_integer_
    }, 1L, USE.NAMES = FALSE)
    if (anyNA(ti)) {
      mt_abort(sprintf("`x` has no forecast of all causes for %s.", gsub("\r", ".", sc$key[which(is.na(ti))[1]])),
               class = "argument")
    }
  }
  periods <- intersect(x$periods, cs$periods)
  if (!length(periods)) mt_abort("`x` and `causes` have no forecast year in common.", class = "argument")
  nd <- dim(x$rates)[4]
  if (dim(cs$rates)[4] != nd) mt_abort("`x` and `causes` need the same number of draws.", class = "argument")
  ord <- seq_len(nd)
  if (!is.null(causes) && pairing == "independent" && nd > 1) ord <- with_seed(seed, sample.int(nd))
  px <- match(periods, x$periods)
  pc <- match(periods, cs$periods)
  na <- length(x$age)
  np <- length(periods)
  ng <- length(ti)
  nc <- length(sc$cause_levels)
  tot <- array(x$rates[, px, ti, , drop = FALSE], c(na, np, ng, nd))
  cau <- array(NA_real_, c(na, np, ng, nc, nd))
  for (k in seq_len(nc)) cau[, , , k, ] <- cs$rates[, pc, sc$causes[, k], ord, drop = FALSE]
  # observed rates of the last training year: the causes should add up to the total
  obs_c <- pess_last_observed(cs)
  obs_t <- pess_last_observed(x)
  if (!is.null(obs_c) && !is.null(obs_t) && identical(obs_c$period, obs_t$period)) {
    pess_check_sum(obs_t$rates[, ti, drop = FALSE],
                   matrix(rowSums(array(obs_c$rates[, as.vector(sc$causes), drop = FALSE], c(na, ng, nc)),
                                  dims = 2), na),
                   sprintf("in the last observed year (%s)", obs_t$period))
  }
  shares <- NULL
  if (need_shares) {
    if (is.null(obs_c)) {
      mt_abort("`by = \"cause\"` needs the observed rates of the training years, which the draws of the causes lack.",
               class = "pessimism")
    }
    s <- array(obs_c$rates[, as.vector(sc$causes), drop = FALSE], c(na, ng, nc))
    s <- s / as.vector(rowSums(s, dims = 2))
    shares <- array(s, c(na, ng, nc, nd))
  }
  list(age = x$age, periods = periods, groups = sc$groups, total = tot, causes = cau, cause_levels = sc$cause_levels,
       shares = shares, var = var)
}

# Observed rates [age, stratum] of the last training year held in draws
pess_last_observed <- function(d) {
  o <- d$observed
  if (is.null(o) || length(dim(o)) != 3) return(NULL)
  nt <- dim(o)[2]
  p <- if (!is.null(d$observed_periods)) d$observed_periods[nt] else dimnames(o)[[2]][nt]
  list(rates = matrix(o[, nt, ], dim(o)[1]), period = as.character(p))
}

# Warn where the causes do not add up to the total (matrices [age, group])
pess_check_sum <- function(tot, sum_causes, when) {
  ok <- is.finite(tot) & is.finite(sum_causes) & tot > 0
  if (!any(ok)) return(invisible())
  gap <- max(abs(sum_causes[ok] / tot[ok] - 1))
  if (gap > 1e-6) {
    mt_warn(sprintf(paste("The causes do not add up to all causes %s: by up to %.3g (relative). The m-ratio then",
                          "also measures the part of all causes that the causes leave out; add a residual cause",
                          "(see mt_group_levels(other = ))."), when, gap), class = "pessimism")
  }
  invisible()
}

# Wilmoth's approximation: forecasts of the total and the causes from their
# rates in two periods (first-to-last method applied to periods, equations
# 8-12), and the shares of the causes in the second period
pess_from_data <- function(x, var, total_level, periods, h, to, n_draws, seed) {
  check_number(n_draws, "n_draws", lower = 1, integerish = TRUE)
  obs <- mt_observed(x, n_draws = n_draws, seed = seed)
  obs$level_totals <- x$level_totals
  obs$hierarchy <- x$hierarchy
  sc <- pess_strata(obs, var, total_level, "x")
  yrs <- obs$periods
  if (is.null(periods)) {
    k <- floor(length(yrs) / 2)
    if (k < 1) mt_abort("The approximation needs at least two observed years.", class = "argument")
    periods <- list(yrs[seq_len(k)], yrs[length(yrs) - k + seq_len(k)])
  }
  if (!is.list(periods) || length(periods) != 2 || !all(vapply(periods, is.numeric, TRUE)) ||
      any(lengths(periods) == 0) || anyNA(unlist(periods))) {
    mt_abort("`periods` must be a list of two vectors of observed years.", class = "argument")
  }
  miss <- setdiff(unlist(periods), yrs)
  if (length(miss)) mt_abort(sprintf("Years of `periods` not observed: %s.", paste(miss, collapse = ", ")),
                             class = "argument")
  if (max(periods[[1]]) >= min(periods[[2]])) {
    mt_abort("The first period of `periods` must end before the second begins.", class = "argument")
  }
  last <- max(yrs)
  if (!is.null(to)) h <- to - last
  check_number(h, "h", lower = 1, integerish = TRUE)
  fy <- last + seq_len(h)
  t1 <- mean(periods[[1]])
  t2 <- mean(periods[[2]])
  E <- obs$exposure
  D <- obs$rates * as.vector(E)
  na <- length(obs$age)
  ng <- length(sc$key)
  nc <- length(sc$cause_levels)
  nd <- dim(D)[4]
  if (!anyNA(sc$total)) {
    sum_c <- rowSums(array(obs$train_D[, , as.vector(sc$causes), drop = FALSE], c(na * dim(D)[2], ng, nc)), dims = 2)
    pess_check_sum(matrix(obs$train_D[, , sc$total, drop = FALSE], na * dim(D)[2]), sum_c, "in the observed data")
  }
  ids <- as.vector(sc$causes)
  period_sum <- function(yy) {
    Dk <- array(0, c(na, ng, nc, nd))
    Ek <- array(0, c(na, ng, nc))
    for (i in match(yy, yrs)) {
      Dk <- Dk + array(D[, i, ids, , drop = FALSE], dim(Dk))
      Ek <- Ek + array(E[, i, ids, drop = FALSE], dim(Ek))
    }
    list(D = Dk, E = Ek)
  }
  p1 <- period_sum(periods[[1]])
  p2 <- period_sum(periods[[2]])
  zero <- p1$D == 0 & p2$D > 0
  if (any(zero)) {
    mt_warn(sprintf(paste("%d cells of age, cause and draw have no deaths in the first period and some in the",
                          "second: their first-period deaths are set to 0.5."), sum(zero)), class = "pessimism")
    p1$D[zero] <- 0.5
  }
  m1 <- p1$D / as.vector(p1$E)
  m2 <- p2$D / as.vector(p2$E)
  r <- m2 / m1
  r[m2 == 0] <- 1
  m1t <- rowSums(aperm(m1, c(1, 2, 4, 3)), dims = 3)
  m2t <- rowSums(aperm(m2, c(1, 2, 4, 3)), dims = 3)
  rt <- m2t / m1t
  rt[m2t == 0] <- 1
  tot <- array(NA_real_, c(na, h, ng, nd))
  cau <- array(NA_real_, c(na, h, ng, nc, nd))
  for (j in seq_len(h)) {
    g <- (fy[j] - t2) / (t2 - t1)
    tot[, j, , ] <- m2t * rt^g
    cau[, j, , , ] <- m2 * r^g
  }
  shares <- m2 / as.vector(aperm(array(m2t, c(na, ng, nd, nc)), c(1, 2, 4, 3)))
  list(age = obs$age, periods = fy, groups = sc$groups, total = tot, causes = cau, cause_levels = sc$cause_levels,
       shares = shares, var = var)
}

# The measures of `by` from the total [age, period, group, draw] and cause
# [age, period, group, cause, draw] forecasts, path by path
pess_measures <- function(src, by, at, probs, paths, ...) {
  tot <- src$total
  d <- dim(src$causes)
  na <- d[1]
  np <- d[2]
  ng <- d[3]
  nc <- d[4]
  nd <- d[5]
  summed <- array(0, c(na, np, ng, nd))
  for (k in seq_len(nc)) summed <- summed + array(src$causes[, , , k, , drop = FALSE], c(na, np, ng, nd))
  age <- src$age
  gi <- function(each, times = 1) rep(rep(seq_len(ng), each = each), times)
  if (by == "total") {
    e_tot <- life_expectancy(tot, age, at = at, ...)
    e_sum <- life_expectancy(summed, age, at = at, ...)
    vals <- list(total = e_tot, causes = e_sum, difference = e_tot - e_sum)
    key <- data.table(.g = gi(np), period = rep(src$periods, ng))
  } else if (by == "age") {
    vals <- list(total = tot, causes = summed, m_ratio = summed / tot)
    key <- data.table(.g = gi(na * np), period = rep(rep(src$periods, each = na), ng),
                      age_lo = rep(age$lower, np * ng), age = rep(age$labels, np * ng))
  } else {
    # [age, cause, period, group, draw]: m_i / m - p_i (equation 13)
    rel <- aperm(src$causes, c(1, 4, 2, 3, 5)) / as.vector(aperm(array(tot, c(na, np, ng, nd, nc)),
                                                                   c(1, 5, 2, 3, 4)))
    sh <- aperm(array(src$shares, c(na, ng, nc, nd, np)), c(1, 3, 5, 2, 4))
    vals <- list(contribution = rel - sh)
    key <- data.table(.g = gi(na * nc * np), period = rep(rep(src$periods, each = na * nc), ng),
                      age_lo = rep(age$lower, nc * np * ng), age = rep(age$labels, nc * np * ng),
                      cause = rep(rep(src$cause_levels, each = na), np * ng))
    setnames(key, "cause", src$var)
  }
  n <- nrow(key)
  blocks <- lapply(seq_along(vals), function(i) {
    M <- matrix(vals[[i]], n)
    b <- cbind(key, .m = i, measure = names(vals)[i])
    if (paths) {
      b <- b[rep(seq_len(n), nd)]
      b[, `:=`(draw = rep(seq_len(nd), each = n), value = as.vector(M))]
      return(b)
    }
    set(b, j = "mean", value = rowMeans(M))
    ok <- rowSums(is.na(M)) == 0
    Q <- matrix(NA_real_, n, length(probs))
    if (any(ok)) Q[ok, ] <- row_quantiles(M[ok, , drop = FALSE], probs)
    for (k in seq_along(probs)) set(b, j = qnames(probs)[k], value = Q[, k])
    b
  })
  out <- rbindlist(blocks)
  if (by == "total") set(out, j = "at", value = at %||% age$lower[1])
  setcolorder(out, c(".g", "period", setdiff(names(key), c(".g", "period")), ".m", "measure"))
  setorderv(out, c(".g", "period", ".m", if (paths) "draw"))
  if (by == "total") setcolorder(out, c(setdiff(names(out), c("at", "mean", "draw", "value", qnames(probs))),
                                        "at"))
  if (length(src$groups)) out <- cbind(src$groups[out$.g], out)
  out[, c(".g", ".m") := NULL]
  out[]
}
