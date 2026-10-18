#' Excess mortality
#'
#' Excess deaths are the observed deaths minus the deaths expected without
#' the shock: the predictive deaths of a baseline model fitted to the years
#' before it, given the observed exposures of the later years. The WHO
#' estimates of the excess mortality of the COVID-19 pandemic follow this
#' approach with their own baseline model (Msemburi et al. 2023); Acosta
#' and Irizarry (2022) describe a model of the excess itself. Here the expected deaths are draws of
#' [mt_forecast()] with `predictive = TRUE`: the uncertainty of the
#' baseline trend plus the Poisson and extra-Poisson variation of the death
#' counts. The observed deaths minus each draw give draws of the excess.
#' When the baseline model holds, the interval of the excess covers the true
#' excess with its nominal probability, and an interval that contains zero
#' means deaths within the range of normal variation.
#'
#' The result gives, for each row, the observed deaths and the mean and
#' quantiles of
#'
#' * `expected`: the expected deaths;
#' * `excess`: observed minus expected deaths;
#' * `p_score`: the excess in percent of the expected deaths,
#'   `100 * excess / expected`, computed draw by draw;
#'
#' and `p_positive`, the share of draws with a positive excess.
#'
#' Sums over ages (`ages`), strata (`by`) and years (`cumulative`) are formed
#' draw by draw before the quantiles are taken. With `quantity =
#' "life_expectancy"`, the excess is the observed life expectancy at age `at`
#' minus that of each draw of expected deaths, so a loss of life expectancy
#' is negative; strata are combined by summing their deaths and exposures.
#'
#' The deaths are compared at their observed resolution: the draws on the
#' grid of the model are summed to each observed cell, so a model fitted on
#' single years of age (for example with `prep = mt_prep(step_ungroup())`)
#' can be compared with deaths in 5-year groups. As in [mt_backtest()],
#' cells whose deaths were estimated (provenance `"estimated"`, for example
#' suppressed cells filled by [mt_fill_suppressed()]) and missing cells are
#' left out unless `include_estimated = TRUE`; columns `cells` and
#' `skipped` count the cells used and left out. For life expectancy, years
#' with a cell left out are skipped with a warning.
#'
#' The baseline should end with the last year before the shock and leave
#' out years of other shocks (wars, earlier epidemics). Models need
#' consecutive years, so the baseline is one span of years. The expected
#' deaths of a long gap between the baseline and the shock (for example
#' 1914-1917 for 1918 with a baseline 1900-1913) are forecasts over the
#' whole gap, with correspondingly wide intervals.
#'
#' `mt_plot_excess()` draws the result: observed and expected values over
#' the years (`type = "expected"`), or the excess or P-score with its
#' interval by age group (when there are several) or year.
#'
#' @param model An [mt_fit] of the baseline model, fitted to the years
#'   before the shock, or a model specification, which is then fitted to
#'   the `baseline` years of `data`.
#' @param data An [mt_data] object with the observed deaths and exposures of
#'   the years of interest (it may also contain the baseline years).
#' @param baseline For a model specification: the first and last year of
#'   the baseline (`c(1990, 2002)`).
#' @param periods Years for which the excess is estimated (default: every
#'   year of `data` after the baseline).
#' @param quantity `"deaths"` or `"life_expectancy"`.
#' @param by Stratum variables kept in the result (default all, one row per
#'   stratum); the deaths of strata with the same values of `by` are summed.
#'   `character()` sums all strata. Strata at the total or group levels of
#'   a hierarchy (see [mt_aggregate_strata()]) are left out of the sums
#'   over their variable.
#' @param ages Age groups of the result: `NULL` sums over all ages, a
#'   numeric vector gives the start of each group (`c(0, 65, 85)`: 0-64,
#'   65-84 and 85+; ages below the first value are left out), and
#'   `"observed"` keeps the age groups of the observed deaths. Observed
#'   cells must lie within the groups. Deaths only.
#' @param cumulative Sum over the years of `periods` up to each year (deaths
#'   only). Column `from` gives the first year of the sum.
#' @param at Age for life expectancy (default the first age).
#' @param n_draws Number of draws of the expected deaths.
#' @param seed Random seed of the forecast draws (`NULL` leaves the random
#'   number generator as is).
#' @param probs Quantile levels.
#' @param include_estimated Use cells whose deaths were estimated.
#' @param prep Preprocessing pipeline ([mt_prep()]) for a model
#'   specification.
#' @param ... Passed to the life-table calculation for life expectancy
#'   (`ax`, `open`, `k`).
#' @param x An `mt_excess()` result.
#' @param type `"expected"` (observed and expected values), `"excess"` or
#'   `"p_score"`.
#' @param level Level of the interval shown (its quantiles must be in `x`).
#' @param colour Stratum variable that gives the colours (default: the one
#'   with the most levels).
#' @return `mt_excess()`: a `data.table` with the `by` variables, `age`,
#'   `age_lo` and `age_hi` (deaths), `period`, `observed`, the mean and
#'   quantiles of `expected` and `excess` (`expected_mean`,
#'   `expected_q0.025`, ...), for deaths those of `p_score`, then
#'   `p_positive`, `cells` and `skipped`; for life expectancy also `at`.
#'   `mt_plot_excess()`: a ggplot object.
#' @references
#' Acosta RJ, Irizarry RA (2022). A flexible statistical framework for
#' estimating excess mortality. *Epidemiology* 33(3): 346-353.
#' \doi{10.1097/EDE.0000000000001445}
#'
#' Msemburi W, Karlinsky A, Knutson V, Aleshin-Guendel S, Chatterji S,
#' Wakefield J (2023). The WHO estimates of excess mortality associated with
#' the COVID-19 pandemic. *Nature* 613: 130-137.
#' \doi{10.1038/s41586-022-05522-2}
#' @examples
#' x <- mt_example(years = 2000:2019, deaths = "single")
#' # a shock: 10% more deaths at ages 70+ in 2018 and 2019
#' x$obs$deaths[period_lo >= 2018 & age_lo >= 70, value := round(1.1 * value)]
#' ex <- mt_excess(lee_carter(), x, baseline = c(2000, 2015), ages = c(30, 70), n_draws = 200, seed = 1)
#' ex[, list(sex, age, period, observed, excess_q0.5, p_score_q0.025, p_score_q0.975)]
#' if (requireNamespace("ggplot2", quietly = TRUE)) mt_plot_excess(ex)
#' @export
mt_excess <- function(model, data, baseline = NULL, periods = NULL, quantity = c("deaths", "life_expectancy"),
                      by = NULL, ages = NULL, cumulative = FALSE, at = NULL, n_draws = 500, seed = NULL,
                      probs = c(0.025, 0.5, 0.975), include_estimated = FALSE, prep = NULL, ...) {
  quantity <- match.arg(quantity)
  check_class(data, "mt_data", "data")
  check_flag(cumulative, "cumulative")
  check_flag(include_estimated, "include_estimated")
  check_number(n_draws, "n_draws", lower = 2, integerish = TRUE)
  if (!is.numeric(probs) || !length(probs) || anyNA(probs) || any(probs <= 0 | probs >= 1)) {
    mt_abort("`probs` must be numbers between 0 and 1.", class = "argument")
  }
  if (quantity == "life_expectancy" && (cumulative || !is.null(ages))) {
    mt_abort("`ages` and `cumulative` apply to deaths, not to life expectancy.", class = "argument")
  }
  st <- data$strata
  vars <- setdiff(names(st), c("stratum", "label"))
  by <- if (is.null(by)) vars else as.character(by)
  if (!all(by %in% vars)) {
    mt_abort(sprintf("`by` must name stratum variables (%s).", paste(vars, collapse = ", ")), class = "argument")
  }
  fit <- excess_fit(model, data, baseline, prep)
  last <- max(fit$periods)
  yrs <- mt_grid(data)$period$lower
  periods <- sort(unique(periods %||% yrs[yrs > last]))
  if (!length(periods)) {
    mt_abort(sprintf("`data` has no years after the baseline (last baseline year %s).", last), class = "argument")
  }
  if (any(periods <= last)) {
    mt_abort(sprintf("`periods` must be years after the baseline (last baseline year %s).", last), class = "argument")
  }
  miss <- setdiff(periods, yrs)
  if (length(miss)) mt_abort(sprintf("`data` has no deaths for %s.", paste(miss, collapse = ", ")), class = "argument")
  sim <- excess_cells(fit, data, periods, n_draws, seed, include_estimated)
  # strata at the total or group levels of a variable that is summed over
  # would be counted twice
  bottom <- st$stratum[excess_bottom(st, object_hierarchy(data), by, quantity)]
  keep <- sim$cells$stratum %in% bottom
  if (!any(keep)) mt_abort("No strata left to sum.", class = "argument")
  cell_tab <- cbind(st[match(sim$cells$stratum[keep], st$stratum), by, with = FALSE], sim$cells[keep])
  sim$D <- sim$D[keep, , drop = FALSE]
  sim$skipped <- sim$skipped[sim$skipped$stratum %in% bottom]
  left_tab <- cbind(st[match(sim$skipped$stratum, st$stratum), by, with = FALSE], sim$skipped)
  out <- if (quantity == "deaths") {
    excess_deaths(cell_tab, sim$D, left_tab, by, ages, cumulative, probs)
  } else {
    excess_life_expectancy(cell_tab, sim$D, left_tab, by, at, probs, ...)
  }
  setattr(out, "excess", list(quantity = quantity, probs = probs, model = fit$spec$label,
                              baseline = range(fit$periods), by = by, cumulative = cumulative))
  out[]
}

# Strata to sum: those not at the total or a group level of a hierarchy
# variable left out of `by`. Life expectancy sums exposures, which is wrong
# for variables whose levels share exposures (causes of death)
excess_bottom <- function(st, h, by, quantity) {
  keep <- rep(TRUE, nrow(st))
  if (is.null(h)) return(keep)
  for (v in intersect(setdiff(names(h$vars), by), names(st))) {
    hv <- h$vars[[v]]
    if (quantity == "life_expectancy" && identical(hv$additive, "rates")) {
      mt_abort(sprintf("The levels of `%s` share exposures; keep it in `by` for life expectancy.", v),
               class = "argument")
    }
    keep <- keep & !st[[v]] %in% c(hv$total, names(hv$groups))
  }
  keep
}

# The baseline fit: a fit as given, or a specification fitted to the
# baseline years of the data
excess_fit <- function(model, data, baseline, prep) {
  if (inherits(model, "mt_fit")) {
    if (!is.null(prep)) mt_abort("`prep` applies to a model specification; a fit has its own.", class = "argument")
    if (!is.null(baseline) && !isTRUE(all.equal(range(baseline), range(model$periods)))) {
      mt_abort(sprintf("`baseline` differs from the years of the fit (%s-%s); leave it out for a fit.",
                       min(model$periods), max(model$periods)), class = "argument")
    }
    return(model)
  }
  if (!inherits(model, "mt_spec")) mt_abort("`model` must be a fit (mt_fit) or a model specification.",
                                            class = "argument")
  if (!is.numeric(baseline) || !length(baseline) || anyNA(baseline)) {
    mt_abort("Give the `baseline` years for a model specification, for example baseline = c(1990, 2002).",
             class = "argument")
  }
  mt_fit(model, mt_window(data, from = min(baseline), to = max(baseline)), prep = prep)
}

# Observed deaths and draws of expected deaths at the observed cells of
# `periods`: cells (stratum, age_lo, age_hi, period, observed, exposure), D
# (cells x draws) and the cells left out
excess_cells <- function(fit, data, periods, n_draws, seed, include_estimated) {
  last <- max(fit$periods)
  to <- max(periods)
  test <- mt_window(data, from = last + 1, to = to)
  test <- test_on_age_axis(test, fit$age)
  E_obs <- test_exposures(test, fit)
  # exposures of every forecast year; years without data (not reported) take
  # those of the nearest year with data
  fy <- seq(last + 1, to)
  tp <- mt_grid(test)$period$lower
  E <- E_obs[, vapply(fy, function(y) which.min(abs(tp - y)), 1L), , drop = FALSE]
  dimnames(E)[[2]] <- as.character(fy)
  fc <- mt_forecast(fit, to = to, n_draws = n_draws, seed = seed, predictive = TRUE, exposure = E)
  tst <- test$strata
  sf <- match(tst$label, fc$strata$label)
  if (anyNA(sf)) {
    mt_abort(sprintf("The baseline fit lacks strata of `data`: %s.", paste(utils::head(tst$label[is.na(sf)], 5),
                                                                          collapse = ", ")), class = "argument")
  }
  obs <- test$obs$deaths
  obs <- obs[obs$period_lo %in% periods]
  g <- list(age = fit$age, period = single_years(last + 1, to))
  idx <- obs_grid_index(obs, g, "deaths", strict = TRUE)
  inside <- !is.na(idx$a1)
  obs <- obs[inside]
  idx <- idx[inside]
  left_out <- is.na(obs$value) | (!include_estimated & obs$provenance == "estimated")
  na <- length(fit$age)
  rows <- list()
  mats <- list()
  failed <- character()
  for (s in unique(obs$stratum[!left_out])) {
    j <- match(s, tst$stratum)
    Ds <- fc$rates[, , sf[j], , drop = FALSE] * as.vector(E[, , j])
    if (anyNA(Ds)) {
      failed <- c(failed, tst$label[j])
      left_out[obs$stratum == s] <- TRUE
      next
    }
    sel <- which(!left_out & obs$stratum == s)
    C <- build_operator(idx[sel], na, length(fy))
    mats[[length(mats) + 1L]] <- as.matrix(C %*% matrix(Ds, ncol = n_draws))
    rows[[length(rows) + 1L]] <- data.table(stratum = s, age_lo = obs$age_lo[sel], age_hi = obs$age_hi[sel],
                                            period = obs$period_lo[sel], observed = obs$value[sel],
                                            exposure = as.numeric(C %*% as.vector(E[, , j])))
  }
  if (length(failed)) {
    mt_warn(sprintf("The baseline model has no forecast for %d stratum/strata (%s); they are left out.",
                    length(failed), paste(utils::head(failed, 3), collapse = ", ")), class = "excess")
  }
  if (!length(rows)) mt_abort("No observed cells to compare with the baseline.", class = "excess")
  list(cells = rbindlist(rows), D = do.call(rbind, mats),
       skipped = data.table(stratum = obs$stratum[left_out], age_lo = obs$age_lo[left_out],
                            age_hi = obs$age_hi[left_out], period = obs$period_lo[left_out]))
}

# Age group of each cell: list(lo, hi), NA for cells outside the groups
excess_age_groups <- function(lo, hi, ages, all_lo, all_hi) {
  if (is.null(ages)) return(list(lo = rep(min(all_lo), length(lo)), hi = rep(max(all_hi), length(lo))))
  if (identical(ages, "observed")) return(list(lo = lo, hi = hi))
  if (!is.numeric(ages) || !length(ages) || anyNA(ages) || anyDuplicated(ages)) {
    mt_abort("`ages` must be NULL, \"observed\" or the starts of age groups (for example c(0, 65, 85)).",
             class = "argument")
  }
  b <- sort(ages)
  top <- max(all_hi)
  i <- findInterval(lo + 1e-9, b)
  i[i == 0L] <- NA_integer_
  glo <- b[i]
  ghi <- c(b[-1], top)[i]
  bad <- which(!is.na(i) & hi > ghi + 1e-9)
  if (length(bad)) {
    mt_abort(sprintf("The observed age group %s crosses the limits of `ages`; choose limits of the observed groups.",
                     interval_label(lo[bad[1]], hi[bad[1]])), class = "argument")
  }
  list(lo = glo, hi = ghi)
}

# Mean and quantiles of the rows of X, columns <prefix>_mean and
# <prefix>_q<p>; rows with a value that is not finite are NA
path_summary <- function(X, prefix, probs) {
  ok <- rowSums(!is.finite(X)) == 0
  Q <- matrix(NA_real_, nrow(X), length(probs))
  if (any(ok)) Q[ok, ] <- row_quantiles(X[ok, , drop = FALSE], probs)
  mn <- rep(NA_real_, nrow(X))
  mn[ok] <- rowMeans(X[ok, , drop = FALSE])
  out <- data.table(mn)
  setnames(out, paste0(prefix, "_mean"))
  for (k in seq_along(probs)) set(out, j = paste0(prefix, "_", qnames(probs[k])), value = Q[, k])
  out
}

excess_deaths <- function(cell_tab, D, left_tab, by, ages, cumulative, probs) {
  all_lo <- cell_tab$age_lo
  all_hi <- cell_tab$age_hi
  grp <- excess_age_groups(all_lo, all_hi, ages, all_lo, all_hi)
  cell_tab[, `:=`(group_lo = grp$lo, group_hi = grp$hi)]
  keep <- !is.na(cell_tab$group_lo)
  if (!any(keep)) mt_abort("No observed cells in the age groups of `ages`.", class = "argument")
  cell_tab <- cell_tab[keep]
  D <- D[keep, , drop = FALSE]
  keys <- c(by, "group_lo", "group_hi", "period")
  cell_tab[, idx := .GRP, by = keys]
  Dg <- rowsum(D, cell_tab$idx, reorder = TRUE)
  tab <- cell_tab[, list(observed = sum(observed), cells = .N), by = c("idx", keys)]
  setorderv(tab, "idx")
  if (nrow(left_tab)) {
    sg <- excess_age_groups(left_tab$age_lo, left_tab$age_hi, ages, all_lo, all_hi)
    left_tab[, `:=`(group_lo = sg$lo, group_hi = sg$hi)]
    sk <- left_tab[!is.na(left_tab$group_lo), list(skipped = .N), by = keys]
    tab <- merge(tab, sk, by = keys, all.x = TRUE, sort = FALSE)
    tab[is.na(skipped), skipped := 0L]
  } else {
    tab[, skipped := 0L]
  }
  setorderv(tab, c(by, "group_lo", "period"))
  Dg <- Dg[tab$idx, , drop = FALSE]
  if (cumulative) {
    tab[, group_ := .GRP, by = c(by, "group_lo", "group_hi")]
    for (sr in unique(tab$group_)) {
      i <- which(tab$group_ == sr)
      if (length(i) > 1) Dg[i, ] <- apply(Dg[i, , drop = FALSE], 2, cumsum)
    }
    tab[, `:=`(from = min(period), observed = cumsum(observed), cells = cumsum(cells), skipped = cumsum(skipped)),
        by = "group_"]
    tab[, group_ := NULL]
  }
  X <- tab$observed - Dg
  out <- tab[, c(by, "group_lo", "group_hi", if (cumulative) "from", "period", "observed"), with = FALSE]
  setnames(out, c("group_lo", "group_hi"), c("age_lo", "age_hi"))
  out[, age := interval_label(age_lo, age_hi)]
  setcolorder(out, c(by, "age", "age_lo", "age_hi"))
  out <- cbind(out, path_summary(Dg, "expected", probs), path_summary(X, "excess", probs),
               path_summary(100 * X / Dg, "p_score", probs))
  out[, `:=`(p_positive = rowMeans(X > 0), cells = tab$cells, skipped = tab$skipped)]
  out[]
}

excess_life_expectancy <- function(cell_tab, D, left_tab, by, at, probs, ...) {
  keys <- c(by, "period")
  cell_tab[, idx := .GRP, by = c(keys, "age_lo", "age_hi")]
  Dg <- rowsum(D, cell_tab$idx, reorder = TRUE)
  tab <- cell_tab[, list(observed = sum(observed), exposure = sum(exposure)), by = c("idx", keys, "age_lo", "age_hi")]
  setorderv(tab, "idx")
  tab[, group_ := .GRP, by = keys]
  sk <- if (nrow(left_tab)) left_tab[, list(skipped = .N), by = keys] else NULL
  rows <- list()
  lost <- character()
  for (k in unique(tab$group_)) {
    i <- which(tab$group_ == k)
    i <- i[order(tab$age_lo[i])]
    t1 <- tab[i[1]]
    n_skip <- if (is.null(sk)) NA_integer_ else sk[t1[, keys, with = FALSE], on = keys][["skipped"]]
    lo <- tab$age_lo[i]
    hi <- tab$age_hi[i]
    lab <- paste(c(unlist(t1[, by, with = FALSE]), t1$period), collapse = " ")
    if (!is.na(n_skip) || length(i) < 2 || any(abs(hi[-length(hi)] - lo[-1]) > 1e-9)) {
      lost <- c(lost, lab)
      next
    }
    ax <- mt_axis(lo, hi, kind = "age")
    a0 <- at %||% lo[1]
    if (!a0 %in% lo) mt_abort("`at` must be the start of an observed age group.", class = "argument")
    E <- tab$exposure[i]
    e_draws <- life_expectancy(Dg[i, , drop = FALSE] / E, ax, at = a0, ...)
    e_obs <- life_expectancy(tab$observed[i] / E, ax, at = a0, ...)
    X <- matrix(e_obs - e_draws, 1)
    r <- cbind(t1[, keys, with = FALSE], data.table(at = a0, observed = e_obs),
               path_summary(matrix(e_draws, 1), "expected", probs), path_summary(X, "excess", probs))
    r[, `:=`(p_positive = mean(X > 0), cells = length(i), skipped = 0L)]
    rows[[length(rows) + 1L]] <- r
  }
  if (length(lost)) {
    mt_warn(sprintf(paste("Life expectancy needs observed deaths at all ages without gaps; %d stratum-year(s) are",
                          "skipped (%s)."), length(lost), paste(utils::head(lost, 3), collapse = ", ")),
            class = "excess")
  }
  if (!length(rows)) mt_abort("No stratum-year has observed deaths at all ages.", class = "excess")
  out <- rbindlist(rows)
  setorderv(out, keys)
  out[]
}

#' @rdname mt_excess
#' @export
mt_plot_excess <- function(x, type = c("expected", "excess", "p_score"), level = 0.95, colour = NULL) {
  need_ggplot()
  type <- match.arg(type)
  info <- attr(x, "excess")
  if (!is.data.frame(x) || !all(c("period", "observed", "expected_mean", "excess_mean") %in% names(x))) {
    mt_abort("`x` must be a result of mt_excess().", class = "argument")
  }
  d <- copy(as.data.table(x))
  quantity <- info$quantity %||% if ("p_score_mean" %in% names(d)) "deaths" else "life_expectancy"
  if (type == "p_score" && quantity != "deaths") {
    mt_abort("The P-score is defined for deaths, not for life expectancy.", class = "argument")
  }
  prefix <- if (type == "expected") "expected" else type
  lo_col <- paste0(prefix, "_", qnames((1 - level) / 2))
  hi_col <- paste0(prefix, "_", qnames(1 - (1 - level) / 2))
  mid <- paste0(prefix, "_q0.5")
  if (!all(c(lo_col, hi_col) %in% names(d))) {
    mt_abort(sprintf("`x` lacks the quantiles of a %s%% interval; give them as `probs` in mt_excess().",
                     level * 100), class = "argument")
  }
  if (!mid %in% names(d)) mid <- paste0(prefix, "_mean")
  by <- info$by %||% intersect(names(d), setdiff(names(d), c("age", "age_lo", "age_hi", "from", "period", "at")))
  by <- intersect(by, names(d))
  nlev <- vapply(by, function(v) length(unique(d[[v]])), 0)
  if (!is.null(colour) && !colour %in% by) {
    mt_abort(sprintf("`colour` must be a stratum variable of `x` (%s).", paste(by, collapse = ", ")),
             class = "argument")
  }
  cand <- by[nlev > 1]
  colour <- colour %||% if (length(cand)) cand[which.max(nlev[cand])] else NULL
  d[, col := if (is.null(colour)) "all" else as.character(get(colour))]
  levs <- unique(d$col)
  d[, col := factor(col, levels = levs)]
  pvars <- setdiff(cand, colour)
  has_age <- "age" %in% names(d) && length(unique(d$age)) > 1
  if (has_age) d[, age := factor(age, levels = unique(age[order(age_lo)]))]
  np <- length(unique(d$period))
  d[, `:=`(med = get(mid), lo = get(lo_col), hi = get(hi_col))]
  vals <- stats::setNames(mt_palette(length(levs), if (looks_ordered(levs)) "ordered" else "categorical"), levs)
  ylab <- switch(type,
                 expected = if (quantity == "deaths") "Deaths" else sprintf("Life expectancy at age %s", d$at[1]),
                 excess = if (quantity == "deaths") "Excess deaths" else "Excess life expectancy (years)",
                 p_score = "P-score (% of expected deaths)")
  lvl <- paste0(level * 100, "%")
  base <- if (!is.null(info)) sprintf("Baseline %s, %s-%s", info$model, info$baseline[1], info$baseline[2]) else NULL
  by_age <- has_age && (type != "expected" || np == 1)
  panel_vars <- c(pvars, if (has_age && !by_age) "age", if (by_age && np > 1) "period")
  if (length(panel_vars)) {
    d[, panel := do.call(paste, c(lapply(panel_vars, function(v) as.character(get(v))), sep = ", "))]
    d[, panel := factor(panel, levels = unique(panel))]
  }
  dodge <- ggplot2::position_dodge(width = 0.5)
  if (type == "expected" && np > 1) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = period, group = col)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = col), alpha = 0.2) +
      ggplot2::geom_line(ggplot2::aes(y = med, colour = col), linetype = "dashed", linewidth = 0.6) +
      ggplot2::geom_line(ggplot2::aes(y = observed, colour = col), linewidth = 0.5) +
      ggplot2::geom_point(ggplot2::aes(y = observed, colour = col), size = 1.2)
    sub <- sprintf("Points and solid lines: observed; dashed: expected (median) with %s interval", lvl)
    xlab <- NULL
  } else {
    xv <- if (by_age) "age" else "period"
    p <- ggplot2::ggplot(d, ggplot2::aes(x = if (xv == "period") factor(period) else age, colour = col))
    if (type == "expected") {
      p <- p + ggplot2::geom_pointrange(ggplot2::aes(y = med, ymin = lo, ymax = hi), shape = 1, position = dodge) +
        ggplot2::geom_point(ggplot2::aes(y = observed), size = 2, position = dodge)
      sub <- sprintf("Filled points: observed; open points: expected (median) with %s interval", lvl)
    } else {
      p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
        ggplot2::geom_pointrange(ggplot2::aes(y = med, ymin = lo, ymax = hi), position = dodge)
      sub <- sprintf("Median and %s interval%s", lvl,
                     if (isTRUE(info$cumulative)) ", cumulative over the years" else "")
    }
    xlab <- if (xv == "age") "Age" else NULL
  }
  if (length(panel_vars)) p <- p + ggplot2::facet_wrap(ggplot2::vars(panel), scales = "free_y")
  if (np == 1) sub <- paste0(sub, "; ", d$period[1])
  ribbon <- type == "expected" && np > 1
  p <- p + ggplot2::scale_colour_manual(values = vals, breaks = levs, aesthetics = c("colour", if (ribbon) "fill")) +
    ggplot2::labs(x = xlab, y = ylab, colour = colour, title = base, subtitle = sub)
  if (ribbon) p <- p + ggplot2::labs(fill = colour)
  p + theme_mortensor(legend = if (is.null(colour)) "none" else "bottom")
}
