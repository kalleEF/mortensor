# Compositional data (CoDa) models of life-table deaths: coda() and
# coda_causes(). Written from Bergeron-Boucher, Canudas-Romo, Oeppen and
# Vaupel (2017) and Oeppen (2008); the numerical core is in
# model-coda-core.R. Earlier versions adapted MortalityForecast 0.8.1 (GPL-3,
# Pascariu and Bergeron-Boucher); its random walk with drift for the index
# and AR(1) deviations remain available through `order` and `order_dev`.

#' Compositional data (CoDa) model of life-table deaths
#'
#' `coda()` is the Lee-Carter model in compositional form (Oeppen 2008) as
#' published by Bergeron-Boucher, Canudas-Romo, Oeppen and Vaupel (2017). It
#' models the life-table distribution of deaths by age, which sums to 1. With
#' `coherent = TRUE` or `common_by`, it is their coherent model for several
#' populations (CoDa-coherent). Status: published, with the deviations listed
#' below. Page and equation numbers refer to Bergeron-Boucher et al. (2017).
#'
#' **Model (Section 2.3, Appendix A).** For one stratum:
#'
#' 1. the death rates of each year give life-table deaths `d[x, t]` with
#'    radix 1 at the first age ([life_table()] conventions);
#' 2. `a[x]` is the geometric mean of `d[x, t]` over the years, closed to 1;
#' 3. the deaths are centred by perturbation (`d[x, t] / a[x]`, closed) and
#'    transformed by the centred log-ratio (clr), eq. (5);
#' 4. the first component of a singular value decomposition gives
#'    `clr(d[x, t] / a[x]) = b[x] k[t]`, eq. (6); `b` and `k` sum to zero, and
#'    their signs make `k` increase over the fitting period (p. 533);
#' 5. `k[t]` follows an ARIMA(0,1,1) model with drift, fitted by maximum
#'    likelihood (`order`, `drift`; p. 536 and Appendix A, p. 555, where this
#'    order had the best BIC for most of 15 countries);
#' 6. forecasts are transformed back (inverse clr, perturbation by `a`,
#'    closure, eq. 7) to life-table deaths and then to death rates.
#'
#' **Jump-off (Appendix G, p. 566).** With `jump_off = "observed"` (the
#' default), forecasts are the life-table deaths of the last year `T`
#' perturbed by the change of the model (eq. 17): `d[x, T + j]` is
#' proportional to `d[x, T] exp(b[x] (k[T + j] - k[T]))`. With `"fitted"`,
#' forecasts start from the fitted deaths of year `T`. The point forecast of
#' the ARIMA(0,1,1) index with drift is shifted so that its first step equals
#' the drift (eq. 20), which removes the break that the moving-average term
#' causes in the first forecast year; the package applies the same
#' correction to any ARIMA(0,1,q) index with drift.
#'
#' **Coherent model (Section 2.4, Appendix C).** The common component
#' `B[x] K[t]` is the model above fitted to the average mortality of a group:
#' the arithmetic mean of the death rates of its strata (equal weights, p.
#' 557), turned into life-table deaths. Each stratum `i` has its own centre
#' `a[x, i]` and deviation `clr(d[x, t, i] / a[x, i] / C[exp(B K)]) = b[x, i]
#' k[t, i]` (eq. 8a), the first component of the SVD after removing the
#' common component. `K` follows the model of `order` and `drift` (with eq.
#' 20), each `k[t, i]` an ARIMA(1,1,0) model without drift (`order_dev`,
#' `drift_dev`; p. 537), whose point forecast levels off, so that the point
#' forecasts of the strata do not diverge; their intervals still widen with
#' the horizon. `order_dev = c(1, 0, 0)` gives an AR(1) process with a mean,
#' which Li and Lee (2005) suggest for their deviations (as cited on p. 531)
#' and whose deviations stay bounded in distribution. Groups are all strata,
#' or each level of `common_by`; if the data carry a hierarchy (see
#' [mt_aggregate_strata()]), the average uses the bottom strata only, as
#' [li_lee()] does.
#'
#' **Draws (Appendix B, pp. 556-557).** Draws carry the uncertainty of the
#' estimates and of the index forecasts: the residuals of eq. (5) are
#' resampled cell by cell and added to the fitted component, `b` and `k` are
#' estimated again (`a` is kept, as in the paper), the ARIMA model is fitted
#' again to the new index, and future paths of the index are simulated with
#' innovations resampled from its residuals. `ceiling(sqrt(n_draws))` refits
#' share the draws (the paper uses 100 refits with 100 paths each). In the
#' coherent model, the common and the stratum components are simulated
#' independently and added (Appendix B); the draws of the common component
#' are the same for all strata of a group. The draws are then centred on the
#' point forecast: each draw is perturbed by the ratio of the point forecast
#' of the first forecast year to the median of the draws in that year, part
#' by part. This follows the authors' code (CoDaForecast.R, supplementary
#' material of the paper); eq. (19) of the paper shifts the bounds of the
#' intervals of life expectancy so that their median is the point forecast.
#' Draws are expected rates; predictive draws (`mt_forecast(predictive =
#' TRUE)`) add Poisson variation of the deaths given the rates and the
#' exposures, and no other observation noise. The model weights all parts
#' equally on the clr scale, so the noise of ages with few deaths enters `k`;
#' its intervals can be wider than those of [lee_carter()] (p. 543).
#'
#' **Old ages (Appendix C, pp. 557-558).** With single-year age groups from
#' 80 (at least to 89) and `open = "auto"` (the default) or `"kannisto"`, the
#' rates from age 80 are replaced in every year by a Kannisto model (Thatcher,
#' Kannisto and Vaupel 1998) fitted by Poisson maximum likelihood to the
#' deaths and exposures of the single ages 80 to 100 of the grid (the
#' authors' code fits ages 80-100 of the Human Mortality Database, open
#' interval 110+) and applied up to age 120, where the life table closes; the
#' single ages 80-120 are parts of the composition. The rate of an open age interval is the forecast deaths over
#' the forecast person-years of its parts, and the forecasts of the ages from
#' 80 are those of the smoothed rates. Other age groups (for example 5-year
#' groups) do not allow this extension. Then, with `open = "auto"` or
#' `"regression"`, the deaths of an open interval are the survivors at its
#' lower bound, and its log death rate gets loadings on the same indexes,
#' estimated by least squares: `log m[open, t] = a_o + b_o k[t]` (coherent
#' model: `a_o[i] + B_o K[t] + b_o[i] k[t, i]`, `B_o` from the log of the
#' average rate of the group); a closed age range adds the survivors at its
#' upper bound as a part. Both are package-specific (reference unknown).
#'
#' **Zero deaths (Appendix C, p. 558).** Zero counts of the ages that keep
#' their own rates (below 80 with the Kannisto extension, all ages otherwise)
#' are replaced multiplicatively in each year: each zero gets the share
#' `delta` of the deaths of the year, the other shares are multiplied by
#' `1 - z delta` (`z` zeros), and `delta` is half the smallest positive count
#' over all ages and years divided by the deaths of the year. A warning is
#' given when more than 10% of the cells are zero. Package safeguards: the
#' imputed share of a year is at most 1/2, and a year whose deaths at ages
#' from 80 allow no Kannisto fit takes the fit to all years together (with a
#' warning).
#'
#' **Deviations from the paper.** Besides the package-specific treatment of
#' other age groups: the life tables use `a_x = n / 2` at every age (the
#' authors' code uses 0.06 at age 0); the residuals of the bootstrap are
#' those of eq. (5), as the paper describes (the authors' code takes them
#' around the in-sample one-step forecasts of the ARIMA model); the number of
#' refits follows `n_draws`; the draws are centred as in the authors' code
#' (see above).
#'
#' @param jump_off `"observed"` (forecasts start from the observed deaths of
#'   the last year) or `"fitted"` (from the fitted deaths).
#' @param coherent Fit the coherent model with a common component per group.
#' @param common_by Stratum variable(s) defining groups with their own common
#'   component (implies `coherent = TRUE`); `NULL` with `coherent = TRUE` puts
#'   all strata in one group.
#' @param order ARIMA order `c(p, d, q)` of the index (of the common index in
#'   the coherent model); `c(0, 1, 1)` by default (Bergeron-Boucher et al.
#'   2017). `c(0, 1, 0)` with `drift = TRUE` is the random walk with drift of
#'   Lee and Carter (1992), the default of the R package MortalityForecast.
#' @param drift Include a drift (a constant of the differenced index) when
#'   `d = 1`; ignored otherwise (`d = 0` includes a mean, `d = 2` no
#'   constant).
#' @param order_dev,drift_dev ARIMA order and drift of the indexes of the
#'   deviations of each stratum in the coherent model; `c(1, 1, 0)` without
#'   drift by default (Bergeron-Boucher et al. 2017).
#' @param open Treatment of the old ages: `"auto"` (the Kannisto extension
#'   where the age groups allow it, otherwise the regression of the open
#'   interval), `"kannisto"` or `"regression"`.
#' @return A model specification.
#' @references
#' Oeppen J (2008). Coherent forecasting of multiple-decrement life tables:
#' a test using Japanese cause of death data. CODAWORK'08, Universitat de
#' Girona. \url{https://hdl.handle.net/10256/742}
#'
#' Bergeron-Boucher M-P, Canudas-Romo V, Oeppen J, Vaupel JW (2017).
#' Coherent forecasts of mortality with compositional data analysis.
#' *Demographic Research* 37(17): 527-566. \doi{10.4054/DemRes.2017.37.17}
#'
#' Thatcher AR, Kannisto V, Vaupel JW (1998). *The Force of Mortality at Ages
#' 80 to 120*. Odense University Press.
#'
#' Li N, Lee R (2005). Coherent mortality forecasts for a group of
#' populations: an extension of the Lee-Carter method. *Demography* 42(3):
#' 575-594. \doi{10.1353/dem.2005.0021}
#'
#' Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
#' *Journal of the American Statistical Association* 87(419): 659-671.
#' @seealso [coda_causes()] for causes of death, [lee_carter()] for the model
#'   of death rates, [li_lee()] for its coherent version.
#' @examples
#' x <- mt_example(deaths = "single", years = 1995:2012)
#' fc <- mt_forecast(mt_fit(coda(), x), h = 5, n_draws = 50, seed = 1)
#' fc2 <- mt_forecast(mt_fit(coda(coherent = TRUE), x), h = 5, n_draws = 50, seed = 1)
#' # random walk with drift for the index (Lee and Carter 1992)
#' fc3 <- mt_forecast(mt_fit(coda(order = c(0, 1, 0)), x), h = 5, n_draws = 50, seed = 1)
#' @export
coda <- function(jump_off = c("observed", "fitted"), coherent = FALSE, common_by = NULL, order = c(0, 1, 1),
                 drift = TRUE, order_dev = c(1, 1, 0), drift_dev = FALSE, open = c("auto", "kannisto", "regression")) {
  jump_off <- match.arg(jump_off)
  open <- match.arg(open)
  check_flag(coherent, "coherent")
  check_flag(drift, "drift")
  check_flag(drift_dev, "drift_dev")
  order <- check_arima_order(order, "order")
  order_dev <- check_arima_order(order_dev, "order_dev")
  if (!is.null(common_by)) {
    if (!is.character(common_by) || !length(common_by)) {
      mt_abort("`common_by` must name stratum variables.", class = "argument")
    }
    if (!missing(coherent) && !coherent) {
      mt_abort("`common_by` defines groups of the coherent model; it needs `coherent = TRUE`.", class = "argument")
    }
    coherent <- TRUE
  }
  label <- if (!coherent) "CoDa" else if (is.null(common_by)) "Coherent CoDa" else
    sprintf("Coherent CoDa (common by %s)", paste(common_by, collapse = ", "))
  ts_label <- c(if (!identical(order, c(0, 1, 1)) || !drift) arima_label(order, drift),
                if (coherent && (!identical(order_dev, c(1, 1, 0)) || drift_dev))
                  paste("deviations", arima_label(order_dev, drift_dev)))
  new_mt_spec("coda", list(jump_off = jump_off, coherent = coherent, common_by = common_by, order = order,
                           drift = drift, order_dev = order_dev, drift_dev = drift_dev, open = open),
              joint = coherent, min_periods = if (coherent) 8 else 6,
              label = paste0(label, if (length(ts_label)) paste0(", ", paste(ts_label, collapse = ", ")),
                             if (open == "regression") ", open regression" else "",
                             if (jump_off == "fitted") ", fitted jump-off" else ""))
}

#' @export
fit_model.mt_spec_coda <- function(spec, D, E, ages, periods, strata = NULL, hierarchy = NULL, ...) {
  p <- spec$params
  if (isTRUE(p$coherent)) return(coda_fit_coherent(spec, D, E, ages, strata, hierarchy))
  layout <- coda_layout(ages, p$open %||% "auto")
  sc <- coda_schedule(D, E, layout, "the stratum")
  fit <- coda_core(coda_clr(sc$m, layout), 1, p$order %||% c(0, 1, 1), p$drift %||% TRUE)
  if (!is.null(sc$open_log)) fit$open <- coda_open_fit(sc$open_log, fit$k)
  fit$kannisto <- sc$kannisto
  fit$layout <- layout
  fit$parts <- layout$parts
  fit$jump_off <- p$jump_off
  fit$fitted <- coda_layout_rates(fit$alpha + fit$b %*% t(fit$k), layout, coda_open_fitted(fit$open, fit$k))
  fit
}

#' @export
simulate_model.mt_spec_coda <- function(spec, fit, h, n_draws, ...) {
  if (isTRUE(spec$params$coherent)) return(coda_simulate_coherent(fit, h, n_draws))
  pt <- coda_point(fit, h)
  bt <- coda_boot(fit, h, n_draws, open_y = fit$open$y)
  sims <- coda_centre(fit$alpha + bt$Z, pt$Z[, 1], bt$open, pt$open[1])
  array(coda_layout_rates(sims$Z, fit$layout, sims$open), c(length(fit$layout$ages), h, n_draws))
}

# Point forecast of a CoDa fit of one stratum: clr life-table deaths
# [part, h] and the log rate of the open interval (method "regression"),
# from the observed or fitted deaths of the last year
coda_point <- function(fit, h) {
  kf <- coda_index_mean(fit, h)
  kT <- fit$k[nrow(fit$k), ]
  obs <- identical(fit$jump_off, "observed")
  base <- if (obs) fit$Z_last - drop(fit$b %*% kT) else fit$alpha
  open <- NULL
  if (!is.null(fit$open)) {
    open <- (if (obs) fit$open$last - sum(fit$open$b * kT) else fit$open$a) + drop(kf %*% fit$open$b)
  }
  list(Z = base + fit$b %*% t(kf), open = open, k = kf)
}

coda_open_fitted <- function(open, k) if (is.null(open)) NULL else open$a + drop(as.matrix(k) %*% open$b)

# Coherent CoDa model: per group, a common component from the average death
# rates of its (bottom) strata, and per stratum the first component of the
# deviations from it
coda_fit_coherent <- function(spec, D, E, ages, strata, hierarchy) {
  p <- spec$params
  na <- dim(D)[1]
  nt <- dim(D)[2]
  ns <- dim(D)[3]
  group <- rep(1L, ns)
  if (!is.null(p$common_by)) {
    if (is.null(strata) || !all(p$common_by %in% names(strata))) {
      mt_abort("`common_by` must name stratum variables.", class = "argument")
    }
    key <- do.call(paste, c(lapply(p$common_by, function(v) strata[[v]]), sep = "."))
    group <- match(key, unique(key))
  }
  pool <- rep(TRUE, ns)
  if (!is.null(hierarchy) && !is.null(strata)) pool <- hierarchy_structure(strata, hierarchy)$is_bottom
  labels <- if (is.null(strata)) as.character(seq_len(ns)) else strata$label
  layout <- coda_layout(ages, p$open %||% "auto")
  sched <- lapply(seq_len(ns), function(s) {
    coda_schedule(D[, , s], E[, , s], layout, sprintf("stratum %s", labels[s]))
  })
  Z <- lapply(sched, function(sc) coda_clr(sc$m, layout))
  reg <- layout$method == "regression"
  commons <- lapply(seq_len(max(group)), function(g) {
    use <- group == g & pool
    if (!any(use)) use <- group == g
    # average mortality of the group: arithmetic mean of the rates
    mbar <- Reduce(`+`, lapply(sched[use], `[[`, "m")) / sum(use)
    out <- coda_core(coda_clr(mbar, layout), 1, p$order %||% c(0, 1, 1), p$drift %||% TRUE)
    if (reg) {
      obar <- Reduce(`+`, lapply(sched[use], function(sc) exp(sc$open_log))) / sum(use)
      out$open <- coda_open_fit(log(obar), out$k)
    }
    out$members <- which(use)
    out
  })
  fitted <- array(NA_real_, c(na, nt, ns))
  own <- lapply(seq_len(ns), function(s) {
    cm <- commons[[group[s]]]
    common <- cm$b %*% t(cm$k)
    out <- coda_core(Z[[s]], 1, p$order_dev %||% c(1, 1, 0), p$drift_dev %||% FALSE, common = common,
                     adjust = FALSE)
    out$kannisto <- sched[[s]]$kannisto
    lo_fit <- NULL
    if (reg) {
      lo <- sched[[s]]$open_log
      out$open <- coda_open_fit(lo - drop(cm$k %*% cm$open$b), out$k)
      out$open$last <- lo[nt]
      lo_fit <- out$open$a + drop(cm$k %*% cm$open$b) + drop(out$k %*% out$open$b)
    }
    fitted[, , s] <<- coda_layout_rates(out$alpha + common + out$b %*% t(out$k), layout, lo_fit)
    out
  })
  list(commons = commons, group = group, strata = own, fitted = fitted, layout = layout, jump_off = p$jump_off,
       parts = layout$parts)
}

# Point forecasts of the coherent model: clr deaths [part, h] and log open
# rates of stratum s
coda_point_coherent <- function(fit, s, h, Kf = NULL) {
  cm <- fit$commons[[fit$group[s]]]
  st <- fit$strata[[s]]
  Kf <- Kf %||% coda_index_mean(cm, h)
  kf <- coda_index_mean(st, h)
  KT <- cm$k[nrow(cm$k), ]
  kT <- st$k[nrow(st$k), ]
  obs <- identical(fit$jump_off, "observed")
  base <- if (obs) st$Z_last - drop(cm$b %*% KT) - drop(st$b %*% kT) else st$alpha
  open <- NULL
  if (!is.null(st$open)) {
    base_o <- if (obs) st$open$last - sum(cm$open$b * KT) - sum(st$open$b * kT) else st$open$a
    open <- base_o + drop(Kf %*% cm$open$b) + drop(kf %*% st$open$b)
  }
  list(Z = base + cm$b %*% t(Kf) + st$b %*% t(kf), open = open, K = Kf, k = kf)
}

coda_simulate_coherent <- function(fit, h, n_draws) {
  layout <- fit$layout
  na <- length(layout$ages)
  ns <- length(fit$strata)
  common <- lapply(fit$commons, function(cm) {
    list(boot = coda_boot(cm, h, n_draws, open_y = cm$open$y), K = coda_index_mean(cm, h))
  })
  out <- array(NA_real_, c(na, h, ns, n_draws))
  for (s in seq_len(ns)) {
    g <- fit$group[s]
    cm <- fit$commons[[g]]
    st <- fit$strata[[s]]
    pt <- coda_point_coherent(fit, s, h, common[[g]]$K)
    y_dev <- if (!is.null(st$open)) st$open$y else NULL
    dev <- coda_boot(st, h, n_draws, open_y = y_dev)
    Zs <- st$alpha + common[[g]]$boot$Z + dev$Z
    open <- NULL
    if (!is.null(st$open)) open <- common[[g]]$boot$open - cm$open$a + dev$open
    sims <- coda_centre(Zs, pt$Z[, 1], open, pt$open[1])
    out[, , s, ] <- coda_layout_rates(sims$Z, layout, sims$open)
  }
  out
}

#' Multiple-decrement CoDa model of causes of death (experimental)
#'
#' `coda_causes()` forecasts death rates by cause with the multiple-decrement
#' form of the compositional Lee-Carter model, which Oeppen (2008) proposed
#' and tested on Japanese causes of death. Causes add up to all causes in
#' every draw. The model is experimental. Status: published, with the
#' deviations listed below; section numbers refer to Oeppen (2008).
#'
#' **Model (Sections 3, 6 and 8.3).** In each group (the strata that differ
#' only in `var`, for example one group per sex):
#'
#' 1. the cause-specific rates `m[x, c]` add up to all-cause rates `m[x]`,
#'    which give life-table deaths `d[x]` with radix 1 ([life_table()]
#'    conventions); the deaths are split by cause, `d[x, c] = d[x] m[x, c] /
#'    m[x]` (relations of the multiple-decrement life table, eq. 1 and 7);
#' 2. the array of deaths by age and cause of each year is unfolded into one
#'    composition (Section 6), plus the survivors at the upper bound of a
#'    closed age range;
#' 3. the compositions are centred by their closed geometric mean over the
#'    years and transformed by the centred log-ratio; the first `rank`
#'    components of a singular value decomposition give `sum_j b_j[x, c]
#'    k_j[t]` (Section 3, steps 1-7);
#' 4. forecast deaths `d[x, c]` give all-cause rates from the all-cause
#'    life-table deaths `d[x] = sum_c d[x, c]`, and cause-specific rates
#'    `m[x, c] = m[x] d[x, c] / d[x]` (eq. 7, Section 8.3).
#'
#' The causes add up to the all-cause rates in every draw, and the strata of
#' `total_level` get their sum. If the listed causes leave part of the
#' all-cause deaths unexplained (`residual`), that part is a residual cause
#' inside the model and counts towards the all-cause rates; it is not a
#' stratum of the output (package-specific).
#'
#' **Indexes and rank.** For Japanese women 1951-1990, Oeppen (2008, Sections
#' 6 and 7) uses two components and selects ARIMA(0,2,2) models for the
#' indexes by AICc; Section 3 notes, for the Lee-Carter model, that one
#' component is usually adequate for advanced economies and that its index
#' can be forecast by a simple ARIMA model such as the random walk with drift.
#' These choices belong to one data set; the defaults here are one component
#' and the index model of [coda()], ARIMA(0,1,1) with drift (Bergeron-Boucher
#' et al. 2017, with Oeppen as co-author). `rank = 2, order = c(0, 2, 2)`
#' gives Oeppen's specification for Japan. With several components, each
#' index has its own ARIMA model, independent of the others (the paper does
#' not say how the second index was forecast).
#'
#' **Deviations from the paper.** Oeppen (2008) gives point forecasts only:
#' the jump-off (forecasts start from the observed deaths of the last year,
#' `jump_off = "observed"`) and the draws (bootstrap of the residuals and of
#' the index paths, centred on the point forecast) are those of [coda()]
#' (Bergeron-Boucher et al. 2017, Appendices B and G). Oeppen (2008, Section
#' 3) replaces zero counts by lightly smoothing the deaths of each year with
#' a penalised histogram smoother whose smoothing parameter is not given;
#' here zero counts are replaced multiplicatively in the composition of each
#' year (Bergeron-Boucher et al. 2017, Appendix C; see [coda()]). How the
#' death rate of the open age interval follows from the forecast deaths is
#' not stated; here the all-cause log rate of an open interval is regressed
#' on the indexes and split by cause through the composition
#' (package-specific, as `open = "regression"` of [coda()]).
#'
#' Draws are expected rates. Predictive draws (`mt_forecast(predictive =
#' TRUE)`) draw Poisson all-cause deaths given the all-cause rates and
#' exposures and split them multinomially between the causes by the cause
#' shares of the draw, so that the causes also add up to all causes in
#' predictive draws; no other observation noise is added (package-specific).
#'
#' The extensions of Kjærgaard, Ergemen, Kallestrup-Lamb, Oeppen and
#' Lindahl-Jacobsen (2019) are not implemented. Compare the model with
#' [reconciled()] separate models in a backtest before using it.
#'
#' @param var Stratum variable with the causes.
#' @param total_level Level of `var` for all causes.
#' @param jump_off `"observed"` (forecasts start from the observed deaths of
#'   the last year) or `"fitted"` (from the fitted deaths).
#' @param residual Add a residual cause, the part of all-cause mortality not
#'   covered by the listed causes: `TRUE`, `FALSE`, or `"auto"` (added when
#'   the listed causes leave more than 0.1% of all-cause mortality
#'   unexplained in the median cell).
#' @param rank Number of components (1 by default; Oeppen 2008 used 2 for
#'   Japan).
#' @param order,drift ARIMA order and drift of each index, as in [coda()].
#' @return A model specification.
#' @references
#' Oeppen J (2008). Coherent forecasting of multiple-decrement life tables:
#' a test using Japanese cause of death data. CODAWORK'08, Universitat de
#' Girona. \url{https://hdl.handle.net/10256/742}
#'
#' Bergeron-Boucher M-P, Canudas-Romo V, Oeppen J, Vaupel JW (2017).
#' Coherent forecasts of mortality with compositional data analysis.
#' *Demographic Research* 37(17): 527-566. \doi{10.4054/DemRes.2017.37.17}
#'
#' Kjærgaard S, Ergemen YE, Kallestrup-Lamb M, Oeppen J, Lindahl-Jacobsen R
#' (2019). Forecasting causes of death by using compositional data analysis:
#' the case of cancer deaths. *Journal of the Royal Statistical Society
#' Series C* 68(5): 1351-1370. \doi{10.1111/rssc.12357}
#' @seealso [reconciled()] with a hierarchy of causes, [coda()].
#' @examples
#' x <- mt_example_causes(years = 2005:2019)
#' fit <- mt_fit(coda_causes(), x)
#' fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
#' @export
coda_causes <- function(var = "cause", total_level = "all", jump_off = c("observed", "fitted"), residual = "auto",
                        rank = 1, order = c(0, 1, 1), drift = TRUE) {
  jump_off <- match.arg(jump_off)
  if (!is.character(var) || length(var) != 1L) mt_abort("`var` must be one stratum variable.", class = "argument")
  if (length(total_level) != 1L) mt_abort("`total_level` must be one level.", class = "argument")
  if (!identical(residual, "auto") && !(is.logical(residual) && length(residual) == 1L && !is.na(residual))) {
    mt_abort("`residual` must be TRUE, FALSE or \"auto\".", class = "argument")
  }
  check_number(rank, "rank", lower = 1, upper = 5, integerish = TRUE)
  check_flag(drift, "drift")
  order <- check_arima_order(order, "order")
  ts_label <- c(if (rank > 1) sprintf("rank %d", as.integer(rank)),
                if (!identical(order, c(0, 1, 1)) || !drift) arima_label(order, drift))
  new_mt_spec("coda_causes", list(var = var, total_level = total_level, jump_off = jump_off, residual = residual,
                                  rank = as.integer(rank), order = order, drift = drift),
              joint = TRUE, min_periods = 6,
              label = paste0("CoDa causes (multiple decrement)",
                             if (length(ts_label)) paste0(", ", paste(ts_label, collapse = ", ")),
                             if (jump_off == "fitted") ", fitted jump-off" else ""))
}

#' @export
mt_fit.mt_spec_coda_causes <- function(spec, data, prep = NULL, cores = 1, ...) {
  fit <- NextMethod()
  class(fit) <- c("mt_fit_coda_causes", class(fit))
  fit
}

#' @export
fit_model.mt_spec_coda_causes <- function(spec, D, E, ages, periods, strata = NULL, hierarchy = NULL, ...) {
  p <- spec$params
  if (is.null(strata) || !p$var %in% names(strata)) {
    mt_abort(sprintf("`%s` is not a stratum variable.", p$var), class = "argument")
  }
  na <- dim(D)[1]
  nt <- dim(D)[2]
  ns <- dim(D)[3]
  lev <- as.character(strata[[p$var]])
  others <- setdiff(names(strata), c("stratum", "label", p$var))
  gkey <- if (length(others)) do.call(paste, c(lapply(others, function(o) strata[[o]]), sep = ".")) else
    rep("all", ns)
  groups <- unique(gkey)
  causes <- setdiff(unique(lev), p$total_level)
  if (!length(causes)) mt_abort("No cause levels besides the total.", class = "argument")
  find <- function(g, level, what) {
    i <- which(gkey == g & lev == level)
    if (length(i) != 1) mt_abort(sprintf("Group '%s' %s.", g, what), class = "model_data")
    i
  }
  tot_id <- vapply(groups, function(g) find(g, p$total_level, "has no all-cause stratum"), 1L, USE.NAMES = FALSE)
  part_id <- matrix(vapply(causes, function(cc) {
    vapply(groups, function(g) find(g, cc, sprintf("lacks cause '%s'", cc)), 1L, USE.NAMES = FALSE)
  }, integer(length(groups))), nrow = length(groups))
  K <- length(causes)
  layout <- coda_layout(ages, "regression")
  open <- layout$method == "regression"
  # rates of the causes [age, period, cause] and of all causes, per group
  mc <- lapply(seq_along(groups), function(g) {
    array(D[, , part_id[g, ], drop = FALSE] / E[, , part_id[g, ], drop = FALSE], c(na, nt, K))
  })
  mt <- lapply(seq_along(groups), function(g) matrix(D[, , tot_id[g]] / E[, , tot_id[g]], na))
  short <- unlist(lapply(seq_along(groups), function(g) (1 - rowSums(mc[[g]], dims = 2) / mt[[g]])[mt[[g]] > 0]))
  use_residual <- if (identical(p$residual, "auto")) stats::median(short) > 1e-3 else isTRUE(p$residual)
  nc <- K + use_residual
  rank <- p$rank %||% 1L
  fitted <- array(NA_real_, c(na, nt, ns))
  fits <- lapply(seq_along(groups), function(g) {
    m <- mc[[g]]
    if (use_residual) m <- array(c(m, pmax(mt[[g]] - rowSums(m, dims = 2), 0)), c(na, nt, nc))
    m[!is.finite(m)] <- 0
    # zero deaths: multiplicative replacement in the composition by age and
    # cause of each year (deaths on the all-cause exposures)
    Et <- matrix(E[, , tot_id[g]], na)
    Dg <- matrix(aperm(m * as.vector(Et), c(1, 3, 2)), na * nc, nt)
    Dg <- coda_replace_zeros(Dg, sprintf("the causes of group %s", groups[g]))
    m <- aperm(array(Dg, c(na, nc, nt)), c(1, 3, 2)) / as.vector(Et)
    out <- coda_core(md_forward(m, ages), rank, p$order %||% c(0, 1, 1), p$drift %||% TRUE)
    if (open) out$open <- coda_open_fit(log(rowSums(m, dims = 2)[na, ]), out$k)
    out$jump_off <- p$jump_off
    R <- md_rates(out$alpha + out$b %*% t(out$k), layout, nc, coda_open_fitted(out$open, out$k))
    for (j in seq_len(K)) fitted[, , part_id[g, j]] <<- R[, , j]
    fitted[, , tot_id[g]] <<- md_total(R)
    out
  })
  parts <- c(paste(rep(ages$labels, nc), rep(c(causes, if (use_residual) "residual"), each = na)),
             if (!open) "survivors")
  list(groups = fits, group_labels = groups, tot_id = tot_id, part_id = part_id, causes = causes, K = K,
       residual = use_residual, parts = parts, ages = ages, layout = layout, fitted = fitted)
}

#' @export
simulate_model.mt_spec_coda_causes <- function(spec, fit, h, n_draws, ...) {
  na <- length(fit$ages)
  ns <- length(c(fit$tot_id, fit$part_id))
  out <- array(NA_real_, c(na, h, ns, n_draws))
  nc <- fit$K + fit$residual
  for (g in seq_along(fit$groups)) {
    gf <- fit$groups[[g]]
    pt <- coda_point(gf, h)
    bt <- coda_boot(gf, h, n_draws, open_y = gf$open$y)
    sims <- coda_centre(gf$alpha + bt$Z, pt$Z[, 1], bt$open, pt$open[1])
    R <- md_rates(sims$Z, fit$layout, nc, sims$open)
    for (j in seq_len(fit$K)) out[, , fit$part_id[g, j], ] <- R[, , j]
    out[, , fit$tot_id[g], ] <- md_total(R)
  }
  out
}

#' @export
mt_forecast.mt_fit_coda_causes <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                                           exposure = NULL, ...) {
  if (!is.null(to)) {
    if (!missing(h)) mt_warn("Both `h` and `to` are given; `to` is used.", class = "argument")
    h <- to - max(fit$periods)
  }
  if (!predictive) return(mt_forecast.mt_fit(fit, h = h, n_draws = n_draws, seed = seed, ...))
  with_seed(seed, {
    fc <- mt_forecast.mt_fit(fit, h = h, n_draws = n_draws, ...)
    fc$rates[] <- coda_causes_predictive(fc$rates, fit, exposure)
    fc$predictive <- TRUE
    fc$seed <- seed
    fc
  })
}

# Predictive draws: Poisson all-cause deaths given the all-cause rates and
# exposures, split multinomially by the cause shares of the draw (successive
# binomial draws; without a residual cause the last cause takes the rest)
coda_causes_predictive <- function(arr, fit, exposure) {
  d <- dim(arr)
  E <- expand_exposure(exposure %||% fit$exposure_last, d, dimnames(arr))
  f <- fit$fits
  out <- arr
  for (g in seq_along(f$tot_id)) {
    i <- f$tot_id[g]
    tot <- arr[, , i, , drop = FALSE]
    mu <- tot * as.vector(E[, , i])
    ok <- !is.na(mu)
    Dt <- array(NA_real_, dim(tot))
    Dt[ok] <- stats::rpois(sum(ok), mu[ok])
    out[, , i, ] <- Dt / pmax(as.vector(E[, , i]), 1e-12)
    left <- Dt
    rest <- array(1, dim(tot))
    for (k in seq_len(f$K)) {
      j <- f$part_id[g, k]
      q <- arr[, , j, , drop = FALSE] / tot
      dk <- left
      if (k < f$K || isTRUE(f$residual)) {
        pr <- q / rest
        pr[!is.finite(pr)] <- 1
        okk <- !is.na(left)
        dk[okk] <- stats::rbinom(sum(okk), left[okk], pmin(pmax(pr[okk], 0), 1))
      }
      out[, , j, ] <- dk / pmax(as.vector(E[, , j]), 1e-12)
      left <- left - dk
      rest <- pmax(rest - q, 0)
    }
  }
  out
}

# Multiple-decrement composition: cause rates [age, period, cause] -> clr of
# the life-table deaths by age and cause [age x cause (+ survivors), period]
md_forward <- function(mc, ages) {
  d <- dim(mc)
  na <- d[1]
  m <- rowSums(mc, dims = 2)
  lt <- lt_columns(m, ages)
  dxc <- mc / as.vector(m) * as.vector(lt$dx)
  comp <- matrix(aperm(dxc, c(1, 3, 2)), na * d[3], d[2])
  if (!is.infinite(ages$upper[na])) comp <- rbind(comp, lt$lx[na, ] - lt$dx[na, ])
  z <- log(pmax(comp, 1e-12))
  sweep(z, 2, colMeans(z))
}

# Inverse of md_forward(): cause rates [age, column, cause] from clr deaths
# [part, column]; all-cause rates from the all-cause life-table deaths
# (layout "regression" or "closed"), the open interval from `open_log` (log
# all-cause rates); cause rates by eq. (7) of Oeppen (2008)
md_rates <- function(Z, layout, nc, open_log = NULL) {
  na <- length(layout$ages)
  Zm <- matrix(Z, nrow(Z))
  P <- exp(Zm)
  P <- sweep(P, 2, colSums(P), "/")
  dxc <- array(P[seq_len(na * nc), , drop = FALSE], c(na, nc, ncol(P)))
  dx <- matrix(dxc[, 1, ], na)
  for (j in seq_len(nc)[-1]) dx <- dx + dxc[, j, ]
  lz <- log(if (nrow(P) > na * nc) rbind(dx, P[nrow(P), ]) else dx)
  m <- matrix(coda_layout_rates(lz, layout, open_log), na)
  out <- array(NA_real_, c(na, ncol(P), nc))
  for (j in seq_len(nc)) out[, , j] <- m * matrix(dxc[, j, ], na) / dx
  out
}

# All-cause rates as the sum of the cause rates [age, column, cause]
md_total <- function(R) {
  tot <- R[, , 1]
  for (j in seq_len(dim(R)[3])[-1]) tot <- tot + R[, , j]
  tot
}

# Parameters of CoDa fits, for mt_parameters(): the closed centre `a`, the
# loadings `b` (parts of the composition: ages, ages by cause, survivors),
# the indexes with the coefficients and innovation standard deviation of
# their ARIMA models, the share of variation explained, the loadings of the
# open interval and the Kannisto parameters of each year
coda_parameters <- function(fit, sp = fit$spec) {
  lab <- fit$strata$label
  yrs <- as.character(fit$periods)
  f <- fit$fits
  row <- function(stratum, parameter, index, value) {
    data.table(stratum = stratum, parameter = parameter, index = as.character(index), value = as.numeric(value))
  }
  centre <- function(alpha) exp(alpha) / sum(exp(alpha))
  one <- function(who, x, parts, idx = "k") {
    b_name <- if (idx == "K") "B" else "b"
    r <- ncol(x$k)
    suffix <- if (r > 1) seq_len(r) else ""
    out <- list(row(who, "a", parts, centre(x$alpha)))
    for (j in seq_len(r)) {
      tsj <- x$ts[[j]]
      cf <- tsj$coef
      out <- c(out, list(row(who, paste0(b_name, suffix[j]), parts, x$b[, j]),
                         row(who, paste0(idx, suffix[j]), yrs, x$k[, j]),
                         row(who, paste0("order_", idx, suffix[j]), c("p", "d", "q"), tsj$order),
                         if (length(cf)) row(who, paste0("arima_", idx, suffix[j]), names(cf), cf),
                         row(who, paste0("sd_", idx, suffix[j]), "", tsj$sigma),
                         row(who, paste0("explained", suffix[j]), "", x$explained[j])))
    }
    if (!is.null(x$open)) {
      out <- c(out, list(row(who, "a_open", "", x$open$a), row(who, paste0(b_name, "_open"), suffix, x$open$b),
                         row(who, "sd_open", "", x$open$sd)))
    }
    if (!is.null(x$kannisto)) {
      out <- c(out, list(row(who, "kannisto_a", yrs, exp(x$kannisto[1, ])), row(who, "kannisto_b", yrs,
                                                                                   x$kannisto[2, ])))
    }
    rbindlist(out)
  }
  if (inherits(sp, "mt_spec_coda_causes")) {
    return(rbindlist(lapply(seq_along(f$groups), function(g) {
      one(paste0("causes: ", f$group_labels[g]), f$groups[[g]], f$parts)
    })))
  }
  if (!isTRUE(sp$params$coherent)) {
    return(rbindlist(lapply(seq_along(f), function(s) if (!is.null(f[[s]])) one(lab[s], f[[s]], f[[s]]$parts))))
  }
  common <- rbindlist(lapply(seq_along(f$commons), function(g) {
    cm <- f$commons[[g]]
    one(sprintf("common %d (%s)", g, paste(lab[f$group == g], collapse = ", ")), cm, f$parts, "K")
  }))
  own <- rbindlist(lapply(seq_along(f$strata), function(s) one(lab[s], f$strata[[s]], f$parts)))
  rbind(common, own)
}

# Statement of the model and its uncertainty, for summary()
coda_assumptions <- function(spec) {
  p <- spec$params
  draws <- paste("Draws (Bergeron-Boucher et al. 2017, Appendix B): residuals resampled and the components and",
                 "ARIMA models estimated again, index paths with resampled innovations; centred on the point",
                 "forecast in the first forecast year.")
  ts <- arima_label(p$order %||% c(0, 1, 1), p$drift %||% TRUE)
  if (inherits(spec, "mt_spec_coda_causes")) {
    return(paste("Multiple-decrement CoDa model (Oeppen 2008; experimental): per group, clr of the life-table",
                 sprintf("deaths by age and cause = a + sum of %d component(s) b k(t) of the SVD; k(t) %s.",
                         p$rank %||% 1L, ts),
                 "Causes add up to all causes in every draw; predictive draws split Poisson all-cause deaths",
                 "multinomially.", draws, "Jump-off:", p$jump_off, "deaths. The log rate of an open age interval",
                 "is regressed on the indexes."))
  }
  old <- switch(p$open %||% "auto",
                regression = "The log rate of an open age interval is regressed on the indexes.",
                paste("Ages from 80 by a Kannisto model extended to 120 where the age groups are single years from 80",
                      "(otherwise the log rate of an open interval is regressed on the indexes)."))
  if (!isTRUE(p$coherent)) {
    return(paste("CoDa Lee-Carter model (Oeppen 2008; Bergeron-Boucher et al. 2017): clr of the life-table deaths",
                 sprintf("= clr(a) + b k(t), first SVD component; k(t) %s.", ts), draws, "Jump-off:", p$jump_off,
                 "deaths.", old))
  }
  paste("Coherent CoDa model (Bergeron-Boucher et al. 2017): clr of the life-table deaths = clr(a_i) + B K(t) +",
        "b_i k_i(t), with B K from the average death rates of the strata of each group",
        if (!is.null(p$common_by)) sprintf("(groups by %s)", paste(p$common_by, collapse = ", ")) else "",
        sprintf("; K(t) %s, k_i(t) %s.", ts, arima_label(p$order_dev %||% c(1, 1, 0), p$drift_dev %||% FALSE)),
        draws, "Jump-off:", p$jump_off, "deaths.", old)
}
