# Multiway compositional (CoDa) model of life-table deaths for several
# populations: coda_multiway(). Written from Bergeron-Boucher, Simonacci,
# Oeppen and Gallo (2018); no code or supplementary material of the paper
# was found. The Tucker3 fit uses hooi() and mode_product() of
# model-tucker3.R; life tables, old ages, zero deaths and the ARIMA models of
# the indexes are those of coda() (model-coda-core.R).

#' Multiway compositional (CoDa) model of life-table deaths for subpopulations
#'
#' `coda_multiway()` is the three-way compositional model of Bergeron-Boucher,
#' Simonacci, Oeppen and Gallo (2018): a Tucker3 decomposition of the centred
#' log-ratios of the life-table deaths of several populations, with time
#' indexes shared by all of them. The strata are the populations (for example
#' the regions of a country). Status: published, with the choices listed
#' below where the paper leaves the method open. Section and equation numbers
#' refer to that paper. The model is experimental: compare it with
#' [coda()] (`coherent = TRUE`) and [li_lee()] in a backtest before using it.
#'
#' **Model (Section 2.2, eq. 5).** For population `i`:
#'
#' 1. the death rates of each year give life-table deaths `d[x, t, i]` with
#'    radix 1 at the first age ([life_table()] conventions), as in [coda()];
#' 2. `a[x, i]` is the geometric mean of `d[x, t, i]` over the years, closed
#'    to 1;
#' 3. the deaths are centred by perturbation (`d[x, t, i] / a[x, i]`, closed)
#'    and transformed by the centred log-ratio (clr), eq. (3);
#' 4. the array of centred log-ratios `[age, year, population]` follows the
#'    Tucker3 model `clr(d[x, t, i] / a[x, i]) = sum_pqr g[p, q, r] b[x, p]
#'    k[t, q] gamma[i, r]`, with age patterns `b`, time indexes `k` and
#'    population markers `gamma` (the paper's beta, kappa and gamma), each of
#'    unit length, and the core array `g`; `ranks = c(P, Q, R)` gives the
#'    numbers of components of age, time and population. The default
#'    `c(1, 1, 1)` is the model of the paper's forecasts (Sections 4 and 5);
#' 5. the time index follows a random walk with drift (`order = c(0, 1, 0)`,
#'    `drift = TRUE`), the ARIMA model that the paper selected by AIC for the
#'    index of Canadian women (Section 5);
#' 6. forecasts are transformed back (inverse clr, perturbation by `a[x, i]`,
#'    closure) to life-table deaths and then to death rates.
#'
#' With one component per mode, population `i` changes by `g gamma[i] b[x]
#' k[t]`: all populations follow one index and one age pattern, at a pace set
#' by their markers. The forecasts are coherent in this sense; populations
#' with different markers keep changing at different paces (Section 5.1).
#' With one population, the fit and the point forecasts are those of
#' [coda()] with the same index model and jump-off; the draws of [coda()]
#' also carry the uncertainty of the estimates.
#'
#' **Estimation (Sections 2.2 and 4).** Least squares by alternating least
#' squares: higher-order orthogonal iteration (De Lathauwer, De Moor and
#' Vandewalle 2000), started from the higher-order singular value
#' decomposition. The paper uses the Tucker3 model for compositional data of
#' Gallo (2015), an alternating least-squares algorithm for the same
#' criterion. The clr of a composition sums to zero and the centring removes
#' the mean over the years, so `b` and `k` sum to zero. Signs: each time
#' index increases over the fitting period, each set of population markers
#' has a positive sum, and the largest core element of each age pattern is
#' positive; with one component, life-table deaths then move from ages with
#' negative `b` to ages with positive `b` as `k` rises. `fit$fits$fit_pct`
#' is the share of the sum of squares of the centred array that the model
#' explains (the explained variance of Section 4.1), `fit$fits$core_pct` the
#' share of each core element.
#'
#' **Draws (Section 5).** The paper draws its prediction intervals from
#' simulations of the time index with resampled errors. Here each time index
#' is simulated from its ARIMA model with innovations resampled from its
#' residuals; `a`, `b`, `gamma`, `g` and the drift are fixed, and the draws of
#' the index are shared by all populations. Draws are expected rates;
#' predictive draws (`mt_forecast(predictive = TRUE)`) add Poisson variation
#' of the deaths given the rates and the exposures, and no other observation
#' noise.
#'
#' **Choices where the paper is silent.**
#'
#' * Jump-off: the paper states no adjustment, so forecasts start from the
#'   fitted deaths (`jump_off = "fitted"`). `jump_off = "observed"` starts
#'   each population from its observed deaths of the last year, perturbed by
#'   the change of the model, as in [coda()] (Bergeron-Boucher, Canudas-Romo,
#'   Oeppen and Vaupel 2017, Appendix G).
#' * Several time indexes (`ranks[2] > 1`): the paper forecasts with one
#'   index. Each index then has its own ARIMA model, independent of the
#'   others, as for several components of a CoDa model in Shang and Haberman
#'   (2020).
#' * Other ARIMA orders (`order`, `drift`): the point forecast of an
#'   ARIMA(0,1,q) index with drift gets the first-step correction of [coda()]
#'   (Bergeron-Boucher et al. 2017, eq. 20); its draws are shifted by the
#'   same amount.
#' * Old ages and zero deaths, as in [coda()]: with single-year age groups
#'   from 80 (`open = "auto"` or `"kannisto"`), the rates from age 80 come
#'   from a Kannisto model fitted per year and population and extended to age
#'   120 (Bergeron-Boucher et al. 2017, Appendix C); with other age groups,
#'   the log death rate of the open interval of each population is regressed
#'   on the time indexes (package-specific). Zero deaths are replaced
#'   multiplicatively (Bergeron-Boucher et al. 2017, Appendix C). The paper
#'   uses five-year age groups to 100+ and five-year periods and combines the
#'   smallest populations to avoid zero counts; the models of this package
#'   need single calendar years.
#' * Strata of a hierarchy (see [mt_aggregate_strata()]) enter the
#'   decomposition as populations of their own; [reconciled()] makes their
#'   sums coherent.
#'
#' The paper chooses the numbers of components by the DifFit criterion
#' (Timmerman and Kiers 2000), gives bootstrap confidence intervals of the
#' parameters (Kiers 2004) and compares models by the AIC of Poisson
#' deviances (Table 2); these are not implemented. It fits women and men
#' separately, as `stratified("sex", female = coda_multiway(), male =
#' coda_multiway())` does.
#'
#' @param ranks Numbers of components `c(age, time, population)`, by default
#'   `c(1, 1, 1)`. They are reduced where the data allow fewer components
#'   (at most the number of parts of the composition minus 1, of years minus
#'   1 and of populations, and no rank above the product of the other two),
#'   with a warning (class `mt_warning_argument`) that names the requested
#'   and the ranks used; the label of the fit and [mt_parameters()]
#'   (`rank`) give the ranks used.
#' @param order ARIMA order `c(p, d, q)` of each time index; `c(0, 1, 0)` by
#'   default (random walk, Bergeron-Boucher et al. 2018, Section 5).
#' @param drift Include a drift when `d = 1`; ignored otherwise (`d = 0`
#'   includes a mean, `d = 2` no constant).
#' @param jump_off `"fitted"` (forecasts start from the fitted deaths) or
#'   `"observed"` (from the observed deaths of the last year).
#' @param open Treatment of the old ages: `"auto"` (the Kannisto extension
#'   where the age groups allow it, otherwise the regression of the open
#'   interval), `"kannisto"` or `"regression"`, as in [coda()].
#' @return A model specification.
#' @references
#' Bergeron-Boucher M-P, Simonacci V, Oeppen J, Gallo M (2018). Coherent
#' modeling and forecasting of mortality patterns for subpopulations using
#' multiway analysis of compositions: an application to Canadian provinces
#' and territories. *North American Actuarial Journal* 22(1): 92-118.
#' \doi{10.1080/10920277.2017.1377620}
#'
#' Bergeron-Boucher M-P, Canudas-Romo V, Oeppen J, Vaupel JW (2017).
#' Coherent forecasts of mortality with compositional data analysis.
#' *Demographic Research* 37(17): 527-566. \doi{10.4054/DemRes.2017.37.17}
#'
#' Gallo M (2015). Tucker3 model for compositional data. *Communications in
#' Statistics - Theory and Methods* 44(21): 4441-4453.
#' \doi{10.1080/03610926.2013.798664}
#'
#' Shang HL, Haberman S (2020). Forecasting age distribution of death counts:
#' an application to annuity pricing. *Annals of Actuarial Science* 14(1):
#' 150-169. \doi{10.1017/S1748499519000101}
#'
#' De Lathauwer L, De Moor B, Vandewalle J (2000). On the best rank-1 and
#' rank-(R1, R2, ..., RN) approximation of higher-order tensors. *SIAM Journal
#' on Matrix Analysis and Applications* 21(4): 1324-1342.
#' \doi{10.1137/S0895479898346995}
#'
#' Timmerman ME, Kiers HAL (2000). Three-mode principal components
#' analysis: choosing the numbers of components and sensitivity to local
#' optima. *British Journal of Mathematical and Statistical Psychology*
#' 53(1): 1-16. \doi{10.1348/000711000159132}
#'
#' Kiers HAL (2004). Bootstrap confidence intervals for three-way methods.
#' *Journal of Chemometrics* 18(1): 22-36. \doi{10.1002/cem.841}
#' @seealso [coda()] for one population and for the coherent model with a
#'   common component and deviations, [tucker3()] for the three-way model of
#'   log death rates.
#' @examples
#' x <- mt_example_regions(years = 2000:2015, regions = 3)
#' fit <- mt_fit(coda_multiway(), x)
#' fit$fits$fit_pct
#' fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
#' @export
coda_multiway <- function(ranks = c(1, 1, 1), order = c(0, 1, 0), drift = TRUE, jump_off = c("fitted", "observed"),
                          open = c("auto", "kannisto", "regression")) {
  jump_off <- match.arg(jump_off)
  open <- match.arg(open)
  check_flag(drift, "drift")
  order <- check_arima_order(order, "order")
  if (!is.numeric(ranks) || length(ranks) != 3 || anyNA(ranks) || any(ranks < 1) || any(ranks != round(ranks))) {
    mt_abort("`ranks` must be three whole numbers of at least 1 (age, time, population).", class = "argument",
             call = NULL)
  }
  ranks <- as.numeric(ranks)
  p <- list(ranks = ranks, order = order, drift = drift, jump_off = jump_off, open = open)
  new_mt_spec("coda_multiway", p, joint = TRUE, min_periods = 6, label = coda_multiway_label(p, ranks))
}

# Label of coda_multiway() with the ranks `ranks` (those requested, or those
# of the fit where the data allow fewer)
coda_multiway_label <- function(p, ranks) {
  ts_label <- if (!identical(p$order, c(0, 1, 0)) || !p$drift) arima_label(p$order, p$drift)
  paste0("Multiway CoDa",
         if (!identical(as.numeric(ranks), c(1, 1, 1))) sprintf(" (ranks %s)", paste(ranks, collapse = "/")),
         if (length(ts_label)) paste0(", ", ts_label),
         if (p$open == "regression") ", open regression" else "",
         if (p$jump_off == "observed") ", observed jump-off" else "")
}

#' @export
fit_model.mt_spec_coda_multiway <- function(spec, D, E, ages, periods, strata = NULL, ...) {
  p <- spec$params
  na <- length(ages)
  nt <- length(periods)
  ns <- length(D) %/% (na * nt)
  D <- array(D, c(na, nt, ns))
  E <- array(E, c(na, nt, ns))
  labels <- if (is.null(strata)) as.character(seq_len(ns)) else strata$label
  layout <- coda_layout(ages, p$open %||% "auto")
  sched <- lapply(seq_len(ns), function(s) {
    coda_schedule(D[, , s], E[, , s], layout, sprintf("stratum %s", labels[s]))
  })
  np <- length(layout$parts)
  Z <- array(unlist(lapply(sched, function(sc) coda_clr(sc$m, layout))), c(np, nt, ns))
  requested <- p$ranks %||% c(1, 1, 1)
  tk <- coda_tucker(Z, requested)
  capped <- any(tk$ranks != requested)
  if (capped) {
    mt_warn(sprintf(paste("The ranks c(%s) exceed what the data allow (the number of parts of the",
                          "composition minus 1, of years minus 1 and of populations, and no rank above the product",
                          "of the other two); ranks c(%s) are used."), paste(requested, collapse = ", "),
                    paste(tk$ranks, collapse = ", ")), class = "argument")
  }
  k <- tk$U[[2]]
  order <- p$order %||% c(0, 1, 0)
  drift <- p$drift %||% TRUE
  ts <- lapply(seq_len(ncol(k)), function(q) index_arima_fit(k[, q], order, drift))
  fit <- c(tk, list(k = k, rank = ncol(k), ts = ts, order = order, drift = drift, adjust = TRUE, layout = layout,
                    parts = layout$parts,
                    jump_off = p$jump_off %||% "fitted", Z_last = matrix(Z[, nt, ], np), labels = labels))
  # the label of the fit names the ranks used (mt_fit())
  if (capped) {
    fit$effective_label <- coda_multiway_label(list(order = order, drift = drift, open = p$open %||% "auto",
                                                    jump_off = fit$jump_off), tk$ranks)
  }
  # loadings of each population on the time indexes [part, index]
  fit$load <- lapply(seq_len(ns), function(s) coda_tucker_loadings(tk, s))
  fit$kannisto <- lapply(sched, `[[`, "kannisto")
  fit$open <- NULL
  if (layout$method == "regression") {
    fit$open <- lapply(sched, function(sc) coda_open_fit(sc$open_log, k))
  }
  fitted <- array(NA_real_, c(na, nt, ns))
  for (s in seq_len(ns)) {
    lo <- if (!is.null(fit$open)) coda_open_fitted(fit$open[[s]], k)
    fitted[, , s] <- coda_layout_rates(fit$alpha[, s] + fit$load[[s]] %*% t(k), layout, lo)
  }
  fit$fitted <- fitted
  fit
}

#' @export
simulate_model.mt_spec_coda_multiway <- function(spec, fit, h, n_draws, ...) {
  layout <- fit$layout
  na <- length(layout$ages)
  np <- length(fit$parts)
  ns <- length(fit$load)
  K <- coda_multiway_paths(fit, h, n_draws)
  out <- array(NA_real_, c(na, h, ns, n_draws))
  for (s in seq_len(ns)) {
    base <- coda_multiway_base(fit, s)
    Zs <- array(base$Z, c(np, h, n_draws))
    open <- NULL
    if (!is.null(base$open)) open <- matrix(base$open, h, n_draws)
    for (q in seq_len(fit$rank)) {
      Zs <- Zs + outer(fit$load[[s]][, q], K[, , q])
      if (!is.null(open)) open <- open + fit$open[[s]]$b[q] * K[, , q]
    }
    out[, , s, ] <- coda_layout_rates(Zs, layout, open)
  }
  out
}

# Tucker3 model of clr life-table deaths Z [part, year, population]
# (Bergeron-Boucher et al. 2018, eq. 5): the centre `alpha` [part,
# population] (clr of the closed geometric mean over the years), and the
# Tucker3 decomposition of the centred array with ranks c(age, time,
# population): factor matrices U (unit-length columns) and core G. Ranks are
# reduced to the dimensions of the data (the parts span at most np - 1
# dimensions after the clr, the years nt - 1 after centring) and so that no
# rank exceeds the product of the other two. Signs: time indexes increase,
# population markers have a positive sum, and the largest core element of
# each age pattern is positive.
coda_tucker <- function(Z, ranks) {
  dims <- dim(Z)
  alpha <- matrix(apply(Z, c(1, 3), mean), dims[1])
  X <- Z - array(alpha[, rep(seq_len(dims[3]), each = dims[2])], dims)
  ss <- sum(X^2)
  if (!is.finite(ss) || ss <= .Machine$double.eps) {
    mt_abort("The life-table deaths do not change over the years; the multiway model needs variation.",
             class = "model_fit")
  }
  r <- pmax(pmin(ranks, dims - c(1, 1, 0)), 1)
  repeat {
    lim <- c(r[2] * r[3], r[1] * r[3], r[1] * r[2])
    if (all(r <= lim)) break
    r <- pmin(r, lim)
  }
  dec <- hooi(X, r, maxit = 1000, tol = 1e-14)
  # signs of the components of mode n (sg = -1 or 1); the core absorbs them
  flip <- function(dec, n, sg) {
    dec$U[[n]] <- sweep(dec$U[[n]], 2, sg, "*")
    dec$G <- mode_product(dec$G, diag(sg, length(sg)), n)
    dec
  }
  B <- dec$U[[2]]
  dec <- flip(dec, 2, ifelse(B[nrow(B), ] < B[1, ], -1, 1))
  dec <- flip(dec, 3, ifelse(colSums(dec$U[[3]]) < 0, -1, 1))
  dec <- flip(dec, 1, vapply(seq_len(r[1]), function(p) {
    g <- dec$G[p, , ]
    if (g[which.max(abs(g))] < 0) -1 else 1
  }, 1))
  G <- dec$G
  list(alpha = alpha, U = dec$U, G = G, ranks = r, fit_pct = 100 * sum(G^2) / ss, core_pct = 100 * G^2 / ss,
       resid = X - tucker_reconstruct(G, dec$U), iterations = dec$iterations)
}

# Loadings [part, time index] of population s: sum over p and r of
# g[p, q, r] b[x, p] gamma[s, r]
coda_tucker_loadings <- function(tk, s) {
  A <- tk$U[[1]]
  C <- tk$U[[3]]
  G <- tk$G
  W <- matrix(0, ncol(A), dim(G)[2])
  for (r in seq_len(ncol(C))) W <- W + C[s, r] * matrix(G[, , r], ncol(A))
  A %*% W
}

# Start of the forecast of population s: clr deaths [part] and log rate of
# the open interval at index 0, from the fitted deaths (centre) or the
# observed deaths of the last year
coda_multiway_base <- function(fit, s) {
  obs <- identical(fit$jump_off, "observed")
  kT <- fit$k[nrow(fit$k), ]
  Z <- if (obs) fit$Z_last[, s] - drop(fit$load[[s]] %*% kT) else fit$alpha[, s]
  open <- NULL
  if (!is.null(fit$open)) {
    o <- fit$open[[s]]
    open <- if (obs) o$last - sum(o$b * kT) else o$a
  }
  list(Z = Z, open = open)
}

# Point forecast of population s: clr deaths [part, h], log open rates [h]
# and the indexes [h, rank]
coda_multiway_point <- function(fit, s, h) {
  kf <- coda_index_mean(fit, h)
  base <- coda_multiway_base(fit, s)
  open <- if (!is.null(base$open)) base$open + drop(kf %*% fit$open[[s]]$b)
  list(Z = base$Z + fit$load[[s]] %*% t(kf), open = open, k = kf)
}

# Simulated paths of the time indexes [h, draw, index]: each from its ARIMA
# model with resampled innovations, shifted by the first-step correction of
# its point forecast (coda_index_mean())
coda_multiway_paths <- function(fit, h, n_draws) {
  kf <- coda_index_mean(fit, h)
  out <- array(NA_real_, c(h, n_draws, fit$rank))
  for (q in seq_len(fit$rank)) {
    f <- fit$ts[[q]]
    shift <- kf[, q] - as.vector(index_arima_paths(f, h, 1, innovations = "none"))
    out[, , q] <- index_arima_paths(f, h, n_draws) + shift
  }
  out
}

# Parameters of a multiway CoDa fit, for mt_parameters()
coda_multiway_parameters <- function(fit) {
  f <- fit$fits
  lab <- fit$strata$label
  yrs <- as.character(fit$periods)
  row <- param_rows
  suffix <- function(n) if (n > 1) seq_len(n) else rep("", n)
  r <- f$ranks
  G <- f$G
  common <- c(
    list(row("common", "rank", c("age", "time", "population"), r)),
    lapply(seq_len(r[1]), function(p) row("common", paste0("b", suffix(r[1])[p]), f$parts, f$U[[1]][, p])),
    unlist(lapply(seq_len(r[2]), function(q) {
      tsq <- f$ts[[q]]
      sq <- suffix(r[2])[q]
      list(row("common", paste0("k", sq), yrs, f$U[[2]][, q]),
           row("common", paste0("order_k", sq), c("p", "d", "q"), tsq$order),
           if (length(tsq$coef)) row("common", paste0("arima_k", sq), names(tsq$coef), tsq$coef),
           row("common", paste0("sd_k", sq), "", tsq$sigma))
    }), recursive = FALSE),
    lapply(seq_len(r[3]), function(j) row("common", paste0("gamma", suffix(r[3])[j]), lab, f$U[[3]][, j])),
    list(row("common", "g", apply(arrayInd(seq_along(G), dim(G)), 1, paste, collapse = ","), G),
         row("common", "explained", "", f$fit_pct / 100)))
  own <- lapply(seq_along(lab), function(s) {
    a <- exp(f$alpha[, s])
    o <- f$open[[s]]
    kan <- f$kannisto[[s]]
    rbindlist(list(row(lab[s], "a", f$parts, a / sum(a)),
                   if (!is.null(o)) row(lab[s], "a_open", "", o$a),
                   if (!is.null(o)) row(lab[s], "b_open", suffix(r[2]), o$b),
                   if (!is.null(o)) row(lab[s], "sd_open", "", o$sd),
                   if (!is.null(kan)) row(lab[s], "kannisto_a", yrs, exp(kan[1, ])),
                   if (!is.null(kan)) row(lab[s], "kannisto_b", yrs, kan[2, ])))
  })
  rbindlist(c(common, own))
}
