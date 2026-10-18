# Compositional data (CoDa) models of life-table deaths.
# coda() adapts model.Oeppen(), predict.Oeppen() and get_dx_values(), and
# coda(coherent = TRUE) adapts model.OeppenC() and predict.OeppenC(), of
# MortalityForecast 0.8.1 (https://github.com/mpascariu/MortalityForecast,
# files R/model_Oeppen.R and R/model_OeppenC.R), GPL-3, by Marius D. Pascariu
# and Marie-Pier Bergeron-Boucher. Zero rates are replaced as by the "min"
# method of replace.zeros() (R/utils.R of the same package). Changes: life-table
# deaths come from the death rates through life_table() conventions
# (lt_columns()); the clr transform is computed directly (no compositions
# package); the benchmark of the coherent model is the closed geometric mean of
# the bottom strata of each group (the reference takes a benchmark population
# from the user); the forecast of year T + j uses the index j steps ahead (the
# reference uses j + 1 steps and corrects the first step by its jump-off, which
# cancels for the random walk with drift and the observed jump-off); draws are
# simulated paths of the indexes, with drift uncertainty; the log rate of an
# open age interval is regressed on the indexes (the reference forecasts the
# deaths only). coda_causes() extends the model to the multiple-decrement form
# described by Oeppen (2008).

#' Compositional data (CoDa) model of life-table deaths
#'
#' `coda()` is the Lee-Carter model in compositional form (Oeppen 2008), as
#' formulated by Bergeron-Boucher, Canudas-Romo, Oeppen and Vaupel (2017). It
#' models the life-table distribution of deaths by age (which sums to 1). With
#' `coherent = TRUE` or `common_by`, it is their coherent model for several
#' populations. The code follows the reference implementation, `model.Oeppen()`
#' and `model.OeppenC()` of the R package MortalityForecast 0.8.1 (Pascariu
#' and Bergeron-Boucher), and reproduces its fitted values and point forecasts.
#'
#' The model of one stratum has six steps:
#'
#' 1. the death rates of each year give life-table deaths `d[x, t]` with
#'    radix 1 at the first age ([life_table()] conventions); for a closed age
#'    range, the survivors at its upper bound are an extra part, so that the
#'    parts sum to 1;
#' 2. `a[x]` is the geometric mean of `d[x, t]` over the years, closed to 1;
#' 3. the deaths are centred by perturbation (`d[x, t] / a[x]`, closed) and
#'    transformed by the centred log-ratio (clr);
#' 4. the first component of a singular value decomposition gives
#'    `clr(d[x, t] / a[x]) = b[x] k[t]` (`b` sums to zero);
#' 5. `k[t]` follows a random walk with drift (ARIMA(0,1,0) with drift, the
#'    default of the reference implementation);
#' 6. forecasts are transformed back (inverse clr, perturbation by `a`,
#'    closure) to life-table deaths and then to death rates.
#'
#' With `jump_off = "observed"` (the default, as `jumpchoice = "actual"` of
#' the reference implementation), forecasts are the observed deaths of the
#' last year `T` perturbed by the change of the model:
#' `d[x, T + j]` is proportional to `d_obs[x, T] exp(b[x] (k[T + j] - k[T]))`.
#' With `"fitted"`, forecasts start from the fitted deaths of the last year.
#' The forecast of year `T + j` uses the index `j` steps ahead; the reference
#' implementation uses `j + 1` steps, which makes no difference for the
#' observed jump-off with a random walk with drift.
#'
#' **Coherent model.** With `coherent = TRUE` or `common_by` (for example
#' `"sex"`), the strata of each group (all strata, or each level of
#' `common_by`) share a common component (Bergeron-Boucher et al. 2017):
#' `clr(d[x, t, i] / a[x, i]) = B[x] K[t] + b[x, i] k[t, i]`. The common
#' component `B K` is the first component of the CoDa model of the average
#' distribution of the group, the closed geometric mean of the life-table
#' deaths of its strata (the centre of the compositions). The deviations of
#' each stratum are the first component of the remaining variation. `K`
#' follows a random walk with drift, and `k[t, i]` an AR(1) process with mean
#' (estimated by maximum likelihood), as in the reference implementation
#' (`model.OeppenC()`, defaults `order.B = c(0, 1, 0)` with drift and
#' `order.D = c(1, 0, 0)` without drift), so the deviations stay bounded. If
#' the data carry a hierarchy (see [mt_aggregate_strata()]), the average uses
#' the bottom strata only, as [li_lee()] does.
#'
#' **Open age interval.** Life-table deaths of an open interval are the
#' survivors at its lower bound; they do not determine its death rate, and
#' the reference implementation forecasts the deaths only. Here the log death
#' rate of the open interval gets loadings on the same indexes, estimated by
#' least squares: `log m[open, t] = a_o + b_o k[t]`, and in the coherent model
#' `log m[open, t, i] = a_o[i] + B_o K[t] + b_o[i] k[t, i]`, with `B_o` from the
#' mean log rate of the group. The open interval thus follows the rest of the
#' age pattern (and the jump-off rule). An auxiliary random walk on the
#' observed log rate, as in [tucker3()], accumulates the noise of the few
#' deaths at the highest ages: on French data (ages 100+, forecasts of 10
#' years) a quarter to a half of its paths exceeded a rate of 1.
#'
#' **Zero deaths.** Cells with zero deaths have no logarithm. As
#' `replace.zeros()` of the reference implementation (method `"min"`), zero
#' rates are replaced by half the smallest positive rate of the stratum in
#' the training years; a warning is given when more than 10% of the cells
#' are affected.
#'
#' Draws simulate the indexes: the random walk with drift includes the
#' uncertainty of the drift (normal, with the variance of the mean change),
#' the AR(1) processes their innovations; `a` and `b` are fixed. Draws are
#' expected rates; predictive draws (`mt_forecast(predictive = TRUE)`) add
#' Poisson variation of the deaths given the rates and the exposures, and no
#' other observation noise. If intervals are too narrow on real data, that
#' is a property of the model (see [recalibrated()], [ensemble()], or
#' overdispersed models such as `stanmomo(family = "nb")`).
#' The model weights all ages equally on the clr scale, so the noise of ages
#' with few deaths enters `k`; its intervals can be wider than those of
#' [lee_carter()].
#'
#' @param jump_off `"observed"` (forecasts start from the observed deaths of
#'   the last year) or `"fitted"` (from the fitted deaths).
#' @param coherent Fit the coherent model with a common component per group.
#' @param common_by Stratum variable(s) defining groups with their own common
#'   component (implies `coherent = TRUE`); `NULL` with `coherent = TRUE` puts
#'   all strata in one group.
#' @return A model specification.
#' @references
#' Oeppen J (2008). Coherent forecasting of multiple-decrement life tables:
#' a test using Japanese cause of death data. CODAWORK'08, Universitat de
#' Girona.
#'
#' Bergeron-Boucher M-P, Canudas-Romo V, Oeppen J, Vaupel JW (2017).
#' Coherent forecasts of mortality with compositional data analysis.
#' *Demographic Research* 37(17): 527-566. \doi{10.4054/DemRes.2017.37.17}
#' @seealso [coda_causes()] for causes of death, [lee_carter()] for the model
#'   of death rates, [li_lee()] for its coherent version.
#' @examples
#' x <- mt_example(deaths = "single", years = 1995:2012)
#' fc <- mt_forecast(mt_fit(coda(), x), h = 5, n_draws = 50, seed = 1)
#' fc2 <- mt_forecast(mt_fit(coda(coherent = TRUE), x), h = 5, n_draws = 50, seed = 1)
#' @export
coda <- function(jump_off = c("observed", "fitted"), coherent = FALSE, common_by = NULL) {
  jump_off <- match.arg(jump_off)
  check_flag(coherent, "coherent")
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
  new_mt_spec("coda", list(jump_off = jump_off, coherent = coherent, common_by = common_by), joint = coherent,
              min_periods = if (coherent) 8 else 5,
              label = paste0(label, if (jump_off == "fitted") ", fitted jump-off" else ""))
}

#' @export
fit_model.mt_spec_coda <- function(spec, D, E, ages, periods, strata = NULL, hierarchy = NULL, ...) {
  if (isTRUE(spec$params$coherent)) return(coda_fit_coherent(spec, D, E, ages, strata, hierarchy))
  fit <- coda_fit_one(coda_positive_rates(D / E), ages)
  fit$jump_off <- spec$params$jump_off
  fit$ages <- ages
  fit$fitted <- coda_rates(fit$alpha + outer(fit$b, fit$k), ages, coda_open_fitted(fit$open, fit$k))
  fit
}

#' @export
simulate_model.mt_spec_coda <- function(spec, fit, h, n_draws, ...) {
  if (isTRUE(spec$params$coherent)) return(coda_simulate_coherent(fit, h, n_draws))
  kv <- as.vector(coda_rw_paths(fit, h, n_draws))
  obs <- identical(fit$jump_off, "observed")
  kT <- fit$k[length(fit$k)]
  Zf <- (if (obs) fit$Z_last - fit$b * kT else fit$alpha) + outer(fit$b, kv)
  open <- NULL
  if (!is.null(fit$open)) open <- (if (obs) fit$open$last - fit$open$b * kT else fit$open$a) + fit$open$b * kv
  array(coda_rates(Zf, fit$ages, open), c(length(fit$ages), h, n_draws))
}

# CoDa Lee-Carter fit of one stratum (or group of causes) from its clr
# life-table deaths Z [part, period]; `open_log`: log rates of the open
# interval (NULL for a closed age range)
coda_fit_clr <- function(Z, open_log = NULL) {
  alpha <- rowMeans(Z)
  sv <- coda_svd(Z - alpha)
  dk <- diff(sv$k)
  drift <- mean(dk)
  list(alpha = alpha, b = sv$b, k = sv$k, drift = drift,
       sigma_k = sqrt(sum((dk - drift)^2) / max(length(dk) - 1, 1)), explained = sv$explained,
       Z_last = Z[, ncol(Z)], open = if (!is.null(open_log)) coda_open_fit(open_log, sv$k))
}

# CoDa Lee-Carter fit of one stratum from positive rates [age, period]
coda_fit_one <- function(m, ages) {
  na <- length(ages)
  open <- is.infinite(ages$upper[na])
  if (na + as.integer(!open) < 2) mt_abort("The CoDa model needs at least two age groups.", class = "model_fit")
  out <- coda_fit_clr(matrix(coda_forward(m, ages), ncol = ncol(m)), if (open) log(m[na, ]))
  out$parts <- c(ages$labels, if (!open) "survivors")
  out
}

# First component of the SVD of centred clr data [part, period]: k = d u
# (periods), b = v (parts); the largest entry of b is positive
coda_svd <- function(C) {
  s <- svd(t(C), nu = 1, nv = 1)
  b <- s$v[, 1]
  k <- s$d[1] * s$u[, 1]
  if (b[which.max(abs(b))] < 0) {
    b <- -b
    k <- -k
  }
  list(b = b, k = k, explained = s$d[1]^2 / max(sum(s$d^2), .Machine$double.xmin))
}

# Least-squares loadings of the log rate of the open interval on the period
# index(es) (columns of X), the open-interval analogue of b[x]
coda_open_fit <- function(lo, X) {
  X <- cbind(1, X)
  cf <- qr.coef(qr(X), lo)
  cf[is.na(cf)] <- 0
  res <- lo - X %*% cf
  list(a = cf[[1]], b = unname(cf[-1]), last = lo[length(lo)],
       sd = sqrt(sum(res^2) / max(length(lo) - ncol(X), 1)))
}

coda_open_fitted <- function(open, k) if (is.null(open)) NULL else open$a + open$b * k

# Index paths [h, draw]: random walk with drift and normal innovations,
# including the uncertainty of the drift (as lee_carter())
coda_rw_paths <- function(f, h, n_draws) {
  dk <- diff(f$k)
  matrix(rw_paths(f$k[length(f$k)], f$drift, matrix(dk - f$drift, nrow = 1), h, n_draws, "gaussian",
                  sigma = f$sigma_k, exp_out = FALSE), h, n_draws)
}

# Rates from clr life-table deaths [part, ...]; `open_log`: log rates of the
# open interval (one per column), ignored for a closed age range
coda_rates <- function(Z, ages, open_log = NULL) {
  open <- NULL
  if (is.infinite(ages$upper[length(ages)])) {
    if (is.null(open_log)) mt_abort("The open interval needs its rates.", class = "model")
    open <- exp(open_log)
  }
  coda_backward(Z, ages, open_rates = open)
}

# Zero rates replaced by half the smallest positive rate (replace.zeros(),
# method "min", of MortalityForecast)
coda_positive_rates <- function(m, what = "the stratum") {
  pos <- is.finite(m) & m > 0
  if (!any(pos)) mt_abort(sprintf("No deaths in %s in the training years.", what), class = "model_fit")
  if (!all(pos)) {
    share <- mean(!pos)
    if (share > 0.1) {
      mt_warn(sprintf(paste("%.0f%% of the cells of %s have no deaths; their rates were set to half the smallest",
                            "positive rate. Wider age groups (mt_regroup()) or smoothing (step_smooth()) avoid",
                            "this."), 100 * share, what), class = "sparse")
    }
    m[!pos] <- min(m[pos]) / 2
  }
  m
}

# Coherent CoDa model: a common component per group from the closed
# geometric mean of its (bottom) strata, and AR(1) deviations per stratum
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
  open <- is.infinite(ages$upper[na])
  if (na + as.integer(!open) < 2) mt_abort("The CoDa model needs at least two age groups.", class = "model_fit")
  m <- lapply(seq_len(ns), function(s) {
    coda_positive_rates(matrix(D[, , s] / E[, , s], na), sprintf("stratum %s", labels[s]))
  })
  # clr of the life-table deaths [part, period] and log open rates of each stratum
  Z <- lapply(m, function(ms) matrix(coda_forward(ms, ages), ncol = nt))
  lo <- if (open) matrix(t(vapply(m, function(ms) log(ms[na, ]), numeric(nt))), ns, nt) else NULL
  commons <- lapply(seq_len(max(group)), function(g) {
    use <- group == g & pool
    if (!any(use)) use <- group == g
    # closed geometric mean of the strata: mean of their clr
    out <- coda_fit_clr(Reduce(`+`, Z[use]) / sum(use), if (open) colMeans(lo[use, , drop = FALSE]))
    out$members <- which(use)
    out
  })
  fitted <- array(NA_real_, c(na, nt, ns))
  own <- lapply(seq_len(ns), function(s) {
    cm <- commons[[group[s]]]
    alpha <- rowMeans(Z[[s]])
    common <- outer(cm$b, cm$k)
    sv <- coda_svd(Z[[s]] - alpha - common)
    out <- list(alpha = alpha, b = sv$b, k = sv$k, ar = ar1_mean_fit(sv$k), explained = sv$explained,
                Z_last = Z[[s]][, nt])
    lo_fit <- NULL
    if (open) {
      out$open <- coda_open_fit(lo[s, ] - cm$open$b * cm$k, sv$k)
      out$open$last <- lo[s, nt]
      lo_fit <- out$open$a + cm$open$b * cm$k + out$open$b * sv$k
    }
    fitted[, , s] <<- coda_rates(alpha + common + outer(sv$b, sv$k), ages, lo_fit)
    out
  })
  list(commons = commons, group = group, strata = own, fitted = fitted, ages = ages, jump_off = p$jump_off,
       parts = c(ages$labels, if (!open) "survivors"))
}

coda_simulate_coherent <- function(fit, h, n_draws) {
  na <- length(fit$ages)
  ns <- length(fit$strata)
  obs <- identical(fit$jump_off, "observed")
  Kp <- lapply(fit$commons, function(cm) as.vector(coda_rw_paths(cm, h, n_draws)))
  out <- array(NA_real_, c(na, h, ns, n_draws))
  for (s in seq_len(ns)) {
    g <- fit$group[s]
    cm <- fit$commons[[g]]
    st <- fit$strata[[s]]
    kv <- as.vector(ar1_mean_paths(st$ar, h, n_draws))
    KT <- cm$k[length(cm$k)]
    kT <- st$k[length(st$k)]
    base <- if (obs) st$Z_last - cm$b * KT - st$b * kT else st$alpha
    Zf <- base + outer(cm$b, Kp[[g]]) + outer(st$b, kv)
    open <- NULL
    if (!is.null(st$open)) {
      base_o <- if (obs) st$open$last - cm$open$b * KT - st$open$b * kT else st$open$a
      open <- base_o + cm$open$b * Kp[[g]] + st$open$b * kv
    }
    out[, , s, ] <- coda_rates(Zf, fit$ages, open)
  }
  out
}

# AR(1) with mean by maximum likelihood (stats::arima(), as
# forecast::Arima(order = c(1, 0, 0), method = "ML") in the reference)
ar1_mean_fit <- function(k) {
  n <- length(k)
  last <- k[n]
  if (n < 4 || stats::sd(k) < 1e-10 * (1 + abs(mean(k)))) {
    return(list(phi = 0, mean = mean(k), sigma = if (n > 1) stats::sd(k) else 0, last = last))
  }
  fit <- NULL
  for (method in c("ML", "CSS-ML")) {
    fit <- tryCatch(suppressWarnings(stats::arima(k, order = c(1, 0, 0), method = method)),
                    error = function(e) NULL)
    if (!is.null(fit)) break
  }
  if (is.null(fit)) {
    a <- ar1_fit(k - mean(k))
    return(list(phi = a$phi, mean = mean(k), sigma = a$sigma, last = last))
  }
  cf <- stats::coef(fit)
  # innovation variance with two degrees of freedom less, as forecast::Arima()
  list(phi = unname(cf[["ar1"]]), mean = unname(cf[["intercept"]]),
       sigma = sqrt(sum(stats::residuals(fit)^2) / (n - 2)), last = last)
}

# Paths [h, draw] of an AR(1) process with mean from its last value
ar1_mean_paths <- function(a, h, n_draws) {
  out <- matrix(NA_real_, h, n_draws)
  prev <- rep(a$last, n_draws)
  for (t in seq_len(h)) {
    prev <- a$mean + a$phi * (prev - a$mean) + stats::rnorm(n_draws, 0, a$sigma)
    out[t, ] <- prev
  }
  out
}

#' Multiple-decrement CoDa model of causes of death (experimental)
#'
#' `coda_causes()` forecasts death rates by cause with the multiple-decrement
#' form of the compositional Lee-Carter model, which Oeppen (2008) proposed
#' and tested on Japanese causes of death. Causes add up to all causes in
#' every draw. The model is experimental.
#'
#' Only the abstract of Oeppen (2008) is published. This function implements
#' the natural multiple-decrement form of [coda()]:
#'
#' 1. in each group (the strata that differ only in `var`, for example one
#'    group per sex), the cause-specific rates `m[x, c]` add up to all-cause
#'    rates `m[x]`, which give life-table deaths `d[x]` with radix 1
#'    ([life_table()] conventions); the deaths are split by cause,
#'    `d[x, c] = d[x] m[x, c] / m[x]` (multiple-decrement life table);
#' 2. the composition is the joint distribution of the deaths `d[x, c]` over
#'    age and cause, plus the survivors at the upper bound of a closed age
#'    range;
#' 3. as in [coda()], the composition is centred by its closed geometric
#'    mean over the years and transformed by the centred log-ratio; the first
#'    component of a singular value decomposition gives `b[x, c] k[t]`, and
#'    `k[t]` follows a random walk with drift;
#' 4. forecast deaths `d[x, c]` give all-cause rates from the all-cause
#'    life-table deaths `d[x] = sum_c d[x, c]`, and cause-specific rates
#'    `m[x, c] = m[x] d[x, c] / d[x]`.
#'
#' The causes add up to the all-cause rates in every draw, and the strata of
#' `total_level` get their sum. If the listed causes leave part of the
#' all-cause deaths unexplained (`residual`), that part is a residual cause
#' inside the model and counts towards the all-cause rates; it is not a
#' stratum of the output.
#'
#' One index `k[t]` drives all ages and causes, so the causes cannot follow
#' different trends beyond the pattern `b[x, c]`. The extensions of
#' Kjærgaard, Ergemen, Kallestrup-Lamb, Oeppen and Lindahl-Jacobsen (2019)
#' are not implemented. In a backtest on [mt_example_causes()] the model was
#' less accurate than the deprecated [cause_shares()]; compare it with
#' [reconciled()] separate models in a backtest before using it.
#'
#' Jump-off, the open age interval and zero deaths are treated as in
#' [coda()]: the all-cause log rate of an open interval is regressed on
#' `k[t]` (its split by cause comes from the composition), and zero rates of
#' a cause are replaced by half its smallest positive rate in the group (as
#' `replace.zeros()` of the R package MortalityForecast). Draws are expected
#' rates. Predictive draws (`mt_forecast(predictive = TRUE)`) draw Poisson
#' all-cause deaths given the all-cause rates and exposures and split them
#' multinomially between the causes by the cause shares of the draw, so that
#' the causes also add up to all causes in predictive draws; no other
#' observation noise is added. If intervals are too narrow on real data,
#' that is a property of the model (see [recalibrated()], [ensemble()], or
#' overdispersed models such as `stanmomo(family = "nb")`).
#'
#' @param var Stratum variable with the causes.
#' @param total_level Level of `var` for all causes.
#' @param jump_off `"observed"` (forecasts start from the observed deaths of
#'   the last year) or `"fitted"` (from the fitted deaths).
#' @param residual Add a residual cause, the part of all-cause mortality not
#'   covered by the listed causes: `TRUE`, `FALSE`, or `"auto"` (added when
#'   the listed causes leave more than 0.1% of all-cause mortality
#'   unexplained in the median cell).
#' @return A model specification.
#' @references
#' Oeppen J (2008). Coherent forecasting of multiple-decrement life tables:
#' a test using Japanese cause of death data. CODAWORK'08, Universitat de
#' Girona.
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
coda_causes <- function(var = "cause", total_level = "all", jump_off = c("observed", "fitted"), residual = "auto") {
  jump_off <- match.arg(jump_off)
  if (!is.character(var) || length(var) != 1L) mt_abort("`var` must be one stratum variable.", class = "argument")
  if (length(total_level) != 1L) mt_abort("`total_level` must be one level.", class = "argument")
  if (!identical(residual, "auto") && !(is.logical(residual) && length(residual) == 1L && !is.na(residual))) {
    mt_abort("`residual` must be TRUE, FALSE or \"auto\".", class = "argument")
  }
  new_mt_spec("coda_causes", list(var = var, total_level = total_level, jump_off = jump_off, residual = residual),
              joint = TRUE, min_periods = 5,
              label = paste0("CoDa causes (multiple decrement)", if (jump_off == "fitted") ", fitted jump-off" else ""))
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
  open <- is.infinite(ages$upper[na])
  # rates of the causes [age, period, cause] and of all causes, per group
  mc <- lapply(seq_along(groups), function(g) {
    array(D[, , part_id[g, ], drop = FALSE] / E[, , part_id[g, ], drop = FALSE], c(na, nt, K))
  })
  mt <- lapply(seq_along(groups), function(g) matrix(D[, , tot_id[g]] / E[, , tot_id[g]], na))
  short <- unlist(lapply(seq_along(groups), function(g) (1 - rowSums(mc[[g]], dims = 2) / mt[[g]])[mt[[g]] > 0]))
  use_residual <- if (identical(p$residual, "auto")) stats::median(short) > 1e-3 else isTRUE(p$residual)
  nc <- K + use_residual
  fitted <- array(NA_real_, c(na, nt, ns))
  fits <- lapply(seq_along(groups), function(g) {
    m <- mc[[g]]
    if (use_residual) m <- array(c(m, pmax(mt[[g]] - rowSums(m, dims = 2), 0)), c(na, nt, nc))
    for (j in seq_len(nc)) {
      what <- sprintf("cause %s of group %s", if (j <= K) causes[j] else "residual", groups[g])
      m[, , j] <- coda_positive_rates(matrix(m[, , j], na), what)
    }
    out <- coda_fit_clr(md_forward(m, ages), if (open) log(rowSums(m, dims = 2)[na, ]))
    out$jump_off <- p$jump_off
    R <- md_rates(out$alpha + outer(out$b, out$k), ages, nc, coda_open_fitted(out$open, out$k))
    for (j in seq_len(K)) fitted[, , part_id[g, j]] <<- R[, , j]
    fitted[, , tot_id[g]] <<- md_total(R)
    out
  })
  parts <- c(paste(rep(ages$labels, nc), rep(c(causes, if (use_residual) "residual"), each = na)),
             if (!open) "survivors")
  list(groups = fits, group_labels = groups, tot_id = tot_id, part_id = part_id, causes = causes, K = K,
       residual = use_residual, parts = parts, ages = ages, fitted = fitted)
}

#' @export
simulate_model.mt_spec_coda_causes <- function(spec, fit, h, n_draws, ...) {
  na <- length(fit$ages)
  ns <- length(c(fit$tot_id, fit$part_id))
  out <- array(NA_real_, c(na, h, ns, n_draws))
  nc <- fit$K + fit$residual
  for (g in seq_along(fit$groups)) {
    gf <- fit$groups[[g]]
    kv <- as.vector(coda_rw_paths(gf, h, n_draws))
    obs <- identical(gf$jump_off, "observed")
    kT <- gf$k[length(gf$k)]
    Zf <- (if (obs) gf$Z_last - gf$b * kT else gf$alpha) + outer(gf$b, kv)
    open <- NULL
    if (!is.null(gf$open)) open <- (if (obs) gf$open$last - gf$open$b * kT else gf$open$a) + gf$open$b * kv
    R <- md_rates(Zf, fit$ages, nc, open)
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
# [part, column]; all-cause rates from the all-cause life-table deaths, the
# open interval from `open_log` (log all-cause rates)
md_rates <- function(Z, ages, nc, open_log = NULL) {
  na <- length(ages)
  Zm <- matrix(Z, nrow(Z))
  P <- exp(Zm)
  P <- sweep(P, 2, colSums(P), "/")
  dxc <- array(P[seq_len(na * nc), , drop = FALSE], c(na, nc, ncol(P)))
  dx <- matrix(dxc[, 1, ], na)
  for (j in seq_len(nc)[-1]) dx <- dx + dxc[, j, ]
  lz <- log(if (nrow(P) > na * nc) rbind(dx, P[nrow(P), ]) else dx)
  m <- matrix(coda_rates(lz, ages, open_log), na)
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
# loadings `b` (parts of the composition: ages, ages by cause, survivors) and
# the indexes, with the drift and standard deviation of the random walks, the
# AR(1) coefficients of the coherent model and the loadings of the open interval
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
    out <- rbind(row(who, "a", parts, centre(x$alpha)), row(who, b_name, parts, x$b), row(who, idx, yrs, x$k),
                 row(who, paste0("drift_", idx), "", x$drift), row(who, paste0("sd_", idx), "", x$sigma_k),
                 row(who, "explained", "", x$explained))
    if (!is.null(x$open)) out <- rbind(out, coda_open_rows(row, who, x$open, b_name))
    out
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
  own <- rbindlist(lapply(seq_along(f$strata), function(s) {
    st <- f$strata[[s]]
    out <- rbind(row(lab[s], "a", f$parts, centre(st$alpha)), row(lab[s], "b", f$parts, st$b),
                 row(lab[s], "k", yrs, st$k), row(lab[s], "phi", "", st$ar$phi),
                 row(lab[s], "mean_k", "", st$ar$mean), row(lab[s], "sd_k", "", st$ar$sigma),
                 row(lab[s], "explained", "", st$explained))
    if (!is.null(st$open)) out <- rbind(out, coda_open_rows(row, lab[s], st$open, "b"))
    out
  }))
  rbind(common, own)
}

coda_open_rows <- function(row, who, open, b_name) {
  rbind(row(who, "a_open", "", open$a), row(who, paste0(b_name, "_open"), "", open$b),
        row(who, "sd_open", "", open$sd))
}

# Statement of the model and its uncertainty, for summary()
coda_assumptions <- function(spec) {
  p <- spec$params
  open <- paste("The log rate of an open age interval is regressed on the same indexes.")
  if (inherits(spec, "mt_spec_coda_causes")) {
    return(paste("Multiple-decrement CoDa model (Oeppen 2008; experimental): per group, clr of the life-table",
                 "deaths by age and cause = a + b k(t), first SVD component; k(t) a random walk with drift (draws",
                 "include innovations and drift uncertainty; a, b fixed). Causes add up to all causes in every",
                 "draw; predictive draws split Poisson all-cause deaths multinomially. Jump-off:", p$jump_off,
                 "deaths.", open))
  }
  if (!isTRUE(p$coherent)) {
    return(paste("CoDa Lee-Carter model (Oeppen 2008; Bergeron-Boucher et al. 2017): clr of the life-table deaths",
                 "= clr(a) + b k(t), first SVD component; k(t) a random walk with drift (draws include",
                 "innovations and drift uncertainty; a, b fixed). Jump-off:", p$jump_off, "deaths.", open))
  }
  paste("Coherent CoDa model (Bergeron-Boucher et al. 2017): clr of the life-table deaths = clr(a_i) + B K(t) +",
        "b_i k_i(t), with B K from the closed geometric mean of the strata of each group",
        if (!is.null(p$common_by)) sprintf("(groups by %s)", paste(p$common_by, collapse = ", ")) else "",
        "; K(t) a random walk with drift, k_i(t) AR(1) with mean, so deviations stay bounded. Draws include the",
        "innovations, the drift uncertainty of K; a, B, b and the AR coefficients are fixed. Jump-off:",
        p$jump_off, "deaths.", open)
}
