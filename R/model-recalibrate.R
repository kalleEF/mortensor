#' Recalibrate forecast distributions on recent years
#'
#' Wraps any model specification so that the spread (and optionally the
#' trend) of its forecast draws is recalibrated on the last `years` observed
#' years of the training data. Each time the model is fitted (also at every
#' origin of an outer [mt_backtest()]):
#'
#' 1. the base model is refitted at `n_origins` origins inside the
#'    calibration window (training data up to the origin) and forecasts the
#'    remaining years of the window;
#' 2. the draws are transformed per cell on the log scale,
#'    \deqn{\log m^* = \tilde m + b\,h + s(h)\,(\log m - \tilde m),}
#'    where \eqn{\tilde m} is the median of the draws in the cell, \eqn{h} the
#'    horizon, \eqn{s(h) = s_1 h^{\gamma}} a spread factor (`scale =
#'    "horizon"`; \eqn{\gamma = 0} for `scale = "constant"`), constant
#'    beyond the longest horizon of the window, and \eqn{b} a correction of
#'    the drift per year (`bias = TRUE`);
#' 3. \eqn{s_1}, \eqn{\gamma} and \eqn{b} minimise the mean CRPS (or the
#'    interval score, `target`) of the log death rates in the window, with
#'    Poisson (and the model's extra-Poisson) variation of the observed
#'    deaths, as in [mt_backtest()]. Draws of an ensemble get the observation
#'    noise of the component they come from.
#'
#' The final model is fitted on all training data and its forecasts are
#' transformed with the estimated parameters. This is minimum-CRPS
#' estimation of a location-scale recalibration (Gneiting et al. 2005)
#' applied to simulated paths: ranks within each cell are kept, so paths stay
#' paths, and derived quantities such as life expectancy can be computed
#' from the calibrated draws.
#'
#' The wrapper works with any specification, including [ensemble()],
#' [bootstrapped()] and joint models. For [reconciled()] models the bottom
#' strata are calibrated and the aggregates recomputed from them, so the
#' draws stay coherent. With `by`, separate parameters are estimated for
#' groups of strata (for example `by = "sex"`); strata outside these groups
#' (aggregates) use the pooled parameters. Groups with different factors
#' scale the common variation of joint models (such as the common trend of
#' [li_lee()]) differently, so it no longer cancels in contrasts between the
#' groups, and intervals of gaps between them widen; use pooled parameters
#' when contrasts between the groups are the target.
#'
#' Calibration needs `years` more training years than the base model. A
#' spread factor above 1 means that the base model was too narrow in the
#' window; it widens future intervals in proportion. The parameters are
#' estimated from few overlapping forecasts (with the defaults, 5 origins and
#' 30 forecast years, horizons 1 to 10), so they carry their own uncertainty,
#' and they cannot anticipate changes after the window. `n_boot` refits them
#' on forecasts (stratum x origin) resampled with replacement and reports
#' 95% percentile intervals. Beyond the longest horizon of the window, the
#' spread factor and the drift correction stay constant. The inner models are
#' fitted to shorter series than the final model, which matters for models
#' whose uncertainty depends on the length of the series (the drift of a
#' random walk). Compare calibrated and uncalibrated models in an outer
#' backtest.
#'
#' **Age bands.** With `ages` (break points, for example `c(50, 75)`), the
#' factors are estimated separately for each age band, so that a model that
#' is too narrow at old ages and too wide at young ages is corrected at
#' both. Observed cells belong to the band of their lower age; use break
#' points that are bounds of the observed age groups.
#'
#' **Quantity.** By default the log death rates of single cells are scored.
#' Intervals of life expectancy and ASR also depend on how the ages move
#' together, which the transformation keeps from the base model;
#' `quantity = "life_expectancy"` (at birth) or `"asr"` scores these instead,
#' per test year and stratum (fewer scored values, so the factors are less
#' certain).
#'
#' **Tuning and recalibration.** `recalibrated()` keeps the model and its
#' settings and rescales the width of its forecast distribution. Choosing
#' settings (smoothing, number of components, fitting window) by forecast
#' accuracy is tuning, done by [tuned()] and [mt_tune()]; the two combine
#' as `recalibrated(tuned(...))`.
#'
#' Calibration widens the intervals of models that were too narrow in the
#' recent past and narrows those that were too wide, for example of
#' ensembles whose intervals already covered. With few training years it
#' rests on few forecasts and can make intervals worse; check it in a
#' backtest, as in the vignette "Tuning and recalibration".
#'
#' `target = "crps"` (default) balances calibration and sharpness;
#' `target = "interval"` minimises the interval score of the central `level`
#' interval (Gneiting and Raftery 2007), whose optimum is the correct
#' coverage of that interval.
#'
#' @param spec A model specification.
#' @param years Length of the calibration window (the last observed years).
#' @param n_origins Number of origins in the window (evenly spaced from
#'   `years` to 2 years before the end of the training data).
#' @param scale `"horizon"` (spread factor \eqn{s_1 h^{\gamma}}) or
#'   `"constant"`.
#' @param bias Estimate a correction of the drift (log rate per year).
#' @param target `"crps"` or `"interval"` (interval score at `level`).
#' @param level Level of the central interval for `target = "interval"`.
#' @param by Stratum variables with separate calibration parameters.
#' @param n_draws Draws per inner forecast.
#' @param seed Random seed of the inner forecasts and of the simulated
#'   observation noise.
#' @param ages Break points between age bands with their own factors
#'   (`NULL`: one set for all ages).
#' @param quantity What is scored: `"rate"` (log death rates of the observed
#'   cells), `"life_expectancy"` (at birth) or `"asr"`.
#' @param n_boot Number of bootstrap refits for intervals of the factors (0:
#'   none; each refit repeats the optimisation, not the inner forecasts).
#' @return A model specification. The fitted object has a `calibration`
#'   element: a `data.table` with one row per group and age band
#'   (`scale_h1`, `slope`, `drift`, the number of observed cells and of
#'   scored values, the number of inner forecasts (stratum x origin), CRPS
#'   and coverage of the inner forecasts before and after calibration, and
#'   with `n_boot` the intervals `scale_h1_lo`, `scale_h1_hi`, ...).
#' @references
#' Gneiting T, Raftery AE, Westveld AH, Goldman T (2005). Calibrated
#' probabilistic forecasting using ensemble model output statistics and
#' minimum CRPS estimation. *Monthly Weather Review* 133(5): 1098-1118.
#' \doi{10.1175/MWR2904.1}
#'
#' Gneiting T, Raftery AE (2007). Strictly proper scoring rules, prediction,
#' and estimation. *Journal of the American Statistical Association*
#' 102(477): 359-378. \doi{10.1198/016214506000001437}
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2019)
#' fit <- mt_fit(recalibrated(lee_carter(), years = 8, n_origins = 3, n_draws = 50), x)
#' fit$calibration
#' fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
#' @export
recalibrated <- function(spec, years = 10, n_origins = 5, scale = c("horizon", "constant"), bias = FALSE,
                       target = c("crps", "interval"), level = 0.95, by = NULL, n_draws = 100, seed = 1,
                       ages = NULL, quantity = c("rate", "life_expectancy", "asr"), n_boot = 0) {
  check_class(spec, "mt_spec", "spec")
  scale <- match.arg(scale)
  target <- match.arg(target)
  quantity <- match.arg(quantity)
  check_number(level, "level", lower = 0.5, upper = 0.99)
  check_number(years, "years", lower = 2, integerish = TRUE)
  check_number(n_origins, "n_origins", lower = 1, integerish = TRUE)
  check_number(n_draws, "n_draws", lower = 10, integerish = TRUE)
  check_number(n_boot, "n_boot", lower = 0, integerish = TRUE)
  check_flag(bias, "bias")
  if (!is.null(ages)) {
    if (!is.numeric(ages) || anyNA(ages) || is.unsorted(ages, strictly = TRUE)) {
      mt_abort("`ages` must be increasing break points between age bands, for example c(50, 75).",
               class = "argument")
    }
    if (quantity != "rate") {
      mt_abort("Age bands (`ages`) need `quantity = \"rate\"`: life expectancy and ASR combine all ages.",
               class = "argument")
    }
  }
  if (n_origins > years - 1) n_origins <- years - 1
  new_mt_spec("recalibrated", list(base = spec, years = years, n_origins = n_origins, scale = scale, bias = bias,
                                 target = target, level = level, by = by, n_draws = n_draws, seed = seed,
                                 ages = ages, quantity = quantity, n_boot = n_boot),
              joint = spec$joint, min_periods = spec$min_periods + years, smooth = spec$smooth,
              label = paste0(spec$label, ", recalibrated", if (length(by)) paste0(" by ", paste(by, collapse = " x ")),
                             if (length(ages)) " by age band",
                             if (quantity != "rate") sprintf(" on %s", gsub("_", " ", quantity)),
                             if (target == "interval") sprintf(" (%g%% interval)", 100 * level)))
}

#' @export
mt_fit.mt_spec_recalibrated <- function(spec, data, prep = NULL, cores = 1, ...) {
  p <- spec$params
  g <- mt_grid(data)
  last <- max(g$period$lower)
  first_ok <- min(g$period$lower) + p$base$min_periods - 1
  origins <- calibration_origins(last, p$years, p$n_origins)
  if (any(origins < first_ok)) {
    origins <- origins[origins >= first_ok]
    if (!length(origins)) {
      mt_abort(sprintf(paste("Calibration on the last %d years needs at least %d training years; give a shorter",
                             "`years`."), p$years, p$base$min_periods + p$years), class = "model_data")
    }
    mt_warn(sprintf("Only %d calibration origin(s) leave enough training years.", length(origins)),
            class = "calibration")
  }
  fit <- mt_fit(p$base, data, prep = prep, cores = cores, ...)
  bottom <- calibration_targets(fit)
  pieces <- par_lapply(origins, function(o) {
    tryCatch(calibration_origin(p$base, data, o, last, prep, p$n_draws, if (!is.null(p$seed)) p$seed + o,
                                fit$strata$label[bottom], p$quantity %||% "rate"), error = function(e) e)
  }, cores, streams = TRUE)
  bad <- vapply(pieces, inherits, logical(1), what = "error")
  if (any(bad)) {
    mt_warn(sprintf("Calibration origin(s) %s failed: %s", paste(origins[bad], collapse = ", "),
                    conditionMessage(pieces[[which(bad)[1]]])), class = "calibration")
  }
  pieces <- unlist(pieces[!bad], recursive = FALSE)
  if (!length(pieces)) mt_abort("Every calibration origin failed.", class = "calibration")
  if (sum(!bad) == 1) {
    mt_warn(paste("The calibration rests on the forecasts of one origin: its horizons are not replicated and the",
                  "factors are uncertain. Give a longer `years` or more `n_origins`."), class = "calibration")
  }
  est <- estimate_calibration(pieces, fit$strata, p)
  wrap_fit(fit, spec, "mt_fit_recalibrated", list(calibration = est$table, calibration_theta = est$theta,
                                                  calibration_origins = origins[!bad],
                                                  calibration_hmax = last - min(origins[!bad])))
}

calibration_origins <- function(last, years, n_origins) {
  unique(last - years + floor((seq_len(n_origins) - 1) * years / n_origins))
}

# Strata whose draws are transformed: the bottom strata of reconciled fits,
# otherwise all
calibration_targets <- function(fit) {
  hs <- fit$reconcile$structure
  if (is.null(hs)) seq_len(nrow(fit$strata)) else hs$bottom
}

# Does the predictive distribution add Poisson and extra-Poisson noise to the
# rate draws (models of smooth rates), or are the draws observed rates?
spec_adds_noise <- function(spec) {
  if (inherits(spec, "mt_spec_ensemble")) return(any(vapply(spec$params$components, spec_adds_noise, TRUE)))
  if (!is.null(spec$params$base) && inherits(spec$params$base, "mt_spec")) return(spec_adds_noise(spec$params$base))
  isTRUE(spec$smooth)
}

# The same for a fit of a wrapper or a single model: the innermost of the
# fitted specs decides (tuned() fits the spec it chose, which can differ from
# the spec the wrapper was given)
fit_adds_noise <- function(fit) {
  specs <- c(list(fit$spec), fit$inner_specs)
  # any model fitted to smoothed rates (step_smooth()) forecasts smooth rates
  spec_adds_noise(specs[[length(specs)]]) || isTRUE(fit$smoothed_input)
}

# Observation noise of each draw: whether Poisson noise is added
# (`noise(s)`: one flag per draw for stratum s) and the extra-Poisson sd by
# age (`sd(s)`: ages x draws). Draws of an ensemble follow the component
# they come from (`fc$component`: one per draw, or [stratum, draw] after
# reordering; without it, draws are allocated to the components in
# proportion to the weights); strata of a fit split into groups
# (stratified(), tuned(by = )) follow the fit of their group.
draw_noise <- function(fit, fc, n_draws) {
  na <- length(fit$age)
  comps <- fit$components
  if (!is.null(comps)) {
    parts <- lapply(comps, draw_noise, fc = NULL, n_draws = n_draws)
    k <- fc$component
    if (is.null(k)) {
      w <- fit$weights %||% rep(1, length(comps))
      k <- findInterval((seq_len(n_draws) - 0.5) / n_draws, cumsum(w) / sum(w)) + 1L
    }
    of <- function(s) if (is.matrix(k)) k[s, ] else k
    return(list(noise = function(s) {
      ks <- of(s)
      out <- logical(n_draws)
      for (j in unique(ks)) out[ks == j] <- parts[[j]]$noise(s)[ks == j]
      out
    }, sd = function(s) {
      ks <- of(s)
      out <- matrix(0, na, n_draws)
      for (j in unique(ks)) out[, ks == j] <- parts[[j]]$sd(s)[, ks == j, drop = FALSE]
      out
    }))
  }
  if (!is.null(fit$group_fits)) {
    grp <- pos <- integer(nrow(fit$strata))
    for (g in seq_along(fit$groups)) {
      grp[fit$groups[[g]]] <- g
      pos[fit$groups[[g]]] <- seq_along(fit$groups[[g]])
    }
    parts <- lapply(fit$group_fits, draw_noise, fc = NULL, n_draws = n_draws)
    return(list(noise = function(s) parts[[grp[s]]]$noise(pos[s]), sd = function(s) parts[[grp[s]]]$sd(pos[s])))
  }
  flag <- fit_adds_noise(fit)
  sd <- fit$noise_sd
  list(noise = function(s) rep(flag, n_draws), sd = function(s) {
    matrix(if (flag && !is.null(sd)) sd[, s] else 0, na, n_draws)
  })
}

# One inner origin: forecast draws of the rest of the window and the observed
# cells they are scored against, per stratum
calibration_origin <- function(spec, data, origin, last, prep, n_draws, seed, targets, quantity = "rate") {
  train <- mt_window(data, to = origin)
  fit <- mt_fit(spec, train, prep = prep)
  last_fit <- max(fit$periods)
  test <- mt_window(data, from = last_fit + 1, to = last)
  test <- mt_set_grid(test, age = fit$age)
  E_grid <- test_exposures(test, fit)
  fc <- suppressWarnings(mt_forecast(fit, to = last, n_draws = n_draws, seed = seed, exposure = E_grid))
  obs <- test$obs$deaths
  obs <- obs[period_lo > origin & provenance != "estimated"]
  gt <- mt_grid(test)
  idx <- obs_grid_index(obs, gt, "deaths", strict = TRUE)
  ok <- !is.na(idx$a1)
  obs <- obs[ok]
  idx <- idx[ok]
  na <- length(gt$age)
  np <- length(gt$period)
  nz <- draw_noise(fit, fc, n_draws)
  out <- list()
  for (s in unique(obs$stratum)) {
    lab <- data$strata$label[s]
    sf <- match(lab, fit$strata$label)
    if (is.na(sf) || !lab %in% targets) next
    R <- matrix(fc$rates[, , sf, ], ncol = n_draws)
    if (anyNA(R)) next
    sel <- obs$stratum == s
    C <- build_operator(idx[sel], na, np)
    used <- which(Matrix::colSums(C) > 0)
    C <- C[, used, drop = FALSE]
    # exposures of the test data (strata of `data`; a reconciled fit of
    # bottom strata also forecasts the aggregates)
    Eg <- as.vector(E_grid[, , s])[used]
    Ecell <- as.numeric(C %*% Eg)
    o <- obs[sel]
    L <- log(pmax(R[used, , drop = FALSE], 1e-12))
    out[[length(out) + 1L]] <- list(
      label = lab, origin = origin, L = L, m = row_quantiles(L, 0.5)[, 1],
      h = rep(fc$periods - origin, each = na)[used],
      E = Eg, C = C, Ecell = Ecell, y = log(pmax(o$value, 0.5) / Ecell), noise = nz$noise(sf),
      sd = nz$sd(sf)[rep(seq_len(na), np)[used], , drop = FALSE], obs_age = o$age_lo,
      qg = if (quantity != "rate") quantity_groups(o, Ecell, quantity))
  }
  out
}

# Observed cells of complete age axes per test year, grouped by axis: column
# j of K holds the cells of one year (sorted by age), qo its observed life
# expectancy or log ASR
quantity_groups <- function(o, Ecell, quantity) {
  out <- list()
  for (pp in sort(unique(o$period_lo))) {
    k <- which(o$period_lo == pp)
    k <- k[order(o$age_lo[k])]
    lo <- o$age_lo[k]
    hi <- o$age_hi[k]
    if (length(k) < 2 || any(abs(hi[-length(hi)] - lo[-1]) > 1e-9)) next
    key <- paste(lo, hi, collapse = ",")
    if (is.null(out[[key]])) {
      ax <- mt_axis(lo, hi, kind = "age")
      out[[key]] <- list(K = NULL, qo = numeric(), ax = ax, quantity = quantity,
                         w = if (quantity == "asr") standard_weights(ax, standard_population("esp2013")))
    }
    r_obs <- o$value[k] / Ecell[k]
    g <- out[[key]]
    qo <- if (quantity == "life_expectancy") life_expectancy(r_obs, g$ax) else log(sum(r_obs * g$w) / sum(g$w))
    out[[key]]$K <- cbind(g$K, k)
    out[[key]]$qo <- c(g$qo, qo)
  }
  unname(out)
}

# The cells of a piece in each age band (by the lower age of the observed
# cell; grid cells follow the observed cell they belong to)
split_bands <- function(piece, breaks) {
  bo <- findInterval(piece$obs_age, breaks) + 1L
  bg <- as.vector(Matrix::t(piece$C) %*% bo)
  lapply(stats::setNames(sort(unique(bo)), band_labels(breaks)[sort(unique(bo))]), function(b) {
    rows <- which(bo == b)
    cols <- which(bg == b)
    list(label = piece$label, origin = piece$origin, L = piece$L[cols, , drop = FALSE], m = piece$m[cols],
         h = piece$h[cols], E = piece$E[cols], C = piece$C[rows, cols, drop = FALSE], Ecell = piece$Ecell[rows],
         y = piece$y[rows], noise = piece$noise, sd = piece$sd[cols, , drop = FALSE], obs_age = piece$obs_age[rows])
  })
}

band_labels <- function(breaks) {
  if (!length(breaks)) return("all ages")
  c(paste0("<", breaks[1]), if (length(breaks) > 1) paste0(breaks[-length(breaks)], "-", breaks[-1]),
    paste0(breaks[length(breaks)], "+"))
}

# Stack the pieces of a group and draw the fixed random numbers of the
# observation noise (common random numbers across parameter values)
stack_pieces <- function(pieces, n_draws) {
  N <- vapply(pieces, function(p) length(p$y), 0)
  P <- list(L = do.call(rbind, lapply(pieces, `[[`, "L")), m = unlist(lapply(pieces, `[[`, "m")),
            h = unlist(lapply(pieces, `[[`, "h")), E = unlist(lapply(pieces, `[[`, "E")),
            C = Matrix::bdiag(lapply(pieces, `[[`, "C")), Ecell = unlist(lapply(pieces, `[[`, "Ecell")),
            y = unlist(lapply(pieces, `[[`, "y")))
  # Poisson noise per observed cell and draw (draws of an ensemble differ)
  P$pois <- do.call(rbind, lapply(seq_along(pieces), function(i) {
    matrix(pieces[[i]]$noise, N[i], n_draws, byrow = TRUE)
  }))
  sd <- do.call(rbind, lapply(pieces, `[[`, "sd"))
  P$extra <- if (any(sd > 0)) matrix(exp(stats::rnorm(length(sd), 0, sd)), nrow(sd)) else 1
  P$U <- matrix(stats::runif(sum(N) * n_draws), sum(N))
  P$Z <- stats::qnorm(P$U)
  # life expectancy or ASR: cell indices shifted to the stacked rows
  off <- cumsum(c(0, N))
  qg <- unlist(lapply(seq_along(pieces), function(i) {
    lapply(pieces[[i]]$qg, function(g) {
      g$K <- g$K + off[i]
      g
    })
  }), recursive = FALSE)
  if (length(qg)) P$qg <- qg
  P
}

# Forecasts and observations that are scored: log rates of the observed
# cells, or life expectancy / log ASR per test year and stratum
calibration_scored <- function(X, P) {
  if (is.null(P$qg)) return(list(S = X, obs = P$y))
  S <- do.call(rbind, lapply(P$qg, function(g) {
    M <- matrix(exp(X[as.vector(g$K), , drop = FALSE]), nrow = nrow(g$K))
    q <- if (g$quantity == "life_expectancy") lt_columns(M, g$ax)$ex[1, ] else log(colSums(M * g$w) / sum(g$w))
    matrix(q, ncol(g$K))
  }))
  list(S = S, obs = unlist(lapply(P$qg, `[[`, "qo")))
}

# Spread factor and drift correction per grid cell; beyond the longest
# horizon of the calibration window both stay constant (the drift correction
# is a level shift there, not extrapolated further)
calibration_shift <- function(theta, h, scale, bias, h_max = Inf) {
  hc <- pmin(h, h_max)
  s <- exp(theta[1] + if (scale == "horizon") theta[2] * log(hc) else 0)
  b <- if (bias) theta[length(theta)] * hc else 0
  list(s = s, b = b)
}

# Log rates of the observed cells under parameters theta, with observation noise
calibration_sample <- function(theta, P, scale, bias) {
  sh <- calibration_shift(theta, P$h, scale, bias)
  Mu <- exp(P$m + sh$b + sh$s * (P$L - P$m)) * P$E * P$extra
  Lam <- as.matrix(P$C %*% Mu)
  D <- Lam
  small <- P$pois & Lam < 50
  large <- P$pois & !small
  D[small] <- stats::qpois(P$U[small], Lam[small])
  D[large] <- pmax(Lam[large] + sqrt(Lam[large]) * P$Z[large], 0)
  log(pmax(D, 0.5) / P$Ecell)
}

crps_rows <- function(X, y) {
  n <- ncol(X)
  S <- matrix(X[order(rep(seq_len(nrow(X)), n), X)], ncol = n, byrow = TRUE)
  rowMeans(abs(X - y)) - as.numeric(S %*% (2 * seq_len(n) - n - 1)) / n^2
}

# Mean interval score of the central interval at `level` (Gneiting and Raftery 2007)
interval_score <- function(X, y, level) {
  a <- 1 - level
  Q <- row_quantiles(X, c(a / 2, 1 - a / 2))
  mean(Q[, 2] - Q[, 1] + 2 / a * pmax(Q[, 1] - y, 0) + 2 / a * pmax(y - Q[, 2], 0))
}

calibration_fit <- function(P, scale, bias, target = "crps", level = 0.95, summary = TRUE) {
  k <- 1 + (scale == "horizon") + bias
  lower <- c(log(0.1), if (scale == "horizon") -1, if (bias) -0.05)
  upper <- c(log(20), if (scale == "horizon") 1, if (bias) 0.05)
  f <- function(theta) {
    if (any(theta < lower | theta > upper)) return(1e10)
    sc <- calibration_scored(calibration_sample(theta, P, scale, bias), P)
    if (target == "crps") mean(crps_rows(sc$S, sc$obs)) else interval_score(sc$S, sc$obs, level)
  }
  theta <- if (k == 1) {
    stats::optimize(f, c(lower, upper))$minimum
  } else {
    stats::optim(rep(0, k), f, method = "Nelder-Mead", control = list(reltol = 1e-6, maxit = 300))$par
  }
  if (!summary) return(list(theta = theta))
  summ <- function(th) {
    sc <- calibration_scored(calibration_sample(th, P, scale, bias), P)
    Q <- row_quantiles(sc$S, c(0.025, 0.1, 0.9, 0.975))
    c(crps = mean(crps_rows(sc$S, sc$obs)), cov80 = mean(sc$obs >= Q[, 2] & sc$obs <= Q[, 3]),
      cov95 = mean(sc$obs >= Q[, 1] & sc$obs <= Q[, 4]), n = length(sc$obs))
  }
  b0 <- summ(rep(0, k))
  b1 <- summ(theta)
  list(theta = theta, table = data.table(
    scale_h1 = exp(theta[1]), slope = if (scale == "horizon") theta[2] else 0,
    drift = if (bias) theta[k] else 0, cells = length(P$y), scored = b0[["n"]], crps_before = b0[["crps"]],
    crps_after = b1[["crps"]], cov80_before = b0[["cov80"]], cov80_after = b1[["cov80"]],
    cov95_before = b0[["cov95"]], cov95_after = b1[["cov95"]]))
}

# Percentile intervals of the factors from refits on pieces (stratum x
# origin) resampled with replacement
calibration_boot <- function(pieces, p, n_draws) {
  th <- t(vapply(seq_len(p$n_boot), function(r) {
    gp <- pieces[sample.int(length(pieces), replace = TRUE)]
    calibration_fit(stack_pieces(gp, n_draws), p$scale, p$bias, p$target %||% "crps", p$level %||% 0.95,
                    summary = FALSE)$theta
  }, numeric(1 + (p$scale == "horizon") + p$bias)))
  if (nrow(th) == 1L && length(th) > 1L) th <- t(th)
  q <- function(v) stats::quantile(v, c(0.025, 0.975), names = FALSE)
  k <- ncol(th)
  out <- data.table(scale_h1_lo = exp(q(th[, 1]))[1], scale_h1_hi = exp(q(th[, 1]))[2])
  if (p$scale == "horizon") out[, `:=`(slope_lo = q(th[, 2])[1], slope_hi = q(th[, 2])[2])]
  if (p$bias) out[, `:=`(drift_lo = q(th[, k])[1], drift_hi = q(th[, k])[2])]
  out
}

estimate_calibration <- function(pieces, strata, p) {
  labels <- vapply(pieces, `[[`, "", "label")
  groups <- list(all = pieces)
  if (!is.null(p$by)) {
    miss <- setdiff(p$by, names(strata))
    if (length(miss)) mt_abort(sprintf("`by` names unknown stratum variables: %s.", paste(miss, collapse = ", ")),
                               class = "argument")
    key <- do.call(paste, c(as.list(strata[match(labels, strata$label), p$by, with = FALSE]), sep = "."))
    groups <- c(groups, split(pieces, key))
  }
  n_draws <- ncol(pieces[[1]]$L)
  res <- with_seed(p$seed, lapply(groups, function(gp) {
    bands <- if (length(p$ages)) {
      sub <- lapply(gp, split_bands, breaks = p$ages)
      labs <- intersect(band_labels(p$ages), unlist(lapply(sub, names)))
      stats::setNames(lapply(labs, function(b) Filter(Negate(is.null), lapply(sub, `[[`, b))), labs)
    } else {
      list("all ages" = gp)
    }
    lapply(bands, function(bp) {
      r <- calibration_fit(stack_pieces(bp, n_draws), p$scale, p$bias, p$target %||% "crps", p$level %||% 0.95)
      set(r$table, j = "forecasts", value = length(bp))
      if ((p$n_boot %||% 0) > 0) r$table <- cbind(r$table, calibration_boot(bp, p, n_draws))
      r
    })
  }))
  tab <- rbindlist(lapply(names(res), function(g) {
    rbindlist(lapply(names(res[[g]]), function(b) cbind(data.table(group = g, band = b), res[[g]][[b]]$table)))
  }), fill = TRUE)
  list(table = tab, theta = lapply(res, function(r) lapply(r, `[[`, "theta")))
}

# Parameters of a group of strata by age band (the pooled ones for strata
# outside the groups; fits from before age bands hold one vector per group)
calibration_theta_of <- function(theta, key) {
  th <- theta[[key]] %||% theta$all
  if (is.numeric(th)) list("all ages" = th) else th
}

#' @export
mt_forecast.mt_fit_recalibrated <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                                          exposure = NULL, ...) {
  p <- fit$spec$params
  key <- if (is.null(p$by)) rep("all", nrow(fit$strata)) else
    do.call(paste, c(as.list(fit$strata[, p$by, with = FALSE]), sep = "."))
  band <- band_labels(p$ages)[findInterval(fit$age$lower, p$ages) + 1L]
  transform <- function(arr, targets, fc) {
    d <- dim(arr)
    hh <- rep(seq_len(d[2]), each = d[1])
    bb <- rep(band, d[2])
    for (s in targets) {
      th <- calibration_theta_of(fit$calibration_theta, key[s])
      L <- matrix(log(pmax(arr[, , s, ], 1e-12)), ncol = d[4])
      if (anyNA(L)) next
      m <- row_quantiles(L, 0.5)[, 1]
      sv <- bv <- numeric(length(hh))
      for (b in unique(bb)) {
        i <- which(bb == b)
        theta <- th[[b]] %||% th[[1]]
        sh <- calibration_shift(theta, hh[i], p$scale, p$bias, fit$calibration_hmax %||% Inf)
        sv[i] <- sh$s
        bv[i] <- sh$b
      }
      arr[, , s, ] <- exp(m + bv + sv * (L - m))
    }
    arr
  }
  res <- wrapped_forecast(fit, "mt_fit_recalibrated", h, to, n_draws, seed, predictive, exposure, transform,
                          h_missing = missing(h))
  res$calibration <- fit$calibration
  res
}

# Forecast of a wrapper whose fit is the base fit with an extra class `cls`:
# draws of the base model (underlying rates) are transformed per stratum by
# `transform(arr, targets, fc)`, observation noise is added for predictive
# draws (per draw of an ensemble from its component), and the aggregates of
# reconciled models are recomputed from the bottom strata
wrapped_forecast <- function(fit, cls, h, to, n_draws, seed, predictive, exposure, transform, h_missing = TRUE) {
  base_fit <- unwrap_fit(fit, cls)
  if (!is.null(to)) {
    if (!h_missing) mt_warn("Both `h` and `to` are given; `to` is used.", class = "argument")
    h <- to - max(fit$periods)
  }
  out <- with_seed(seed, {
    fc <- suppressWarnings(mt_forecast(base_fit, h = h, n_draws = n_draws, exposure = exposure))
    arr <- fc$rates
    d <- dim(arr)
    targets <- calibration_targets(fit)
    arr <- transform(arr, targets, fc)
    hs <- fit$reconcile$structure
    E <- expand_exposure(exposure %||% fit$exposure_last, d, dimnames(arr), hs)
    if (predictive) arr <- add_observation_noise(arr, base_fit, fc, targets, E)
    if (!is.null(hs)) arr <- aggregate_bottom(arr, E, hs, fit$reconcile$hierarchy$additive)
    list(fc = fc, rates = arr)
  })
  res <- out$fc
  res$rates <- out$rates
  dimnames(res$rates) <- dimnames(out$fc$rates)
  check_draws_plausible(res$rates, fit$age, fit$strata, fit$spec$label)
  res$model <- fit$spec$label
  res$spec <- fit$spec
  res$predictive <- predictive
  res
}

# A wrapper's fit is the fit of its base model under the wrapper's spec and
# class, with the wrapper's own `fields` (such as `reconcile` or `anchor`).
# The specs fitted underneath are kept as a stack: the fitted base spec can
# differ from spec$params$base (tuned() fits the spec it chose). So are the
# values that the wrapper's fields had in the base fit, which unwrap_fit()
# restores: an inner wrapper must not see the fields of an outer one (the
# hierarchy of reconciled(anchored(...)) would make the anchored draws
# coherent before they are reconciled)
wrap_fit <- function(fit, spec, cls, fields = list()) {
  fit$inner_specs <- c(list(fit$spec), fit$inner_specs)
  fit$inner_fields <- c(list(lapply(stats::setNames(nm = names(fields)), function(n) fit[[n]])), fit$inner_fields)
  for (n in names(fields)) fit[[n]] <- fields[[n]]
  fit$spec <- spec
  class(fit) <- c(cls, class(fit))
  fit
}

unwrap_fit <- function(fit, cls) {
  base <- fit
  class(base) <- class(fit)[-match(cls, class(fit))]
  base$spec <- if (length(fit$inner_specs)) fit$inner_specs[[1]] else fit$spec$params$base
  base$inner_specs <- fit$inner_specs[-1]
  if (length(fit$inner_fields)) {
    prev <- fit$inner_fields[[1]]
    for (n in names(prev)) base[[n]] <- prev[[n]]
    base$inner_fields <- fit$inner_fields[-1]
  }
  base
}

# Poisson and extra-Poisson noise of the base model added to rate draws.
# Cells can have zero deaths (rate 0), as in predictive_rates()
add_observation_noise <- function(arr, base_fit, fc, targets, E) {
  d <- dim(arr)
  nz <- draw_noise(base_fit, fc, d[4])
  for (s in targets) {
    R <- arr[, , s, ]
    if (anyNA(R) || !any(nz$noise(s))) next
    Es <- as.vector(E[, , s])
    sd <- array(nz$sd(s)[rep(seq_len(d[1]), d[2]), ], d[c(1, 2, 4)])
    mu <- R * exp(stats::rnorm(length(R), 0, sd)) * Es
    D <- stats::rpois(length(mu), mu) / pmax(Es, 1e-12)
    k <- rep(nz$noise(s), each = d[1] * d[2])
    R[k] <- D[k]
    arr[, , s, ] <- R
  }
  arr
}

# Recompute the aggregate strata of a hierarchy from the bottom strata
aggregate_bottom <- function(arr, E, hs, additive) {
  d <- dim(arr)
  b <- hs$bottom
  by_deaths <- !identical(additive, "rates")
  Eb <- E[, , b, drop = FALSE]
  vb <- if (by_deaths) arr[, , b, , drop = FALSE] * as.vector(Eb) else arr[, , b, , drop = FALSE]
  En <- array(NA_real_, d[1:3])
  for (t in seq_len(d[2])) {
    En[, t, ] <- matrix(Eb[, t, ], d[1]) %*% t(hs$Sx %||% hs$S)
    for (j in seq_len(d[4])) arr[, t, , j] <- matrix(vb[, t, , j], d[1]) %*% t(hs$S)
  }
  if (by_deaths) arr / as.vector(En) else arr
}
