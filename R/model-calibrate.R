#' Calibrate forecast distributions on recent years
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
#' and they cannot anticipate changes after the window. Compare calibrated
#' and uncalibrated models in an outer backtest.
#'
#' In backtests on England & Wales men (HMD via StMoMo, origins 1996-2001,
#' h = 10), calibration raised the 95% coverage of [lee_carter()] from 0.75
#' to 0.88 and lowered its CRPS by 4% (0.0687 to 0.0657); `bias = TRUE`
#' gave the same. On the ONS deprivation data (6-year window), the coverage
#' of Lee-Carter and Li-Lee rose with `target = "interval"` (0.83 to 0.90 and
#' 0.84 to 0.89) at unchanged CRPS, and `bias = TRUE` made the CRPS worse
#' (0.069 to 0.082). An equal-weight [ensemble()] of Lee-Carter and [rwd()]
#' was already calibrated; calibration then narrowed it (95% coverage 0.93
#' to 0.85 on ONS, 0.95 to 0.90 on England & Wales) at the same CRPS
#' (see the vignette "Forecasting mortality by deprivation quintile").
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
#' @return A model specification. The fitted object has a `calibration`
#'   element: a `data.table` with one row per group (`scale_h1`, `slope`,
#'   `drift`, the number of scored cells, and CRPS and coverage of the inner
#'   forecasts before and after calibration).
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
#' fit <- mt_fit(calibrated(lee_carter(), years = 8, n_origins = 3, n_draws = 50), x)
#' fit$calibration
#' fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
#' @export
calibrated <- function(spec, years = 10, n_origins = 5, scale = c("horizon", "constant"), bias = FALSE,
                       target = c("crps", "interval"), level = 0.95, by = NULL, n_draws = 100, seed = 1) {
  check_class(spec, "mt_spec", "spec")
  scale <- match.arg(scale)
  target <- match.arg(target)
  check_number(level, "level", lower = 0.5, upper = 0.99)
  check_number(years, "years", lower = 2, integerish = TRUE)
  check_number(n_origins, "n_origins", lower = 1, integerish = TRUE)
  check_number(n_draws, "n_draws", lower = 10, integerish = TRUE)
  check_flag(bias, "bias")
  if (n_origins > years - 1) n_origins <- years - 1
  new_mt_spec("calibrated", list(base = spec, years = years, n_origins = n_origins, scale = scale, bias = bias,
                                 target = target, level = level, by = by, n_draws = n_draws, seed = seed),
              joint = spec$joint, min_periods = spec$min_periods + years, smooth = spec$smooth,
              label = paste0(spec$label, ", calibrated", if (length(by)) paste0(" by ", paste(by, collapse = " x ")),
                             if (target == "interval") sprintf(" (%g%% interval)", 100 * level)))
}

#' @export
mt_fit.mt_spec_calibrated <- function(spec, data, prep = NULL, cores = 1, ...) {
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
    tryCatch(calibration_origin(p$base, data, o, last, prep, p$n_draws, p$seed + o, fit$strata$label[bottom]),
             error = function(e) e)
  }, cores)
  bad <- vapply(pieces, inherits, logical(1), what = "error")
  if (any(bad)) {
    mt_warn(sprintf("Calibration origin(s) %s failed: %s", paste(origins[bad], collapse = ", "),
                    conditionMessage(pieces[[which(bad)[1]]])), class = "calibration")
  }
  pieces <- unlist(pieces[!bad], recursive = FALSE)
  if (!length(pieces)) mt_abort("Every calibration origin failed.", class = "calibration")
  est <- estimate_calibration(pieces, fit$strata, p)
  fit$calibration <- est$table
  fit$calibration_theta <- est$theta
  fit$calibration_origins <- origins[!bad]
  fit$calibration_hmax <- last - min(origins[!bad])
  fit$spec <- spec
  class(fit) <- c("mt_fit_calibrated", class(fit))
  fit
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

# Observation noise of each draw: whether Poisson noise is added (`noise`, one
# per draw) and the extra-Poisson sd by age (`sd(s)`: ages x draws for stratum
# s). Draws of an ensemble follow the component they come from.
draw_noise <- function(fit, fc, n_draws) {
  na <- length(fit$age)
  comps <- fit$components
  if (!is.null(comps) && !is.null(fc$component)) {
    flags <- vapply(comps, function(f) spec_adds_noise(f$spec), TRUE)
    k <- fc$component
    return(list(noise = flags[k], sd = function(s) {
      matrix(vapply(k, function(j) {
        if (!flags[j] || is.null(comps[[j]]$noise_sd)) rep(0, na) else comps[[j]]$noise_sd[, s]
      }, numeric(na)), na)
    }))
  }
  flag <- spec_adds_noise(fit$spec)
  sd <- fit$noise_sd
  list(noise = rep(flag, n_draws), sd = function(s) {
    matrix(if (flag && !is.null(sd)) sd[, s] else 0, na, n_draws)
  })
}

# One inner origin: forecast draws of the rest of the window and the observed
# cells they are scored against, per stratum
calibration_origin <- function(spec, data, origin, last, prep, n_draws, seed, targets) {
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
    Eg <- as.vector(E_grid[, , sf])[used]
    Ecell <- as.numeric(C %*% Eg)
    o <- obs[sel]
    L <- log(pmax(R[used, , drop = FALSE], 1e-12))
    out[[length(out) + 1L]] <- list(
      label = lab, L = L, m = row_quantiles(L, 0.5)[, 1], h = rep(fc$periods - origin, each = na)[used],
      E = Eg, C = C, Ecell = Ecell, y = log(pmax(o$value, 0.5) / Ecell), noise = nz$noise,
      sd = nz$sd(sf)[rep(seq_len(na), np)[used], , drop = FALSE])
  }
  out
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
  P
}

# Spread factor and drift correction per grid cell; beyond the longest
# horizon of the calibration window the spread factor stays constant
calibration_shift <- function(theta, h, scale, bias, h_max = Inf) {
  s <- exp(theta[1] + if (scale == "horizon") theta[2] * log(pmin(h, h_max)) else 0)
  b <- if (bias) theta[length(theta)] * h else 0
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

calibration_fit <- function(P, scale, bias, target = "crps", level = 0.95) {
  k <- 1 + (scale == "horizon") + bias
  lower <- c(log(0.1), if (scale == "horizon") -1, if (bias) -0.05)
  upper <- c(log(20), if (scale == "horizon") 1, if (bias) 0.05)
  f <- function(theta) {
    if (any(theta < lower | theta > upper)) return(1e10)
    X <- calibration_sample(theta, P, scale, bias)
    if (target == "crps") mean(crps_rows(X, P$y)) else interval_score(X, P$y, level)
  }
  theta <- if (k == 1) {
    stats::optimize(f, c(lower, upper))$minimum
  } else {
    stats::optim(rep(0, k), f, method = "Nelder-Mead", control = list(reltol = 1e-6, maxit = 300))$par
  }
  summ <- function(th) {
    X <- calibration_sample(th, P, scale, bias)
    Q <- row_quantiles(X, c(0.025, 0.1, 0.9, 0.975))
    c(crps = mean(crps_rows(X, P$y)), cov80 = mean(P$y >= Q[, 2] & P$y <= Q[, 3]),
      cov95 = mean(P$y >= Q[, 1] & P$y <= Q[, 4]))
  }
  b0 <- summ(rep(0, k))
  b1 <- summ(theta)
  list(theta = theta, table = data.table(
    scale_h1 = exp(theta[1]), slope = if (scale == "horizon") theta[2] else 0,
    drift = if (bias) theta[k] else 0, cells = length(P$y), crps_before = b0[["crps"]],
    crps_after = b1[["crps"]], cov80_before = b0[["cov80"]], cov80_after = b1[["cov80"]],
    cov95_before = b0[["cov95"]], cov95_after = b1[["cov95"]]))
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
  res <- with_seed(p$seed, lapply(groups, function(gp) {
    calibration_fit(stack_pieces(gp, ncol(gp[[1]]$L)), p$scale, p$bias, p$target %||% "crps", p$level %||% 0.95)
  }))
  tab <- rbindlist(lapply(names(res), function(n) cbind(data.table(group = n), res[[n]]$table)))
  list(table = tab, theta = lapply(res, `[[`, "theta"))
}

#' @export
mt_forecast.mt_fit_calibrated <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                                          exposure = NULL, ...) {
  base_fit <- fit
  class(base_fit) <- setdiff(class(fit), "mt_fit_calibrated")
  base_fit$spec <- fit$spec$params$base
  p <- fit$spec$params
  if (!is.null(to)) h <- to - max(fit$periods)
  out <- with_seed(seed, {
    fc <- suppressWarnings(mt_forecast(base_fit, h = h, n_draws = n_draws, exposure = exposure))
    arr <- fc$rates
    d <- dim(arr)
    targets <- calibration_targets(fit)
    key <- if (is.null(p$by)) rep("all", d[3]) else
      do.call(paste, c(as.list(fit$strata[, p$by, with = FALSE]), sep = "."))
    hh <- rep(seq_len(d[2]), each = d[1])
    for (s in targets) {
      theta <- fit$calibration_theta[[key[s]]] %||% fit$calibration_theta$all
      L <- matrix(log(pmax(arr[, , s, ], 1e-12)), ncol = d[4])
      if (anyNA(L)) next
      m <- row_quantiles(L, 0.5)[, 1]
      sh <- calibration_shift(theta, hh, p$scale, p$bias, fit$calibration_hmax %||% Inf)
      arr[, , s, ] <- exp(m + sh$b + sh$s * (L - m))
    }
    E <- expand_exposure(exposure %||% fit$exposure_last, d, dimnames(arr))
    if (predictive) {
      nz <- draw_noise(base_fit, fc, d[4])
      if (any(nz$noise)) {
        for (s in targets) {
          R <- arr[, , s, ]
          if (anyNA(R)) next
          Es <- as.vector(E[, , s])
          sd <- array(nz$sd(s)[rep(seq_len(d[1]), d[2]), ], d[c(1, 2, 4)])
          mu <- R * exp(stats::rnorm(length(R), 0, sd)) * Es
          D <- pmax(stats::rpois(length(mu), mu), 0.5) / Es
          k <- rep(nz$noise, each = d[1] * d[2])
          R[k] <- D[k]
          arr[, , s, ] <- R
        }
      }
    }
    hs <- fit$reconcile$structure
    if (!is.null(hs)) arr <- aggregate_bottom(arr, E, hs, fit$reconcile$hierarchy$additive)
    list(fc = fc, rates = arr)
  })
  res <- out$fc
  res$rates <- out$rates
  dimnames(res$rates) <- dimnames(out$fc$rates)
  check_draws_plausible(res$rates, fit$age, fit$strata, fit$spec$label)
  res$model <- fit$spec$label
  res$predictive <- predictive
  res$calibration <- fit$calibration
  res
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
