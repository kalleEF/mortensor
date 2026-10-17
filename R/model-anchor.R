#' Control the jump-off of any model
#'
#' Forecasts of most models start from their fitted rates in the last
#' observed year. Where the fit misses the level of the last years (for
#' example after a recent change of trend), the forecast starts with a step
#' away from the observations; where the observations are noisy, starting
#' from them carries the noise into every forecast year. `anchored()` wraps
#' any model specification and moves the start of its forecasts towards the
#' observed level, per age and stratum:
#'
#' 1. the offset is the log ratio of observed to fitted deaths, summed over
#'    the last `years` training years (`to = "observed"`: the rates before
#'    smoothing; `to = "smoothed"`: the rates after preprocessing, which
#'    differs when the pipeline smooths);
#' 2. offsets of cells with few deaths are shrunk towards the deaths-weighted
#'    mean offset of their stratum (empirical Bayes, with the Poisson
#'    variance of each offset; `shrink = TRUE`);
#' 3. the forecast draws are multiplied by `exp(offset * w(h))`, where the
#'    weight `w(h)` falls linearly from 1 at the first forecast year to 0
#'    after `decay` years (`decay = Inf` keeps the offset, `decay = 0`
#'    leaves the model's own jump-off).
#'
#' With the default `decay = 10`, forecasts more than ten years ahead are
#' those of the model itself; during the fade the offsets add a small
#' apparent trend. `years = 1, decay = Inf, shrink = FALSE` is the observed jump-off of
#' Lee-Carter models (`lee_carter(jump_off = "observed")`). Models that
#' have no fitted rates (random walks, time series per age) start from the
#' observed rates already and are left unchanged; in an [ensemble()] each
#' draw is anchored with the offset of the component it comes from. For
#' [reconciled()] models the bottom strata are anchored and the aggregates
#' recomputed, so the draws stay coherent. [mt_jump_off()] reports the step
#' between the last observed or fitted year and the first forecast year.
#'
#' Anchoring helps where the fitted rates of the last year are far from the
#' observed level. With large counts the last observed year is a good
#' anchor; where the last observed rates are noisy, several years are
#' better. Where the model's own jump-off is already close to the observed
#' level, anchoring changes little and can narrow the intervals of the
#' first years too much. The vignette "Where a forecast starts: the
#' jump-off" shows the effect in a backtest.
#'
#' @param spec A model specification.
#' @param to `"observed"` or `"smoothed"`.
#' @param years Number of last training years whose deaths define the
#'   offset.
#' @param decay Years over which the offset fades out (`Inf`: constant).
#' @param shrink Shrink offsets of cells with few deaths.
#' @return A model specification. The fitted object has the offsets in
#'   `$anchor$offset` (`[age, stratum]`, log scale).
#' @examples
#' x <- mt_example(deaths = "single", years = 1995:2019)
#' fit <- mt_fit(anchored(lee_carter(), years = 3, decay = 10), x)
#' summary(as.vector(fit$anchor$offset))
#' fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
#' mt_jump_off(fc)
#' @export
anchored <- function(spec, to = c("observed", "smoothed"), years = 3, decay = 10, shrink = TRUE) {
  check_class(spec, "mt_spec", "spec")
  to <- match.arg(to)
  check_number(years, "years", lower = 1, integerish = TRUE)
  check_number(decay, "decay", lower = 0)
  check_flag(shrink, "shrink")
  new_mt_spec("anchored", list(base = spec, to = to, years = years, decay = decay, shrink = shrink),
              joint = spec$joint, min_periods = max(spec$min_periods, years), smooth = spec$smooth,
              label = sprintf("%s, anchored (%s, %d year%s, decay %s)", spec$label, to, years,
                              if (years > 1) "s" else "", format(decay)))
}

#' @export
mt_fit.mt_spec_anchored <- function(spec, data, prep = NULL, cores = 1, ...) {
  p <- spec$params
  fit <- mt_fit(p$base, data, prep = prep, cores = cores, ...)
  D <- if (p$to == "observed") fit$raw_D %||% fit$train_D else fit$train_D
  offset_of <- function(f) anchor_offsets(f, D, fit$train_E, p$years, p$shrink)
  wrap_fit(fit, spec, "mt_fit_anchored", list(anchor = list(
    offset = offset_of(fit), components = if (!is.null(fit$components)) lapply(fit$components, offset_of))))
}

# Log ratio of observed to fitted deaths over the last `years` years per age
# and stratum, shrunk towards the stratum mean; 0 where the model has no
# fitted rates
anchor_offsets <- function(fit, D, E, years, shrink) {
  f <- tryCatch(stats::fitted(fit), error = function(e) NULL)
  d <- dim(D)
  out <- matrix(0, d[1], d[3])
  if (is.null(f) || !identical(dim(f), d)) return(out)
  last <- seq(d[2] - years + 1, d[2])
  Dl <- apply(D[, last, , drop = FALSE], c(1, 3), sum)
  Fl <- apply(f[, last, , drop = FALSE] * E[, last, , drop = FALSE], c(1, 3), sum)
  r <- log(pmax(Dl, 0.5) / Fl)
  r[!is.finite(r)] <- 0
  if (shrink) {
    for (s in seq_len(d[3])) {
      v <- 1 / pmax(Dl[, s], 0.5)
      w <- pmax(Dl[, s], 0.5)
      rbar <- sum(w * r[, s]) / sum(w)
      # weighted variance of the offsets minus their weighted sampling
      # variance (sum(w * v) / sum(w) with w = 1 / v)
      tau2 <- max(sum(w * (r[, s] - rbar)^2) / sum(w) - sum(w * v) / sum(w), 0)
      r[, s] <- rbar + tau2 / (tau2 + v) * (r[, s] - rbar)
    }
  }
  out[] <- r
  out
}

#' @export
mt_forecast.mt_fit_anchored <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                                        exposure = NULL, ...) {
  decay <- fit$spec$params$decay
  transform <- function(arr, targets, fc) {
    d <- dim(arr)
    w <- if (decay == 0) rep(0, d[2]) else pmax(0, 1 - (seq_len(d[2]) - 1) / decay)
    comp <- fc$component
    offs <- fit$anchor$components
    for (s in targets) {
      if (!is.null(offs) && !is.null(comp)) {
        # one component per draw, or per stratum and draw after reordering (reconciled())
        cs <- if (is.matrix(comp)) comp[s, ] else comp
        for (k in unique(cs)) {
          j <- which(cs == k)
          arr[, , s, j] <- arr[, , s, j] * as.vector(exp(outer(offs[[k]][, s], w)))
        }
      } else {
        arr[, , s, ] <- arr[, , s, ] * as.vector(exp(outer(fit$anchor$offset[, s], w)))
      }
    }
    arr
  }
  wrapped_forecast(fit, "mt_fit_anchored", h, to, n_draws, seed, predictive, exposure, transform,
                   h_missing = missing(h))
}

#' Step between the last observed year and the first forecast year
#'
#' Reports, per stratum, how far the first forecast year lies from the last
#' observed and the last fitted year: for each age the log ratio of the
#' forecast median to the observed and to the fitted rate, and for life
#' expectancy the difference in years. Large steps against the observations
#' point to a jump-off problem (see [anchored()]); steps against the fitted
#' rates are the forecast of one year.
#'
#' @param draws An [mt_draws] object.
#' @param quantity `"life_expectancy"` (one row per stratum) or `"rate"`
#'   (one row per stratum and age).
#' @param at Age for life expectancy (default: the first grid age).
#' @return A `data.table` with the stratum columns, the observed and fitted
#'   values of the last training year, the median of the first forecast
#'   year and the steps (`step_observed`, `step_fitted`; log ratios for
#'   rates, years for life expectancy).
#' @examples
#' x <- mt_example(deaths = "single", years = 1995:2019)
#' fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 100, seed = 1)
#' mt_jump_off(fc)
#' head(mt_jump_off(fc, quantity = "rate"))
#' @export
mt_jump_off <- function(draws, quantity = c("life_expectancy", "rate"), at = NULL) {
  check_class(draws, "mt_draws", "draws")
  quantity <- match.arg(quantity)
  o <- draws$observed
  if (is.null(o)) mt_abort("The draws hold no observed rates.", class = "argument")
  nt <- dim(o)[2]
  obs_T <- o[, nt, , drop = FALSE]
  fit_T <- if (!is.null(draws$fitted)) draws$fitted[, nt, , drop = FALSE] else obs_T * NA
  first <- draws$rates[, 1, , , drop = FALSE]
  st <- draws$strata
  if (quantity == "rate") {
    med <- apply(first, c(1, 3), stats::median)
    out <- data.table(stratum = rep(st$stratum, each = length(draws$age)),
                      age = rep(draws$age$labels, nrow(st)), observed = as.vector(obs_T[, 1, ]),
                      fitted = as.vector(fit_T[, 1, ]), forecast = as.vector(med))
    out[, `:=`(step_observed = log(forecast / observed), step_fitted = log(forecast / fitted))]
  } else {
    e_obs <- life_expectancy(obs_T, draws$age, at = at)
    e_fit <- if (all(is.na(fit_T))) rep(NA_real_, nrow(st)) else life_expectancy(fit_T, draws$age, at = at)
    e_first <- apply(life_expectancy(first, draws$age, at = at), 2, stats::median)
    out <- data.table(stratum = st$stratum, observed = as.vector(e_obs), fitted = as.vector(e_fit),
                      forecast = as.vector(e_first))
    out[, `:=`(step_observed = forecast - observed, step_fitted = forecast - fitted)]
  }
  out <- merge(st, out, by = "stratum", sort = FALSE)
  out[, stratum := NULL]
  out[]
}
