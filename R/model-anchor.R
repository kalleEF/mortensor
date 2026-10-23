#' Control the jump-off of any model
#'
#' Forecasts of most models start from their fitted rates in the last
#' observed year. Where the fit misses the level of the last years (for
#' example after a recent change of trend), the forecast starts with a step
#' away from the observations. `anchored()` wraps any model specification
#' and moves the start of its forecasts to the observed level, per age and
#' stratum. With the defaults this is the jump-off adjustment of Lee and
#' Miller (2001): the forecasts start from the observed rates of the last
#' training year and keep this offset in every forecast year.
#'
#' 1. the offset is the log ratio of observed to fitted deaths, summed over
#'    the last `years` training years (`to = "observed"`: the rates before
#'    smoothing; `to = "smoothed"`: the rates after preprocessing, which
#'    differs when the pipeline smooths);
#' 2. with `shrink = TRUE`, offsets of ages with few deaths are shrunk
#'    towards the mean offset of their stratum (see below);
#' 3. the forecast draws are multiplied by `exp(offset * w(h))`, where the
#'    weight `w(h)` falls linearly from 1 at the first forecast year to 0
#'    after `decay` years (`decay = Inf` keeps the offset; `decay = 0`
#'    leaves the model's own jump-off).
#'
#' The defaults `years = 1, decay = Inf, shrink = FALSE` give, for any
#' model, the observed jump-off of Lee-Carter models
#' (`lee_carter(jump_off = "observed")`). The offset multiplies every draw
#' of a cell by the same factor, so it moves the forecast distribution of
#' the log rates without changing its width. Models that have no fitted
#' rates (random walks, time series per age) start from the observed rates
#' already and are left unchanged; in an [ensemble()] each draw is anchored
#' with the offset of the component it comes from, which moves the
#' components relative to each other and can change the width of the
#' mixture. For [reconciled()]
#' models the bottom strata are anchored and the aggregates recomputed, so
#' the draws stay coherent. [mt_jump_off()] reports the step between the
#' last observed or fitted year and the first forecast year.
#'
#' **Jump-off options of Stoeldraijer et al. (2018).** The jump-off rates
#' are the rates of the last year of the fitting period, the jump-off year.
#' Stoeldraijer et al. compare six options for the Lee-Carter model (their
#' names in quotes):
#'
#' * "model values", the fitted rates of the jump-off year (Lee and Carter
#'   1992): the model itself, without `anchored()` (or `decay = 0`);
#' * "last observed values" (Lee and Miller 2001): `anchored()` with its
#'   defaults, for any model; for Lee-Carter also `lee_carter(jump_off =
#'   "observed")`;
#' * "average n years observed" (`n` = 2 to 5): `lee_carter(jump_off =
#'   "observed", jump_off_years = n)`, as the paper defines it: the mean of
#'   the observed log rates of the last `n` years, with the index 0 at the
#'   midpoint of these years. `anchored(years = n)` is a package-specific
#'   variant for any model (reference unknown): its offset is the log ratio
#'   of the observed to the fitted deaths of these years, so the forecast
#'   also keeps the change of the fitted rates from these years to the last
#'   year.
#'
#' They call a forecast accurate if its out-of-sample errors in historical
#' data are small (mean absolute error of the forecast life expectancy at
#' 65 one and five years ahead), and robust if it changes little when one
#' more year of data is added (standard deviation of the changes of the
#' forecast life expectancy at 65 in 2020 between successive fitting
#' periods). In eight Western European countries the last observed values
#' were the most accurate for most countries and an average of observed
#' years the most robust; the more years averaged, the more robust and the
#' less accurate the forecast. The model values were the most accurate where
#' the model fitted the last years well. They recommend moving from the most
#' accurate jump-off in the first forecast year to the most robust one by
#' the fifth (Appendix 6, without a formula for the weights). Where the
#' last observed values are the most accurate and the model values the most
#' robust, `anchored(decay = 4)` does this: the weight of the observed
#' offset falls linearly from 1 in the first forecast year to 0 in the fifth
#' (the linear weights are package-specific). The filtered jump-off of De Jong
#' and Tickle (2006) is `lee_carter(jump_off = "filtered")`.
#'
#' **Package-specific options** (reference unknown). Starting from the
#' observed rates of one year carries their random variation into every
#' forecast year, which matters where deaths are few. The other settings
#' are package-specific: averaging the offset over the last `years` years
#' (see above), fading it out over `decay` years (forecasts more than
#' `decay` years ahead are then those of the model itself, and during the
#' fade the offsets add a small apparent trend) and shrinking it (`shrink =
#' TRUE`). The shrinkage
#' is empirical Bayes within each stratum: the offset \eqn{r_x} of age
#' \eqn{x} has the approximate Poisson sampling variance \eqn{v_x = 1/D_x},
#' where \eqn{D_x} is the number of deaths in the last `years` (at least
#' 0.5); the variance \eqn{\tau^2} of the true offsets between the ages of
#' the stratum is estimated with the moment estimator of DerSimonian and
#' Laird (1986); and each offset becomes
#' \eqn{\hat\mu + \tau^2 / (\tau^2 + v_x) (r_x - \hat\mu)}, where
#' \eqn{\hat\mu} is the mean offset weighted by \eqn{1 / (\tau^2 + v_x)}.
#'
#' Anchoring helps where the fitted rates of the last year are far from the
#' observed level; where the model's own jump-off is already close to it,
#' anchoring changes little. The vignette "Where a forecast starts: the
#' jump-off" shows the effect in a backtest.
#'
#' @param spec A model specification.
#' @param to `"observed"` or `"smoothed"`.
#' @param years Number of last training years whose deaths define the
#'   offset (default 1, the last training year; see above for the averages
#'   of Stoeldraijer et al. 2018).
#' @param decay Years over which the offset fades out (default `Inf`:
#'   constant).
#' @param shrink Shrink offsets of ages with few deaths towards the mean of
#'   their stratum (default `FALSE`).
#' @return A model specification. The fitted object has the offsets in
#'   `$anchor$offset` (`[age, stratum]`, log scale).
#' @references
#' Stoeldraijer L, van Duin C, van Wissen L, Janssen F (2018). Comparing
#' strategies for matching mortality forecasts to the most recently observed
#' data: exploring the trade-off between accuracy and robustness. *Genus*
#' 74: 16. \doi{10.1186/s41118-018-0041-y}
#'
#' Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
#' *Journal of the American Statistical Association* 87(419): 659-671.
#'
#' Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
#' method for forecasting mortality. *Demography* 38(4): 537-549.
#'
#' De Jong P, Tickle L (2006). Extending Lee-Carter mortality forecasting.
#' *Mathematical Population Studies* 13(1): 1-18.
#' \doi{10.1080/08898480500452109}
#'
#' DerSimonian R, Laird N (1986). Meta-analysis in clinical trials.
#' *Controlled Clinical Trials* 7(3): 177-188.
#' @examples
#' x <- mt_example(deaths = "single", years = 1995:2019)
#' fit <- mt_fit(anchored(lee_carter()), x)
#' summary(as.vector(fit$anchor$offset))
#' fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
#' mt_jump_off(fc)
#' # package-specific: the mean offset of three years, fading out over ten years
#' fit3 <- mt_fit(anchored(lee_carter(), years = 3, decay = 10, shrink = TRUE), x)
#' @export
anchored <- function(spec, to = c("observed", "smoothed"), years = 1, decay = Inf, shrink = FALSE) {
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
# and stratum, optionally shrunk towards the stratum mean; 0 where the model
# has no fitted rates
anchor_offsets <- function(fit, D, E, years, shrink) {
  f <- tryCatch(stats::fitted(fit), error = function(e) NULL)
  d <- dim(D)
  out <- matrix(0, d[1], d[3])
  if (is.null(f) || !identical(dim(f), d)) return(out)
  last <- seq(d[2] - years + 1, d[2])
  Dl <- apply(D[, last, , drop = FALSE], c(1, 3), sum)
  Fl <- apply(f[, last, , drop = FALSE] * E[, last, , drop = FALSE], c(1, 3), sum)
  r <- log(pmax(Dl, 0.5) / Fl)
  ok <- is.finite(r)
  r[!ok] <- 0
  if (shrink) {
    for (s in seq_len(d[3])) {
      k <- which(ok[, s])
      if (length(k) > 1) r[k, s] <- shrink_dl(r[k, s], 1 / pmax(Dl[k, s], 0.5))$shrunk
    }
  }
  out[] <- r
  out
}

# Empirical Bayes shrinkage of estimates y with sampling variances v towards
# their mean, with the between-estimate variance tau2 from the moment
# estimator of DerSimonian and Laird (1986):
#   tau2 = max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))), w = 1 / v,
# where Q is the weighted sum of squares around the mean weighted by w. The
# estimates are shrunk to the mean weighted by 1 / (tau2 + v)
shrink_dl <- function(y, v) {
  k <- length(y)
  w <- 1 / v
  ybar <- sum(w * y) / sum(w)
  q <- sum(w * (y - ybar)^2)
  denom <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (k > 1 && denom > 0) max(0, (q - (k - 1)) / denom) else 0
  ws <- 1 / (tau2 + v)
  mu <- sum(ws * y) / sum(ws)
  list(shrunk = mu + tau2 / (tau2 + v) * (y - mu), tau2 = tau2, mu = mu)
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
