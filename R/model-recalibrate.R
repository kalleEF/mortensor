#' Recalibrate forecast distributions on recent years
#'
#' Wraps any model specification so that its forecast distribution is
#' recalibrated on the last `years` observed years of the training data:
#' the spread (and optionally the trend) of its draws is rescaled (`method =
#' "spread"`, default), or the draws are replaced by the point forecast of
#' the model plus resampled errors of its past forecasts (`method =
#' "errors"`). Each time the model is fitted (also at every origin of an
#' outer [mt_backtest()]):
#'
#' 1. the base model is refitted at `n_origins` origins inside the
#'    calibration window (training data up to the origin) and forecasts the
#'    remaining years of the window;
#' 2. the draws are transformed per cell on the log scale,
#'    \deqn{\log m^* = \tilde m + b\,h + s(h)\,(\log m - \tilde m),}
#'    where \eqn{\tilde m} is the median of the draws in the cell, \eqn{h} the
#'    horizon, \eqn{s(h)} a spread factor, one value \eqn{s_1} for all
#'    horizons (`scale = "constant"`, default) or \eqn{s_1 h^{\gamma}}
#'    (`scale = "horizon"`, deprecated), held constant beyond the longest horizon of the
#'    window, or one value \eqn{s_h} for each horizon (`scale =
#'    "per_horizon"`), and \eqn{b} a correction of the drift per year (`bias =
#'    TRUE`; with `scale = "per_horizon"` one \eqn{b_h} for each horizon);
#' 3. the parameters minimise the mean CRPS (or the interval score,
#'    `target`) of the log death rates in the window, with the Poisson (or
#'    the model's negative binomial) variation of the observed deaths, as in
#'    [mt_backtest()]. Draws of an ensemble get the observation noise of the
#'    component they come from.
#'
#' The final model is fitted on all training data and its forecasts are
#' transformed with the estimated parameters. Ranks within each cell are
#' kept, so paths stay paths, and derived quantities such as life expectancy
#' can be computed from the calibrated draws.
#'
#' **Published methods.** The closest published method is the
#' member-by-member post-processing of ensemble forecasts (Van Schaeybroeck
#' and Vannitsem 2015), which shifts every member of an ensemble and
#' rescales its deviation from the ensemble mean, with parameters for each
#' lead time estimated by maximum likelihood or minimum CRPS. Here each
#' simulated path is rescaled around the median of the draws of each cell,
#' on the log scale, with one spread factor for all horizons (`scale =
#' "constant"`) estimated by minimum CRPS on inner rolling-origin forecasts.
#' Minimum-CRPS estimation follows ensemble model output statistics (EMOS,
#' Gneiting et al. 2005), which fits a parametric predictive distribution
#' instead of transforming the draws.
#'
#' **Factors per horizon** (`scale = "per_horizon"`). The published
#' post-processing methods fit their parameters to the forecasts of one
#' lead time: Gneiting et al. (2005) fit EMOS to 48-hour forecasts (section
#' 3) and note that distinct forecast horizons may need distinct training
#' periods (section 3a), and member-by-member post-processing is done for
#' each lead time, location and variable separately (Van Schaeybroeck and
#' Vannitsem 2015). The factor \eqn{s_h} minimises the mean CRPS (or the
#' interval score) of the inner forecasts of horizon \eqn{h} alone, by the
#' minimum-CRPS estimation of the constant factor (Gneiting et al. 2005,
#' section 2c); with `bias = TRUE` each horizon also gets its own shift
#' \eqn{b_h h}. Where the spread of the model is wrong by the same factor at
#' every horizon, the factors of all horizons estimate the constant factor,
#' each from fewer forecasts.
#'
#' The factor of horizon \eqn{h} rests on the inner origins at least \eqn{h}
#' years before the last training year. A horizon gets a factor only where
#' at least two origins reach it: the forecast of one origin has the errors
#' of a single year, which all its ages and strata share. When every origin
#' has enough training years, the longest horizon with a factor is `years -
#' floor(years / n_origins)` (8 with the defaults). The fit reports it as
#' `calibration_hmax`, and the number of origins of each horizon in column
#' `origins`. A forecast further ahead stops with an error of class
#' `mt_error_recalibrate`: no factor is extrapolated. For forecasts of
#' \eqn{h} years, use for example `years = h + 1` and `n_origins = h` (an
#' origin in every year from \eqn{h + 1} to 2 years before the last
#' training year); horizon \eqn{h} then rests on two origins and horizon 1
#' on \eqn{h}.
#'
#' **Resampled past errors** (`method = "errors"`). The draws are the point
#' forecast of the base model plus errors of its own past forecasts. The
#' inner forecasts are those of step 1, so each inner fit sees the data up
#' to its origin only. The error vector \eqn{e_{o,h}} of origin \eqn{o} and
#' horizon \eqn{h} holds the observed log rates of year \eqn{o + h} minus the
#' median of the inner draws, for every age of the grid and every stratum.
#' Draw \eqn{j} of the forecast at horizon \eqn{h} is
#' \deqn{\log m^*_j = \tilde m + \pi_h\, e_{o_j,h},}
#' where \eqn{\tilde m} is the median of the draws of the base model in the
#' cell (its point forecast) and the origin \eqn{o_j} is drawn with
#' replacement from the inner origins that reach horizon \eqn{h}, for each
#' draw and horizon.
#'
#' * Keilman, Pham and Hetland (2002, section 2) describe the extrapolation
#'   of errors of past forecasts: a deterministic forecast gives the
#'   expected value, and the errors of historical forecasts give the
#'   distribution around it, which they summarise for each forecast
#'   duration (section 3.1.1.2).
#' * Shang and Hyndman (2017, section 4.4 of the preprint, pointwise
#'   intervals) build intervals from h-step errors, observed minus forecast
#'   values, resampled with replacement: lower and upper bounds from the
#'   resampled errors of each cell are scaled by a tuning factor
#'   \eqn{\pi_h}, the one with which a share `level` of all errors of the
#'   horizon lies within the scaled bounds. Here the bounds are the
#'   quantiles at \eqn{(1 - level) / 2} and \eqn{(1 + level) / 2} of the
#'   bootstrap distribution itself (the empirical distribution of the
#'   errors of the cell, the limit of many resamples), and \eqn{\pi_h} is
#'   the smallest such factor. These bounds contain a
#'   share `level` of the errors of each cell, so \eqn{\pi_h \le 1}; with
#'   fewer than \eqn{1 / (1 - level)} origins at a horizon (20 at 95%), they
#'   are the smallest and the largest error and \eqn{\pi_h = 1}.
#' * One origin per draw and horizon gives all ages and strata the errors of
#'   the same past forecast, which keeps their dependence: the joint
#'   bootstrap of Panagiotelis et al. (2023).
#' * The errors are used as they are, as by Shang and Hyndman (2017), so a
#'   mean error of the inner forecasts shifts the forecast distribution away
#'   from the point forecast. Keilman, Pham and Hetland (2002, section 2)
#'   assume an expected error of zero.
#'
#' Deviations and package-specific parts (reference unknown): Shang and
#' Hyndman (2017) forecast with principal components estimated from all
#' years; here the model is refitted at each origin, as in the historical
#' forecasts of Keilman, Pham and Hetland (2002). The origins are those of
#' `years` and `n_origins` (Shang and Hyndman use every year). The point
#' forecast is the median of the base draws in each cell. The origins are
#' drawn independently for each horizon, so a draw combines the errors of
#' different origins in different years: quantities that combine years,
#' such as cohort measures or deaths summed over years, miss the dependence
#' of the errors over time. Panagiotelis et al. (2023) note that
#' multi-step forecasts need block-wise bootstraps to keep such serial
#' dependence. The errors are those of the observed log rates, with half a
#' death in empty cells as in [mt_backtest()], so the draws describe
#' observed rates and predictive draws add nothing. The method needs
#' observed deaths on the age groups of the grid in every cell of the
#' window, and at least two origins at each forecast horizon, as `scale =
#' "per_horizon"`; otherwise the fit or the forecast stops with an error.
#' The method gives no coverage guarantee: its intervals cover future
#' errors as far as these resemble the errors of the inner forecasts.
#' Keilman, Pham and Hetland (2002) name two limits of past errors: their
#' series are short, and errors may shrink over successive forecast rounds.
#' The draws of a cell at horizon \eqn{h} lie within \eqn{\pi_h} times the
#' range of its \eqn{M} errors. For exchangeable errors, the range of
#' \eqn{M} errors contains a further one with probability
#' \eqn{(M - 1) / (M + 1)}: 0.82 for 10 origins, 0.33 for two. Horizons that
#' few origins reach therefore get narrow intervals; a long `years` with
#' `n_origins = years - 1` gives each horizon more origins.
#'
#' **Package-specific options** (reference unknown): the power law of the
#' spread factor over the horizon (`scale = "horizon"`), the correction of
#' the drift (`bias = TRUE`), separate factors for age bands (`ages`) and
#' calibration on life expectancy or the age-standardised rate (`quantity`).
#' The drift correction and the age bands have published analogues of
#' another form. EMOS corrects the bias of the forecasts in its predictive
#' mean, a linear function of the forecasts with an intercept, estimated
#' together with the spread by minimum CRPS (Gneiting et al. 2005, eq. 10
#' and section 2c); `bias = TRUE` shifts the
#' log rates by a drift per year of the horizon, one drift for all horizons
#' (with `scale = "per_horizon"` one shift per horizon, a bias correction
#' per lead time). Member-by-member post-processing is done for each
#' location separately (Van Schaeybroeck and Vannitsem 2015), and Gneiting
#' et al. (2005, section 3a) note that distinct regions may need distinct
#' training periods; `ages` estimates separate factors for bands of ages
#' that the user chooses.
#'
#' **Estimation.** The parameters are bounded: \eqn{s_1} (and \eqn{s_h})
#' between 0.1 and 20, \eqn{\gamma} between -1 and 1, and \eqn{b} between
#' -0.05 and 0.05 (log rate per year). With one parameter (`scale =
#' "constant"` or, for each horizon, `"per_horizon"`, with `bias = FALSE`)
#' the log factor is found by [stats::optimize()] within its bounds;
#' otherwise the parameters are found by Nelder-Mead ([stats::optim()],
#' relative tolerance 1e-6, at most 300 iterations), started at
#' \eqn{s_1 = 1}, \eqn{\gamma = 0}, \eqn{b = 0}, with values outside the
#' bounds rejected. The observation noise of the scored forecasts is
#' simulated once (common random numbers), so the score is a deterministic
#' function of the parameters.
#'
#' **Long horizons.** The parameters are estimated on horizons 1 to
#' \eqn{h_{max}}, the distance from the first origin to the last training
#' year (10 years with the defaults). With `scale = "constant"` or
#' `"horizon"`, the spread factor and the drift correction are held at
#' their values at \eqn{h_{max}} beyond it. A forecast more than twice as
#' long as \eqn{h_{max}} gives a warning (class `mt_warning_recalibrate`),
#' because nothing in the training data shows how the spread of the model
#' is wrong at such horizons. With the defaults this concerns forecasts of
#' more than 20 years (`mt_forecast()` forecasts 10 by default); give a
#' longer `years` for long projections, or accept the constant factor. The
#' per-horizon methods stop instead (see above).
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
#' 95% percentile intervals (with `scale = "per_horizon"`, for each horizon
#' from its forecasts). The inner models are fitted to shorter series
#' than the final model, which matters for models whose uncertainty depends
#' on the length of the series (the drift of a random walk). Compare
#' calibrated and uncalibrated models in an outer backtest.
#'
#' **Age bands** (package-specific). With `ages` (break points, for example
#' `c(50, 75)`), the factors are estimated separately for each age band, so
#' that a model that is too narrow at old ages and too wide at young ages is
#' corrected at both. Observed cells belong to the band of their lower age;
#' use break points that are bounds of the observed age groups.
#'
#' **Quantity** (package-specific). By default the log death rates of single
#' cells are scored. Intervals of life expectancy and ASR also depend on how
#' the ages move together, which the transformation keeps from the base model;
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
#' `mt_score_table(by = c("model", "horizon"))` gives the coverage of each
#' horizon.
#'
#' `target = "crps"` (default) balances calibration and sharpness;
#' `target = "interval"` minimises the interval score of the central `level`
#' interval (Gneiting and Raftery 2007), whose optimum is the correct
#' coverage of that interval.
#'
#' @param spec A model specification.
#' @param years Length of the calibration window (the last observed years).
#' @param n_origins Number of origins in the window (evenly spaced from
#'   `years` to 2 years before the end of the training data; at most
#'   `years - 1`). The per-horizon methods need two origins at each
#'   forecast horizon (see "Factors per horizon").
#' @param scale `"constant"` (default: one spread factor \eqn{s_1} for all
#'   horizons), `"per_horizon"` (one factor for each horizon) or the
#'   deprecated `"horizon"` (spread factor \eqn{s_1 h^{\gamma}},
#'   package-specific; use `"per_horizon"`).
#' @param bias Estimate a correction of the drift (log rate per year;
#'   package-specific; with `scale = "per_horizon"` one for each horizon).
#' @param target `"crps"` or `"interval"` (interval score at `level`).
#' @param level Level of the central interval for `target = "interval"`,
#'   and of the interval whose coverage of the inner errors is tuned with
#'   `method = "errors"`.
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
#' @param method `"spread"` (default: rescale the spread of the draws) or
#'   `"errors"` (point forecast plus resampled errors of the inner
#'   forecasts; see "Resampled past errors"). `method = "errors"` takes
#'   `years`, `n_origins`, `level`, `n_draws` and `seed`; the other
#'   arguments keep their defaults.
#' @return A model specification. The fitted object has a `calibration`
#'   element: a `data.table` with one row per group and age band
#'   (`scale_h1`, `slope`, `drift`, the number of observed cells and of
#'   scored values, the number of inner forecasts (stratum x origin), CRPS
#'   and coverage of the inner forecasts before and after calibration, and
#'   with `n_boot` the intervals `scale_h1_lo`, `scale_h1_hi`, ...). With
#'   `scale = "per_horizon"` it has one row per group, age band and horizon,
#'   with `horizon`, `origins` (inner origins of the horizon), `scale` and
#'   `drift` in place of `scale_h1`, `slope` and `drift` (and `scale_lo`,
#'   `scale_hi`, ... with `n_boot`). With `method = "errors"` it has one row
#'   per horizon: `horizon`, `origins` (error vectors), `cells` (ages x
#'   strata of each vector), `mean_error`, `tuning` (\eqn{\pi_h}) and
#'   `coverage` (share of the inner errors within the tuned interval). For
#'   the per-horizon methods, `calibration_hmax` is the longest horizon that
#'   forecasts can reach.
#' @references
#' Van Schaeybroeck B, Vannitsem S (2015). Ensemble post-processing using
#' member-by-member approaches: theoretical aspects. *Quarterly Journal of
#' the Royal Meteorological Society* 141(688): 807-818.
#' \doi{10.1002/qj.2397}
#'
#' Gneiting T, Raftery AE, Westveld AH, Goldman T (2005). Calibrated
#' probabilistic forecasting using ensemble model output statistics and
#' minimum CRPS estimation. *Monthly Weather Review* 133(5): 1098-1118.
#' \doi{10.1175/MWR2904.1}
#'
#' Gneiting T, Raftery AE (2007). Strictly proper scoring rules, prediction,
#' and estimation. *Journal of the American Statistical Association*
#' 102(477): 359-378. \doi{10.1198/016214506000001437}
#'
#' Keilman N, Pham DQ, Hetland A (2002). Why population forecasts should be
#' probabilistic: illustrated by the case of Norway. *Demographic Research*
#' 6(15): 409-454. \doi{10.4054/DemRes.2002.6.15}
#'
#' Panagiotelis A, Gamakumara P, Athanasopoulos G, Hyndman RJ (2023).
#' Probabilistic forecast reconciliation: properties, evaluation and score
#' optimisation. *European Journal of Operational Research* 306(2):
#' 693-706. \doi{10.1016/j.ejor.2022.07.040}
#'
#' Shang HL, Hyndman RJ (2017). Grouped functional time series forecasting:
#' an application to age-specific mortality rates. *Journal of
#' Computational and Graphical Statistics* 26(2): 330-343.
#' \doi{10.1080/10618600.2016.1237877}
#' @examples
#' x <- mt_example(deaths = "single", years = 1990:2019)
#' fit <- mt_fit(recalibrated(lee_carter(), years = 8, n_origins = 3, n_draws = 50), x)
#' fit$calibration
#' fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
#' # one factor per horizon: horizons 1 to 5 rest on at least two origins
#' fh <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 5, n_draws = 50,
#'                           scale = "per_horizon"), x)
#' fh$calibration[, c("horizon", "origins", "scale")]
#' # point forecast plus resampled errors of the inner forecasts
#' fe <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 5, n_draws = 50,
#'                           method = "errors"), x)
#' fe$calibration
#' fce <- mt_forecast(fe, h = 5, n_draws = 100, seed = 1)
#' @export
recalibrated <- function(spec, years = 10, n_origins = 5, scale = c("constant", "horizon", "per_horizon"),
                         bias = FALSE, target = c("crps", "interval"), level = 0.95, by = NULL, n_draws = 100,
                         seed = 1, ages = NULL, quantity = c("rate", "life_expectancy", "asr"), n_boot = 0,
                         method = c("spread", "errors")) {
  check_class(spec, "mt_spec", "spec")
  scale <- match.arg(scale)
  if (identical(scale, "horizon")) {
    mt_warn(paste("`scale = \"horizon\"` (a power law of the factor over the horizon, package-specific) is",
                  "deprecated; use `scale = \"per_horizon\"`, one factor for each horizon (Gneiting et al. 2005)."),
            class = "deprecated")
  }
  target <- match.arg(target)
  quantity <- match.arg(quantity)
  method <- match.arg(method)
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
  if (method == "errors") {
    unused <- c(scale = scale != "constant", bias = bias, target = target != "crps", by = !is.null(by),
                ages = !is.null(ages), quantity = quantity != "rate", n_boot = n_boot > 0)
    if (any(unused)) {
      mt_abort(sprintf(paste("`method = \"errors\"` resamples the errors of past forecasts and has no spread factors:",
                             "leave out %s (`level` sets the interval whose inner coverage is tuned)."),
                       paste0("`", names(unused)[unused], "`", collapse = ", ")), class = "argument")
    }
  }
  if (n_origins > years - 1) n_origins <- years - 1
  if ((scale == "per_horizon" || method == "errors") && n_origins < 2) {
    mt_abort(paste("Factors or errors per horizon need inner forecasts of at least two origins at each horizon:",
                   "give `n_origins` of 2 or more (and `years` of 3 or more)."), class = "argument")
  }
  label <- if (method == "errors") {
    paste0(spec$label, ", recalibrated by past errors")
  } else {
    paste0(spec$label, ", recalibrated", if (scale == "per_horizon") " per horizon",
           if (length(by)) paste0(" by ", paste(by, collapse = " x ")), if (length(ages)) " by age band",
           if (quantity != "rate") sprintf(" on %s", gsub("_", " ", quantity)),
           if (target == "interval") sprintf(" (%g%% interval)", 100 * level))
  }
  # the draws of the error method are point forecasts plus errors of observed
  # rates: they describe observed rates, and predictive draws add nothing
  new_mt_spec("recalibrated", list(base = spec, years = years, n_origins = n_origins, scale = scale, bias = bias,
                                 target = target, level = level, by = by, n_draws = n_draws, seed = seed,
                                 ages = ages, quantity = quantity, n_boot = n_boot, method = method),
              joint = spec$joint, min_periods = spec$min_periods + years,
              smooth = if (method == "errors") FALSE else spec$smooth, label = label)
}

# Does a spec describe draws of the error method (observed rates)?
is_error_recalibrated <- function(spec) {
  inherits(spec, "mt_spec_recalibrated") && identical(spec$params$method, "errors")
}

# Per-horizon methods: the factor (or the errors) of a horizon need inner
# forecasts of at least two origins
per_horizon_method <- function(p) identical(p$scale, "per_horizon") || identical(p$method, "errors")

# Largest horizon that the inner forecasts of at least two origins reach
horizon_reach <- function(origins, last) {
  r <- sort(last - origins, decreasing = TRUE)
  if (length(r) < 2) 0 else r[2]
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
  ok <- origins[!bad]
  if (per_horizon_method(p)) {
    p$h_reach <- horizon_reach(ok, last)
    if (p$h_reach < 1) {
      mt_abort(sprintf(paste("%s needs inner forecasts of at least two origins at each horizon, but %d origin(s)",
                             "gave forecasts. Give more `n_origins` (and a longer `years`)."),
                       if (identical(p$method, "errors")) "`method = \"errors\"`" else "`scale = \"per_horizon\"`",
                       length(ok)), class = "calibration")
    }
    if (identical(p$method, "errors")) {
      est <- estimate_errors(pieces, fit, ok, last, p)
      return(wrap_fit(fit, spec, "mt_fit_recalibrated", list(calibration = est$table, calibration_errors = est,
                                                             calibration_origins = ok,
                                                             calibration_hmax = p$h_reach)))
    }
  } else if (sum(!bad) == 1) {
    mt_warn(paste("The calibration rests on the forecasts of one origin: its horizons are not replicated and the",
                  "factors are uncertain. Give a longer `years` or more `n_origins`."), class = "calibration")
  }
  est <- estimate_calibration(pieces, fit$strata, p)
  wrap_fit(fit, spec, "mt_fit_recalibrated", list(calibration = est$table, calibration_theta = est$theta,
                                                  calibration_origins = ok,
                                                  calibration_hmax = p$h_reach %||% (last - min(ok))))
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

# Does the predictive distribution add Poisson (or negative binomial) deaths
# to the rate draws (models of expected rates), or are the draws observed
# rates?
spec_adds_noise <- function(spec) {
  # draws of recalibrated(method = "errors") contain the observation noise
  if (is_error_recalibrated(spec)) return(FALSE)
  if (inherits(spec, "mt_spec_ensemble")) return(any(vapply(spec$params$components, spec_adds_noise, TRUE)))
  if (!is.null(spec$params$base) && inherits(spec$params$base, "mt_spec")) return(spec_adds_noise(spec$params$base))
  isTRUE(spec$smooth)
}

# The same for a fit of a wrapper or a single model: the innermost of the
# fitted specs decides (tuned() fits the spec it chose, which can differ from
# the spec the wrapper was given)
fit_adds_noise <- function(fit) {
  specs <- c(list(fit$spec), fit$inner_specs)
  # draws of recalibrated(method = "errors"), also inside other wrappers,
  # describe observed rates
  if (any(vapply(specs, is_error_recalibrated, logical(1)))) return(FALSE)
  # any model fitted to smoothed rates (step_smooth()) forecasts smooth rates
  spec_adds_noise(specs[[length(specs)]]) || isTRUE(fit$smoothed_input)
}

# Observation noise of each draw: whether deaths are drawn given the rates
# (`noise(s)`: one flag per draw for stratum s) and the size of their
# negative binomial distribution (`size(s)`: one per draw; Inf for Poisson).
# The sizes are those of the posterior draws of `fc` (stanmomo(family =
# "nb")), paired by position, or else drawn from the posterior of the fit.
# Draws of an ensemble follow the component they come from (`fc$component`:
# one per draw, or [stratum, draw] after reordering; without it, draws are
# allocated to the components in proportion to the weights); strata of a fit
# split into groups (stratified(), tuned(by = )) follow the fit of their
# group.
draw_noise <- function(fit, fc, n_draws) {
  # a reconciled bottom-up model fitted to the bottom strata (also below
  # other wrappers): the noise of its base fit, stratum s of all nodes being
  # bottom stratum pos[s] of the base fit
  base <- bottom_base_fit(fit)
  if (!is.null(base)) {
    b <- fit$reconcile$structure$bottom
    fb <- fc
    if (is.matrix(fc$component)) fb$component <- fc$component[b, , drop = FALSE]
    if (is.matrix(fc$nb_size)) fb$nb_size <- fc$nb_size[b, , drop = FALSE]
    parts <- draw_noise(base, fb, n_draws)
    pos <- match(seq_len(nrow(fit$strata)), b)
    return(list(noise = function(s) parts$noise(pos[s]), size = function(s) parts$size(pos[s])))
  }
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
    }, size = function(s) {
      ks <- of(s)
      out <- rep(Inf, n_draws)
      for (j in unique(ks)) out[ks == j] <- parts[[j]]$size(s)[ks == j]
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
    return(list(noise = function(s) parts[[grp[s]]]$noise(pos[s]),
                size = function(s) parts[[grp[s]]]$size(pos[s])))
  }
  flag <- fit_adds_noise(fit)
  size <- if (!is.null(fc$nb_size) && ncol(fc$nb_size) == n_draws) fc$nb_size else sample_nb_size(fit, n_draws)
  list(noise = function(s) rep(flag, n_draws), size = function(s) {
    if (flag && !is.null(size)) size[s, ] else rep(Inf, n_draws)
  })
}

# Sizes of the negative binomial distribution of the deaths [stratum, draw]
# drawn from the posterior of fits with that likelihood (stanmomo(family =
# "nb")); NULL for other models
sample_nb_size <- function(fit, n_draws) {
  if (isTRUE(fit$spec$joint) || !is.list(fit$fits)) return(NULL)
  post <- lapply(fit$fits, function(f) if (is.list(f) && !inherits(f, "mt_fit")) f$nb_size)
  if (all(vapply(post, is.null, logical(1)))) return(NULL)
  do.call(rbind, lapply(post, function(p) {
    if (is.null(p)) rep(Inf, n_draws) else p[sample.int(length(p), n_draws, replace = TRUE)]
  }))
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
      h = rep(fc$periods - origin, each = na)[used], grid_age = rep(gt$age$lower, np)[used],
      E = Eg, C = C, Ecell = Ecell, y = log(pmax(o$value, 0.5) / Ecell), noise = nz$noise(sf),
      size = nz$size(sf), obs_age = o$age_lo,
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
      out[[key]] <- list(K = NULL, E = NULL, qo = numeric(), ax = ax, quantity = quantity,
                         standard = if (quantity == "asr") standard_population("esp2013"))
    }
    r_obs <- o$value[k] / Ecell[k]
    g <- out[[key]]
    # ASR by direct standardisation on the standard's age groups with the
    # exposures of the test year, as in mt_backtest()
    qo <- if (quantity == "life_expectancy") life_expectancy(r_obs, g$ax) else
      log(asr_columns(matrix(r_obs), g$ax, g$standard, Ecell[k], per = 1))
    out[[key]]$K <- cbind(g$K, k)
    out[[key]]$E <- cbind(g$E, Ecell[k])
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
         y = piece$y[rows], noise = piece$noise, size = piece$size, obs_age = piece$obs_age[rows])
  })
}

band_labels <- function(breaks) {
  if (!length(breaks)) return("all ages")
  c(paste0("<", breaks[1]), if (length(breaks) > 1) paste0(breaks[-length(breaks)], "-", breaks[-1]),
    paste0(breaks[length(breaks)], "+"))
}

# The cells of a piece at each horizon, named by the horizon (an observed
# cell lies in one year, so its grid cells share their horizon); life
# expectancy or ASR groups keep the years of the horizon
split_horizons <- function(piece) {
  ho <- as.vector(piece$C %*% piece$h) / Matrix::rowSums(piece$C)
  hs <- sort(unique(ho))
  out <- lapply(hs, function(hv) {
    rows <- which(ho == hv)
    cols <- which(piece$h == hv)
    sub <- piece
    sub$L <- piece$L[cols, , drop = FALSE]
    for (f in c("m", "h", "E", "grid_age")) if (!is.null(piece[[f]])) sub[[f]] <- piece[[f]][cols]
    for (f in c("Ecell", "y", "obs_age")) sub[[f]] <- piece[[f]][rows]
    sub$C <- piece$C[rows, cols, drop = FALSE]
    if (!is.null(piece$qg)) sub$qg <- subset_quantity_groups(piece$qg, rows)
    sub
  })
  stats::setNames(out, hs)
}

# Life expectancy or ASR groups restricted to the years whose cells are all
# among `rows`, with the cell indices renumbered to positions in `rows`
subset_quantity_groups <- function(qg, rows) {
  out <- lapply(qg, function(g) {
    inside <- matrix(g$K %in% rows, nrow(g$K))
    keep <- which(colSums(inside) == nrow(g$K))
    if (!length(keep)) return(NULL)
    g$K <- matrix(match(g$K[, keep], rows), nrow(g$K))
    g$E <- g$E[, keep, drop = FALSE]
    g$qo <- g$qo[keep]
    g
  })
  Filter(Negate(is.null), out)
}

# One spread factor (and with `bias` one drift) per horizon 1 to p$h_reach,
# each estimated like the constant factor from the inner forecasts of that
# horizon only
calibration_fit_per_horizon <- function(pieces, p, n_draws, where = "") {
  sub <- lapply(pieces, split_horizons)
  pc <- p
  pc$scale <- "constant"
  res <- lapply(seq_len(p$h_reach), function(hv) {
    ph <- Filter(Negate(is.null), lapply(sub, `[[`, as.character(hv)))
    n_o <- length(unique(vapply(ph, function(x) x$origin, 0)))
    if (n_o < 2) {
      mt_abort(sprintf(paste("Horizon %d has observed cells from %d inner origin(s)%s; a factor per horizon needs",
                             "two. Check for missing observations, or give more `n_origins`."), hv, n_o, where),
               class = "calibration")
    }
    r <- calibration_fit(stack_pieces(ph, n_draws), "constant", p$bias, p$target %||% "crps", p$level %||% 0.95)
    tab <- r$table
    set(tab, j = "slope", value = NULL)
    setnames(tab, "scale_h1", "scale")
    tab <- cbind(data.table(horizon = hv, origins = n_o), tab, data.table(forecasts = length(ph)))
    if ((p$n_boot %||% 0) > 0) {
      bt <- calibration_boot(ph, pc, n_draws)
      setnames(bt, c("scale_h1_lo", "scale_h1_hi"), c("scale_lo", "scale_hi"))
      tab <- cbind(tab, bt)
    }
    list(theta = r$theta, table = tab)
  })
  list(theta = do.call(rbind, lapply(res, `[[`, "theta")), table = rbindlist(lapply(res, `[[`, "table")))
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
  # negative binomial size per observed cell and draw (Inf: Poisson)
  P$size <- do.call(rbind, lapply(seq_along(pieces), function(i) {
    matrix(pieces[[i]]$size %||% Inf, N[i], n_draws, byrow = TRUE)
  }))
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
    q <- if (g$quantity == "life_expectancy") lt_columns(M, g$ax)$ex[1, ] else
      log(asr_columns(M, g$ax, g$standard, g$E, per = 1))
    matrix(q, ncol(g$K))
  }))
  list(S = S, obs = unlist(lapply(P$qg, `[[`, "qo")))
}

# Spread factor and drift correction per grid cell; beyond the longest
# horizon of the calibration window both stay constant (the drift correction
# is a level shift there, not extrapolated further)
calibration_shift <- function(theta, h, scale, bias, h_max = Inf) {
  if (identical(scale, "per_horizon")) {
    # one row of parameters per horizon; no factor beyond the last row
    th <- matrix(theta, ncol = 1 + bias)
    if (any(h > nrow(th))) mt_abort(sprintf("No spread factor for horizon %d.", max(h)), class = "recalibrate")
    return(list(s = exp(th[h, 1]), b = if (bias) th[h, 2] * h else 0))
  }
  hc <- pmin(h, h_max)
  s <- exp(theta[1] + if (scale == "horizon") theta[2] * log(hc) else 0)
  b <- if (bias) theta[length(theta)] * hc else 0
  list(s = s, b = b)
}

# Log rates of the observed cells under parameters theta, with observation noise
calibration_sample <- function(theta, P, scale, bias) {
  sh <- calibration_shift(theta, P$h, scale, bias)
  Mu <- exp(P$m + sh$b + sh$s * (P$L - P$m)) * P$E
  Lam <- as.matrix(P$C %*% Mu)
  D <- Lam
  # negative binomial deaths (finite size): an observed cell that sums grid
  # cells gets the size that matches the variance of the sum
  V2 <- if (any(is.finite(P$size))) as.matrix(P$C %*% Mu^2) / P$size else 0
  nb <- P$pois & V2 > 0
  small <- P$pois & !nb & Lam < 50
  large <- P$pois & !nb & !small
  D[small] <- stats::qpois(P$U[small], Lam[small])
  D[large] <- pmax(Lam[large] + sqrt(Lam[large]) * P$Z[large], 0)
  if (any(nb)) D[nb] <- stats::qnbinom(P$U[nb], size = Lam[nb]^2 / V2[nb], mu = Lam[nb])
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
  res <- with_seed(p$seed, lapply(stats::setNames(nm = names(groups)), function(g) {
    gp <- groups[[g]]
    bands <- if (length(p$ages)) {
      sub <- lapply(gp, split_bands, breaks = p$ages)
      labs <- intersect(band_labels(p$ages), unlist(lapply(sub, names)))
      stats::setNames(lapply(labs, function(b) Filter(Negate(is.null), lapply(sub, `[[`, b))), labs)
    } else {
      list("all ages" = gp)
    }
    lapply(stats::setNames(nm = names(bands)), function(b) {
      bp <- bands[[b]]
      if (identical(p$scale, "per_horizon")) {
        return(calibration_fit_per_horizon(bp, p, n_draws, where = sprintf(" (group %s, %s)", g, b)))
      }
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

# Errors of the inner forecasts for recalibrated(method = "errors"): for
# each horizon 1 to p$h_reach an array V [age, stratum, origin] of observed
# minus forecast log rates (the median of the inner draws), with every age
# of the grid and every calibrated stratum, and the tuning factor of
# Shang and Hyndman (2017) for the central `level` interval
estimate_errors <- function(pieces, fit, origins, last, p) {
  labs <- fit$strata$label[calibration_targets(fit)]
  ages <- fit$age$lower
  na <- length(ages)
  ns <- length(labs)
  for (pc in pieces) {
    if (any(Matrix::rowSums(pc$C) != 1)) {
      mt_abort(paste("`method = \"errors\"` needs the observed deaths on the age groups of the model grid (one grid",
                     "cell per observed cell); the observed cells are coarser than the grid."), class = "calibration")
    }
  }
  per_h <- lapply(seq_len(p$h_reach), function(hv) {
    os <- sort(origins[origins + hv <= last])
    V <- array(NA_real_, c(na, ns, length(os)), dimnames = list(age = fit$age$labels, stratum = labs, origin = os))
    for (pc in pieces) {
      k <- match(pc$origin, os)
      s <- match(pc$label, labs)
      if (is.na(k) || is.na(s)) next
      # the grid cell of each observed cell
      j <- as.vector(pc$C %*% seq_len(ncol(pc$C)))
      rows <- which(pc$h[j] == hv)
      a <- match(round(pc$grid_age[j[rows]], 8), round(ages, 8))
      V[a, s, k] <- pc$y[rows] - pc$m[j[rows]]
    }
    if (anyNA(V)) {
      i <- which(is.na(V), arr.ind = TRUE)[1, ]
      mt_abort(sprintf(paste("`method = \"errors\"` needs observed deaths in every cell of the calibration window;",
                             "the error of horizon %d is missing for age %s, stratum %s, origin %s (year %s)."),
                       hv, fit$age$labels[i[1]], labs[i[2]], os[i[3]], os[i[3]] + hv), class = "calibration")
    }
    tu <- error_tuning(matrix(V, na * ns), p$level)
    list(V = V, origins = os, tuning = tu$tuning, coverage = tu$coverage, mean_error = mean(V))
  })
  table <- data.table(horizon = seq_len(p$h_reach), origins = vapply(per_h, function(e) length(e$origins), 0),
                      cells = na * ns, mean_error = vapply(per_h, `[[`, 0, "mean_error"),
                      tuning = vapply(per_h, `[[`, 0, "tuning"), coverage = vapply(per_h, `[[`, 0, "coverage"))
  list(table = table, V = lapply(per_h, `[[`, "V"), origins = lapply(per_h, `[[`, "origins"),
       tuning = vapply(per_h, `[[`, 0, "tuning"), strata = labs, level = p$level)
}

# Tuning factor of Shang and Hyndman (2017, section 4.4, pointwise
# intervals): `E` [cell, origin] holds the errors of one horizon. The bounds
# of each cell are the quantiles at (1 - level) / 2 and (1 + level) / 2 of
# its bootstrap distribution, the empirical distribution of its errors
# (inverse of the empirical distribution function, the limit of many
# resamples). The factor is the smallest pi > 0 with which a share `level`
# of all errors lies in [pi * lower, pi * upper]; it is at most 1, because
# these bounds contain at least that share of the errors of each cell.
error_tuning <- function(E, level) {
  M <- ncol(E)
  S <- if (M == 1) E else t(apply(E, 1, sort))
  lo <- rep(S[, ceiling((1 - level) / 2 * M)], M)
  hi <- rep(S[, ceiling((1 + level) / 2 * M)], M)
  e <- as.vector(E)
  # the factors pi > 0 that cover each error form an interval [a, b]:
  # pi * lo <= e and e <= pi * hi
  a <- rep(0, length(e))
  b <- rep(Inf, length(e))
  never <- rep(FALSE, length(e))
  a <- ifelse(lo < 0, pmax(a, e / lo), a)
  b <- ifelse(lo > 0, pmin(b, e / lo), b)
  never <- never | (lo == 0 & e < 0)
  a <- ifelse(hi > 0, pmax(a, e / hi), a)
  b <- ifelse(hi < 0, pmin(b, e / hi), b)
  never <- never | (hi == 0 & e > 0)
  ok <- !never & a <= b
  covered <- function(x) {
    (findInterval(x, sort(a[ok])) - findInterval(x, sort(b[ok]), left.open = TRUE)) / length(e)
  }
  cand <- sort(unique(c(a[ok & a > 0], 1)))
  cov <- covered(cand)
  tuning <- cand[which(cov >= level - 1e-12)[1]]
  list(tuning = tuning, coverage = covered(tuning))
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
  hmax <- fit$calibration_hmax %||% Inf
  h_fc <- if (is.null(to)) h else to - max(fit$periods)
  if (per_horizon_method(p) && is.numeric(h_fc) && length(h_fc) == 1 && h_fc > hmax) {
    mt_abort(sprintf(paste("The forecast reaches %d years ahead, but %s exist for horizons 1 to %d only (horizons",
                           "with inner forecasts of at least two origins). Give a longer `years` or more",
                           "`n_origins`: for h years, for example years = h + 1 and n_origins = h."),
                     as.integer(h_fc),
                     if (identical(p$method, "errors")) "past errors" else "spread factors per horizon",
                     as.integer(hmax)), class = "recalibrate")
  }
  if (identical(p$method, "errors")) {
    return(forecast_past_errors(fit, h, to, n_draws, seed, predictive, exposure, h_missing = missing(h)))
  }
  if (is.numeric(h_fc) && length(h_fc) == 1 && is.finite(hmax) && h_fc > 2 * hmax) {
    mt_warn(sprintf(paste("The forecast reaches %d years ahead, more than twice the longest horizon of the",
                          "calibration (%d years); the %s of horizon %d %s used beyond it. Give a longer `years`",
                          "to calibrate longer horizons."), as.integer(h_fc), as.integer(hmax),
                    if (isTRUE(p$bias)) "spread factor and the drift correction" else "spread factor",
                    as.integer(hmax), if (isTRUE(p$bias)) "are" else "is"),
            class = "recalibrate")
  }
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

# recalibrated(method = "errors"): at horizon h, draw j is the median of the
# base draws in each cell plus the tuned error vector of one inner origin,
# drawn with replacement for each horizon and draw; the same origin for all
# ages and strata keeps their dependence. The draws describe observed rates,
# so predictive draws add no observation noise.
forecast_past_errors <- function(fit, h, to, n_draws, seed, predictive, exposure, h_missing) {
  er <- fit$calibration_errors
  transform <- function(arr, targets, fc) {
    d <- dim(arr)
    pick <- lapply(seq_len(d[2]), function(t) sample.int(length(er$origins[[t]]), d[4], replace = TRUE))
    for (s in targets) {
      k <- match(fit$strata$label[s], er$strata)
      L <- matrix(log(pmax(arr[, , s, ], 1e-12)), ncol = d[4])
      if (anyNA(L) || is.na(k)) next
      m <- matrix(row_quantiles(L, 0.5)[, 1], d[1])
      for (t in seq_len(d[2])) {
        arr[, t, s, ] <- exp(m[, t] + er$tuning[t] * er$V[[t]][, k, pick[[t]], drop = FALSE])
      }
    }
    arr
  }
  res <- wrapped_forecast(fit, "mt_fit_recalibrated", h, to, n_draws, seed, predictive = FALSE, exposure,
                          transform, h_missing = h_missing)
  res$predictive <- predictive
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

# Observation noise of the base model added to rate draws: Poisson (or
# negative binomial) deaths given the rates and exposures, divided by the
# exposures. Cells can have zero deaths (rate 0), as in predictive_rates().
# `draw` [age, period, stratum, draw]: the draw of `fc` each value comes
# from, where values were reordered cell by cell (coupled()); by default
# the draw itself
add_observation_noise <- function(arr, base_fit, fc, targets, E, draw = NULL) {
  d <- dim(arr)
  nz <- draw_noise(base_fit, fc, d[4])
  for (s in targets) {
    R <- arr[, , s, ]
    if (anyNA(R) || !any(nz$noise(s))) next
    j <- if (is.null(draw)) rep(seq_len(d[4]), each = d[1] * d[2]) else as.vector(draw[, , s, ])
    Es <- as.vector(E[, , s])
    mu <- R * Es
    D <- r_deaths(mu, nz$size(s)[j]) / pmax(Es, 1e-12)
    k <- nz$noise(s)[j]
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
