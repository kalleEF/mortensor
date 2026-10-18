# Recalibrate forecast distributions on recent years

Wraps any model specification so that the spread (and optionally the
trend) of its forecast draws is recalibrated on the last `years`
observed years of the training data. Each time the model is fitted (also
at every origin of an outer [`mt_backtest()`](mt_backtest.md)):

## Usage

``` r
recalibrated(
  spec,
  years = 10,
  n_origins = 5,
  scale = c("constant", "horizon"),
  bias = FALSE,
  target = c("crps", "interval"),
  level = 0.95,
  by = NULL,
  n_draws = 100,
  seed = 1,
  ages = NULL,
  quantity = c("rate", "life_expectancy", "asr"),
  n_boot = 0
)
```

## Arguments

- spec:

  A model specification.

- years:

  Length of the calibration window (the last observed years).

- n_origins:

  Number of origins in the window (evenly spaced from `years` to 2 years
  before the end of the training data).

- scale:

  `"constant"` (default: one spread factor \\s_1\\ for all horizons) or
  `"horizon"` (spread factor \\s_1 h^{\gamma}\\, package-specific).

- bias:

  Estimate a correction of the drift (log rate per year;
  package-specific).

- target:

  `"crps"` or `"interval"` (interval score at `level`).

- level:

  Level of the central interval for `target = "interval"`.

- by:

  Stratum variables with separate calibration parameters.

- n_draws:

  Draws per inner forecast.

- seed:

  Random seed of the inner forecasts and of the simulated observation
  noise.

- ages:

  Break points between age bands with their own factors (`NULL`: one set
  for all ages).

- quantity:

  What is scored: `"rate"` (log death rates of the observed cells),
  `"life_expectancy"` (at birth) or `"asr"`.

- n_boot:

  Number of bootstrap refits for intervals of the factors (0: none; each
  refit repeats the optimisation, not the inner forecasts).

## Value

A model specification. The fitted object has a `calibration` element: a
`data.table` with one row per group and age band (`scale_h1`, `slope`,
`drift`, the number of observed cells and of scored values, the number
of inner forecasts (stratum x origin), CRPS and coverage of the inner
forecasts before and after calibration, and with `n_boot` the intervals
`scale_h1_lo`, `scale_h1_hi`, ...).

## Details

1.  the base model is refitted at `n_origins` origins inside the
    calibration window (training data up to the origin) and forecasts
    the remaining years of the window;

2.  the draws are transformed per cell on the log scale, \$\$\log m^\* =
    \tilde m + b\\h + s(h)\\(\log m - \tilde m),\$\$ where \\\tilde m\\
    is the median of the draws in the cell, \\h\\ the horizon, \\s(h)\\
    a spread factor, one value \\s_1\\ for all horizons
    (`scale = "constant"`, default) or \\s_1 h^{\gamma}\\
    (`scale = "horizon"`), held constant beyond the longest horizon of
    the window, and \\b\\ a correction of the drift per year
    (`bias = TRUE`);

3.  \\s_1\\, \\\gamma\\ and \\b\\ minimise the mean CRPS (or the
    interval score, `target`) of the log death rates in the window, with
    the Poisson (or the model's negative binomial) variation of the
    observed deaths, as in [`mt_backtest()`](mt_backtest.md). Draws of
    an ensemble get the observation noise of the component they come
    from.

The final model is fitted on all training data and its forecasts are
transformed with the estimated parameters. Ranks within each cell are
kept, so paths stay paths, and derived quantities such as life
expectancy can be computed from the calibrated draws.

**Published methods.** The closest published method is the
member-by-member post-processing of ensemble forecasts (Van Schaeybroeck
and Vannitsem 2015), which shifts every member of an ensemble and
rescales its deviation from the ensemble mean, with parameters for each
lead time estimated by maximum likelihood or minimum CRPS. Here each
simulated path is rescaled around the median of the draws of each cell,
on the log scale, with one spread factor for all horizons
(`scale = "constant"`) estimated by minimum CRPS on inner rolling-origin
forecasts. Minimum-CRPS estimation follows ensemble model output
statistics (EMOS, Gneiting et al. 2005), which fits a parametric
predictive distribution instead of transforming the draws.

**Package-specific options** (reference unknown): the power law of the
spread factor over the horizon (`scale = "horizon"`), the correction of
the drift (`bias = TRUE`), separate factors for age bands (`ages`) and
calibration on life expectancy or the age-standardised rate
(`quantity`).

**Estimation.** The parameters are bounded: \\s_1\\ between 0.1 and 20,
\\\gamma\\ between -1 and 1, and \\b\\ between -0.05 and 0.05 (log rate
per year). With one parameter (`scale = "constant"`, `bias = FALSE`)
\\\log s_1\\ is found by
[`stats::optimize()`](https://rdrr.io/r/stats/optimize.html) within its
bounds; otherwise the parameters are found by Nelder-Mead
([`stats::optim()`](https://rdrr.io/r/stats/optim.html), relative
tolerance 1e-6, at most 300 iterations), started at \\s_1 = 1\\,
\\\gamma = 0\\, \\b = 0\\, with values outside the bounds rejected. The
observation noise of the scored forecasts is simulated once (common
random numbers), so the score is a deterministic function of the
parameters.

**Long horizons.** The parameters are estimated on horizons 1 to
\\h\_{max}\\, the distance from the first origin to the last training
year (10 years with the defaults). Beyond \\h\_{max}\\ the spread factor
and the drift correction are held at their values at \\h\_{max}\\. A
forecast more than twice as long as \\h\_{max}\\ gives a warning (class
`mt_warning_recalibrate`), because nothing in the training data shows
how the spread of the model is wrong at such horizons. With the defaults
this concerns forecasts of more than 20 years
([`mt_forecast()`](mt_forecast.md) forecasts 10 by default); give a
longer `years` for long projections, or accept the constant factor.

The wrapper works with any specification, including
[`ensemble()`](ensemble.md), [`bootstrapped()`](bootstrapped.md) and
joint models. For [`reconciled()`](reconciled.md) models the bottom
strata are calibrated and the aggregates recomputed from them, so the
draws stay coherent. With `by`, separate parameters are estimated for
groups of strata (for example `by = "sex"`); strata outside these groups
(aggregates) use the pooled parameters. Groups with different factors
scale the common variation of joint models (such as the common trend of
[`li_lee()`](li_lee.md)) differently, so it no longer cancels in
contrasts between the groups, and intervals of gaps between them widen;
use pooled parameters when contrasts between the groups are the target.

Calibration needs `years` more training years than the base model. A
spread factor above 1 means that the base model was too narrow in the
window; it widens future intervals in proportion. The parameters are
estimated from few overlapping forecasts (with the defaults, 5 origins
and 30 forecast years, horizons 1 to 10), so they carry their own
uncertainty, and they cannot anticipate changes after the window.
`n_boot` refits them on forecasts (stratum x origin) resampled with
replacement and reports 95% percentile intervals. The inner models are
fitted to shorter series than the final model, which matters for models
whose uncertainty depends on the length of the series (the drift of a
random walk). Compare calibrated and uncalibrated models in an outer
backtest.

**Age bands** (package-specific). With `ages` (break points, for example
`c(50, 75)`), the factors are estimated separately for each age band, so
that a model that is too narrow at old ages and too wide at young ages
is corrected at both. Observed cells belong to the band of their lower
age; use break points that are bounds of the observed age groups.

**Quantity** (package-specific). By default the log death rates of
single cells are scored. Intervals of life expectancy and ASR also
depend on how the ages move together, which the transformation keeps
from the base model; `quantity = "life_expectancy"` (at birth) or
`"asr"` scores these instead, per test year and stratum (fewer scored
values, so the factors are less certain).

**Tuning and recalibration.** `recalibrated()` keeps the model and its
settings and rescales the width of its forecast distribution. Choosing
settings (smoothing, number of components, fitting window) by forecast
accuracy is tuning, done by [`tuned()`](mt_tune.md) and
[`mt_tune()`](mt_tune.md); the two combine as
`recalibrated(tuned(...))`.

Calibration widens the intervals of models that were too narrow in the
recent past and narrows those that were too wide, for example of
ensembles whose intervals already covered. With few training years it
rests on few forecasts and can make intervals worse; check it in a
backtest, as in the vignette "Tuning and recalibration".

`target = "crps"` (default) balances calibration and sharpness;
`target = "interval"` minimises the interval score of the central
`level` interval (Gneiting and Raftery 2007), whose optimum is the
correct coverage of that interval.

## References

Van Schaeybroeck B, Vannitsem S (2015). Ensemble post-processing using
member-by-member approaches: theoretical aspects. *Quarterly Journal of
the Royal Meteorological Society* 141(688): 807-818.
[doi:10.1002/qj.2397](https://doi.org/10.1002/qj.2397)

Gneiting T, Raftery AE, Westveld AH, Goldman T (2005). Calibrated
probabilistic forecasting using ensemble model output statistics and
minimum CRPS estimation. *Monthly Weather Review* 133(5): 1098-1118.
[doi:10.1175/MWR2904.1](https://doi.org/10.1175/MWR2904.1)

Gneiting T, Raftery AE (2007). Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical
Association* 102(477): 359-378.
[doi:10.1198/016214506000001437](https://doi.org/10.1198/016214506000001437)

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2019)
fit <- mt_fit(recalibrated(lee_carter(), years = 8, n_origins = 3, n_draws = 50), x)
fit$calibration
#>     group     band  scale_h1 slope drift cells scored crps_before crps_after
#>    <char>   <char>     <num> <num> <num> <int>  <num>       <num>      <num>
#> 1:    all all ages 0.7630612     0     0  2414   2414  0.03388302 0.03386637
#>    cov80_before cov80_after cov95_before cov95_after forecasts
#>           <num>       <num>        <num>       <num>     <int>
#> 1:    0.7531069   0.7311516    0.8956089   0.8860812         6
fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
```
