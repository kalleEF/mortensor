# Recalibrate forecast distributions on recent years

Wraps any model specification so that its forecast distribution is
recalibrated on the last `years` observed years of the training data:
the spread (and optionally the trend) of its draws is rescaled
(`method = "spread"`, default), or the draws are replaced by the point
forecast of the model plus resampled errors of its past forecasts
(`method = "errors"`). Each time the model is fitted (also at every
origin of an outer [`mt_backtest()`](mt_backtest.md)):

## Usage

``` r
recalibrated(
  spec,
  years = 10,
  n_origins = 5,
  scale = c("constant", "per_horizon"),
  bias = FALSE,
  target = c("crps", "interval"),
  level = 0.95,
  by = NULL,
  n_draws = 100,
  seed = 1,
  ages = NULL,
  quantity = c("rate", "life_expectancy", "asr"),
  n_boot = 0,
  method = c("spread", "errors")
)
```

## Arguments

- spec:

  A model specification.

- years:

  Length of the calibration window (the last observed years).

- n_origins:

  Number of origins in the window (evenly spaced from `years` to 2 years
  before the end of the training data; at most `years - 1`). The
  per-horizon methods need two origins at each forecast horizon (see
  "Factors per horizon").

- scale:

  `"constant"` (default: one spread factor \\s_1\\ for all horizons) or
  `"per_horizon"` (one factor for each horizon).

- bias:

  Estimate a correction of the drift (log rate per year;
  package-specific; with `scale = "per_horizon"` one for each horizon).

- target:

  `"crps"` or `"interval"` (interval score at `level`).

- level:

  Level of the central interval for `target = "interval"`, and of the
  interval whose coverage of the inner errors is tuned with
  `method = "errors"`.

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

- method:

  `"spread"` (default: rescale the spread of the draws) or `"errors"`
  (point forecast plus resampled errors of the inner forecasts; see
  "Resampled past errors"). `method = "errors"` takes `years`,
  `n_origins`, `level`, `n_draws` and `seed`; the other arguments keep
  their defaults.

## Value

A model specification. The fitted object has a `calibration` element: a
`data.table` with one row per group and age band (`scale_h1`, `drift`,
the number of observed cells and of scored values, the number of inner
forecasts (stratum x origin), CRPS and coverage of the inner forecasts
before and after calibration, and with `n_boot` the intervals
`scale_h1_lo`, `scale_h1_hi`, ...). With `scale = "per_horizon"` it has
one row per group, age band and horizon, with `horizon`, `origins`
(inner origins of the horizon), `scale` and `drift` in place of
`scale_h1` and `drift` (and `scale_lo`, `scale_hi`, ... with `n_boot`).
With `method = "errors"` it has one row per horizon: `horizon`,
`origins` (error vectors), `cells` (ages x strata of each vector),
`mean_error`, `tuning` (\\\pi_h\\) and `coverage` (share of the inner
errors within the tuned interval). For the per-horizon methods,
`calibration_hmax` is the longest horizon that forecasts can reach.

## Details

1.  the base model is refitted at `n_origins` origins inside the
    calibration window (training data up to the origin) and forecasts
    the remaining years of the window;

2.  the draws are transformed per cell on the log scale, \$\$\log m^\* =
    \tilde m + b\\h + s(h)\\(\log m - \tilde m),\$\$ where \\\tilde m\\
    is the median of the draws in the cell, \\h\\ the horizon, \\s(h)\\
    a spread factor, one value \\s_1\\ for all horizons
    (`scale = "constant"`, default) or one value \\s_h\\ for each
    horizon (`scale = "per_horizon"`), and \\b\\ a correction of the
    drift per year (`bias = TRUE`; with `scale = "per_horizon"` one
    \\b_h\\ for each horizon);

3.  the parameters minimise the mean CRPS (or the interval score,
    `target`) of the log death rates in the window, with the Poisson (or
    the model's negative binomial) variation of the observed deaths, as
    in [`mt_backtest()`](mt_backtest.md). Draws of an ensemble get the
    observation noise of the component they come from.

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

**Factors per horizon** (`scale = "per_horizon"`). The published
post-processing methods fit their parameters to the forecasts of one
lead time: Gneiting et al. (2005) fit EMOS to 48-hour forecasts (section
3) and note that distinct forecast horizons may need distinct training
periods (section 3a), and member-by-member post-processing is done for
each lead time, location and variable separately (Van Schaeybroeck and
Vannitsem 2015). The factor \\s_h\\ minimises the mean CRPS (or the
interval score) of the inner forecasts of horizon \\h\\ alone, by the
minimum-CRPS estimation of the constant factor (Gneiting et al. 2005,
section 2c); with `bias = TRUE` each horizon also gets its own shift
\\b_h h\\. Where the spread of the model is wrong by the same factor at
every horizon, the factors of all horizons estimate the constant factor,
each from fewer forecasts.

The factor of horizon \\h\\ rests on the inner origins at least \\h\\
years before the last training year. A horizon gets a factor only where
at least two origins reach it: the forecast of one origin has the errors
of a single year, which all its ages and strata share. When every origin
has enough training years, the longest horizon with a factor is
`years - floor(years / n_origins)` (8 with the defaults). The fit
reports it as `calibration_hmax`, and the number of origins of each
horizon in column `origins`. A forecast further ahead stops with an
error of class `mt_error_recalibrate`: no factor is extrapolated. For
forecasts of \\h\\ years, use for example `years = h + 1` and
`n_origins = h` (an origin in every year from \\h + 1\\ to 2 years
before the last training year); horizon \\h\\ then rests on two origins
and horizon 1 on \\h\\.

**Resampled past errors** (`method = "errors"`). The draws are the point
forecast of the base model plus errors of its own past forecasts. The
inner forecasts are those of step 1, so each inner fit sees the data up
to its origin only. The error vector \\e\_{o,h}\\ of origin \\o\\ and
horizon \\h\\ holds the observed log rates of year \\o + h\\ minus the
median of the inner draws, for every age of the grid and every stratum.
Draw \\j\\ of the forecast at horizon \\h\\ is \$\$\log m^\*\_j = \tilde
m + \pi_h\\ e\_{o_j,h},\$\$ where \\\tilde m\\ is the median of the
draws of the base model in the cell (its point forecast) and the origin
\\o_j\\ is drawn with replacement from the inner origins that reach
horizon \\h\\, for each draw and horizon.

- Keilman, Pham and Hetland (2002, section 2) describe the extrapolation
  of errors of past forecasts: a deterministic forecast gives the
  expected value, and the errors of historical forecasts give the
  distribution around it, which they summarise for each forecast
  duration (section 3.1.1.2).

- Shang and Hyndman (2017, section 4.4 of the preprint, pointwise
  intervals) build intervals from h-step errors, observed minus forecast
  values, resampled with replacement: lower and upper bounds from the
  resampled errors of each cell are scaled by a tuning factor \\\pi_h\\,
  the one with which a share `level` of all errors of the horizon lies
  within the scaled bounds. Here the bounds are the quantiles at \\(1 -
  level) / 2\\ and \\(1 + level) / 2\\ of the bootstrap distribution
  itself (the empirical distribution of the errors of the cell, the
  limit of many resamples), and \\\pi_h\\ is the smallest such factor.
  These bounds contain a share `level` of the errors of each cell, so
  \\\pi_h \le 1\\; with fewer than \\1 / (1 - level)\\ origins at a
  horizon (20 at 95%), they are the smallest and the largest error and
  \\\pi_h = 1\\.

- One origin per draw and horizon gives all ages and strata the errors
  of the same past forecast, which keeps their dependence: the joint
  bootstrap of Panagiotelis et al. (2023).

- The errors are used as they are, as by Shang and Hyndman (2017), so a
  mean error of the inner forecasts shifts the forecast distribution
  away from the point forecast. Keilman, Pham and Hetland (2002,
  section 2) assume an expected error of zero.

Deviations and package-specific parts (reference unknown): Shang and
Hyndman (2017) forecast with principal components estimated from all
years; here the model is refitted at each origin, as in the historical
forecasts of Keilman, Pham and Hetland (2002). The origins are those of
`years` and `n_origins` (Shang and Hyndman use every year). The point
forecast is the median of the base draws in each cell. The origins are
drawn independently for each horizon, so a draw combines the errors of
different origins in different years: quantities that combine years,
such as cohort measures or deaths summed over years, miss the dependence
of the errors over time. Panagiotelis et al. (2023) note that multi-step
forecasts need block-wise bootstraps to keep such serial dependence. The
errors are those of the observed log rates, with half a death in empty
cells as in [`mt_backtest()`](mt_backtest.md), so the draws describe
observed rates and predictive draws add nothing. The method needs
observed deaths on the age groups of the grid in every cell of the
window, and at least two origins at each forecast horizon, as
`scale = "per_horizon"`; otherwise the fit or the forecast stops with an
error. The method gives no coverage guarantee: its intervals cover
future errors as far as these resemble the errors of the inner
forecasts. Keilman, Pham and Hetland (2002) name two limits of past
errors: their series are short, and errors may shrink over successive
forecast rounds. The draws of a cell at horizon \\h\\ lie within
\\\pi_h\\ times the range of its \\M\\ errors. For exchangeable errors,
the range of \\M\\ errors contains a further one with probability
\\(M - 1) / (M + 1)\\: 0.82 for 10 origins, 0.33 for two. Horizons that
few origins reach therefore get narrow intervals; a long `years` with
`n_origins = years - 1` gives each horizon more origins.

**Package-specific options** (reference unknown): the correction of the
drift (`bias = TRUE`), separate factors for age bands (`ages`) and
calibration on life expectancy or the age-standardised rate
(`quantity`). The drift correction and the age bands have published
analogues of another form. EMOS corrects the bias of the forecasts in
its predictive mean, a linear function of the forecasts with an
intercept, estimated together with the spread by minimum CRPS (Gneiting
et al. 2005, eq. 10 and section 2c); `bias = TRUE` shifts the log rates
by a drift per year of the horizon, one drift for all horizons (with
`scale = "per_horizon"` one shift per horizon, a bias correction per
lead time). Member-by-member post-processing is done for each location
separately (Van Schaeybroeck and Vannitsem 2015), and Gneiting et al.
(2005, section 3a) note that distinct regions may need distinct training
periods; `ages` estimates separate factors for bands of ages that the
user chooses.

**Estimation.** The parameters are bounded: \\s_1\\ (and \\s_h\\)
between 0.1 and 20, and \\b\\ between -0.05 and 0.05 (log rate per
year). With one parameter (`scale = "constant"` or, for each horizon,
`"per_horizon"`, with `bias = FALSE`) the log factor is found by
[`stats::optimize()`](https://rdrr.io/r/stats/optimize.html) within its
bounds; otherwise the parameters are found by Nelder-Mead
([`stats::optim()`](https://rdrr.io/r/stats/optim.html), relative
tolerance 1e-6, at most 300 iterations), started at \\s_1 = 1\\ and \\b
= 0\\, with values outside the bounds rejected. The observation noise of
the scored forecasts is simulated once (common random numbers), so the
score is a deterministic function of the parameters. In this calibration
sample, cells with 50 or more expected deaths get Poisson deaths from
the normal approximation (mean and variance equal to the expected
deaths, truncated at zero); the other cells get quantiles of the Poisson
or negative binomial distribution (package-specific, reference unknown).

**Long horizons.** The parameters are estimated on horizons 1 to
\\h\_{max}\\, the distance from the first origin to the last training
year (10 years with the defaults). With `scale = "constant"`, the spread
factor and the drift correction of \\h\_{max}\\ are used beyond it. A
forecast more than twice as long as \\h\_{max}\\ gives a warning (class
`mt_warning_recalibrate`), because nothing in the training data shows
how the spread of the model is wrong at such horizons. With the defaults
this concerns forecasts of more than 20 years
([`mt_forecast()`](mt_forecast.md) forecasts 10 by default); give a
longer `years` for long projections, or accept the constant factor. The
per-horizon methods stop instead (see above).

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
replacement and reports 95% percentile intervals (with
`scale = "per_horizon"`, for each horizon from its forecasts). The inner
models are fitted to shorter series than the final model, which matters
for models whose uncertainty depends on the length of the series (the
drift of a random walk). Compare calibrated and uncalibrated models in
an outer backtest.

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
`mt_score_table(by = c("model", "horizon"))` gives the coverage of each
horizon.

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

Keilman N, Pham DQ, Hetland A (2002). Why population forecasts should be
probabilistic: illustrated by the case of Norway. *Demographic Research*
6(15): 409-454.
[doi:10.4054/DemRes.2002.6.15](https://doi.org/10.4054/DemRes.2002.6.15)

Panagiotelis A, Gamakumara P, Athanasopoulos G, Hyndman RJ (2023).
Probabilistic forecast reconciliation: properties, evaluation and score
optimisation. *European Journal of Operational Research* 306(2):
693-706.
[doi:10.1016/j.ejor.2022.07.040](https://doi.org/10.1016/j.ejor.2022.07.040)

Shang HL, Hyndman RJ (2017). Grouped functional time series forecasting:
an application to age-specific mortality rates. *Journal of
Computational and Graphical Statistics* 26(2): 330-343.
[doi:10.1080/10618600.2016.1237877](https://doi.org/10.1080/10618600.2016.1237877)

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2019)
fit <- mt_fit(recalibrated(lee_carter(), years = 8, n_origins = 3, n_draws = 50), x)
fit$calibration
#>     group     band  scale_h1 drift cells scored crps_before crps_after
#>    <char>   <char>     <num> <num> <int>  <num>       <num>      <num>
#> 1:    all all ages 0.7630612     0  2414   2414  0.03388302 0.03386637
#>    cov80_before cov80_after cov95_before cov95_after forecasts
#>           <num>       <num>        <num>       <num>     <int>
#> 1:    0.7531069   0.7311516    0.8956089   0.8860812         6
fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
# one factor per horizon: horizons 1 to 5 rest on at least two origins
fh <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 5, n_draws = 50,
                          scale = "per_horizon"), x)
fh$calibration[, c("horizon", "origins", "scale")]
#>    horizon origins     scale
#>      <int>   <int>     <num>
#> 1:       1       5 0.1241023
#> 2:       2       5 0.8355160
#> 3:       3       4 0.1000042
#> 4:       4       3 1.3174286
#> 5:       5       2 0.1005425
# point forecast plus resampled errors of the inner forecasts
fe <- mt_fit(recalibrated(lee_carter(), years = 6, n_origins = 5, n_draws = 50,
                          method = "errors"), x)
fe$calibration
#>    horizon origins cells    mean_error tuning coverage
#>      <int>   <num> <int>         <num>  <num>    <num>
#> 1:       1       5   142 -0.0061264434      1        1
#> 2:       2       5   142 -0.0052819251      1        1
#> 3:       3       4   142 -0.0044615082      1        1
#> 4:       4       3   142 -0.0029353580      1        1
#> 5:       5       2   142  0.0003896597      1        1
fce <- mt_forecast(fe, h = 5, n_draws = 100, seed = 1)
```
