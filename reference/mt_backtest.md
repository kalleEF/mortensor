# Rolling-origin backtest

For each forecast origin, the data are cut at the origin, the
preprocessing pipeline and the model are fitted on the training window
only, and predictive draws for the next `h` years are compared with the
observations (rolling-origin evaluation; Tashman 2000). Nothing after
the origin enters the fit, including ungrouping and smoothing (they are
part of `prep`).

## Usage

``` r
mt_backtest(
  spec,
  data,
  origins = NULL,
  h = 10,
  prep = NULL,
  n_draws = 200,
  n_origins = 5,
  seed = 1,
  cores = 1,
  score_estimated = FALSE,
  quantities = NULL,
  contrast = NULL
)
```

## Arguments

- spec:

  A model specification.

- data:

  An [mt_data](mt_data.md) object.

- origins:

  Last training years. Default: the last `n_origins` years whose full
  horizon `h` is observed, so that every horizon is scored equally
  often.

- h:

  Forecast horizon.

- prep:

  Preprocessing pipeline ([`mt_prep()`](mt_prep.md)).

- n_draws:

  Number of draws per origin.

- n_origins:

  Number of origins when `origins` is `NULL`.

- seed:

  Random seed (forecast draws and the randomised PIT); `NULL` for none.

- cores:

  Cores used over origins (forked processes; 1 on Windows).

- score_estimated:

  Score cells whose deaths were estimated (for example suppressed cells
  filled by [`mt_fill_suppressed()`](mt_fill_suppressed.md)); by default
  only observed cells are scored.

- quantities:

  Derived quantities scored in addition to the rates:
  `"life_expectancy"` (at the first age of the grid, in years) and/or
  `"asr"` (age-standardised rate, European Standard Population 2013, on
  the log scale; direct standardisation on the age groups of the
  standard with the exposures of the test year, see
  [`asr()`](derived.md)), computed from the predictive draws of deaths
  aggregated to the observed age groups of each stratum and test year,
  like the observed values. Years whose observed cells do not cover the
  full age range without gaps (missing, suppressed or estimated cells)
  are skipped. A zero ASR (no deaths at any age) is replaced by the
  smallest value that half a death gives. They are added as rows with
  `target` `"life_expectancy"` or `"asr"`; summarise them with
  `mt_score(..., target = )`.

- contrast:

  Optional `c(var = , ref = )`: a stratum variable and its reference
  level. With `quantities`, the differences between every other level
  and the reference level (with the same values of the other stratum
  variables) are scored too, path by path: targets
  `"life_expectancy_gap"` (difference in years) and `"asr_gap"`
  (difference of log ASR, the log rate ratio). This checks forecasts of
  inequality, whose uncertainty depends on how the paths of the strata
  move together; the rows carry the stratum of the compared level.

## Value

An object of class `mt_backtest`: a `data.table` with one row per
observed cell, origin and horizon: `model`, `prep` (the preprocessing
steps, `"none"` without), `origin`, `horizon`, stratum columns,
`age_lo`, `age_hi`, `period`, `observed` (log rate), the MASE scale
`scale` (mean absolute year-to-year change in the training data),
`target` (`"rate"`, or the derived quantity), the forecast `median` and
`mean`, and the scores `crps`, `wis`, `ae` (absolute error of the
median), `cov50`, `cov80`, `cov95` and `pit` (randomised for ties), and
for the opt-in summaries of [`mt_score()`](mt_score.md) (`measures`) the
standard deviation of the draws `sd`, `cov90` and the half-widths
`hw80`, `hw90` and `hw95` of the central intervals. Rows of derived
quantities have no `age_lo`, `age_hi` and `scale`. Summarise with
[`mt_score()`](mt_score.md).

## Details

Evaluation happens at the resolution of the observed deaths: grid draws
of deaths (Poisson given the exposures of the test years, which are
treated as known) are summed to each observed cell, so models on single
years of age can be scored against data in 5-year groups. Exposures of
the test years are ungrouped if needed (no deaths of the test years are
used for that).

Data whose deaths were ungrouped beforehand with
[`mt_ungroup()`](mt_ungroup.md), and data smoothed with
[`mt_smooth()`](mt_smooth.md), are refused: the forecasts would be
scored against estimates instead of observations, and an ungrouping or
smoothing over time has seen the test years. Ungroup or smooth inside
each fit instead (`prep = mt_prep(step_ungroup())`,
[`step_smooth()`](mt_prep.md)). With [`step_regroup()`](mt_prep.md), the
test data are summed to the age groups of the fit.

The log rate of an observed cell is `log(max(D, 0.5) / E)`: a cell
without deaths counts as half a death, so that its log rate is finite (a
convention of this package; reference unknown). The predictive draws of
deaths are first summed from the grid to the observed cell and then
floored in the same way, so that draws and observations are treated
alike; a floor per grid cell would put the draws of a grouped cell with
few deaths above the observation.

Strata without forecast draws at an origin (the model failed in them or
gave missing rates) are not scored, with a warning of class
`mt_warning_backtest`. Compare models on the cells all of them scored
([`mt_score_table()`](mt_compare.md) does so for
[`mt_compare()`](mt_compare.md) results).

Scores (see [`mt_score()`](mt_score.md) for their summaries): the CRPS
is the empirical-CDF estimator from the draws (Grimit et al. 2006), the
mean absolute difference between the draws and the observation minus
half the mean absolute difference between pairs of draws. It is the CRPS
of the empirical distribution of the draws, which exceeds that of the
predictive distribution by about `E|X - X'| / (2 n_draws)` on average,
so compare models backtested with the same `n_draws`. The weighted
interval score uses the 23 quantiles of the COVID-19 forecast hubs
(Bracher et al. 2021); the PIT is randomised where draws equal the
observation (Czado, Gneiting and Held 2009), as with zero deaths.

## References

Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
forecasts in an interval format. *PLoS Computational Biology* 17(2):
e1008618.
[doi:10.1371/journal.pcbi.1008618](https://doi.org/10.1371/journal.pcbi.1008618)

Czado C, Gneiting T, Held L (2009). Predictive model assessment for
count data. *Biometrics* 65(4): 1254-1261.

Gneiting T, Raftery AE (2007). Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical
Association* 102(477): 359-378.
[doi:10.1198/016214506000001437](https://doi.org/10.1198/016214506000001437)

Grimit EP, Gneiting T, Berrocal VJ, Johnson NA (2006). The continuous
ranked probability score for circular variables and its application to
mesoscale forecast ensemble verification. *Quarterly Journal of the
Royal Meteorological Society* 132: 2925-2942.

Tashman LJ (2000). Out-of-sample tests of forecasting accuracy: an
analysis and review. *International Journal of Forecasting* 16(4):
437-450.

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2012)
bt <- mt_backtest(rwd(), x, origins = c(2005, 2007), h = 5, n_draws = 100)
mt_score(bt)
#>                     model       crps        wis        mae      rmse     mase
#>                    <char>      <num>      <num>      <num>     <num>    <num>
#> 1: Random walk with drift 0.05507087 0.04923878 0.07405228 0.1188856 1.265727
#>        cov50    cov80     cov95  pit_mean     n     unit
#>        <num>    <num>     <num>     <num> <int>   <char>
#> 1: 0.6415493 0.878169 0.9591549 0.4867606  1420 log rate
```
