# Rolling-origin backtest

For each forecast origin, the data are cut at the origin, the
preprocessing pipeline and the model are fitted on the training window
only, and predictive draws for the next `h` years are compared with the
observations. Nothing after the origin enters the fit, including
ungrouping and smoothing (they are part of `prep`).

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
  quantities = NULL
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

  Random seed (forecast draws and the randomised PIT).

- cores:

  Cores used over origins (forked processes; 1 on Windows).

- score_estimated:

  Score cells whose deaths were estimated (for example suppressed cells
  filled by [`mt_fill_suppressed()`](mt_fill_suppressed.md)); by default
  only observed cells are scored.

- quantities:

  Derived quantities scored in addition to the rates:
  `"life_expectancy"` (at the first observed age, in years) and/or
  `"asr"` (age-standardised rate, European Standard Population 2013, on
  the log scale), computed from the predictive draws aggregated to the
  observed age groups of each stratum and test year (years with missing
  or suppressed cells are skipped). They are added as rows with `target`
  `"life_expectancy"` or `"asr"`; summarise them with
  `mt_score(..., target = )`.

## Value

An object of class `mt_backtest`: a `data.table` with one row per
observed cell, origin and horizon: `model`, `prep` (the preprocessing
steps, `"none"` without), `origin`, `horizon`, stratum columns,
`age_lo`, `age_hi`, `period`, `observed` (log rate), `median`, `mean`,
quantiles, and the scores `crps`, `wis`, `ae` (absolute error of the
median), `cov50`, `cov80`, `cov95`, `pit` and the MASE scale `scale`
(mean absolute year-to-year change in the training data), and `target`
(`"rate"`, or the derived quantity). Summarise with
[`mt_score()`](mt_score.md).

## Details

Evaluation happens at the resolution of the observed deaths: grid draws
of deaths (Poisson given the exposures of the test years, which are
treated as known) are summed to each observed cell, so models on single
years of age can be scored against data in 5-year groups. Exposures of
the test years are ungrouped if needed (no deaths of the test years are
used for that).

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2012)
bt <- mt_backtest(rwd(), x, origins = c(2005, 2007), h = 5, n_draws = 100)
mt_score(bt)
#>                     model       crps        wis        mae     mase     cov50
#>                    <char>      <num>      <num>      <num>    <num>     <num>
#> 1: Random walk with drift 0.05507087 0.04923878 0.07405228 1.265727 0.6415493
#>       cov80     cov95  pit_mean     n     unit
#>       <num>     <num>     <num> <int>   <char>
#> 1: 0.878169 0.9591549 0.4867606  1420 log rate
```
