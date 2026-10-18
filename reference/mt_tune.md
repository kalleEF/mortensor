# Tune model and preprocessing settings by forecast accuracy

`mt_tune()` evaluates every combination of settings in `grid` by how
well it forecasts the last years of `data` (an inner backtest) and
returns the scores and the best specification. `tuned()` wraps the same
procedure into a model specification: each time it is fitted (for
example at every origin of an outer [`mt_backtest()`](mt_backtest.md)),
the tuning uses only the training data it receives, so the selection
never sees the years being forecast. The chosen model is then fitted to
all training data.

## Usage

``` r
mt_tune(
  spec_fun,
  grid,
  data,
  prep = NULL,
  inner_origins = 3,
  inner_h = 5,
  metric = "crps",
  n_draws = 100,
  strategy = c("full", "two_stage"),
  target = c("rate", "life_expectancy", "asr"),
  select = c("min", "one_se"),
  by = NULL
)

tuned(
  spec_fun,
  grid,
  inner_origins = 3,
  inner_h = 5,
  metric = "crps",
  n_draws = 100,
  label = "Tuned model",
  strategy = c("full", "two_stage"),
  target = c("rate", "life_expectancy", "asr"),
  select = c("min", "one_se"),
  by = NULL
)
```

## Arguments

- spec_fun:

  A function that returns a model specification (or
  `list(spec = , prep = )`) from named arguments, for example
  `function(order) hyndman_ullah(order = order)`.

- grid:

  A data frame of settings (columns match the arguments of `spec_fun`),
  ordered from the simplest to the most complex setting for
  `select = "one_se"`.

- data:

  An [mt_data](mt_data.md) object.

- prep:

  Preprocessing pipeline used by every candidate (when `spec_fun`
  returns specifications only).

- inner_origins:

  Number of inner origins (`1` for a single hold-out).

- inner_h:

  Horizon of the inner backtest (the number of years held out from each
  origin).

- metric:

  `"crps"`, `"wis"`, `"mae"`, `"rmse"` or `"mase"`.

- n_draws:

  Draws per inner forecast.

- strategy:

  `"full"` evaluates every row of `grid`; `"two_stage"` first evaluates
  a coarse grid (every second value of each column, including the first
  and the last) and then the rows of `grid` within one step of the best
  coarse combination in every column. For large grids, such as Tucker
  ranks, this needs far fewer fits. The two-stage search is a package
  heuristic: it can miss the best row when the score has several minima
  along a column.

- target:

  `"rate"`, `"life_expectancy"` or `"asr"`.

- select:

  `"min"` or `"one_se"`.

- by:

  Stratum variables: tune each group of strata separately.

- label:

  Display name.

## Value

`mt_tune()`: a list with `scores` (one row per setting and group: the
settings, the score, its standard error `se` and whether it was
`selected`), `by_horizon` (scores by forecast horizon), `best` (the
chosen specification; with `by`, a list per group), `best_prep` (the
preprocessing of the chosen setting), `chosen` and `origins`. `tuned()`:
a model specification; its fit carries `tuning` (the scores) and
`tuning_choice`.

## Details

**Design.** Each candidate is refitted at `inner_origins` consecutive
origins and forecasts `inner_h` years from each; the last origin is
`inner_h` years before the end of the data. `inner_origins = 1` is a
single hold-out: fit to the years up to `T - inner_h` and forecast the
last `inner_h` years. Several origins average over more periods and make
the choice depend less on the events of one period (Tashman 2000).
`tuned()` asks for at least `inner_origins + inner_h + 5` training
years, so that the first inner fit has six years; this minimum is a
package heuristic. Candidates whose models need more years fail at the
first inner origins and are then scored on fewer origins than the
others, or not at all.

**What can be tuned.** `spec_fun` receives one row of `grid` as named
arguments and returns a model specification, or `list(spec = , prep = )`
to tune preprocessing together with the model: for example the smoothing
parameter of [`step_smooth()`](mt_prep.md), the fitting window of
[`step_window()`](mt_prep.md) and the number of components of
[`hyndman_ullah()`](hyndman_ullah.md). Preprocessing is given either by
`spec_fun` or as `prep`.

**Criterion.** `metric` is computed for `target`: the log death rates of
the observed cells (`"rate"`), life expectancy at birth
(`"life_expectancy"`, years) or the log age-standardised rate (`"asr"`),
per test year and stratum. `"crps"` and `"wis"` score the whole forecast
distribution; `"mae"` and `"rmse"` the error of the median (point
accuracy); `"mase"` scales the absolute error by the mean year-to-year
change of the training data (rates only).

**Selection.** `select = "min"` takes the lowest score.
`select = "one_se"` takes the first row of `grid` whose score is within
one standard error of the lowest (the one-standard-error rule of Hastie,
Tibshirani and Friedman 2009, Section 7.10), so order `grid` from the
simplest to the most complex setting. The standard error comes from the
mean scores per origin and stratum (per test year for a hold-out of a
single stratum); these blocks are correlated, so it is approximate and
tends to be too small. With few test years, differences between settings
are often within this error: the selection is then uncertain, and
`scores` shows it.

**Groups.** With `by`, the strata are split by these stratum variables
and each group is tuned and fitted on its own (for example women and men
with their own number of components). A joint model such as
[`li_lee()`](li_lee.md) then pools only within each group.

## References

Tashman LJ (2000). Out-of-sample tests of forecasting accuracy: an
analysis and review. *International Journal of Forecasting* 16(4):
437-450.

Hastie T, Tibshirani R, Friedman J (2009). *The Elements of Statistical
Learning*, 2nd edition. Springer, New York.

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2012)
tn <- mt_tune(function(drift_window) rwd(drift_window = drift_window),
              data.frame(drift_window = c(5, 10, 20)), x, inner_origins = 2, inner_h = 3)
tn$scores
#>    drift_window       crps          se selected
#>           <num>      <num>       <num>   <lgcl>
#> 1:            5 0.06098906 0.003459354    FALSE
#> 2:           10 0.05295919 0.003492327    FALSE
#> 3:           20 0.04738829 0.003376065     TRUE
# a single hold-out of the last 5 years, tuning the smoothing too
smoothed_lc <- function(lambda) {
  list(spec = lee_carter(), prep = mt_prep(step_smooth(lambda = lambda)))
}
tn2 <- mt_tune(smoothed_lc, data.frame(lambda = c(1, 100)), x, inner_origins = 1, inner_h = 5,
               metric = "rmse", n_draws = 50)
tn2$scores
#>    lambda       rmse          se selected
#>     <num>      <num>       <num>   <lgcl>
#> 1:      1 0.06993299 0.008076431    FALSE
#> 2:    100 0.06808111 0.008657170     TRUE
```
