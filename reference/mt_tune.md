# Tune model parameters on training data

`mt_tune()` evaluates every combination of parameters in `grid` with an
inner rolling-origin backtest on `data` and returns the scores and the
best specification. `tuned()` wraps the same procedure into a model
specification: each time it is fitted (for example at every origin of an
outer [`mt_backtest()`](mt_backtest.md)), the tuning uses only the
training data it receives, so the selection never sees the years being
forecast.

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
  strategy = c("full", "two_stage")
)

tuned(
  spec_fun,
  grid,
  inner_origins = 3,
  inner_h = 5,
  metric = "crps",
  n_draws = 100,
  label = "Tuned model",
  strategy = c("full", "two_stage")
)
```

## Arguments

- spec_fun:

  A function that returns a model specification from named arguments
  (for example `function(drift_window) rwd(drift_window)`).

- grid:

  A data frame of parameter combinations (columns match the arguments of
  `spec_fun`).

- data:

  An [mt_data](mt_data.md) object.

- prep:

  Preprocessing pipeline.

- inner_origins:

  Number of inner origins (the last years that leave `inner_h` years of
  test data).

- inner_h:

  Horizon of the inner backtest.

- metric:

  Column of [`mt_score()`](mt_score.md) to minimise (`"crps"`, `"wis"`,
  `"mae"`, `"mase"`).

- n_draws:

  Draws per inner forecast.

- strategy:

  `"full"` evaluates every row of `grid`; `"two_stage"` first evaluates
  a coarse grid (every second value of each column, including the first
  and the last) and then the rows of `grid` within one step of the best
  coarse combination in every column. For large grids, such as Tucker
  ranks, this needs far fewer fits.

- label:

  Display name.

## Value

`mt_tune()`: a list with `scores` (one row per combination) and `best`
(the chosen specification). `tuned()`: a model specification.

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2012)
tn <- mt_tune(function(drift_window) rwd(drift_window = drift_window),
              data.frame(drift_window = c(5, 10, 20)), x, inner_origins = 2, inner_h = 3)
tn$scores
#>    drift_window       crps
#>           <num>      <num>
#> 1:            5 0.06098906
#> 2:           10 0.05295919
#> 3:           20 0.04738829
```
