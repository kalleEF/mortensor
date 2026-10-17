# Preprocessing pipelines

`mt_prep()` bundles preprocessing steps that are applied to the data
before a model is fitted. Because the steps are part of the fit, a
backtest ([`mt_backtest()`](mt_backtest.md)) re-runs them on every
training window, so that ungrouping or smoothing never sees data after
the forecast origin.

## Usage

``` r
mt_prep(...)

step_ungroup(...)

step_smooth(...)

step_regroup(age = NULL, period = NULL)

step_grid(age = NULL, period = NULL)

step_window(years)

step_shocks(periods, action = c("interpolate", "exclude"))

mt_prep_apply(prep, data)
```

## Arguments

- ...:

  Steps, applied in order.

- age, period:

  Axes for `step_regroup()` and `step_grid()`.

- years:

  Number of most recent training years kept by `step_window()`.

- periods:

  Calendar years affected by a shock.

- action:

  `"interpolate"` or `"exclude"`.

- prep:

  An `mt_prep` object or `NULL` (no preprocessing).

- data:

  An [mt_data](mt_data.md) object.

## Value

An object of class `mt_prep`; `mt_prep_apply()` returns the processed
[mt_data](mt_data.md) object.

## Details

Steps:

- `step_ungroup(...)`: [`mt_ungroup()`](mt_ungroup.md);

- `step_smooth(...)`: [`mt_smooth()`](mt_smooth.md). The model is fitted
  to the smoothed rates; residuals and observation noise refer to the
  raw deaths. Smoothing over time (`along = "period"`, or `dims = 2`)
  narrows the forecast intervals of models with a period index, so check
  coverage in a backtest and consider
  [`recalibrated()`](recalibrated.md);

- `step_regroup(age, period)`: [`mt_regroup()`](mt_regroup.md);

- `step_grid(age, period)`: [`mt_set_grid()`](mt_set_grid.md);

- `step_window(years)`: keeps the last `years` years of the training
  data, a fitting period chosen relative to the forecast origin (for
  example after a change of trend; compare windows with
  [`tuned()`](mt_tune.md) or [`mt_backtest()`](mt_backtest.md));

- `step_shocks(periods, action)`: periods with exceptional mortality
  (for example 2020-2022) are either replaced by log-linear
  interpolation of the rates of the neighbouring periods
  (`"interpolate"`, deaths become exposure times the interpolated rate)
  or removed (`"exclude"`, only at the start or end of the data). With
  `"interpolate"`, shock periods at the start or end of the data, which
  have no neighbour on one side, are removed, as in a backtest whose
  training window ends inside the shock.

## Examples

``` r
prep <- mt_prep(step_ungroup(), step_shocks(2005, action = "interpolate"))
prep
#> <mt_prep> 2 step(s)
#>   step_ungroup()
#>   step_shocks(periods=2005, action=interpolate)
x <- mt_prep_apply(prep, mt_example(years = 2001:2008))
```
