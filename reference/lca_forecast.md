# Lee-Carter forecast via demography package

Fits Lee-Carter models per group using
[`demography::lca()`](https://pkg.robjhyndman.com/demography/reference/lca.html)
and forecasts.

## Usage

``` r
lca_forecast(
  mx_array,
  pop_array,
  years,
  ages,
  groups_dt,
  h,
  return_fitted = FALSE,
  adjust = "dt",
  chooseperiod = FALSE
)
```

## Arguments

- mx_array:

  3D array `[n_a x n_t x n_g]` of mortality rates.

- pop_array:

  3D array of population exposure (same dims).

- years:

  Integer vector of calendar years.

- ages:

  Integer vector of age-group starts.

- groups_dt:

  A data.frame with group definitions (must have `sex`, `quintile`
  columns, or similar identifiers).

- h:

  Integer forecast horizon.

- return_fitted:

  Logical; return fitted values?

- adjust:

  Adjustment method: `"dt"`, `"dxt"`, `"e0"`, or `"none"`.

- chooseperiod:

  Logical; let lca pick the best fitting period?

## Value

A list with `mx_forecast` `[n_a x h x n_g]` and optionally `mx_fitted`.
