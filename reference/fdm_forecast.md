# Hyndman-Ullah Functional Demographic Model (FDM) forecast

Hyndman-Ullah Functional Demographic Model (FDM) forecast

## Usage

``` r
fdm_forecast(
  mx_array,
  pop_array,
  years,
  ages,
  groups_dt,
  h,
  return_fitted = FALSE,
  order = 6,
  method = "arima"
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

- order:

  Number of principal components (default 6).

- method:

  Time-series extrapolation: `"arima"` or `"ets"`.

## Value

A list with `mx_forecast` and optionally `mx_fitted`.
