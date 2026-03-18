# Random Walk with Drift forecast (per-series)

Applies a random walk with drift model to each (age, group) series of
log(mx) independently.

## Usage

``` r
rwd_forecast(mx_array, years, h, drift_window = 0L, return_fitted = FALSE)
```

## Arguments

- mx_array:

  3D array `[n_a x n_t x n_g]` of mortality rates.

- years:

  Integer vector of calendar years.

- h:

  Integer forecast horizon.

- drift_window:

  Integer; if \> 0, use only the last `drift_window` years to estimate
  the drift (default 0 = all years).

- return_fitted:

  Logical; return fitted values?

## Value

A list with `mx_forecast` `[n_a x h x n_g]` and optionally `mx_fitted`.
