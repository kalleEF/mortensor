# Tucker4 forecast on ILR cause-shares

Decomposes the 4D ILR tensor, forecasts the B (time) mode, reconstructs,
and inverts ILR to produce share forecasts.

## Usage

``` r
tucker4_coherent_forecast(
  X_ilr,
  years,
  h,
  R_A = 3,
  R_B = 3,
  R_C = 2,
  R_D = 3,
  ts_method = "ets",
  fast = TRUE,
  return_full = FALSE
)
```

## Arguments

- X_ilr:

  4D array \[n_a x n_t x n_g x K\] of ILR coordinates.

- years:

  Integer vector of observed years.

- h:

  Forecast horizon.

- R_A, R_B, R_C, R_D:

  Tucker4 ranks.

- ts_method:

  `"ets"`, `"arima"`, or `"rwd"`.

- fast:

  Use fast convergence settings?

- return_full:

  Return decomposition components for bootstrap?

## Value

List with `shares_forecast` \[n_a x h x n_g x K+1\] and `fit_pct`.
