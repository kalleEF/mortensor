# Tucker3-Direct forecast

Forecasts mortality by applying Tucker3 decomposition directly to
log(mx). Pipeline: log(mx) -\> Tucker3 -\> forecast B-mode -\> exp -\>
mx.

## Usage

``` r
tucker3_direct_forecast(
  mx_array,
  years,
  h,
  R_A = 3,
  R_B = 3,
  R_C = 2,
  ts_method = "ets",
  fast = TRUE,
  return_full = FALSE
)
```

## Arguments

- mx_array:

  3D array `[n_a x n_t x n_g]` of mortality rates.

- years:

  Integer vector of calendar years.

- h:

  Integer forecast horizon.

- R_A, R_B, R_C:

  Integer Tucker3 ranks for age, time, group modes.

- ts_method:

  Time-series method for B-mode: `"ets"`, `"arima"`, or `"rwd"` (random
  walk with drift).

- fast:

  Logical; use faster but approximate fitting?

- return_full:

  Logical; if TRUE, return decomposition details.

## Value

A list with `mx_forecast` and `fit_pct`.
