# Tucker3-CoDa forecast

Forecasts mortality by applying Tucker3 tensor decomposition to the
Centered Log-Ratio (CLR) transform of the death distribution. Pipeline:
mx -\> dx -\> CLR -\> Tucker3 -\> forecast B-mode -\> inv-CLR -\> mx.

## Usage

``` r
tucker3_coda_forecast(
  mx_array,
  years,
  ages,
  h,
  R_A = 3,
  R_B = 3,
  R_C = 2,
  ts_method = "ets",
  fast = TRUE,
  return_full = FALSE,
  width = 5
)
```

## Arguments

- mx_array:

  3D array `[n_a x n_t x n_g]` of mortality rates.

- years:

  Integer vector of calendar years.

- ages:

  Integer vector of age-group start values.

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

- width:

  Integer age-group width (default 5).

## Value

A list with `mx_forecast` `[n_a x h x n_g]` and `fit_pct`. If
`return_full = TRUE`, also includes decomposition factors, residuals,
etc.
