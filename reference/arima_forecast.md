# Per-series ARIMA forecast

Fits an auto.arima model to each (age, group) series of log(mx).

## Usage

``` r
arima_forecast(
  mx_array,
  years,
  h,
  max.p = 2L,
  max.q = 1L,
  max.d = 1L,
  return_fitted = FALSE
)
```

## Arguments

- mx_array:

  3D array of mortality rates.

- years:

  Integer vector of calendar years.

- h:

  Integer forecast horizon.

- max.p, max.q, max.d:

  Maximum ARIMA orders.

- return_fitted:

  Logical; return fitted values?

## Value

A list with `mx_forecast` and optionally `mx_fitted`.
