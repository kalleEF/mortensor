# Damped ETS forecast (per-series)

Fits an ETS (exponential smoothing) model to each (age, group) series of
log(mx) with optional damping.

## Usage

``` r
ets_forecast(
  mx_array,
  years,
  h,
  ets_model = "AAN",
  damped = TRUE,
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

- ets_model:

  Character ETS model specification (default `"AAN"`).

- damped:

  Logical; apply trend damping?

- return_fitted:

  Logical; return fitted values?

## Value

A list with `mx_forecast` and optionally `mx_fitted`.
