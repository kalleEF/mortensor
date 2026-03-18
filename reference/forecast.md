# Forecast mortality from a calibrated model

Generic function for producing mortality forecasts. Dispatches on the
class of `object`. After
[`library(mortensor)`](https://kalleEF.github.io/mortensor) this generic
is always available without needing to load the forecast package
separately.

## Usage

``` r
forecast(object, h = NULL, ...)
```

## Arguments

- object:

  A `mort_calibration` object (from
  [`calibrate()`](https://kalleEF.github.io/mortensor/reference/calibrate.md)).

- h:

  Integer forecast horizon (years). Overrides the value stored in the
  calibration object.

- ...:

  Additional arguments passed to the method.

## Value

A `mort_forecast` object.
