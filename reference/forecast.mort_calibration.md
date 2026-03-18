# Produce mortality forecasts from a calibrated model

Produce mortality forecasts from a calibrated model

## Usage

``` r
# S3 method for class 'mort_calibration'
forecast(object, h = NULL, ...)
```

## Arguments

- object:

  A `mort_calibration` object from
  [`calibrate()`](https://kalleEF.github.io/mortensor/reference/calibrate.md).

- h:

  Forecast horizon in years (default: from config or 10).

- ...:

  Ignored.

## Value

A `mort_forecast` object.
