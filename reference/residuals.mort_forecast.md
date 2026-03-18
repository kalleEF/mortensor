# Extract residuals from a mort_forecast

Extract residuals from a mort_forecast

## Usage

``` r
# S3 method for class 'mort_forecast'
residuals(object, type = c("response", "deviance"), ...)
```

## Arguments

- object:

  A `mort_forecast` object.

- type:

  One of `"response"` (raw scale: observed - fitted) or `"deviance"`
  (log scale: log(observed) - log(fitted)).

- ...:

  Ignored.

## Value

A 3D array of residuals, or `NULL` if fitted values are not available.
