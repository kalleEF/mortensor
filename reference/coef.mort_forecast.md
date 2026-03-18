# Extract model coefficients from a mort_forecast

For Tucker3-based models, returns the decomposition factors A (age), B
(time), C (group), and core tensor G. For other models, returns any
extras stored during fitting.

## Usage

``` r
# S3 method for class 'mort_forecast'
coef(object, ...)
```

## Arguments

- object:

  A `mort_forecast` object.

- ...:

  Ignored.

## Value

A named list of model components, or `NULL`.
