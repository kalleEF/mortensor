# Extract fitted values from a mort_forecast

Extract fitted values from a mort_forecast

## Usage

``` r
# S3 method for class 'mort_forecast'
fitted(object, ...)
```

## Arguments

- object:

  A `mort_forecast` object.

- ...:

  Ignored.

## Value

A 3D array `[n_ages x n_train_years x n_groups]` of fitted mortality
rates, or `NULL` if not available.
