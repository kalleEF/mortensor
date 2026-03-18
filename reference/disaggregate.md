# Disaggregate a mort_forecast to single-year ages

Converts grouped-age forecasts to single-year age forecasts using spline
interpolation with parametric tail closure.

## Usage

``` r
disaggregate(
  fc,
  target_ages = 30L:100L,
  method = c("natural_gompertz", "hyman_kannisto", "gam"),
  width = 5L,
  ...
)
```

## Arguments

- fc:

  A `mort_forecast` object with grouped ages.

- target_ages:

  Integer vector of single-year target ages (default `30:100`).

- method:

  One of `"natural_gompertz"`, `"hyman_kannisto"`, or `"gam"`.

- width:

  Width of closed age groups (default 5).

- ...:

  Additional arguments for the disaggregation function.

## Value

A new `mort_forecast` object with single-year ages.
