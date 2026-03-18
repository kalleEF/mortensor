# Compute life expectancy

Generic function for computing life expectancy from mortality data or
forecasts. Dispatches to methods for `mort_data` and `mort_forecast`
objects.

## Usage

``` r
life_expectancy(object, ...)

# S3 method for class 'mort_data'
life_expectancy(object, width = 5, open_age_method = "reciprocal", ...)

# S3 method for class 'mort_forecast'
life_expectancy(object, width = 5, open_age_method = "reciprocal", ...)
```

## Arguments

- object:

  A `mort_data` or `mort_forecast` object.

- ...:

  Additional arguments passed to methods.

- width:

  Age-group width. Default 5.

- open_age_method:

  Method for open interval; see
  [`life_table`](https://kalleEF.github.io/mortensor/reference/life_table.md).

## Value

A `data.table` with life expectancy values.

A `data.table` with columns: `year`, `group`, `ex` (life expectancy at
the youngest observed age).

A `data.table` with columns: `year`, `group`, `ex`, and optionally
`ex_lower`, `ex_upper`.

## Methods (by class)

- `life_expectancy(mort_data)`: Compute life expectancy for observed
  mortality data

  Computes period life expectancy at the youngest age for all groups and
  years in a `mort_data` object.

- `life_expectancy(mort_forecast)`: Compute life expectancy for
  mortality forecasts

  Computes period life expectancy from forecast mortality rates,
  including prediction interval bounds if available.
