# Disaggregate grouped mortality to single-year ages (Natural + Gompertz)

Alternative disaggregation using natural cubic spline interpolation with
Gompertz virtual extension for old ages.

## Usage

``` r
disaggregate_natural_gompertz(
  mx_grouped,
  age_starts,
  target_ages = 30:100,
  width = 5
)
```

## Arguments

- mx_grouped:

  Numeric vector of grouped mortality rates.

- age_starts:

  Integer vector of group lower bounds (e.g., 30, 35, ..., 90).

- target_ages:

  Integer vector of single-year target ages (default 30:100).

- width:

  Width of closed age groups (default 5).

## Value

Numeric vector of single-year mortality rates.
