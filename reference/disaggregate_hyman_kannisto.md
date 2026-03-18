# Disaggregate grouped mortality to single-year ages (Hyman + Kannisto)

Primary disaggregation method using a Hyman monotone spline for closed
age groups and Kannisto logistic extrapolation for open-ended 90+ group.
Iterative smooth correction ensures grouped means are preserved.

## Usage

``` r
disaggregate_hyman_kannisto(
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
