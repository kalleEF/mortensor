# Disaggregate prediction interval bounds using ratio interpolation

Computes ratio of lower/upper PI to point forecast at grouped level,
interpolates to single-year ages via natural splines, and applies to
single-year point forecast.

## Usage

``` r
disaggregate_pi_ratios(
  mx_1yr_point,
  mx_grp_point,
  mx_grp_lo,
  mx_grp_hi,
  mx_grp_med,
  target_ages,
  age_starts,
  width
)
```

## Arguments

- mx_1yr_point:

  Numeric vector of single-year point forecast.

- mx_grp_point:

  Numeric vector of grouped point forecast.

- mx_grp_lo:

  Numeric vector of grouped lower 95 percent PI.

- mx_grp_hi:

  Numeric vector of grouped upper 95 percent PI.

- mx_grp_med:

  Numeric vector of grouped median.

- target_ages:

  Integer vector of single-year ages.

- age_starts:

  Integer vector of age-group starts.

- width:

  Age-group width.

## Value

Named list with `mx_1yr`, `mx_1yr_lo95`, `mx_1yr_hi95`, `mx_1yr_median`.
