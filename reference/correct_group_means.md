# Smooth iterative group-mean correction

Adjusts single-year rates so that the mean within each 5-year interval
matches the original grouped rate. Uses natural spline interpolation of
log-correction factors to avoid discontinuous jumps at group boundaries.

## Usage

``` r
correct_group_means(
  mx_single,
  mx_grouped,
  target_ages,
  age_starts,
  width,
  max_age
)
```

## Arguments

- mx_single:

  Numeric vector of single-year mx.

- mx_grouped:

  Numeric vector of grouped mx.

- target_ages:

  Integer vector of single-year ages.

- age_starts:

  Integer vector of age-group starts.

- width:

  Age-group width (default 5).

- max_age:

  Maximum target age.

## Value

Corrected single-year mx vector.
