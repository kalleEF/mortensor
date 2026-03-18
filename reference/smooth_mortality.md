# Smooth mortality data

Applies either a 2D tensor-product P-spline (via
[`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html)) or a 1D
Whittaker smoother to a `mort_data` object, smoothing log(mx) over age
and year (or year only for Whittaker).

## Usage

``` r
smooth_mortality(
  x,
  method = c("gam", "whittaker"),
  ages = NULL,
  years = NULL,
  covid_years = integer(0),
  lambda = NULL,
  floor_threshold = 1e-07
)
```

## Arguments

- x:

  A `mort_data` object, or a 3D array `[n_a x n_t x n_g]`.

- method:

  One of `"gam"` (2D P-spline, recommended) or `"whittaker"` (1D
  time-only fallback).

- ages:

  Integer vector of ages (required if `x` is a raw array).

- years:

  Integer vector of years (required if `x` is a raw array).

- covid_years:

  Integer vector of years to exclude (weight = 0).

- lambda:

  Optional smoothing parameter for Whittaker (auto-selected if NULL).

- floor_threshold:

  Mortality rates below this are treated as zero-count.

## Value

If `x` is `mort_data`, a new `mort_data` with `smoothed = TRUE`. If `x`
is a raw array, a list with `arr` (smoothed array) and `method`.
