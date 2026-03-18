# Check smoothing quality

Diagnostic checks comparing raw and smoothed arrays: MAD, max deviation,
trend preservation, and correlation.

## Usage

``` r
check_smoothing(
  mx_raw,
  mx_smooth,
  years,
  covid_years = c(2020, 2021),
  mad_threshold = 0.2,
  floor_threshold = 1e-07
)
```

## Arguments

- mx_raw:

  3D array of raw mortality rates.

- mx_smooth:

  3D array of smoothed mortality rates.

- years:

  Integer vector of calendar years.

- covid_years:

  Years to exclude from comparison.

- mad_threshold:

  MAD threshold for over-smoothing warning.

- floor_threshold:

  Cells below this are excluded.

## Value

A list with `mad`, `max_dev`, `trend_ratio`, `correlation`, `pass`.
