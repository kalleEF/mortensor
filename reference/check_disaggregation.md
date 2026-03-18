# Check disaggregation quality

Verifies group-mean preservation, monotonicity from age 50+, and
smoothness (second differences).

## Usage

``` r
check_disaggregation(mx_grouped, mx_single, age_starts, target_ages)
```

## Arguments

- mx_grouped:

  Grouped rates.

- mx_single:

  Single-year rates.

- age_starts:

  Group start ages.

- target_ages:

  Single-year ages.

## Value

Named list with `max_rel_error`, `monotone`, `max_jump`.
