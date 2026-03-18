# Convert a death distribution back to mortality rates

Inverse of
[`compute_dx()`](https://kalleEF.github.io/mortensor/reference/compute_dx.md).
Uses the life-table relationship and Gompertz extrapolation for the open
age interval.

## Usage

``` r
dx_to_mx(dx, ages, width = 5L)
```

## Arguments

- dx:

  Numeric vector of death distribution values.

- ages:

  Integer vector of age-group start values.

- width:

  Integer age-group width (default 5).

## Value

Numeric vector of age-specific mortality rates (mx).
