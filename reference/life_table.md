# Compute a standard abridged life table

Converts age-specific mortality rates (mx) into a full life table with
columns for qx, px, lx, dx, Lx, Tx, and ex. Handles both grouped and
single-year age data, with special treatment of the open age interval.

## Usage

``` r
life_table(
  mx,
  ages,
  width = 5,
  radix = 1e+05,
  open_age_method = c("reciprocal", "kannisto", "gompertz")
)
```

## Arguments

- mx:

  Numeric vector of age-specific mortality rates.

- ages:

  Numeric vector of age-group start values.

- width:

  Integer or numeric vector of age-group widths. If scalar, applied to
  all groups except the open interval. Default 5.

- radix:

  Numeric radix (initial cohort size). Default 100000.

- open_age_method:

  Character; method for the open age interval's ax value: `"kannisto"`
  (default) fits a Kannisto model to extrapolate, `"gompertz"` uses
  Gompertz extrapolation, or `"reciprocal"` uses `1/mx` for the open
  interval.

## Value

A `data.frame` with columns: `age`, `n` (width), `mx`, `qx`, `px`, `lx`,
`dx`, `Lx`, `Tx`, `ex`.
