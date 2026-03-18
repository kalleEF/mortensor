# Convert mortality rates to a death distribution (dx)

Applies a standard life-table calculation to convert age-specific
mortality rates (mx) to a death distribution (dx). Uses the classic nax
approximation for grouped ages.

## Usage

``` r
compute_dx(mx, ages, width = 5L)
```

## Arguments

- mx:

  Numeric vector of age-specific mortality rates.

- ages:

  Integer vector of age-group start values.

- width:

  Integer age-group width (default 5).

## Value

Numeric vector of dx values (sums to ~1).
