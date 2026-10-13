# Map intervals onto an axis

For each interval `[lo, hi)`, returns the index of the axis cell that
contains it (`mode = "within"`) or that equals it (`mode = "exact"`),
and `NA` otherwise.

## Usage

``` r
axis_locate(axis, lo, hi, mode = c("within", "exact"))
```

## Arguments

- axis:

  An [mt_axis](mt_axis.md).

- lo, hi:

  Interval bounds.

- mode:

  `"within"` or `"exact"`.

## Value

Integer vector.
