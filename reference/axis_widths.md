# Interval widths and midpoints of an axis

Interval widths and midpoints of an axis

## Usage

``` r
axis_widths(axis, omega = 110)

axis_midpoints(axis, omega = 110)

is_open(axis)
```

## Arguments

- axis:

  An [mt_axis](mt_axis.md).

- omega:

  Assumed upper bound used in place of `Inf` for an open last interval
  (default 110 for ages).

## Value

Numeric vector.
