# Convert mort_forecast to a data.table

Convert mort_forecast to a data.table

## Usage

``` r
# S3 method for class 'mort_forecast'
as.data.table(x, keep.rownames = FALSE, ..., keep_sims = FALSE)
```

## Arguments

- x:

  A `mort_forecast` object.

- keep.rownames:

  Ignored (required for generic compatibility).

- ...:

  Ignored.

- keep_sims:

  Logical; include simulation array? Default FALSE.

## Value

A `data.table` with columns: year, age_group, age_start, group columns,
mx_point, and optionally mx_median, mx_lo95, mx_hi95.
