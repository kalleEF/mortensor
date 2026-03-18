# Convert a mort_data or mort_forecast object to a long-format data.table

Package-native alternative to `as.data.table()` that works after
[`library(mortensor)`](https://kalleEF.github.io/mortensor) without
needing [`library(data.table)`](https://r-datatable.com). For users who
have data.table loaded, `as.data.table()` dispatches to the same
methods.

## Usage

``` r
to_data_table(x, ...)

# S3 method for class 'mort_data'
as.data.table(x, ...)
```

## Arguments

- x:

  A `mort_data` object.

- ...:

  Ignored.

## Value

A long-format `data.table`.

A `data.table` with columns: year, age_group, age_start, group columns,
mx, and optionally pop.

## Functions

- `as.data.table(mort_data)`: Convert a mort_data object
