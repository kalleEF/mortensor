# Export mortality objects to files

Writes a mort_forecast or mort_calibration to CSV/RDS files.

## Usage

``` r
export(x, ...)

# S3 method for class 'mort_forecast'
export(x, dir = ".", prefix = NULL, format = c("csv", "rds"), ...)

# S3 method for class 'mort_calibration'
export(x, dir = ".", prefix = NULL, format = c("csv", "rds"), ...)
```

## Arguments

- x:

  An object to export.

- ...:

  Additional arguments passed to methods.

- dir:

  Output directory (default `"."`).

- prefix:

  File name prefix (defaults to cause label or `"forecast"`).

- format:

  One of `"csv"` or `"rds"`.

## Methods (by class)

- `export(mort_forecast)`: Export a mort_forecast

- `export(mort_calibration)`: Export a mort_calibration
