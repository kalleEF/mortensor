# Build a 3D mortality array from a data.table

Constructs a `[n_ages x n_years x n_groups]` array from long-format
data.

## Usage

``` r
build_array(dt, value_col, ages, years, groups, floor = 1e-12)
```

## Arguments

- dt:

  A `data.table` with columns for year, sex, quintile/deprivation,
  age_group, and the value column to extract.

- value_col:

  Character name of the column containing mortality rates.

- ages:

  Integer vector of age-group start values.

- years:

  Integer vector of years.

- groups:

  A data.frame with one row per group and columns matching
  stratification variables in `dt` (e.g. `sex`, `quintile`).

- floor:

  Numeric floor for zero/NA cells (default 1e-12). Set to NULL to skip
  floor replacement.

## Value

A 3D array `[n_ages x n_years x n_groups]`.
