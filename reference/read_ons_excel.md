# Read mortality data from ONS-format Excel files

Reads UK Office for National Statistics (ONS) death registration and
population data from Excel files and returns a
[mort_data](https://kalleEF.github.io/mortensor/reference/mort_data.md)
object. Handles IMD decile-to-quintile mapping, cause column
aggregation, population forward-filling, and array construction.

## Usage

``` r
read_ons_excel(deaths_file, pop_files, config)

# Deprecated alias
```

## Arguments

- deaths_file:

  Path to the Excel file with death registration data (Sheet 1).

- pop_files:

  Character vector of paths to population Excel files (Sheet 2, one file
  per year). Multiple files are stacked.

- config:

  A
  [mort_config](https://kalleEF.github.io/mortensor/reference/mort_config.md)
  object specifying cause columns, age groups, year range, etc.

## Value

A
[mort_data](https://kalleEF.github.io/mortensor/reference/mort_data.md)
object.

## Note

This function is specific to the UK ONS data format. For generic
mortality data, use
[`read_mortality()`](https://kalleEF.github.io/mortensor/reference/read_mortality.md)
instead.
