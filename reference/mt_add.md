# Add observations to an `mt_data` object

Add observations to an `mt_data` object

## Usage

``` r
mt_add(
  data,
  measure,
  x,
  by = NULL,
  age = "age",
  period = "year",
  value = NULL,
  age_open = NULL,
  missing = c(".", "..", "...", ":", "x", "X", "/", "NA", ""),
  zero = "-",
  decimal = "."
)

mt_add_deaths(data, x, ...)

mt_add_exposure(data, x, ...)

mt_add_population(data, x, ...)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- measure:

  One of `"deaths"`, `"exposure"`, `"population"`.

- x:

  A data frame in long format.

- by:

  Character vector of stratum columns (for example `"sex"`).

- age:

  Age column: one column with single ages (numeric) or labels such as
  `"30-34"` and `"90+"` (see
  [`parse_age_labels()`](parse_age_labels.md)), or two columns with the
  lower and upper bound (`Inf` for an open interval).

- period:

  Period column: calendar years (numeric), labels such as `"2010-2014"`,
  or two columns with lower and upper bounds.

- value:

  Name of the value column.

- age_open:

  Lower bound of an open age interval coded as a single numeric age (for
  example 110 for the Human Mortality Database "110+").

- missing:

  Values that mark missing cells in character value columns, for example
  the symbols of official statistics for confidential or unknown counts
  (`"."`, `"x"`). Missing cells are left out; suppressed death counts
  can then be estimated with
  [`mt_fill_suppressed()`](mt_fill_suppressed.md).

- zero:

  Values that mark a count of zero (for example `"-"`).

- decimal:

  Decimal mark of values given as text: `"."` (thousands separators such
  as `"12,345"` are recognised) or `","` (as in many European exports:
  `"1.234,5"`).

- ...:

  Arguments passed to `mt_add()`.

## Value

The updated `mt_data` object. An existing table for the same measure is
replaced.

## Examples

``` r
d <- data.frame(year = 2001, age = 30:31, deaths = c(5, 7), exposure = c(9e3, 9e3))
x <- mt_add_deaths(mt_data(), d)
x <- mt_add_exposure(x, d)
x
#> <mt_data> 1 stratum, omega = 110
#>   deaths     single years 30-31; periods 2001-2001 (width 1); 2 cells
#>   exposure   single years 30-31; periods 2001-2001 (width 1); 2 cells
#>   grid       <mt_axis: age, 2 intervals> 30 31
#>              <mt_axis: period, 1 interval> 2001
```
