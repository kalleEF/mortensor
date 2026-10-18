# Mortality data observed at any granularity

`mt_data()` stores death counts, exposures (person-years at risk) and
populations for one or more strata. Each table may use its own age and
period resolution, for example deaths in 5-year age groups with an open
90+ interval and exposures in single years of age. The container keeps
the observations as they were supplied; harmonisation onto a common grid
happens explicitly with [`mt_regroup()`](mt_regroup.md) and
[`mt_ungroup()`](mt_ungroup.md).

## Usage

``` r
mt_data(
  deaths = NULL,
  exposure = NULL,
  population = NULL,
  by = NULL,
  age = "age",
  period = "year",
  value = NULL,
  age_open = NULL,
  omega = 110,
  missing = c(".", "..", "...", ":", "x", "X", "/", "NA", ""),
  zero = "-",
  decimal = "."
)
```

## Arguments

- deaths, exposure, population:

  Data frames in long format, one row per cell, or `NULL`.

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

  Name of the value column. Either one name used for all tables, or a
  named vector or list such as `c(deaths = "D", exposure = "E")`. By
  default the column named after the measure, otherwise `"value"`. A
  table in wide format, with one value column per level of a stratum
  variable (for example columns `Female` and `Male`, or one column per
  cause of death), is read by giving all its value columns, for example
  `value = list(deaths = c(female = "Female", male = "Male"))`: the
  columns become the levels (named by the names of the vector, else the
  column names) of the one variable in `by` that the table lacks.

- age_open:

  Lower bound of an open age interval coded as a single numeric age (for
  example 110 for the Human Mortality Database "110+").

- omega:

  Assumed maximum age, used wherever an open interval needs a finite
  upper bound, for example [`mt_ungroup()`](mt_ungroup.md), which
  expands an open interval into single years up to `omega`. Methods that
  need one representative age per cell ([`mt_smooth()`](mt_smooth.md)
  with the Whittaker smoother or GAMs) do not use it: they place an open
  interval at its lower bound plus half the width of the interval below
  it (a package convention).

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

## Value

An object of class `mt_data`.

## Details

Age labels of many sources are recognised (see
[`parse_age_labels()`](parse_age_labels.md)); periods are years, ranges
of years or financial years (`"2001/02"`). Rows whose age or period is
labelled as a total (`"All ages"`, `"Total"`, `"Insgesamt"`) are left
out with a warning. Values of unknown age (`"Unknown"`, `"not stated"`,
`"unbekannt"`, `NA`) are distributed over the age groups of the same
stratum and period in proportion to the known values, the usual practice
of the Human Mortality Database, with a warning. Rows of unknown period
are left out with a warning. Leading and trailing spaces in stratum
labels are removed; rows without a stratum label are left out. A warning
lists levels of a stratum variable that occur in one table only (for
example `"F"` in one table and `"female"` in another).

A table may leave out stratum columns over which its values are shared.
For example, deaths by region, sex and cause of death and populations by
region and sex: every cause of death then uses the population of its
region and sex. Tables with all stratum columns are read first and
define the strata.

## See also

[`mt_add_deaths()`](mt_add.md), [`mt_grid()`](mt_strata.md),
[`mt_regroup()`](mt_regroup.md), [`mt_ungroup()`](mt_ungroup.md),
[`as_mt_data()`](as_mt_data.md).

## Examples

``` r
d <- data.frame(sex = "female", year = 2001:2002, age = "90+",
                deaths = c(120, 130), exposure = c(800, 820))
mt_data(deaths = d, exposure = d, by = "sex")
#> <mt_data> 1 stratum by sex, omega = 110
#>   deaths     open interval 90+; periods 2001-2002 (width 1); 2 cells
#>   exposure   open interval 90+; periods 2001-2002 (width 1); 2 cells
#>   grid       <mt_axis: age, 1 interval> 90+
#>              <mt_axis: period, 2 intervals> 2001 2002
```
