# Strata, grid and tables of an `mt_data` object

`mt_strata()` returns the stratum table (with a column `level`,
`"bottom"`, `"aggregate"` or `"total"`, when a hierarchy is given),
`mt_grid()` the modelling grid (the coarsest grid on which every table
is exact, see [`mt_set_grid()`](mt_set_grid.md)), `mt_obs()` one
observation table.

## Usage

``` r
mt_strata(data, hierarchy = NULL)

mt_obs(data, measure = "deaths")

mt_grid(data)

# S3 method for class 'mt_data'
as.data.table(x, keep.rownames = FALSE, ...)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object (for `mt_strata()` also an
  [mt_draws](mt_draws.md) object).

- hierarchy:

  Optional [`mt_hierarchy()`](mt_hierarchy.md) for the `level` column.

- measure:

  One of `"deaths"`, `"exposure"`, `"population"`.

- x:

  An [mt_data](mt_data.md) object.

- keep.rownames, ...:

  Ignored.

## Value

A `data.table` or, for `mt_grid()`, a list with `age` and `period` axes.

## Examples

``` r
x <- mt_example()
mt_strata(x)
#>    stratum    sex  label
#>      <int> <char> <char>
#> 1:       1 female female
#> 2:       2   male   male
mt_grid(x)
#> $age
#> <mt_axis: age, 71 intervals> 30 31 32 33 ... 98 99 100+
#> 
#> $period
#> <mt_axis: period, 30 intervals> 1990 1991 1992 1993 ... 2017 2018 2019
#> 
head(mt_obs(x, "deaths"))
#>    stratum    sex  label age_lo age_hi period_lo period_hi value provenance
#>      <int> <char> <char>  <num>  <num>     <num>     <num> <num>     <char>
#> 1:       1 female female     30     35      1990      1991   486   observed
#> 2:       1 female female     35     40      1990      1991   590   observed
#> 3:       1 female female     40     45      1990      1991   536   observed
#> 4:       1 female female     45     50      1990      1991   729   observed
#> 5:       1 female female     50     55      1990      1991  1449   observed
#> 6:       1 female female     55     60      1990      1991  2704   observed
```
