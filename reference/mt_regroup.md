# Regroup data onto coarser age or period intervals

Aggregates every table exactly onto the target axes. Counts (deaths and
exposures) are summed; populations are summed over age and averaged over
periods (weighted by period width). An observed cell that straddles a
target boundary raises an error, because exact aggregation is not
possible. Target cells that are only partly covered by observations are
dropped with a warning. The grid becomes the target axes, restricted to
the range of the data.

## Usage

``` r
mt_regroup(data, age = NULL, period = NULL, measures = NULL)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- age:

  Target age axis ([mt_axis](mt_axis.md)), numeric lower bounds, the
  name of a table (`"deaths"`, `"exposure"`, `"population"`) whose age
  intervals are the target (for example deaths in 5-year groups and
  populations in single years: `age = "deaths"`), or `"common"`: the
  finest groups that every table can be aggregated to (for example 85+
  when the deaths of some years end with 85+ and others with 85-89 and
  90+); `NULL` keeps the age intervals.

- period:

  Target period axis, numeric lower bounds or the name of a table;
  `NULL` keeps the periods.

- measures:

  Tables to regroup (default: all). Regrouping only some tables, for
  example `measures = "deaths"`, keeps the other tables and the grid as
  they are: deaths in age groups with single-year exposures, as in many
  official tables, for example to check ungrouping against known single
  years.

## Value

An `mt_data` object.

## Examples

``` r
x <- mt_example()
mt_regroup(x, age = c(30, seq(40, 90, 10)))
#> <mt_data> 2 strata by sex, omega = 110
#>   deaths     10-year groups 30-89, 90+; periods 1990-2019 (width 1); 420 cells
#>   exposure   10-year groups 30-89, 90+; periods 1990-2019 (width 1); 420 cells
#>   grid       <mt_axis: age, 7 intervals> 30-39 40-49 50-59 60-69 70-79 80-89 90+
#>              <mt_axis: period, 30 intervals> 1990 1991 1992 1993 ... 2017 2018 2019
mt_regroup(x, age = "deaths")
#> <mt_data> 2 strata by sex, omega = 110
#>   deaths     5-year groups 30-89, 90+; periods 1990-2019 (width 1); 780 cells
#>   exposure   5-year groups 30-89, 90+; periods 1990-2019 (width 1); 780 cells
#>   grid       <mt_axis: age, 13 intervals> 30-34 35-39 40-44 45-49 ... 80-84 85-89 90+
#>              <mt_axis: period, 30 intervals> 1990 1991 1992 1993 ... 2017 2018 2019
xs <- mt_example(deaths = "single")
mt_regroup(xs, age = c(30, seq(40, 90, 10)), measures = "deaths")
#> <mt_data> 2 strata by sex, omega = 110
#>   deaths     10-year groups 30-89, 90+; periods 1990-2019 (width 1); 420 cells
#>   exposure   single years 30-99, 100+; periods 1990-2019 (width 1); 4260 cells
#>   grid       <mt_axis: age, 71 intervals> 30 31 32 33 ... 98 99 100+
#>              <mt_axis: period, 30 intervals> 1990 1991 1992 1993 ... 2017 2018 2019
```
