# Restrict data to ages, periods and strata

Restrict data to ages, periods and strata

## Usage

``` r
mt_subset(data, ..., ages = NULL, periods = NULL, drop = FALSE)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- ...:

  Conditions on the stratum variables, for example `region != "Germany"`
  or `sex == "female"`.

- ages:

  Age range `c(lower, upper)`: cells with lower bound at or above
  `lower` and upper bound at or below `upper` (`Inf` for an open
  interval) are kept.

- periods:

  Period range `c(first, last)` (calendar years, inclusive; `NA` for an
  open end, for example `c(NA, 2019)`).

- drop:

  Remove stratum variables that have a single level left (for example
  `cause` after `cause == "all"`).

## Value

The restricted `mt_data` object; strata are renumbered and the grid, if
set, is restricted to the same ranges.

## Examples

``` r
x <- mt_example(years = 2001:2019)
y <- mt_subset(x, sex == "female", ages = c(60, Inf), periods = c(2010, 2019))
y
#> <mt_data> 1 stratum by sex, omega = 110
#>   deaths     5-year groups 60-89, 90+; periods 2010-2019 (width 1); 70 cells
#>   exposure   single years 60-99, 100+; periods 2010-2019 (width 1); 410 cells
#>   grid       <mt_axis: age, 41 intervals> 60 61 62 63 ... 98 99 100+
#>              <mt_axis: period, 10 intervals> 2010 2011 2012 2013 ... 2017 2018 2019
```
