# Simulated deaths by sex, region and deprivation quintile

A grouped structure as in small-area mortality data: `sex` x `region` x
`quintile` bottom strata. Log rates follow a Lee-Carter-type model with
a period index per sex, region-specific deviations (AR(1)), a
deprivation gradient that widens slowly over time, Gompertz-Makeham age
patterns and regional population sizes that differ by a factor of up to
four. Deaths are Poisson. True rates of the bottom strata are attached
as attribute `"truth"`.

## Usage

``` r
mt_example_regions(
  years = 1995:2019,
  regions = 4,
  ages = c("grouped", "single"),
  size = 20000,
  seed = 1
)
```

## Arguments

- years:

  Calendar years.

- regions:

  Number of regions.

- ages:

  `"grouped"` (5-year groups 40-44 ... 85-89 and 90+, also used as grid)
  or `"single"` (single years 40-99 and 100+).

- size:

  Population at age 40 per bottom stratum in the largest region.

- seed:

  Random seed.

## Value

An [mt_data](mt_data.md) object with strata `sex` x `region` x
`quintile` (bottom level only; add aggregates with
[`mt_aggregate_strata()`](mt_aggregate_strata.md)).

## Examples

``` r
x <- mt_example_regions(years = 2010:2019, regions = 3)
x
#> <mt_data> 30 strata by sex, region, quintile, omega = 110
#>   deaths     5-year groups 40-89, 90+; periods 2010-2019 (width 1); 3300 cells
#>   exposure   5-year groups 40-89, 90+; periods 2010-2019 (width 1); 3300 cells
#>   grid       <mt_axis: age, 11 intervals> 40-44 45-49 50-54 55-59 ... 80-84 85-89 90+
#>              <mt_axis: period, 10 intervals> 2010 2011 2012 2013 ... 2017 2018 2019
```
