# Simulated deaths by cause

Extends [`mt_example()`](mt_example.md) with deaths by cause (`"cvd"`,
`"cancer"`, `"other"`) and the all-cause total (`"all"`), single years
of age. Cause shares change smoothly with age and time; deaths are split
multinomially, so causes add up to all causes exactly.

## Usage

``` r
mt_example_causes(years = 1990:2019, size = 1e+05, seed = 1)
```

## Arguments

- years:

  Calendar years.

- size:

  Population size at age 30 per year and stratum.

- seed:

  Random seed.

## Value

An [mt_data](mt_data.md) object with strata `sex` x `cause`.

## Examples

``` r
mt_example_causes(years = 2010:2019)
#> <mt_data> 8 strata by sex, cause, omega = 110
#>   deaths     single years 30-99, 100+; periods 2010-2019 (width 1); 5680 cells
#>   exposure   single years 30-99, 100+; periods 2010-2019 (width 1); 5680 cells
#>   grid       <mt_axis: age, 71 intervals> 30 31 32 33 ... 98 99 100+
#>              <mt_axis: period, 10 intervals> 2010 2011 2012 2013 ... 2017 2018 2019
```
