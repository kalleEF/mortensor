# Simulated mortality data with known truth

Generates a two-stratum dataset (women and men) from a Lee-Carter-type
model with Gompertz-Makeham age pattern, a random walk with drift for
the period index, and a population with cohort-size waves. Deaths are
Poisson draws. By default deaths are reported in 5-year age groups with
an open 90+ interval and exposures in single years of age with an open
100+ interval, as in many official statistics. A model fits these data
after the deaths are ungrouped ([`mt_ungroup()`](mt_ungroup.md) or
[`step_ungroup()`](mt_prep.md)) or all tables are aggregated to common
age groups ([`mt_regroup()`](mt_regroup.md)); the examples of the other
help pages pass `deaths = "single"` (deaths in single years of age) so
that every model fits them directly.

## Usage

``` r
mt_example(
  years = 1990:2019,
  deaths = c("grouped", "single"),
  size = 1e+05,
  seed = 1
)
```

## Arguments

- years:

  Calendar years.

- deaths:

  `"grouped"` (5-year groups and 90+) or `"single"` (single ages and
  100+).

- size:

  Population size at age 30 per year and stratum.

- seed:

  Random seed.

## Value

An [mt_data](mt_data.md) object with attribute `"truth"`.

## Details

The true rates on the default grid are attached as attribute `"truth"`
(array `[age, period, stratum]`), which makes the dataset suitable for
testing ungrouping and forecast calibration.

## Examples

``` r
x <- mt_example()
x
#> <mt_data> 2 strata by sex, omega = 110
#>   deaths     5-year groups 30-89, 90+; periods 1990-2019 (width 1); 780 cells
#>   exposure   single years 30-99, 100+; periods 1990-2019 (width 1); 4260 cells
#>   grid       <mt_axis: age, 71 intervals> 30 31 32 33 ... 98 99 100+
#>              <mt_axis: period, 30 intervals> 1990 1991 1992 1993 ... 2017 2018 2019
dim(attr(x, "truth"))
#> [1] 71 30  2
```
