# Aggregate forecast draws over strata or ages

Computes draws for coarser groups from draws on the grid, path by path:
over strata along a hierarchy (bottom-up: deaths are summed and divided
by summed exposures) and over age groups (exposure-weighted mean rates).

## Usage

``` r
mt_aggregate_draws(draws, hierarchy = NULL, age = NULL, exposure = NULL)
```

## Arguments

- draws:

  An [mt_draws](mt_draws.md) object.

- hierarchy:

  Optional [`mt_hierarchy()`](mt_hierarchy.md); aggregate strata are
  added bottom-up from the bottom strata.

- age:

  Optional coarser [mt_axis](mt_axis.md) (cells must nest).

- exposure:

  Exposures `[age, period, stratum]` or `[age, stratum]` (default: last
  training year).

## Value

An `mt_draws` object.

## Examples

``` r
x <- mt_example(deaths = "single", years = 2000:2019)
fc <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 20, seed = 1)
fa <- mt_aggregate_draws(fc, hierarchy = mt_hierarchy(sex = "both"),
                         age = age_groups(c(30, seq(40, 100, 10))))
fa
#> <mt_draws> Lee-Carter; origin 2019; 8 ages x 3 periods (2020-2022) x 3 strata x 20 draws
```
