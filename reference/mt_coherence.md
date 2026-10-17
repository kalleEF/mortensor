# Check the coherence of forecast draws

For every aggregate stratum of a hierarchy, the ratio of the sum of its
bottom strata to the aggregate itself, on the deaths scale (rates times
the exposures of the last observed year, or `exposure`) or on the rate
scale for variables whose levels share an exposure (causes of death).
Coherent forecasts give ratios of 1 in every draw; separately forecast
strata show how far they disagree.

## Usage

``` r
mt_coherence(draws, hierarchy, exposure = NULL)

mt_plot_coherence(draws, hierarchy, exposure = NULL)
```

## Arguments

- draws:

  An [mt_draws](mt_draws.md) object with the strata of the hierarchy.

- hierarchy:

  An [`mt_hierarchy()`](mt_hierarchy.md).

- exposure:

  Optional exposures (see [`mt_deaths()`](mt_deaths.md)).

## Value

`mt_coherence()`: a `data.table` with the aggregate strata, `period`,
and the median and the 2.5% and 97.5% quantiles of the ratio over ages
and draws. `mt_plot_coherence()`: a ggplot object with the distribution
of the median ratios of the aggregates per year and the range of the
ratios of single paths.

## Examples

``` r
x <- mt_example_regions(years = 2005:2019, regions = 2)
h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
fc <- mt_forecast(mt_fit(lee_carter(), mt_aggregate_strata(x, h)), h = 3, n_draws = 30, seed = 1)
head(mt_coherence(fc, h))
#>       sex region quintile         label period    q0.025      q0.5   q0.975
#>    <char> <char>   <char>        <char>  <num>     <num>     <num>    <num>
#> 1: female     R1      all female.R1.all   2020 0.9432926 1.0016314 1.077586
#> 2: female     R1      all female.R1.all   2021 0.9328553 1.0017219 1.094624
#> 3: female     R1      all female.R1.all   2022 0.9103443 0.9938173 1.172363
#> 4: female     R2      all female.R2.all   2020 0.9408038 1.0014475 1.060165
#> 5: female     R2      all female.R2.all   2021 0.9102682 1.0096708 1.154123
#> 6: female     R2      all female.R2.all   2022 0.8838777 1.0034204 1.166966
head(mt_coherence(mt_reconcile(fc, h, method = "wls"), h))
#>       sex region quintile         label period q0.025  q0.5 q0.975
#>    <char> <char>   <char>        <char>  <num>  <num> <num>  <num>
#> 1: female     R1      all female.R1.all   2020      1     1      1
#> 2: female     R1      all female.R1.all   2021      1     1      1
#> 3: female     R1      all female.R1.all   2022      1     1      1
#> 4: female     R2      all female.R2.all   2020      1     1      1
#> 5: female     R2      all female.R2.all   2021      1     1      1
#> 6: female     R2      all female.R2.all   2022      1     1      1
```
