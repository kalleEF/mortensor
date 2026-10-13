# Step between the last observed year and the first forecast year

Reports, per stratum, how far the first forecast year lies from the last
observed and the last fitted year: for each age the log ratio of the
forecast median to the observed and to the fitted rate, and for life
expectancy the difference in years. Large steps against the observations
point to a jump-off problem (see [`anchored()`](anchored.md)); steps
against the fitted rates are the forecast of one year.

## Usage

``` r
mt_jump_off(draws, quantity = c("life_expectancy", "rate"), at = NULL)
```

## Arguments

- draws:

  An [mt_draws](mt_draws.md) object.

- quantity:

  `"life_expectancy"` (one row per stratum) or `"rate"` (one row per
  stratum and age).

- at:

  Age for life expectancy (default: the first grid age).

## Value

A `data.table` with the stratum columns, the observed and fitted values
of the last training year, the median of the first forecast year and the
steps (`step_observed`, `step_fitted`; log ratios for rates, years for
life expectancy).

## Examples

``` r
x <- mt_example(deaths = "single", years = 1995:2019)
fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 100, seed = 1)
mt_jump_off(fc)
#>       sex  label observed   fitted forecast step_observed step_fitted
#>    <char> <char>    <num>    <num>    <num>         <num>       <num>
#> 1: female female 55.95597 55.95828  56.0453    0.08932909  0.08701991
#> 2:   male   male 51.16230 51.17137  51.2954    0.13309758  0.12402871
head(mt_jump_off(fc, quantity = "rate"))
#>       sex  label    age     observed       fitted     forecast step_observed
#>    <char> <char> <char>        <num>        <num>        <num>         <num>
#> 1: female female     30 0.0004823472 0.0003838037 0.0003757812   -0.24965696
#> 2: female female     31 0.0005561083 0.0004786193 0.0004728507   -0.16218331
#> 3: female female     32 0.0005226850 0.0005191787 0.0005142912   -0.01618946
#> 4: female female     33 0.0005526062 0.0005258894 0.0005194164   -0.06193969
#> 5: female female     34 0.0005105983 0.0005136473 0.0005053864   -0.01025989
#> 6: female female     35 0.0006219789 0.0005503272 0.0005418455   -0.13792518
#>     step_fitted
#>           <num>
#> 1: -0.021124177
#> 2: -0.012125610
#> 3: -0.009458592
#> 4: -0.012385119
#> 5: -0.016213487
#> 6: -0.015532033
```
