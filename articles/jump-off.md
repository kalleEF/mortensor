# Where a forecast starts: the jump-off

A forecast continues from some level in the last observed year, the
*jump-off*. Most models start from their own fitted rates of that year,
which differ from the observed rates: a model with few parameters does
not reproduce every age, year and stratum. The difference shows up as a
step between the last observed and the first forecast year. This
vignette shows when the step is large, the three ways of placing the
start, and what they change.

``` r

library(mortensor)
library(data.table)
x <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile"), age = "age_group",
             value = c(deaths = "all_causes", population = "population"))
x <- mt_as_exposure(x, type = "midyear")
x <- mt_group_levels(x, "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)), name = "quintile")
```

## When the step is large

The step is the distance between the model’s fitted rates and the
observed rates in the last year. Lee-Carter fitted to each stratum
separately has an age pattern and a period index of its own for every
stratum and follows each year closely. Li-Lee pools the quintiles of
each sex: they share one trend and deviate from it by a stationary
factor, so the fit of a single quintile in a single year is looser.
[`mt_jump_off()`](../reference/mt_jump_off.md) reports life expectancy
at birth in the last observed year (2019), as fitted, and in the first
forecast year:

``` r

fc_lc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 500, seed = 1)
fc_ll <- mt_forecast(mt_fit(li_lee(common_by = "sex"), x), h = 10, n_draws = 500, seed = 1)
merge(mt_jump_off(fc_lc)[, .(label, observed, lee_carter_fitted = fitted)],
      mt_jump_off(fc_ll)[, .(label, li_lee_fitted = fitted, li_lee_2020 = forecast)], by = "label")[
  , lapply(.SD, function(v) if (is.numeric(v)) round(v, 2) else v)]
#> Key: <label>
#>         label observed lee_carter_fitted li_lee_fitted li_lee_2020
#>        <char>    <num>             <num>         <num>       <num>
#>  1: female.Q1    79.83             79.85         80.26       80.92
#>  2: female.Q2    82.50             82.51         82.74       82.90
#>  3: female.Q3    83.93             83.97         84.01       84.17
#>  4: female.Q4    85.02             85.05         84.88       85.05
#>  5: female.Q5    86.30             86.29         85.95       86.02
#>  6:   male.Q1    75.43             75.48         75.84       76.10
#>  7:   male.Q2    78.60             78.60         78.74       78.96
#>  8:   male.Q3    80.61             80.66         80.56       80.78
#>  9:   male.Q4    81.82             81.84         81.75       81.96
#> 10:   male.Q5    83.09             83.15         83.02       83.20
```

Lee-Carter reproduces the observed life expectancy of 2019 to a few
hundredths of a year. Li-Lee puts women in Q1 0.4 years above the
observed value and women in Q5 0.3 years below it: it pulls the
quintiles towards their common trend. For women in Q1 its first forecast
year then adds another 0.7 years. The step carries into every forecast
year.

## Three starting points

- **fitted** (default): the forecast continues from the model’s fitted
  rates of the last year;
- **observed**: it continues from the observed rates of the last year
  (`li_lee(jump_off = "observed")`,
  `lee_carter(jump_off = "observed")`), with all their noise;
- **anchored**: [`anchored()`](../reference/anchored.md) moves the start
  of any model to the observed level of the last years (by default
  three), shrinks the offsets of cells with few deaths towards the mean
  of the stratum, and lets the offsets fade out over `decay` years (by
  default ten), so that the forecast returns to the model’s trend.

``` r

fc <- list(fitted = fc_ll,
           observed = mt_forecast(mt_fit(li_lee(common_by = "sex", jump_off = "observed"), x), h = 10,
                                  n_draws = 500, seed = 1),
           anchored = mt_forecast(mt_fit(anchored(li_lee(common_by = "sex"), years = 3), x), h = 10,
                                  n_draws = 500, seed = 1))
```

``` r

mt_plot_compare(fc, strata = c("female.Q1", "female.Q5"), periods = c(2010, 2029), data = x)
```

![](jump-off_files/figure-html/unnamed-chunk-4-1.png)

The solid line is the fitted Li-Lee model, the points the observed life
expectancy. The fitted forecast continues from the solid line, above the
last points for women in Q1. The observed and the anchored forecasts
start at the observed level. By age, the observed start also carries the
noise of single cells, which the anchored start averages over three
years:

``` r

mt_plot_compare(fc, quantity = "rate", ages = c("45-49", "65-69", "85-89"), strata = c("female.Q1", "female.Q5"),
                periods = c(2010, 2025))
```

![](jump-off_files/figure-html/unnamed-chunk-5-1.png)

## Which start forecasts better

A backtest fits each version to the years up to an origin and scores the
following five years:

``` r

bt <- mt_compare(list(fitted = li_lee(common_by = "sex"), observed = li_lee(common_by = "sex", jump_off = "observed"),
                      anchored = anchored(li_lee(common_by = "sex"), years = 3)),
                 x, origins = 2008:2014, h = 5, n_draws = 200, quantities = "life_expectancy")
mt_score_table(bt, metrics = c("crps", "cov95"))
#>       model   crps cov95 crps_le cov95_le
#>      <fctr>  <num> <num>   <num>    <num>
#> 1:   fitted 0.0666 0.847   0.336    0.811
#> 2: observed 0.0709 0.848   0.315    0.826
#> 3: anchored 0.0654 0.859   0.315    0.863
```

``` r

mt_plot_scores(bt, target = "life_expectancy")
```

![](jump-off_files/figure-html/unnamed-chunk-7-1.png)

`crps` scores the log death rates, `crps_le` life expectancy at birth in
years (lower is better; see “Backtesting and scores”). The anchored
start gives the most accurate forecasts of the rates and of life
expectancy. The observed start forecasts life expectancy about as well
but the rates worse, because it carries the noise of single ages and
years into the forecast. For Lee-Carter fitted to each stratum, whose
step is small here, the fitted and the anchored start give the same
scores (CRPS 0.068) and the observed start a worse one (0.078).

## Summary

| Start | How | When |
|----|----|----|
| fitted | default of every model | the model fits the last years closely |
| observed | `jump_off = "observed"` of [`lee_carter()`](../reference/lee_carter.md) and [`li_lee()`](../reference/li_lee.md) | large populations with little noise |
| anchored | `anchored(model, years = 3, decay = 10)` around any model | pooled or smooth models; noisy last years |

[`mt_jump_off()`](../reference/mt_jump_off.md) shows the step of any
forecast; a backtest shows whether another start forecasts better on the
data at hand.
