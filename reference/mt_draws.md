# Forecast draws

An `mt_draws` object holds simulated death rates as an array
`[age, period, stratum, draw]` on the modelling grid, the observed rates
of the training data, and metadata (model, forecast origin, seed,
diagnostics). Summaries are computed from the paths: `mt_summarise()`
gives the mean and quantiles per cell,
[`life_expectancy()`](life_expectancy.md) gives life expectancy path by
path.

## Usage

``` r
mt_summarise(draws, probs = c(0.025, 0.1, 0.5, 0.9, 0.975))

# S3 method for class 'mt_draws'
as.data.table(x, keep.rownames = FALSE, ...)
```

## Arguments

- probs:

  Quantile levels.

- x, draws:

  An `mt_draws` object.

- keep.rownames, ...:

  Ignored.

## Value

`mt_summarise()` and `as.data.table()` return a `data.table` with the
stratum columns, `age_lo`, `age`, `period`, `mean` and one column per
quantile (`q0.025`, ...).

## Examples

``` r
x <- mt_example(deaths = "single")
fc <- mt_forecast(mt_fit(rwd(), x), h = 3, n_draws = 50, seed = 1)
head(mt_summarise(fc))
#>    stratum    sex  label age_lo    age period         mean       q0.025
#>      <int> <char> <char>  <num> <char>  <num>        <num>        <num>
#> 1:       1 female female     30     30   2020 0.0003554981 0.0002247483
#> 2:       1 female female     31     31   2020 0.0004321882 0.0002515748
#> 3:       1 female female     32     32   2020 0.0005193313 0.0003676016
#> 4:       1 female female     33     33   2020 0.0005252941 0.0003741907
#> 5:       1 female female     34     34   2020 0.0004861726 0.0003352379
#> 6:       1 female female     35     35   2020 0.0005551962 0.0004175027
#>            q0.1         q0.5         q0.9       q0.975
#>           <num>        <num>        <num>        <num>
#> 1: 0.0002690735 0.0003214454 0.0004830297 0.0005559621
#> 2: 0.0002938909 0.0004074861 0.0005915575 0.0007685350
#> 3: 0.0004030857 0.0005291799 0.0006550429 0.0006888401
#> 4: 0.0003881232 0.0005262917 0.0006576370 0.0007132757
#> 5: 0.0003744938 0.0004479962 0.0006363024 0.0008072660
#> 6: 0.0004438184 0.0005546124 0.0006482273 0.0007395214
```
