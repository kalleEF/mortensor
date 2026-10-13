# Observed rates as draws

Turns observed deaths and exposures into an [mt_draws](mt_draws.md)
object whose paths are Poisson draws of the death counts given the
exposures. All functions for draws then also describe the past, with
intervals for the observed rates:
[`life_expectancy()`](life_expectancy.md), [`asr()`](derived.md),
[`prob_dying()`](derived.md), [`mt_contrast()`](mt_contrast.md),
[`mt_inequality()`](mt_inequality.md). The derived quantities accept an
[mt_data](mt_data.md) object directly and call `mt_observed()` with its
defaults.

## Usage

``` r
mt_observed(data, n_draws = 500, seed = 1)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object with deaths and exposures on the grid.

- n_draws:

  Number of draws.

- seed:

  Random seed (default 1, so that results are reproducible).

## Value

An `mt_draws` object with the observed periods; draw 1 holds the
observed rates.

## Examples

``` r
x <- mt_example(deaths = "single", years = 2015:2019)
life_expectancy(x, at = 65)
#>     stratum    sex  label period     mean   q0.025     q0.1     q0.5     q0.9
#>       <int> <char> <char>  <num>    <num>    <num>    <num>    <num>    <num>
#>  1:       1 female female   2015 22.10007 22.03790 22.05694 22.09804 22.14582
#>  2:       1 female female   2016 22.20439 22.13830 22.15788 22.20354 22.24883
#>  3:       1 female female   2017 22.25782 22.18558 22.21392 22.25935 22.30164
#>  4:       1 female female   2018 22.42149 22.35275 22.37736 22.42280 22.46542
#>  5:       1 female female   2019 22.45561 22.38683 22.41480 22.45560 22.49756
#>  6:       2   male   male   2015 18.06814 18.00324 18.02954 18.06775 18.10818
#>  7:       2   male   male   2016 18.17609 18.11630 18.13820 18.17605 18.21950
#>  8:       2   male   male   2017 18.23438 18.17299 18.19496 18.23477 18.27489
#>  9:       2   male   male   2018 18.35115 18.29143 18.31039 18.35179 18.39155
#> 10:       2   male   male   2019 18.43478 18.37806 18.39382 18.43281 18.47447
#>       q0.975    at
#>        <num> <num>
#>  1: 22.16661    65
#>  2: 22.26815    65
#>  3: 22.31990    65
#>  4: 22.48665    65
#>  5: 22.51892    65
#>  6: 18.12787    65
#>  7: 18.24392    65
#>  8: 18.29418    65
#>  9: 18.40842    65
#> 10: 18.50264    65
```
