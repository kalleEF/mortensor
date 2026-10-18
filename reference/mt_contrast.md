# Differences and ratios between strata, path by path

Compares a derived quantity (life expectancy, age-standardised rate,
probability of dying, death rate of one age group, or an index of
lifespan variation) between the levels of a stratum variable and a
reference level, path by path, so that the interval of the difference
includes the joint uncertainty of both strata. Strata are matched on all
other stratum variables (for example the deprivation gap within each
sex).

## Usage

``` r
mt_contrast(
  x,
  var,
  ref,
  quantity = c("life_expectancy", "asr", "prob_dying", "rate", "lifespan_variation"),
  type = c("difference", "ratio"),
  probs = c(0.025, 0.5, 0.975),
  ...
)
```

## Arguments

- x:

  An [mt_draws](mt_draws.md) object (forecast) or an
  [mt_data](mt_data.md) object (observed, see
  [`mt_observed()`](mt_observed.md)).

- var:

  Stratum variable to compare across (for example `"quintile"`).

- ref:

  Reference level of `var`.

- quantity:

  `"life_expectancy"`, `"asr"`, `"prob_dying"`, `"rate"` or
  `"lifespan_variation"`.

- type:

  `"difference"` (level minus reference) or `"ratio"`.

- probs:

  Quantile levels.

- ...:

  Arguments of the quantity: `at` for life expectancy, `standard`,
  `per`, `exposure` and `age_groups` for [`asr()`](derived.md), `from`
  and `to` for [`prob_dying()`](derived.md), `age` (an age label) for
  `"rate"`, `measure` (one index, default `"edagger"`) and `at` for
  [`lifespan_variation()`](lifespan_variation.md).

## Value

A `data.table` with the other stratum variables, `level`, `ref`,
`period`, `mean`, quantiles and `p_positive`, the share of paths with a
positive difference (or a ratio above 1). For forecasts this is the
predictive probability of a positive difference under the model. For
observed data it is the share of the Poisson bootstrap draws (see
[`mt_observed()`](mt_observed.md)) with a positive difference, which is
neither a p-value nor a posterior probability.

## Examples

``` r
x <- mt_example(deaths = "single", years = 2001:2015)
fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 100, seed = 1)
mt_contrast(fc, var = "sex", ref = "female", at = 30)
#>     level    ref period      mean    q0.025      q0.5    q0.975 p_positive
#>    <char> <char>  <num>     <num>     <num>     <num>     <num>      <num>
#> 1:   male female   2016 -5.045323 -5.267886 -5.043577 -4.815912          0
#> 2:   male female   2017 -5.066032 -5.465374 -5.047023 -4.744987          0
#> 3:   male female   2018 -5.048624 -5.553453 -5.048755 -4.584807          0
#> 4:   male female   2019 -5.046785 -5.582633 -5.084452 -4.584923          0
#> 5:   male female   2020 -5.053376 -5.656648 -5.065486 -4.416842          0
#>           quantity       type
#>             <char>     <char>
#> 1: life_expectancy difference
#> 2: life_expectancy difference
#> 3: life_expectancy difference
#> 4: life_expectancy difference
#> 5: life_expectancy difference
```
