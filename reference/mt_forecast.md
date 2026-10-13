# Forecast draws from a fitted model

Simulates future death rates on the grid. The result holds simulated
paths (draws), so that uncertainty propagates correctly to derived
quantities such as life expectancy.

## Usage

``` r
mt_forecast(
  fit,
  h = 10,
  to = NULL,
  n_draws = 500,
  seed = NULL,
  predictive = FALSE,
  exposure = NULL,
  ...
)
```

## Arguments

- fit:

  An `mt_fit` object.

- h:

  Forecast horizon in years.

- to:

  Alternatively, the last year to forecast.

- n_draws:

  Number of simulated paths.

- seed:

  Random seed (`NULL` leaves the random number generator as is).

- predictive:

  If `TRUE`, draws describe observable rates: for models of smooth
  rates, extra-Poisson noise (estimated from in-sample residuals) and
  Poisson variation of deaths given `exposure` are added, and zero
  counts are replaced by 0.5 as in the training data. Used by
  [`mt_backtest()`](mt_backtest.md).

- exposure:

  Exposures for predictive draws: array `[age, period, stratum]` or
  matrix `[age, stratum]` (default: last training year).

- ...:

  Passed to methods.

## Value

An [mt_draws](mt_draws.md) object.

## Examples

``` r
x <- mt_example(deaths = "single")
fc <- mt_forecast(mt_fit(rwd(), x), h = 5, n_draws = 100, seed = 1)
fc
#> <mt_draws> Random walk with drift; origin 2019; 71 ages x 5 periods (2020-2024) x 2 strata x 100 draws
```
