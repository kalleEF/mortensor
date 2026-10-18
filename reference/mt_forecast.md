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

  If `TRUE`, draws describe observable rates: for models of expected
  rates (such as [`lee_carter()`](lee_carter.md)), and for any model
  fitted to smoothed rates ([`step_smooth()`](mt_prep.md)), deaths are
  drawn given the rates of each draw and `exposure` and divided by the
  exposures. The deaths are Poisson, the observation model of the
  Poisson Lee-Carter model (Brouhns, Denuit and Vermunt 2002); for
  models with a negative binomial likelihood (`stanmomo(family = "nb")`)
  they are negative binomial with the model's dispersion. No other
  observation noise is added. Models of observed rates
  ([`rwd()`](rwd.md), [`arima_series()`](arima_series.md),
  [`ets_series()`](arima_series.md)) add nothing. Cells can have zero
  deaths (rate 0). Used by [`mt_backtest()`](mt_backtest.md), which sums
  the deaths to the observed cells before taking logarithms. If
  predictive intervals are too narrow on real data, that is a property
  of the model: see [`recalibrated()`](recalibrated.md),
  [`ensemble()`](ensemble.md) or overdispersed models such as
  `stanmomo(family = "nb")`.

- exposure:

  Exposures for predictive draws: array `[age, period, stratum]` or
  matrix `[age, stratum]` (default: last training year).

- ...:

  Passed to methods.

## Value

An [mt_draws](mt_draws.md) object.

## References

Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear
regression approach to the construction of projected lifetables.
*Insurance: Mathematics and Economics* 31(3): 373-393.
[doi:10.1016/S0167-6687(02)00185-3](https://doi.org/10.1016/S0167-6687%2802%2900185-3)

## Examples

``` r
x <- mt_example(deaths = "single")
fc <- mt_forecast(mt_fit(rwd(), x), h = 5, n_draws = 100, seed = 1)
fc
#> <mt_draws> Random walk with drift; origin 2019; 71 ages x 5 periods (2020-2024) x 2 strata x 100 draws
```
