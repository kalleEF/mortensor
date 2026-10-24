# Model specifications

Every model is created by a constructor ([`rwd()`](rwd.md),
[`lee_carter()`](lee_carter.md), [`tucker3()`](tucker3.md),
[`arima_series()`](arima_series.md), [`ets_series()`](arima_series.md),
[`apc_bayes()`](apc_bayes.md), [`ensemble()`](ensemble.md)) that returns
an `mt_spec` object. Specifications are fitted with
[`mt_fit()`](mt_fit.md) and forecast with
[`mt_forecast()`](mt_forecast.md).

## Usage

``` r
new_mt_spec(
  model,
  params = list(),
  joint = FALSE,
  min_periods = 5,
  label = model,
  smooth = TRUE,
  seed = NULL
)
```

## Arguments

- model:

  Model name (used in the class).

- params:

  List of parameters.

- joint:

  Whether the model is fitted jointly across strata.

- min_periods:

  Minimum number of periods in the training data.

- label:

  Display name.

- smooth:

  Whether the draws describe expected (underlying) rates (`TRUE`, for
  example Lee-Carter: predictive draws add Poisson variation of the
  deaths given these rates and the exposures) or observed rates
  (`FALSE`, for example a random walk on observed log rates: predictive
  draws add nothing).

- seed:

  Random seed of the fit, for models whose fit uses random numbers: the
  random number streams of the strata in [`mt_fit()`](mt_fit.md) start
  from this seed. `NULL`: they start from the current generator.

## Value

An object of class `c("mt_spec_<model>", "mt_spec")`.

## Details

`new_mt_spec()` is the constructor for developers who add models: a
model class `mt_spec_<name>` needs methods for the internal generics
`fit_model()` and `simulate_model()`. Per-stratum models receive
matrices `[age, period]` of deaths and exposures; joint models
(`joint = TRUE`) receive arrays `[age, period, stratum]`.
`simulate_model()` returns rates `[age, horizon, draw]` (per stratum) or
`[age, horizon, stratum, draw]` (joint).
