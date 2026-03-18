# Create a mortality forecast object

Create a mortality forecast object

## Usage

``` r
mort_forecast(
  mx_point,
  train,
  forecast_years,
  model,
  mx_median = NULL,
  mx_lower = NULL,
  mx_upper = NULL,
  sims = NULL,
  mx_fitted = NULL,
  mx_raw = NULL,
  fit_info = list(),
  extras = list()
)
```

## Arguments

- mx_point:

  3D array `[n_ages x h x n_groups]` of point forecasts.

- train:

  A `mort_data` object used for training (provides metadata).

- forecast_years:

  Integer vector of forecast years.

- model:

  Character string naming the forecast model.

- mx_median:

  Optional 3D array of bootstrap median forecasts.

- mx_lower:

  Optional 3D array of lower PI bound.

- mx_upper:

  Optional 3D array of upper PI bound.

- sims:

  Optional 4D array `[n_ages x h x n_groups x n_sim]` of bootstrap
  simulations.

- mx_fitted:

  Optional 3D array `[n_ages x n_train_years x n_groups]` of in-sample
  fitted values.

- mx_raw:

  Optional 3D array of raw (unsmoothed) training mortality rates, used
  for residual computation.

- fit_info:

  Optional named list of model fit information, e.g. `fit_pct`, `mase`,
  `model_comparison`.

- extras:

  Named list of additional model-specific outputs (Tucker decomposition
  components, ARIMA fits, etc.).

## Value

A `mort_forecast` object.
