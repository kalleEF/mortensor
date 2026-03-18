# Generate bootstrap simulations for any model

Routes to the appropriate bootstrap function based on the model's method
type.

## Usage

``` r
dispatch_model_sims(
  model_name,
  mx,
  pop,
  years,
  ages,
  groups,
  h,
  n_sim,
  fc_result,
  width = 5L,
  calib_all = NULL,
  mx_raw = NULL
)
```

## Arguments

- model_name:

  Model name string.

- mx:

  3D smoothed mortality array.

- pop:

  3D population array.

- years, ages:

  Integer vectors.

- groups:

  data.frame of groups.

- h:

  Forecast horizon.

- n_sim:

  Number of simulations.

- fc_result:

  Output of
  [`dispatch_forecast()`](https://kalleEF.github.io/mortensor/reference/dispatch_forecast.md).

- width:

  Age-group width.

- calib_all:

  Calibration result for non-Tucker3 models.

- mx_raw:

  Optional raw (unsmoothed) rates.

## Value

4D array \[n_a x h x n_g x n_sim\].
