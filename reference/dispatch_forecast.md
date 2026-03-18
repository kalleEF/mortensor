# Dispatch a full-horizon forecast for any model

Routes to the appropriate forecasting function based on `model_name` and
returns a standardised result list.

## Usage

``` r
dispatch_forecast(
  model_name,
  mx,
  pop,
  years,
  ages,
  groups,
  h,
  calib_coda,
  calib_direct,
  width = 5L,
  calib_all = NULL,
  return_fitted = TRUE
)
```

## Arguments

- model_name:

  One of `"Tucker3-CoDa"`, `"Tucker3-Direct"`, `"Lee-Carter"`, `"FDM"`,
  `"StMoMo-APC"`, `"RWD"`, `"ARIMA-perSeries"`, `"ETS-Damped"`.

- mx:

  3D array of (smoothed) mortality rates.

- pop:

  3D population array.

- years:

  Integer vector.

- ages:

  Integer vector of age-group starts.

- groups:

  data.frame with `sex` and `quintile`.

- h:

  Forecast horizon.

- calib_coda:

  Tucker3-CoDa calibration result.

- calib_direct:

  Tucker3-Direct calibration result.

- width:

  Age-group width.

- calib_all:

  Optional calibration for non-Tucker3 models.

- return_fitted:

  Return fitted in-sample values?

## Value

List: `mx_forecast`, `mx_fitted`, `full_fc`, `method`.
