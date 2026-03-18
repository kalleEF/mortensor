# Calibrate all non-Tucker3 models

Grid-searches over hyperparameters for Lee-Carter, FDM, StMoMo-APC, RWD,
ARIMA, and ETS using expanding-window CV.

## Usage

``` r
calibrate_all_models(
  mx_array,
  pop_array,
  years,
  ages,
  groups,
  h_val = 4,
  n_origins = 4,
  width = 5,
  cv_metric = c("mase", "mase_combined", "mase_trend")
)
```

## Arguments

- mx_array:

  3D array of mortality rates.

- pop_array:

  3D array of population exposure.

- years:

  Integer vector.

- ages:

  Integer vector.

- groups:

  data.frame of group strata.

- h_val:

  Validation horizon.

- n_origins:

  Number of CV origins.

- width:

  Age-group width.

- cv_metric:

  Character; scoring metric for cross-validation. `"mase"` (default),
  `"mase_combined"`, or `"mase_trend"`. See
  [`resolve_cv_metric()`](https://kalleEF.github.io/mortensor/reference/resolve_cv_metric.md).

## Value

Named list of best parameters per model.
