# Create a mortality calibration object

Create a mortality calibration object

## Usage

``` r
mort_calibration(
  data,
  data_raw = NULL,
  models,
  comparison,
  best_model,
  ensemble_names = character(0),
  ensemble_weights = NULL,
  calib_coda = NULL,
  calib_direct = NULL,
  calib_all = NULL,
  config = NULL
)
```

## Arguments

- data:

  A `mort_data` object (smoothed training data).

- data_raw:

  A `mort_data` object (unsmoothed raw data), or `NULL`.

- models:

  Character vector of model names evaluated.

- comparison:

  data.table of model comparison results.

- best_model:

  Character: name of the best model.

- ensemble_names:

  Character vector of ensemble component names.

- ensemble_weights:

  Named numeric vector of ensemble weights.

- calib_coda:

  Tucker3-CoDa calibration result.

- calib_direct:

  Tucker3-Direct calibration result.

- calib_all:

  List of calibration results for other models.

- config:

  A `mort_config` object (optional pipeline settings).

## Value

A `mort_calibration` object.
