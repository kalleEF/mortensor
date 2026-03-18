# Calibrate and compare mortality forecasting models

Runs cross-validated calibration for Tucker3 (CoDa and Direct) and
optionally other models (Lee-Carter, FDM, StMoMo, ARIMA, ETS, RWD), then
compares all on a holdout test set and builds an ensemble.

## Usage

``` r
calibrate(
  data,
  models = c("tucker3_coda", "tucker3_direct", "lee_carter", "fdm", "ets", "arima",
    "rwd"),
  holdout_years = NULL,
  holdout_frac = 0.2,
  ensemble = TRUE,
  ensemble_top_k = 3L,
  rank_selection = c("cv", "bic", "aic"),
  cv_metric = c("mase", "mase_combined", "mase_trend"),
  config = NULL,
  ...
)
```

## Arguments

- data:

  A `mort_data` object. If not yet smoothed, smoothing is applied
  automatically.

- models:

  Character vector of models to evaluate. Available: `"tucker3_coda"`,
  `"tucker3_direct"`, `"lee_carter"`, `"fdm"`, `"stmomo"`, `"arima"`,
  `"ets"`, `"rwd"`.

- holdout_years:

  Integer vector of years for holdout test. If `NULL`, the last
  `holdout_frac` of years is used.

- holdout_frac:

  Fraction of years for holdout (default 0.2). Ignored if
  `holdout_years` is provided.

- ensemble:

  Logical; build an inverse-MASE ensemble? (default TRUE)

- ensemble_top_k:

  Maximum number of models in the ensemble.

- rank_selection:

  Character; rank selection strategy for Tucker3 models. `"cv"`
  (default) uses rolling-origin cross-validation, `"bic"` or `"aic"`
  selects ranks via information criteria (faster). Passed to
  [`calibrate_tucker3()`](https://kalleEF.github.io/mortensor/reference/calibrate_tucker3.md).

- cv_metric:

  Character; scoring metric for cross-validation. `"mase"` (default)
  uses cell-level MASE, `"mase_combined"` uses 0.7 \* cell + 0.3 \*
  trend MASE, `"mase_trend"` uses trend MASE only. Applied to all models
  during calibration.

- config:

  Optional `mort_config` for advanced settings (smoothing method, COVID
  years, etc.).

- ...:

  Additional arguments passed to calibration functions.

## Value

A `mort_calibration` object.
