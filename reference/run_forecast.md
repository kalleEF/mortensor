# Run the full mortality forecast pipeline for a single cause

Executes 17 stages: data loading, array building, smoothing, Tucker3
calibration, model calibration, holdout comparison, model selection /
ensemble, production forecast, bootstrap PIs, PI re-centering, coherence
enforcement, fitted values, plots, disaggregation, diagnostics, and
output saving.

## Usage

``` r
run_forecast(data, config = NULL)
```

## Arguments

- data:

  Either a
  [mort_data](https://kalleEF.github.io/mortensor/reference/mort_data.md)
  object (generic path) or a
  [mort_config](https://kalleEF.github.io/mortensor/reference/mort_config.md)
  object (ONS path, backward compatible).

- config:

  Optional
  [mort_config](https://kalleEF.github.io/mortensor/reference/mort_config.md)
  object providing pipeline settings (forecast horizon, number of
  bootstrap simulations, etc.). Required when `data` is a `mort_data`
  object. Ignored when `data` is itself a `mort_config`.

## Value

Invisibly, a list with all results (data, forecasts, PIs, tables).
