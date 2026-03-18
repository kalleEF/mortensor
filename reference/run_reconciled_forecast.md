# Run the reconciled multi-disease forecast pipeline

Runs independent forecasts for each disease and all-causes, then
reconciles via variance-weighted (WLS) adjustment so that disease
forecasts sum to the all-cause forecast.

## Usage

``` r
run_reconciled_forecast(disease_configs, all_cause_config, shared)
```

## Arguments

- disease_configs:

  Named list of `mort_config` objects (one per cause).

- all_cause_config:

  `mort_config` for all-causes.

- shared:

  Named list of shared parameters (output_dir, forecast_h, n_sim,
  age_width, seed, etc.).

## Value

Invisibly, a list with all reconciled results.
