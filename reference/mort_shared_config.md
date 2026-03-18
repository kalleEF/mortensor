# Create a shared configuration for multi-disease pipelines

Provides common parameters used by
[`run_reconciled_forecast()`](https://kalleEF.github.io/mortensor/reference/run_reconciled_forecast.md)
and
[`run_coherent_forecast()`](https://kalleEF.github.io/mortensor/reference/run_coherent_forecast.md).

## Usage

``` r
mort_shared_config(
  output_dir = ".",
  forecast_h = 10L,
  n_sim = 500L,
  age_width = 5L,
  seed = 42L,
  target_ages = 30L:100L,
  sex_constraint = FALSE,
  ...
)
```

## Arguments

- output_dir:

  Output directory for CSV and plot files.

- forecast_h:

  Forecast horizon in years.

- n_sim:

  Number of bootstrap simulations.

- age_width:

  Width of closed age groups.

- seed:

  Random seed.

- target_ages:

  Integer vector of single-year target ages for disaggregation.

- sex_constraint:

  Logical; enforce Female \<= Male?

- ...:

  Additional named parameters (e.g. `holdout_train`, `holdout_test`,
  `cv_years`, `cv_h`, `cv_min_train`, `age_groups`).

## Value

A `mort_shared_config` object (named list with class).
