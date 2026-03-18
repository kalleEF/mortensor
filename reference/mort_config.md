# Create a validated configuration object

Create a validated configuration object

## Usage

``` r
mort_config(
  cause_cols = NULL,
  cause_label = "All causes",
  cause_short = "all",
  is_cause_specific = FALSE,
  ...
)
```

## Arguments

- cause_cols:

  Character vector of column names to sum for this cause. Set to `NULL`
  when providing data via
  [mort_data](https://kalleEF.github.io/mortensor/reference/mort_data.md)
  objects directly (i.e., not using the ONS loading path).

- cause_label:

  Human-readable cause label.

- cause_short:

  Short tag for filenames (must be unique across configs).

- is_cause_specific:

  Logical; if TRUE, cap against all-causes is applied.

- ...:

  Additional named parameters that override defaults. Common overrides
  include `forecast_h`, `n_sim`, `smoothing_method`, `disagg_method`,
  `pi_method`, `covid_years`, `holdout_train`, `holdout_test`,
  `output_dir`.

## Value

A `mort_config` object.
