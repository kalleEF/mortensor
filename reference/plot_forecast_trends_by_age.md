# Forecast trends faceted by age group

Forecast trends faceted by age group

## Usage

``` r
plot_forecast_trends_by_age(
  obs_dt,
  fitted_dt,
  fc_dt,
  output_dir,
  cause_label,
  cause_short,
  forecast_years,
  group_cols = character(0),
  smooth_dt = NULL
)
```

## Arguments

- obs_dt:

  Observed data.table with `mx_obs`.

- fitted_dt:

  Fitted data.table with `mx_fitted`.

- fc_dt:

  Forecast data.table with `mx_adjusted`, `mx_lo95_adj`, `mx_hi95_adj`.

- output_dir:

  Output directory path.

- cause_label:

  Full cause name.

- cause_short:

  Short cause key.

- forecast_years:

  Integer vector of forecast years.

- group_cols:

  Character vector of group column names in the data.tables.

- smooth_dt:

  Optional data.table with `mx_lo95`, `mx_hi95` for Bayesian CI ribbon.
