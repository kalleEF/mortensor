# Forecast age profile for a single year

Forecast age profile for a single year

## Usage

``` r
plot_forecast_age(
  fc_dt,
  output_dir,
  cause_label,
  cause_short,
  age_groups,
  group_cols = character(0)
)
```

## Arguments

- fc_dt:

  data.table of forecasts.

- output_dir, cause_label, cause_short:

  Standard plot params.

- age_groups:

  Ordered character vector of age-group labels.

- group_cols:

  Character vector of group column names.
