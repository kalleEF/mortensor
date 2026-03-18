# Temporal consistency diagnostic

Shows mortality rates at selected ages across all years (observed +
forecast) for one population group. Useful for verifying that the GAM
approach produces smooth temporal evolution at individual single-year
ages, including those near the blend boundary.

## Usage

``` r
plot_disagg_temporal(
  disagg,
  years,
  forecast_years,
  group_idx = 1L,
  output_dir,
  cause_short,
  group_label = "",
  probe_ages = c(50, 70, 85, 90, 95, 100)
)
```

## Arguments

- disagg:

  Return value from
  [`disaggregate_forecast_gam`](https://kalleEF.github.io/mortensor/reference/disaggregate_forecast_gam.md).

- years:

  Integer vector of observed years.

- forecast_years:

  Integer vector of forecast years.

- group_idx:

  Integer group index (default 1).

- output_dir:

  Output directory.

- cause_short:

  Short cause label.

- group_label:

  Character label for subtitle.

- probe_ages:

  Ages to plot (default: 50, 70, 85, 90, 95, 100).
