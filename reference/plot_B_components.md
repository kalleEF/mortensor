# Plot B-mode time components

Plot B-mode time components

## Usage

``` r
plot_B_components(
  fc_obj,
  forecast_years,
  sim_years,
  output_dir,
  cause_label,
  cause_short
)
```

## Arguments

- fc_obj:

  Tucker3 forecast object (with `B` and `arima_fits`).

- forecast_years:

  Integer vector of future years.

- sim_years:

  Integer vector of observed years.

- output_dir, cause_label, cause_short:

  Standard plot params.
