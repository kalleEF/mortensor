# Temporally-consistent GAM-based disaggregation with logistic blend

Fits a 2D GAM on observed grouped data (age-midpoint x year), predicts
single-year ages for all target ages, blends with a parametric tail
(Kannisto/Gompertz) via logistic weights in log-space, then applies
group-mean correction.

## Usage

``` r
disaggregate_forecast_gam(
  mx_grouped,
  ages,
  years,
  mx_forecast,
  forecast_years,
  config
)
```

## Arguments

- mx_grouped:

  3D array `[n_a x n_t x n_g]` of grouped mx (observed).

- ages:

  Integer vector of age-group starts.

- years:

  Integer vector of observed calendar years.

- mx_forecast:

  3D array `[n_a x h x n_g]` of grouped mx (forecast).

- forecast_years:

  Integer vector of forecast years.

- config:

  Config list with elements `age_width`, `open_interval_midpoint`,
  `extrapolation_method`, `covid_years`, `output_max_age`,
  `blend_width`.

## Value

Named list with:

- mx_obs_1yr:

  Array `[n_target x n_t x n_g]`

- mx_fc_1yr:

  Array `[n_target x h x n_g]`

- target_ages:

  Integer vector of single-year ages

- gam_fits:

  List of fitted GAM objects

- diag:

  List of diagnostic arrays
