# Group-mean correction diagnostic

Bar chart showing the ratio `mean(single-year) / grouped rate` per age
group, before and after correction. Both bars should be close to 1.0
after correction.

## Usage

``` r
plot_disagg_correction(
  disagg,
  mx_forecast,
  ages,
  fc_year_idx = 1L,
  group_idx = 1L,
  output_dir,
  cause_short,
  max_age = 100
)
```

## Arguments

- disagg:

  Return value from
  [`disaggregate_forecast_gam`](https://kalleEF.github.io/mortensor/reference/disaggregate_forecast_gam.md).

- mx_forecast:

  3D array `[n_a x h x n_g]` of grouped forecast rates.

- ages:

  Integer vector of age-group starts.

- fc_year_idx:

  Integer forecast year index (default 1).

- group_idx:

  Integer group index (default 1).

- output_dir:

  Output directory.

- cause_short:

  Short cause label.

- max_age:

  Maximum age (default 100).
