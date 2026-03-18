# Observed vs smoothed comparison plot

Observed vs smoothed comparison plot

## Usage

``` r
plot_obs_vs_smooth(
  obs_dt,
  smooth_dt,
  output_dir,
  cause_label,
  cause_short,
  group_cols = character(0)
)
```

## Arguments

- obs_dt:

  data.table with `mx_obs`.

- smooth_dt:

  data.table with `mx_smooth` (and optionally `mx_lo95`, `mx_hi95`).

- output_dir, cause_label, cause_short:

  Standard plot params.

- group_cols:

  Character vector of group column names.
