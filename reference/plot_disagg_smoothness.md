# Disaggregation smoothness heatmap

Displays \\\|d^2 \log(mx) / d(\text{age})^2\|\\ across ages and groups
as a heatmap. Bright spots indicate kinks or jumps.

## Usage

``` r
plot_disagg_smoothness(mx_1yr_dt, output_dir, cause_short)
```

## Arguments

- mx_1yr_dt:

  `data.table` with columns `age`, `mx_1yr`, `year`, and any group
  columns.

- output_dir:

  Output directory.

- cause_short:

  Short cause label for filename.
