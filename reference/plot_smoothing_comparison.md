# Scatterplot of raw vs smoothed log(mx)

Scatterplot of raw vs smoothed log(mx)

## Usage

``` r
plot_smoothing_comparison(
  mx_raw,
  mx_smooth,
  years,
  ages,
  output_dir,
  cause_short,
  covid_years,
  floor_threshold = 1e-07
)
```

## Arguments

- mx_raw, mx_smooth:

  3D mortality arrays.

- years:

  Integer vector.

- ages:

  Integer vector.

- output_dir, cause_short:

  Standard plot params.

- covid_years:

  Years to exclude from comparison.

- floor_threshold:

  Cells below this are excluded.
