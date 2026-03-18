# Disaggregation profile diagnostic

For a single group and forecast year, overlays grouped rates as
horizontal step bars, GAM raw, parametric, blended (pre-correction), and
final corrected curves. A dashed vertical line marks the blend centre.

## Usage

``` r
plot_disagg_profile(
  disagg,
  mx_grouped,
  ages,
  width = 5,
  fc_year_idx = 1L,
  group_idx = 1L,
  open_mid = 95,
  output_dir,
  cause_short,
  group_label = ""
)
```

## Arguments

- disagg:

  Return value from
  [`disaggregate_forecast_gam`](https://kalleEF.github.io/mortensor/reference/disaggregate_forecast_gam.md).

- mx_grouped:

  Numeric vector of grouped rates for this group/year.

- ages:

  Integer vector of age-group starts.

- width:

  Age-group width (default 5).

- fc_year_idx:

  Integer index into forecast years (default 1).

- group_idx:

  Integer group index (default 1).

- open_mid:

  Blend centre age (default 95).

- output_dir:

  Output directory.

- cause_short:

  Short cause label for filename.

- group_label:

  Character label for subtitle.
