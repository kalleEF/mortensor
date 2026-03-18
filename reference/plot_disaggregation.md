# Plot disaggregated single-year mortality

Plot disaggregated single-year mortality

## Usage

``` r
plot_disaggregation(
  mx_1yr_dt,
  output_dir,
  cause_label,
  cause_short,
  group_cols = character(0)
)
```

## Arguments

- mx_1yr_dt:

  data.table with `age`, `mx_1yr`, `mx_1yr_lo95`, `mx_1yr_hi95`.

- output_dir, cause_label, cause_short:

  Standard plot params.

- group_cols:

  Character vector of group column names.
