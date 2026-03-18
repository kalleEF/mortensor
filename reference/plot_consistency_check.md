# Consistency check boxplot (sum of diseases / all-cause)

Consistency check boxplot (sum of diseases / all-cause)

## Usage

``` r
plot_consistency_check(
  fc_dt_combined,
  ac_fc_dt,
  output_dir,
  group_cols = character(0)
)
```

## Arguments

- fc_dt_combined:

  Combined forecast data.table (all causes).

- ac_fc_dt:

  All-cause forecast data.table.

- output_dir:

  Output directory.

- group_cols:

  Character vector of group column names.
