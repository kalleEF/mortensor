# Load and prepare ONS data for a single cause

Internal helper used by the pipeline functions when a `mort_config` with
ONS file paths is provided. Reads deaths + population Excel files, then
calls
[`prepare_ons_analysis_dt()`](https://kalleEF.github.io/mortensor/reference/prepare_ons_analysis_dt.md).

## Usage

``` r
load_and_prepare_ons(config)
```

## Arguments

- config:

  A
  [mort_config](https://kalleEF.github.io/mortensor/reference/mort_config.md)
  object with `deaths_file` and `pop_files` fields.

## Value

A `data.table` ready for array construction.

## Note

ONS-format specific.
