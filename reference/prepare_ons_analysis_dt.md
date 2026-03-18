# Prepare an ONS analysis data.table

Aggregates IMD deciles to quintiles, computes cause-specific and
all-cause death counts, merges with population, forward-fills missing
population years, and computes mortality rates.

## Usage

``` r
prepare_ons_analysis_dt(deaths, pop, config)
```

## Note

ONS-format specific.
