# Enforce sex-differential coherence on forecasts

Adjusts forecasted mortality rates to ensure Female \<= Male within each
age-group/year/quintile combination. Uses observed sex ratios as the
adjustment target.

## Usage

``` r
enforce_coherence(fc_dt, obs_dt, sex_constraint = TRUE)
```

## Arguments

- fc_dt:

  data.table with columns: `age_group`, `year`, `sex`, `quintile`,
  `mx_point`, `mx_lo95`, `mx_hi95`, `mx_median`.

- obs_dt:

  data.table with column `mx_obs` for observed ratios.

- sex_constraint:

  Logical; enforce Female \<= Male?

## Value

Modified `fc_dt` with added columns: `mx_adjusted`, `mx_lo95_adj`,
`mx_hi95_adj`, `mx_median_adj`.
