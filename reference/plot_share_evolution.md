# Stacked area chart of cause shares (coherent mode)

Stacked area chart of cause shares (coherent mode)

## Usage

``` r
plot_share_evolution(
  shares_obs_4d,
  shares_fc_4d,
  years_obs,
  years_fc,
  disease_order,
  output_dir,
  group_index = 1,
  groups,
  sex_filter = NULL,
  quintile_filter = NULL
)
```

## Arguments

- shares_obs_4d:

  4D observed shares array.

- shares_fc_4d:

  4D forecast shares array.

- years_obs, years_fc:

  Year vectors.

- disease_order:

  Character vector of disease names.

- output_dir:

  Output directory.

- group_index:

  Integer index into `groups` (default 1).

- groups:

  data.frame of groups.

- sex_filter, quintile_filter:

  Deprecated strata filters (backward compat).
