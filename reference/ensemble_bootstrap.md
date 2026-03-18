# MASE-weighted ensemble bootstrap

Combines per-model simulation arrays via weighted averaging.

## Usage

``` r
ensemble_bootstrap(model_names, weights, sim_arrays, n_a, h, n_g)
```

## Arguments

- model_names:

  Character vector of ensemble models.

- weights:

  Numeric weight vector (sum to 1).

- sim_arrays:

  Named list of 4D arrays \[n_a x h x n_g x n_sim\].

- n_a, h, n_g:

  Dimension sizes.

## Value

List with `mx_lower`, `mx_upper`, `mx_median`, `mx_sims`.
