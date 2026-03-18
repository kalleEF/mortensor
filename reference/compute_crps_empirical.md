# Empirical CRPS from bootstrap simulations (log scale)

Computes CRPS using the energy-score formulation on the empirical CDF of
bootstrap samples.

## Usage

``` r
compute_crps_empirical(obs_arr, sims_arr)
```

## Arguments

- obs_arr:

  3D array `[n_a x n_t x n_g]` of observed rates.

- sims_arr:

  4D array `[n_a x n_t x n_g x n_sim]` of simulated rates.

## Value

Scalar mean CRPS.
