# Interval Score for 3D mortality arrays

Convenience wrapper that computes the mean interval score across all
cells of 3D mortality rate arrays on the log scale.

## Usage

``` r
compute_interval_score_array(obs_arr, lower_arr, upper_arr, alpha = 0.05)
```

## Arguments

- obs_arr:

  3D array `[n_a x n_t x n_g]` of observed rates.

- lower_arr:

  3D array of lower PI bounds (same dims).

- upper_arr:

  3D array of upper PI bounds (same dims).

- alpha:

  Nominal miscoverage rate. Default 0.05.

## Value

Scalar mean interval score.
