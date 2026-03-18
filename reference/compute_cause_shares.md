# Compute cause shares from disease-specific and all-cause mortality

Compute cause shares from disease-specific and all-cause mortality

## Usage

``` r
compute_cause_shares(mx_disease_list, mx_all, disease_order)
```

## Arguments

- mx_disease_list:

  Named list of 3D arrays, one per disease.

- mx_all:

  3D array of all-cause mortality.

- disease_order:

  Character vector of disease names.

## Value

4D array \[n_a x n_t x n_g x (K+1)\] of shares summing to 1.
