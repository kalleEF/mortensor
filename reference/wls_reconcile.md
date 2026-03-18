# WLS reconciliation of independent disease forecasts

Adjusts independently-forecast disease components so they sum to the
all-cause forecast. Adjustment is proportional to each component's
forecast variance (wide PIs absorb more gap).

## Usage

``` r
wls_reconcile(mx_all_fc, mx_indep, bs_indep, disease_order, recon_cap = 0.5)
```

## Arguments

- mx_all_fc:

  3D array \[n_a x h x n_g\] of all-cause point forecast.

- mx_indep:

  Named list of 3D arrays: per-disease point forecasts.

- bs_indep:

  Named list of bootstrap results (each with `mx_lower`, `mx_upper`,
  `mx_median`).

- disease_order:

  Character vector of disease names.

- recon_cap:

  Soft cap: max \|delta_k\| as fraction of point forecast.

## Value

Named list with `mx_recon`, `bs_recon`, `delta`, `sigma_sq`.
