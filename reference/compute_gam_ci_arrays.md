# Compute pointwise Bayesian 95% credible intervals from GAM posteriors

Uses `predict(gam_fit, se.fit = TRUE)` for efficient closed-form CIs.

## Usage

``` r
compute_gam_ci_arrays(gam_fits, ages, years)
```

## Arguments

- gam_fits:

  List of [mgcv::gam](https://rdrr.io/pkg/mgcv/man/gam.html) objects
  (one per group).

- ages:

  Numeric vector of ages.

- years:

  Integer vector of years.

## Value

List with `mx_lo95` and `mx_hi95` (3D arrays \[n_a x n_t x n_g\]).
