# Bootstrap PIs for StMoMo-APC models

Uses StMoMo's built-in simulate() or falls back to residual bootstrap.

## Usage

``` r
bootstrap_stmomo(
  mx_array,
  pop_array,
  years,
  ages,
  groups_dt,
  h,
  n_sim,
  mx_raw = NULL
)
```

## Arguments

- mx_array:

  3D array of mortality rates.

- pop_array:

  3D population array.

- years:

  Integer vector.

- ages:

  Integer vector.

- groups_dt:

  data.frame of groups.

- h:

  Forecast horizon.

- n_sim:

  Number of simulations.

- mx_raw:

  Optional raw rates for residuals.

## Value

4D array \[n_a x h x n_g x n_sim\] or NULL.
