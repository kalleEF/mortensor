# Bootstrap PIs for Tucker4 coherent forecast

Bootstrap PIs for Tucker4 coherent forecast

## Usage

``` r
bootstrap_coherent_forecast(fc_obj, mx_all_fc, n_sim, h)
```

## Arguments

- fc_obj:

  Tucker4 forecast object (return_full = TRUE).

- mx_all_fc:

  3D array \[n_a x h x n_g\] of all-cause point forecast.

- n_sim:

  Number of simulations.

- h:

  Forecast horizon.

## Value

List with 4D arrays: `mx_lower`, `mx_upper`, `mx_median` \[n_a x h x n_g
x K\].
