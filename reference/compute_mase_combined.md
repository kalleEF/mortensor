# Combined MASE (cell + trend)

Balances cell-level and group-trend accuracy with weights 0.7 / 0.3.

## Usage

``` r
compute_mase_combined(arr_test, fc_arr, arr_train)
```

## Arguments

- arr_test:

  3D array of observed test-period rates.

- fc_arr:

  3D array of forecast rates.

- arr_train:

  3D array of training-period rates.

## Value

Scalar combined MASE.
