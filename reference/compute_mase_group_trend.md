# Group-trend MASE

For each group, averages mx across ages first, then computes MASE on the
group-level trend. Measures whether the model captures group-specific
temporal dynamics.

## Usage

``` r
compute_mase_group_trend(arr_test, fc_arr, arr_train)
```

## Arguments

- arr_test:

  3D array of observed test-period rates.

- fc_arr:

  3D array of forecast rates.

- arr_train:

  3D array of training-period rates.

## Value

Scalar mean group-trend MASE.
