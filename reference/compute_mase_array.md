# Cell-level MASE across a 3D array

Computes
[`compute_mase()`](https://kalleEF.github.io/mortensor/reference/compute_mase.md)
for each (age, group) series and returns the mean.

## Usage

``` r
compute_mase_array(arr_test, fc_arr, arr_train)
```

## Arguments

- arr_test:

  3D array of observed test-period rates.

- fc_arr:

  3D array of forecast rates (same dims as test).

- arr_train:

  3D array of training-period rates (same ages/groups).

## Value

Scalar mean MASE.
