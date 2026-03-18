# Rolling-origin CV for Tucker4 coherent forecast

Rolling-origin CV for Tucker4 coherent forecast

## Usage

``` r
cv_rolling_origin_tucker4(
  X_ilr,
  mx_all,
  mx_disease_list,
  disease_order,
  years,
  R_A,
  R_B,
  R_C,
  R_D,
  ts_method,
  h_val = 3,
  min_train = 12,
  metric_fn = compute_mase_array
)
```

## Arguments

- X_ilr:

  4D ILR array.

- mx_all:

  3D all-cause array.

- mx_disease_list:

  Named list of 3D disease arrays.

- disease_order:

  Character vector.

- years:

  Integer vector.

- R_A, R_B, R_C, R_D:

  Tucker4 ranks.

- ts_method:

  Time-series method.

- h_val:

  Validation horizon.

- min_train:

  Minimum training years.

- metric_fn:

  Scoring function with signature
  `function(arr_test, fc_arr, arr_train)`. Default:
  [`compute_mase_array()`](https://kalleEF.github.io/mortensor/reference/compute_mase_array.md).

## Value

Scalar mean CV score.
