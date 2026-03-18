# Cross-validate Tucker3 with rolling origin

Cross-validate Tucker3 with rolling origin

## Usage

``` r
cv_rolling_origin(
  mx_array,
  years,
  ages,
  R_A,
  R_B,
  R_C,
  ts_method,
  method = "coda",
  h_val = 3,
  min_train = 12,
  width = 5,
  metric_fn = compute_mase_array
)
```

## Arguments

- mx_array:

  3D array of mortality rates.

- years:

  Integer vector of years.

- ages:

  Integer vector of age-group starts.

- R_A, R_B, R_C:

  Tucker3 ranks.

- ts_method:

  Time-series method.

- method:

  `"coda"` or `"direct"`.

- h_val:

  Validation horizon.

- min_train:

  Minimum training years.

- width:

  Age-group width.

- metric_fn:

  Scoring function with signature
  `function(arr_test, fc_arr, arr_train)`. Default:
  [`compute_mase_array()`](https://kalleEF.github.io/mortensor/reference/compute_mase_array.md).

## Value

Scalar mean CV score.
