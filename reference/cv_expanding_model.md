# Expanding-window cross-validation for a generic model

Expanding-window cross-validation for a generic model

## Usage

``` r
cv_expanding_model(
  mx_array,
  pop_array,
  years,
  ages,
  groups,
  fit_fn,
  h_val = 4,
  n_origins = 4,
  metric_fn = compute_mase_array,
  ...
)
```

## Arguments

- mx_array:

  3D array of mortality rates.

- pop_array:

  3D array of population exposure.

- years:

  Integer vector.

- ages:

  Integer vector.

- groups:

  data.frame of group strata.

- fit_fn:

  A function with signature
  `function(arr_train, pop_train, yrs_train, ages, groups, h, ...)`
  returning a list with element `mx_forecast`.

- h_val:

  Validation horizon.

- n_origins:

  Number of expanding origins.

- metric_fn:

  Scoring function with signature
  `function(arr_test, fc_arr, arr_train)`. Default:
  [`compute_mase_array()`](https://kalleEF.github.io/mortensor/reference/compute_mase_array.md).

- ...:

  Passed to `fit_fn`.

## Value

Scalar mean CV score.
