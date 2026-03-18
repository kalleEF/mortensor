# Compare all candidate models on a train/test split

Runs Tucker3-CoDa, Tucker3-Direct, Lee-Carter, FDM, StMoMo-APC, RWD,
ARIMA-perSeries, and ETS-Damped on the training data; evaluates on the
holdout test set. Optionally builds an ensemble via inverse-MASE
weighting and NNLS stacking.

## Usage

``` r
compare_models(
  mx_array,
  pop_array,
  years,
  ages,
  groups,
  train_idx,
  test_idx,
  calib_coda,
  calib_direct,
  width = 5,
  calib_all = NULL,
  mx_eval_ref = NULL
)
```

## Arguments

- mx_array:

  3D array \[n_a x n_t x n_g\] of (smoothed) mortality rates.

- pop_array:

  3D population exposure array.

- years:

  Integer vector of years.

- ages:

  Integer vector of age-group starts.

- groups:

  data.frame with `sex` and `quintile` columns.

- train_idx:

  Integer indices of training years.

- test_idx:

  Integer indices of test years.

- calib_coda:

  Output of
  [`calibrate_tucker3()`](https://kalleEF.github.io/mortensor/reference/calibrate_tucker3.md)
  with `method="coda"`.

- calib_direct:

  Output of
  [`calibrate_tucker3()`](https://kalleEF.github.io/mortensor/reference/calibrate_tucker3.md)
  with `method="direct"`.

- width:

  Age-group width (default 5).

- calib_all:

  Output of
  [`calibrate_all_models()`](https://kalleEF.github.io/mortensor/reference/calibrate_all_models.md)
  (NULL for defaults).

- mx_eval_ref:

  Optional raw-scale array for MASE evaluation.

## Value

Named list: `results` (data.table), `fc_arrays`, `ens_best_k`,
`ens_best_names`, `ens_best_weights`.
