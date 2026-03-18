# Calibrate Tucker4 via two-stage grid search

Stage 1: coarse grid. Stage 2: fine grid around best from stage 1.

## Usage

``` r
calibrate_tucker4(
  X_ilr,
  mx_all,
  mx_disease_list,
  disease_order,
  years,
  h_val = 3,
  min_train = 12,
  cv_metric = c("mase", "mase_combined", "mase_trend")
)
```

## Arguments

- X_ilr:

  4D ILR array.

- mx_all:

  3D all-cause array.

- mx_disease_list:

  Named list of disease arrays.

- disease_order:

  Disease name vector.

- years:

  Integer vector.

- h_val:

  Validation horizon.

- min_train:

  Minimum training years.

- cv_metric:

  Character; scoring metric for cross-validation. `"mase"` (default),
  `"mase_combined"`, or `"mase_trend"`. See
  [`resolve_cv_metric()`](https://kalleEF.github.io/mortensor/reference/resolve_cv_metric.md).

## Value

Named list of best ranks and `cv_score`.
