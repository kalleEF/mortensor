# Compare coherent Tucker4 vs independent per-disease Tucker3

Compare coherent Tucker4 vs independent per-disease Tucker3

## Usage

``` r
compare_coherent_vs_independent(
  X_ilr,
  mx_all,
  mx_disease_list,
  mx_disease_smooth_list,
  pop_all,
  disease_order,
  disease_configs,
  years,
  ages,
  groups,
  train_idx,
  test_idx,
  calib_tucker4,
  width = 5
)
```

## Arguments

- X_ilr:

  4D ILR array.

- mx_all:

  3D all-cause array.

- mx_disease_list:

  Named list of disease arrays (raw).

- mx_disease_smooth_list:

  Named list of smoothed disease arrays.

- pop_all:

  3D population array.

- disease_order:

  Character vector.

- disease_configs:

  Named list of mort_config objects.

- years, ages, groups:

  Study dimensions.

- train_idx, test_idx:

  Integer index vectors.

- calib_tucker4:

  Output of
  [`calibrate_tucker4()`](https://kalleEF.github.io/mortensor/reference/calibrate_tucker4.md).

- width:

  Age-group width.

## Value

data.table with per-disease coherent vs independent MASE.
