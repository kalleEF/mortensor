# Calibrate Tucker3 via grid search

Searches over rank combinations and time-series methods to find the best
Tucker3 configuration. Supports two selection strategies: rolling-origin
cross-validation (`"cv"`, slower but more robust) or information
criteria (`"bic"`/`"aic"`, fast screening).

## Usage

``` r
calibrate_tucker3(
  mx_array,
  years,
  ages,
  method = "coda",
  width = 5,
  h_val = 3,
  min_train = 12,
  rank_selection = c("cv", "bic", "aic"),
  cv_metric = c("mase", "mase_combined", "mase_trend")
)
```

## Arguments

- mx_array:

  3D array of mortality rates.

- years:

  Integer vector of years.

- ages:

  Integer vector of age-group starts.

- method:

  `"coda"` or `"direct"`.

- width:

  Age-group width.

- h_val:

  Validation horizon (used only for `rank_selection = "cv"`).

- min_train:

  Minimum training years (used only for `rank_selection = "cv"`).

- rank_selection:

  Character; `"cv"` (default) for rolling-origin cross-validation,
  `"bic"` or `"aic"` for information criterion based selection. When
  using IC, the ts_method is selected via a smaller CV using the best
  ranks.

- cv_metric:

  Character; scoring metric for cross-validation. `"mase"` (default),
  `"mase_combined"`, or `"mase_trend"`. See
  [`resolve_cv_metric()`](https://kalleEF.github.io/mortensor/reference/resolve_cv_metric.md).

## Value

Named list with `R_A`, `R_B`, `R_C`, `ts_method`, `cv_score` (or
`ic_value` for BIC/AIC selection).
