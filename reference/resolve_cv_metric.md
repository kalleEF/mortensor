# Resolve a CV metric name to its scoring function

Resolve a CV metric name to its scoring function

## Usage

``` r
resolve_cv_metric(metric = c("mase", "mase_combined", "mase_trend"))
```

## Arguments

- metric:

  Character; one of `"mase"`, `"mase_combined"`, `"mase_trend"`.

## Value

A function with signature `(arr_test, fc_arr, arr_train)`.
