# Evaluate forecast accuracy

Computes a comprehensive set of scoring rules for a mortality forecast.
When called with observed test data, evaluates out-of-sample accuracy.
When called without test data, evaluates in-sample fit (fitted vs raw
training observations).

## Usage

``` r
accuracy(object, ...)

# S3 method for class 'mort_forecast'
accuracy(object, test = NULL, ...)
```

## Arguments

- object:

  A `mort_forecast` object.

- ...:

  Ignored.

- test:

  Optional `mort_data` object containing observed mortality rates for
  the forecast period. If `NULL`, evaluates in-sample accuracy using
  fitted values and raw training data (both must be present).

## Value

A named list with class `mort_accuracy`, containing:

- mase_cell:

  Cell-level MASE across all (age, group) series.

- mase_trend:

  Group-trend MASE.

- mase_combined:

  Combined MASE (0.7 cell + 0.3 trend).

- rmse_log:

  RMSE on the log scale.

- mae:

  Mean absolute error.

- coverage:

  Empirical coverage of 95 percent PI (if PIs present).

- interval_score:

  Mean interval score (if PIs present).

- crps:

  Empirical CRPS from simulations (if sims present).

- n_ages:

  Number of age groups evaluated.

- n_years:

  Number of years evaluated.

- n_groups:

  Number of population groups evaluated.

- type:

  `"out_of_sample"` or `"in_sample"`.

- model:

  Model name.

## Examples

``` r
if (FALSE) { # \dontrun{
# Out-of-sample: compare forecast against held-out observed data
accuracy(fc, test = holdout_data)

# In-sample: evaluate fitted vs raw training data
accuracy(fc)
} # }
```
