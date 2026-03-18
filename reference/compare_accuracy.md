# Compare accuracy across multiple forecasts

Takes several `mort_forecast` objects (or `mort_accuracy` objects) and
returns a summary comparison table.

## Usage

``` r
compare_accuracy(..., test = NULL)
```

## Arguments

- ...:

  Named `mort_forecast` or `mort_accuracy` objects.

- test:

  A `mort_data` object of observed test data. Required when passing
  `mort_forecast` objects. Ignored for `mort_accuracy` inputs.

## Value

A `data.table` with one row per model and columns for each metric.

## Examples

``` r
if (FALSE) { # \dontrun{
compare_accuracy(
  "Tucker3-CoDa" = fc_coda,
  "Lee-Carter"   = fc_lc,
  "Ensemble"     = fc_ens,
  test = holdout_data
)
} # }
```
