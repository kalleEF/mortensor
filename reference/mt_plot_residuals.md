# Residual heat map of a fitted model

Residual heat map of a fitted model

## Usage

``` r
mt_plot_residuals(fit, strata = NULL, type = c("deviance", "pearson"))
```

## Arguments

- fit:

  An `mt_fit` object.

- strata:

  Labels of the strata to show (default: up to four).

- type:

  `"deviance"` or `"pearson"`.

## Value

A ggplot object: residuals by age and period, one panel per stratum.
Diagonal bands indicate cohort effects, vertical bands period shocks.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  x <- mt_example(deaths = "single")
  mt_plot_residuals(mt_fit(lee_carter(), x))
}
```
