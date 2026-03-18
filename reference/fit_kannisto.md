# Fit Kannisto logistic model for old-age mortality extrapolation

Fits mu(x) = a*exp(b*x) / (1 + a*exp(b*x)) to grouped mortality rates
via least-squares on the log scale.

## Usage

``` r
fit_kannisto(mx, ages)
```

## Arguments

- mx:

  Numeric vector of grouped mortality rates.

- ages:

  Numeric vector of age-group midpoints.

## Value

Named list with `a`, `b`, and `predict` function.
