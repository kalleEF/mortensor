# Weighted Interval Score (WIS)

Computes the weighted interval score across multiple quantile levels. A
generalisation of the interval score that evaluates the full predictive
distribution via nested intervals.

## Usage

``` r
compute_wis(obs, quantile_matrix, probs)
```

## Arguments

- obs:

  Numeric vector of observed values.

- quantile_matrix:

  Numeric matrix `[n x K]` of quantile forecasts, where rows correspond
  to observations and columns to quantile levels.

- probs:

  Numeric vector of length K with the quantile probabilities (e.g.,
  `c(0.025, 0.25, 0.5, 0.75, 0.975)`).

## Value

Scalar mean WIS (lower is better).
