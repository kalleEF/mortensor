# Interval Score

Computes the interval score (IS) for prediction intervals. The IS
penalises interval width and observations falling outside the interval.
Operates on the log scale for mortality rates.

## Usage

``` r
compute_interval_score(obs, lower, upper, alpha = 0.05)
```

## Arguments

- obs:

  Numeric vector of observed values (log scale recommended).

- lower:

  Numeric vector of lower prediction interval bounds.

- upper:

  Numeric vector of upper prediction interval bounds.

- alpha:

  Numeric; nominal miscoverage rate. Default 0.05 for a 95 percent PI.

## Value

Scalar mean interval score (lower is better).

## References

Gneiting, T. and Raftery, A.E. (2007). Strictly proper scoring rules,
prediction, and estimation. JASA 102(477), 359-378.
