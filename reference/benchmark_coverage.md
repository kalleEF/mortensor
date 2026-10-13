# Interval calibration benchmark on simulated data

Fits `spec` on 30 years of [`simulate_lc_dgp()`](simulate_lc_dgp.md)
data and scores predictive draws for the next `h` years, over `n_rep`
replicates.

## Usage

``` r
benchmark_coverage(
  spec,
  n_rep = 10,
  h = 10,
  n_draws = 300,
  n_train = 30,
  prep = NULL
)
```

## Arguments

- spec:

  Model specification.

- n_rep:

  Replicates.

- h:

  Horizon.

- n_draws:

  Draws per forecast.

- n_train:

  Training years.

## Value

A `data.table` with coverage by horizon and a summary list in attribute
`"summary"` (`cov95_mean`, `cov95_min`, `cov80_mean`, `cov50_mean`).
