# Proper scores from forecast samples

CRPS (the empirical-CDF estimator, Grimit et al. 2006: the CRPS of the
empirical distribution of the samples, which depends on their number),
weighted interval score (Bracher et al. 2021, with the 23 quantiles of
the COVID-19 forecast hubs: 11 central intervals plus the median),
interval coverage, the PIT (randomised where samples equal the
observation, Czado, Gneiting and Held 2009) and the absolute error of
the median.

## Usage

``` r
score_samples(X, y, probs = wis_probs())
```

## Arguments

- X:

  Matrix of samples (one row per target).

- y:

  Observed values.

- probs:

  Quantile levels (must contain 0.5 and symmetric pairs).

## Value

A `data.table` with `median`, `mean`, `crps`, `wis`, `ae`, `cov50`,
`cov80`, `cov95`, `pit`.

## References

Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
forecasts in an interval format. *PLoS Computational Biology* 17(2):
e1008618.
[doi:10.1371/journal.pcbi.1008618](https://doi.org/10.1371/journal.pcbi.1008618)

Czado C, Gneiting T, Held L (2009). Predictive model assessment for
count data. *Biometrics* 65(4): 1254-1261.

Gneiting T, Raftery AE (2007). Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical
Association* 102(477): 359-378.
[doi:10.1198/016214506000001437](https://doi.org/10.1198/016214506000001437)

Grimit EP, Gneiting T, Berrocal VJ, Johnson NA (2006). The continuous
ranked probability score for circular variables and its application to
mesoscale forecast ensemble verification. *Quarterly Journal of the
Royal Meteorological Society* 132: 2925-2942.
