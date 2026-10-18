# Random walk with drift on log death rates

Benchmark model: the log death rate of every age follows a random walk
with drift, `y[t] = y[t-1] + mu + e[t]`. The drift is estimated from the
last `drift_window` years. Draws include drift uncertainty and
innovations:

## Usage

``` r
rwd(
  drift_window = NULL,
  innovations = c("bootstrap", "gaussian"),
  observation_noise = FALSE
)
```

## Arguments

- drift_window:

  Number of most recent years used for drift and innovations (`NULL`:
  all years).

- innovations:

  `"bootstrap"` or `"gaussian"`.

- observation_noise:

  Deprecated: model the log rates as a random walk plus noise (see
  Details); use `arima_series(order = c(0, 1, 1), drift = TRUE)`.

## Value

A model specification.

## Details

- `"bootstrap"` (default): innovations are age vectors resampled from
  the observed year-to-year changes, which keeps the correlation between
  ages; the drift of each path is the mean of a bootstrap sample of
  those changes;

- `"gaussian"`: independent normal innovations per age with the
  estimated variance, and normal drift uncertainty.

Zero death counts are replaced by 0.5 before taking logs. The draws
describe observed rates, so predictive draws add nothing.

**Deprecated: `observation_noise = TRUE`.** This variant treats the log
rate as a random walk with drift observed with noise (a local level
model with drift), with variances estimated by a moment estimator of
this package without a published source: from the variance and first
autocovariance of the year-to-year changes, the noise variance as the
Poisson part (one over the deaths) plus an extra part common to all
ages, clamped to at most 0.45 times the variance of the changes, the
level variance to at least 0.1 times it, and shrunk towards the median
of the ages with a prior standard deviation of three times that median;
the level is Kalman filtered and the draws are paths of the level, to
which predictive draws add Poisson variation of the deaths. It still
works and warns (class `mt_warning_deprecated`). The published
equivalent is the local level model with drift fitted by maximum
likelihood (Harvey 1989), whose reduced form for the observed log rates
is ARIMA(0,1,1) with drift:
`arima_series(order = c(0, 1, 1), drift = TRUE)`. With a `drift_window`,
the variant is fitted to the years of the window and fitted rates of
earlier years are missing.

## References

Harvey AC (1989). *Forecasting, Structural Time Series Models and the
Kalman Filter*. Cambridge University Press.

## See also

[`arima_series()`](arima_series.md).

## Examples

``` r
x <- mt_example(deaths = "single")
fc <- mt_forecast(mt_fit(rwd(drift_window = 15), x), h = 5, n_draws = 100, seed = 1)
```
