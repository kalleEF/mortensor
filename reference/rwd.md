# Random walk with drift on log death rates

Benchmark model: the log death rate of every age follows a random walk
with drift, `y[t] = y[t-1] + mu + e[t]`. The drift is estimated from the
last `drift_window` years. Draws include drift uncertainty and
innovations:

## Usage

``` r
rwd(drift_window = NULL, innovations = c("bootstrap", "gaussian"))
```

## Arguments

- drift_window:

  Number of most recent years used for drift and innovations (`NULL`:
  all years).

- innovations:

  `"bootstrap"` or `"gaussian"`.

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

A random walk with drift observed with noise (the local level model with
drift, Harvey 1989) has ARIMA(0,1,1) with drift as its reduced form for
the observed log rates: `arima_series(order = c(0, 1, 1), drift = TRUE)`
fits it by maximum likelihood.

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
