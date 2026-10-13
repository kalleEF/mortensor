# Per-age ARIMA and ETS models

Benchmarks that model the log death rate of every age as a separate time
series with forecast: `arima_series()` selects an ARIMA model with
[`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html)
(drift allowed), `ets_series()` an exponential smoothing state space
model without seasonality with
[`forecast::ets()`](https://pkg.robjhyndman.com/forecast/reference/ets.html)
(damped trends allowed). Draws are simulated future paths of the fitted
models (`simulate(..., future = TRUE)`), so uncertainty grows with the
horizon. Ages are simulated independently; parameter uncertainty is not
included.

## Usage

``` r
arima_series(..., bootstrap = FALSE)

ets_series(..., bootstrap = FALSE)
```

## Arguments

- ...:

  Arguments passed to
  [`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html)
  or
  [`forecast::ets()`](https://pkg.robjhyndman.com/forecast/reference/ets.html).

- bootstrap:

  Simulate with resampled residuals instead of normal errors.

## Value

A model specification.

## References

Hyndman RJ, Khandakar Y (2008). Automatic time series forecasting: the
forecast package for R. *Journal of Statistical Software* 27(3): 1-22.
[doi:10.18637/jss.v027.i03](https://doi.org/10.18637/jss.v027.i03)

## Examples

``` r
if (requireNamespace("forecast", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 2000:2019)
  fc <- mt_forecast(mt_fit(ets_series(), x), h = 3, n_draws = 20, seed = 1)
}
```
