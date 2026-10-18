# Per-age ARIMA and ETS models

Benchmarks that model the log death rate of every age as a separate time
series with forecast: `arima_series()` selects an ARIMA model with
[`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html)
(drift allowed) or fits a given order with
[`forecast::Arima()`](https://pkg.robjhyndman.com/forecast/reference/Arima.html)
(`order`), `ets_series()` an exponential smoothing state space model
without seasonality with
[`forecast::ets()`](https://pkg.robjhyndman.com/forecast/reference/ets.html)
(damped trends allowed). Draws are simulated future paths of the fitted
models (`simulate(..., future = TRUE)`), so uncertainty grows with the
horizon. Ages are simulated independently; parameter uncertainty is not
included. The models describe observed rates, so predictive draws add
nothing.

## Usage

``` r
arima_series(..., order = NULL, drift = NULL, bootstrap = FALSE)

ets_series(..., bootstrap = FALSE)
```

## Arguments

- ...:

  Arguments passed to
  [`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html),
  [`forecast::Arima()`](https://pkg.robjhyndman.com/forecast/reference/Arima.html)
  (with `order`) or
  [`forecast::ets()`](https://pkg.robjhyndman.com/forecast/reference/ets.html).

- order:

  ARIMA order `c(p, d, q)` fitted to every age with
  [`forecast::Arima()`](https://pkg.robjhyndman.com/forecast/reference/Arima.html);
  `NULL` (default) selects the order with
  [`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html).

- drift:

  With `order`: include a drift (a constant in the differenced series,
  `include.drift` of
  [`forecast::Arima()`](https://pkg.robjhyndman.com/forecast/reference/Arima.html)).
  Without `order`: whether
  [`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html)
  may select a drift (`allowdrift`; `NULL` keeps its default).

- bootstrap:

  Simulate with resampled residuals instead of normal errors.

## Value

A model specification.

## Details

`arima_series(order = c(0, 1, 1), drift = TRUE)` is the local level
model with drift (a random walk with drift observed with noise; Harvey
1989) in its reduced form: the year-to-year changes of a random walk
with drift plus independent noise are an MA(1) process with mean equal
to the drift. The parameters are estimated by maximum likelihood
(`method = "ML"` of
[`forecast::Arima()`](https://pkg.robjhyndman.com/forecast/reference/Arima.html),
unless `...` gives another `method`). Where deaths are few, it separates
the noise of the observed rates from the trend, which the plain random
walk ([`rwd()`](rwd.md)) accumulates.

## References

Hyndman RJ, Khandakar Y (2008). Automatic time series forecasting: the
forecast package for R. *Journal of Statistical Software* 27(3): 1-22.
[doi:10.18637/jss.v027.i03](https://doi.org/10.18637/jss.v027.i03)

Harvey AC (1989). *Forecasting, Structural Time Series Models and the
Kalman Filter*. Cambridge University Press.

## Examples

``` r
if (requireNamespace("forecast", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 2000:2019)
  fc <- mt_forecast(mt_fit(ets_series(), x), h = 3, n_draws = 20, seed = 1)
  # local level model with drift (random walk with drift plus noise)
  fc2 <- mt_forecast(mt_fit(arima_series(order = c(0, 1, 1), drift = TRUE), x), h = 3, n_draws = 20,
                     seed = 1)
}
```
