# Functional demographic model of Hyndman and Ullah (demography)

Wrapper around
[`demography::fdm()`](https://pkg.robjhyndman.com/demography/reference/fdm.html)
(Hyndman and Ullah 2007): the rates of each year are smoothed over age
with
[`demography::smooth.demogdata()`](https://pkg.robjhyndman.com/demography/reference/smooth.demogdata.html)
(`smooth = TRUE`, by default monotonic regression splines above age 65,
weighted by the exposures), the log rates are decomposed into a mean
function and principal components, and the scores are forecast with
ARIMA or ETS models. Draws are simulated with `demography`'s
[`simulate()`](https://rdrr.io/r/stats/simulate.html) method for
functional forecasts (uncertainty of the coefficients and resampled
model errors). Predictive draws (`mt_forecast(predictive = TRUE)`) add
Poisson variation of the deaths given the rates and the exposures, as
for other models of expected rates, and no other observation noise. If
intervals are too narrow on real data, that is a property of the model
(see [`recalibrated()`](recalibrated.md), [`ensemble()`](ensemble.md),
or overdispersed models such as `stanmomo(family = "nb")`).

## Usage

``` r
hyndman_ullah(
  order = 6,
  method = "arima",
  smooth = TRUE,
  smooth_args = list(),
  ...
)
```

## Arguments

- order:

  Number of principal components.

- method:

  Forecast method for the coefficients (`"arima"`, `"ets"`, `"rwdrift"`,
  ...), passed to `forecast()`.

- smooth:

  Smooth the rates of each year over age before the decomposition.

- smooth_args:

  List of arguments of
  [`demography::smooth.demogdata()`](https://pkg.robjhyndman.com/demography/reference/smooth.demogdata.html)
  (for example `list(k = 20)` or `list(method = "spline")`).

- ...:

  Further arguments of
  [`demography::fdm()`](https://pkg.robjhyndman.com/demography/reference/fdm.html).

## Value

A model specification.

## Details

Without smoothing, the model errors that
[`simulate()`](https://rdrr.io/r/stats/simulate.html) resamples contain
the observation noise, which predictive draws then add a second time (as
Poisson variation): keep `smooth = TRUE` unless the rates were smoothed
in a preprocessing step. `smooth.demogdata()` works on single years of
age; for age groups the smoothed values at the observed ages are used (a
weaker smoothing). The smoothing settings (`smooth_args`, such as `k`,
the number of knots) and the number of components can be chosen by
forecast accuracy with [`tuned()`](mt_tune.md). When the default
monotonic spline cannot be fitted (for example for age ranges that start
in middle age), the rates are smoothed without the monotonicity
constraint and the fit diagnostics say so.

## References

Hyndman RJ, Ullah MS (2007). Robust forecasting of mortality and
fertility rates: a functional data approach. *Computational Statistics &
Data Analysis* 51(10): 4942-4956.

## Examples

``` r
if (requireNamespace("demography", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 1995:2019)
  fc <- mt_forecast(mt_fit(hyndman_ullah(order = 2), x), h = 5, n_draws = 50, seed = 1)
}
```
