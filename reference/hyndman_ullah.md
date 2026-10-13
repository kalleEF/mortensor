# Functional demographic model of Hyndman and Ullah (demography)

Wrapper around
[`demography::fdm()`](https://pkg.robjhyndman.com/demography/reference/fdm.html)
(Hyndman and Ullah 2007): log death rates are smoothed over age and
decomposed into principal components whose scores are forecast with
ARIMA or ETS models. Draws are simulated with `demography`'s
[`simulate()`](https://rdrr.io/r/stats/simulate.html) method for
functional forecasts (uncertainty of the coefficients and model error);
predictive draws add the extra-Poisson and Poisson variation as for
other models of smooth rates.

## Usage

``` r
hyndman_ullah(order = 6, method = "arima", ...)
```

## Arguments

- order:

  Number of principal components.

- method:

  Forecast method for the coefficients (`"arima"`, `"ets"`, `"rwdrift"`,
  ...), passed to `forecast()`.

- ...:

  Further arguments of
  [`demography::fdm()`](https://pkg.robjhyndman.com/demography/reference/fdm.html).

## Value

A model specification.

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
