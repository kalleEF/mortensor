# Two-dimensional P-spline model

`p_splines()` smooths the log death rates over age and year with a
Poisson P-spline surface and forecasts them by extending the surface
over the forecast years (Currie, Durban and Eilers 2004). The code
follows the reference implementation, `Mort2Dsmooth()` and
`predict.Mort2Dsmooth()` of the R package MortalitySmooth 2.3.4
(Camarda), and reproduces its fitted values, forecasts and standard
errors for given smoothing parameters.

## Usage

``` r
p_splines(
  lambda = NULL,
  criterion = c("bic", "aic"),
  pord = 2,
  deg = 3,
  ndx = NULL,
  overdispersion = FALSE,
  lambda_range = c(1e-04, 1e+06),
  infant = FALSE
)
```

## Arguments

- lambda:

  Smoothing parameters over age and years, `c(age, year)`; `NULL`
  (default) selects them by `criterion`.

- criterion:

  `"bic"` (default) or `"aic"`.

- pord:

  Order of the difference penalties, over age and years (one value for
  both).

- deg:

  Degree of the B-splines, over age and years.

- ndx:

  Numbers of intervals between knots over age and years; `NULL` gives
  about one knot per five years (see Details).

- overdispersion:

  Account for overdispersion in the choice of the smoothing parameters
  and in the draws.

- lambda_range:

  Search range of the smoothing parameters (the default of
  MortalitySmooth). A selected value at an end of the range gives a
  warning (class `mt_warning_lambda_bound`) when the criterion improves
  by more than 1 two decades beyond it.

- infant:

  Give the first year of life its own coefficients (Camarda 2019).

## Value

A model specification.

## Details

**Model.** The deaths of age `x` and year `t` are Poisson with mean
`E[x, t] m[x, t]`. The log rates are a tensor-product B-spline surface,
`log m = Ba A t(By)`, with B-splines of degree `deg` on equally spaced
knots over age (`Ba`) and over years (`By`) and a matrix of coefficients
`A`. Differences of order `pord` between neighbouring coefficients are
penalised, over age with weight `lambda[1]` and over years with weight
`lambda[2]`. The coefficients maximise the penalised Poisson
log-likelihood (penalised scoring algorithm, computed in array form).
The smoothing parameters minimise the BIC (default) or the AIC, with the
effective dimension (the trace of the hat matrix) as the number of
parameters.

**Forecast.** The forecast is part of the fit. The basis over years is
extended by B-splines on the same knot spacing beyond the last forecast
year (and below the first year, as `predict.Mort2Dsmooth()` places the
knots), and the forecast years enter with weight 0: they carry no data,
and their coefficients follow from the penalty. With `pord = 2` the
coefficients continue about linearly, so each age continues the trend of
its last years; the order of the penalty determines the form of the
forecast. The model is fitted again on the extended basis with the
smoothing parameters of the fit, as `predict.Mort2Dsmooth()` does. The
penalty over age also acts on the added coefficients, so the fitted
rates of the data years change slightly, and the forecast depends a
little on the horizon through the last knot (up to 0.02 on the log scale
for England and Wales men, ages 30-100, 10 years ahead).

**Uncertainty.** Draws come from the normal approximation of the
coefficients, with the covariance `(B'WB + P)^-1` of the penalised
likelihood (`W` the fitted deaths, `P` the penalty). MortalitySmooth
computes the standard errors of forecast log rates from this covariance.
Over the forecast years the penalty acts as a prior: the variance grows
with the horizon, the faster the smaller `lambda[2]`. Draws are paths of
the whole surface, correlated over ages and years as in the model. They
describe the expected rates; predictive draws
(`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths.
With `overdispersion = TRUE`, the dispersion `psi2` (deviance divided by
the number of cells minus the effective dimension) enters the choice of
the smoothing parameters and the penalty, as in
`Mort2Dsmooth(overdispersion = TRUE)`, and the covariance is
`psi2 (B'WB + psi2 P)^-1`. MortalitySmooth's
[`predict()`](https://rdrr.io/r/stats/predict.html) instead refits the
forecast without `psi2` in the penalty and multiplies the standard
errors by `sqrt(psi2)`.

**Data.** Ages are placed at the midpoints of their cells, an open
interval at its lower bound plus half the width of the interval below it
(as in [`mt_smooth()`](mt_smooth.md)). Cells with missing deaths or
exposures, and cells without exposure, get weight 0, so the surface
interpolates them. `infant = TRUE` gives the first year of life a row of
coefficients of its own, penalised over years but not connected over age
with the following ages (Camarda 2019, section 3.1); it applies when the
grid starts with age group 0-1 and has at least four age groups.

**Package-specific choices (reference unknown).** By default the knots
are about five years of age and five calendar years apart,
`floor(n / 5)` intervals for `n` ages or years as in MortalitySmooth;
for age groups wider than one year, `n` is the number of single years
that the cell ages span. The smoothing parameters are searched on a grid
with steps of 1 on the log10 scale, followed by a continuous search near
the best grid point (MortalitySmooth: two grid searches with steps of 2
and 0.5).

## References

Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
mortality rates. *Statistical Modelling* 4(4): 279-298.
[doi:10.1191/1471082X04st080oa](https://doi.org/10.1191/1471082X04st080oa)

Camarda CG (2015). MortalitySmooth: Smoothing and forecasting Poisson
counts with P-splines. R package version 2.3.4 (CRAN archive).

Camarda CG (2019). Smooth constrained mortality forecasting.
*Demographic Research* 41(38): 1091-1130.
[doi:10.4054/DemRes.2019.41.38](https://doi.org/10.4054/DemRes.2019.41.38)

## See also

[`cp_splines()`](cp_splines.md) for forecasts constrained to past age
patterns and time trends; [`mt_smooth()`](mt_smooth.md) for smoothing
without forecasting.

## Examples

``` r
x <- mt_example(deaths = "single", years = 1995:2012)
fit <- mt_fit(p_splines(), x)
fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
mt_parameters(fit)
#>     stratum   parameter  index        value
#>      <char>      <char> <char>        <num>
#>  1:  female  lambda_age        8.337861e-01
#>  2:  female lambda_year        1.206932e+04
#>  3:  female          ed        3.398214e+01
#>  4:  female         bic        1.804454e+03
#>  5:  female  dispersion        1.255109e+00
#>  6:    male  lambda_age        9.801742e-01
#>  7:    male lambda_year        3.600860e+03
#>  8:    male          ed        3.669204e+01
#>  9:    male         bic        1.575542e+03
#> 10:    male  dispersion        1.057821e+00
```
