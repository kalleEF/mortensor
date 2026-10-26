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

  `FALSE` (default): Poisson deaths. `TRUE`: one dispersion for all
  cells in the choice of the smoothing parameters and in the draws
  (MortalitySmooth). `"age"`: a dispersion for each age that varies
  smoothly over age (Djeundje and Currie 2011; see Details).

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
errors by `sqrt(psi2)`. For England and Wales men (ages 30-100,
1961-2011), the two forecasts differ by up to 0.08 on the log scale ten
years ahead.

**Dispersion by age.** `overdispersion = "age"` replaces the Poisson
assumption by the first two moments of the deaths, `E[D] = E m` and
`var(D) = phi[x] E[D]`, with a dispersion `phi[x]` for each age that
varies smoothly over age (Djeundje and Currie 2011, eq. 2.1; Model 3,
eq. 5.26). Rates and dispersions are estimated by the bias corrected
scheme of the paper (section 4):

1.  Start with `phi = 1` at every age.

2.  Fit the rates by penalised quasi-likelihood: the scoring algorithm
    with the weights of the cells divided by their dispersion (eq.
    3.13). The smoothing parameters minimise the criterion of the scaled
    deviance: the deviance of each cell divided by its dispersion, plus
    `log(n)` (BIC, eq. 3.16) or 2 (AIC) times the effective dimension.

3.  For each age, compute the bias-corrected mean Pearson statistic
    `d[x]`: the sum of `(D - mu)^2 / mu` over its years, divided by its
    number of cells `n[x]` minus its contribution to the effective
    dimension, the sum of its diagonal elements of the hat matrix (eq.
    4.20 and 4.21).

4.  Smooth the statistics over age by penalised quasi-likelihood with
    `E[d] = phi`, `var(d) = tau phi^2 / n[x]` and the link
    `-1 / phi = Ba beta` (eq. 4.22): the B-splines over age of the rates
    and a difference penalty of order `pord[1]` (eq. 5.28). The
    smoothing parameter minimises the scaled BIC of this fit.

5.  Repeat steps 2 to 4 until no dispersion changes by more than 1e-5
    relative to its value (section 5.2.2), at most 50 times.

The forecast fit keeps the dispersions of the ages in the weights, and
the draws come from the covariance `(B'WB + P)^-1` with `W` the fitted
deaths divided by their dispersion. In the paper's example (CMI assured
lives, ages 40-90, section 5.2.3) the dispersion down-weights the fit to
the data and the selected surface is smoother: the effective dimension
is 46 under the Poisson assumption, 36 with one dispersion and 33 with
the dispersion by age, and the confidence intervals of the forecast are
narrower. With a single age, `"age"` gives one dispersion, the
bias-corrected mean Pearson statistic (eq. 4.20). The fit selects the
smoothing parameters in every round, so it takes several times as long
as the Poisson fit.

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
and 0.5). With `overdispersion = "age"`: `tau` is the mean Pearson
statistic of the smooth of step 4,
`sum(n[x] (d[x] - phi[x])^2 / phi[x]^2)` divided by the number of ages
minus its effective dimension, alternated with the choice of its
smoothing parameter as `psi2` with `overdispersion = TRUE` (at most 5
times; the paper estimates `tau` without stating how); that smoothing
parameter is searched over `lambda_range` as those of the rates; ages
without data get weight 0 in step 4, so the smooth gives their
dispersion; predictive draws add Poisson variation of the deaths, as
with `overdispersion = TRUE`, because the quasi-likelihood defines the
mean and variance of the deaths but no distribution, and the paper gives
no predictive distribution. The paper's full extended quasi-likelihood
scheme (mean deviances by age, eq. 2.6) is not offered: it
underestimates the dispersion (section 4, Figure 4), and the paper
reports its results with the bias corrected scheme.

## References

Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
mortality rates. *Statistical Modelling* 4(4): 279-298.
[doi:10.1191/1471082X04st080oa](https://doi.org/10.1191/1471082X04st080oa)

Camarda CG (2015). MortalitySmooth: Smoothing and forecasting Poisson
counts with P-splines. R package version 2.3.4 (CRAN archive).

Camarda CG (2019). Smooth constrained mortality forecasting.
*Demographic Research* 41(38): 1091-1130.
[doi:10.4054/DemRes.2019.41.38](https://doi.org/10.4054/DemRes.2019.41.38)

Djeundje VAB, Currie ID (2011). Smoothing dispersed counts with
applications to mortality data. *Annals of Actuarial Science* 5(1):
33-52.
[doi:10.1017/S1748499510000047](https://doi.org/10.1017/S1748499510000047)

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
# dispersion that varies smoothly over age (Djeundje and Currie 2011),
# smoothing parameters of the rates given
fit_age <- mt_fit(p_splines(lambda = c(1, 1000), overdispersion = "age"), x)
p <- mt_parameters(fit_age)
p[parameter == "dispersion" & stratum == "female"]
#>     stratum  parameter  index       value
#>      <char>     <char> <char>       <num>
#>  1:  female dispersion     30   0.9671517
#>  2:  female dispersion     31   0.9790779
#>  3:  female dispersion     32   0.9909822
#>  4:  female dispersion     33   1.0029657
#>  5:  female dispersion     34   1.0151403
#>  6:  female dispersion     35   1.0276151
#>  7:  female dispersion     36   1.0403023
#>  8:  female dispersion     37   1.0529406
#>  9:  female dispersion     38   1.0652390
#> 10:  female dispersion     39   1.0768806
#> 11:  female dispersion     40   1.0875171
#> 12:  female dispersion     41   1.0965968
#> 13:  female dispersion     42   1.1033598
#> 14:  female dispersion     43   1.1070284
#> 15:  female dispersion     44   1.1068512
#> 16:  female dispersion     45   1.1021625
#> 17:  female dispersion     46   1.0929401
#> 18:  female dispersion     47   1.0799393
#> 19:  female dispersion     48   1.0639887
#> 20:  female dispersion     49   1.0458979
#> 21:  female dispersion     50   1.0264255
#> 22:  female dispersion     51   1.0063595
#> 23:  female dispersion     52   0.9865562
#> 24:  female dispersion     53   0.9677634
#> 25:  female dispersion     54   0.9506179
#> 26:  female dispersion     55   0.9356685
#> 27:  female dispersion     56   0.9233800
#> 28:  female dispersion     57   0.9141380
#> 29:  female dispersion     58   0.9083311
#> 30:  female dispersion     59   0.9063910
#> 31:  female dispersion     60   0.9088306
#> 32:  female dispersion     61   0.9159990
#> 33:  female dispersion     62   0.9275124
#> 34:  female dispersion     63   0.9428972
#> 35:  female dispersion     64   0.9616961
#> 36:  female dispersion     65   0.9834203
#> 37:  female dispersion     66   1.0074322
#> 38:  female dispersion     67   1.0326511
#> 39:  female dispersion     68   1.0576762
#> 40:  female dispersion     69   1.0808251
#> 41:  female dispersion     70   1.1001516
#> 42:  female dispersion     71   1.1137971
#> 43:  female dispersion     72   1.1213488
#> 44:  female dispersion     73   1.1230756
#> 45:  female dispersion     74   1.1194434
#> 46:  female dispersion     75   1.1110671
#> 47:  female dispersion     76   1.0987179
#> 48:  female dispersion     77   1.0836321
#> 49:  female dispersion     78   1.0671696
#> 50:  female dispersion     79   1.0505653
#> 51:  female dispersion     80   1.0349148
#> 52:  female dispersion     81   1.0211342
#> 53:  female dispersion     82   1.0095068
#> 54:  female dispersion     83   0.9999902
#> 55:  female dispersion     84   0.9925590
#> 56:  female dispersion     85   0.9872098
#> 57:  female dispersion     86   0.9839948
#> 58:  female dispersion     87   0.9834807
#> 59:  female dispersion     88   0.9866502
#> 60:  female dispersion     89   0.9945922
#> 61:  female dispersion     90   1.0085854
#> 62:  female dispersion     91   1.0302970
#> 63:  female dispersion     92   1.0634685
#> 64:  female dispersion     93   1.1149781
#> 65:  female dispersion     94   1.1955604
#> 66:  female dispersion     95   1.3242876
#> 67:  female dispersion     96   1.5403639
#> 68:  female dispersion     97   1.9309526
#> 69:  female dispersion     98   2.7450798
#> 70:  female dispersion     99   5.2136555
#> 71:  female dispersion   100+ 262.7173992
#>     stratum  parameter  index       value
#>      <char>     <char> <char>       <num>
```
