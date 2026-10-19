# Age-period-cohort-improvement (APCI) model

The stochastic implementation of the APCI model of Richards, Currie,
Kleinow and Ritchie (2019), fitted per stratum:
`log m[x, y] = alpha[x] + beta[x] (y - ybar) + kappa[y] + gamma[y - x]`
(eq. 8), with `ybar` the mean of the years of the data. The model needs
at least five single years of age, without an open interval, and at
least ten years (package-specific, for the choice of the ARIMA models).

## Usage

``` r
apci(smooth = TRUE, kappa_order = NULL, gamma_order = NULL)
```

## Arguments

- smooth:

  Smooth `alpha` and `beta` over age with P-splines, the smoothing
  parameters chosen by BIC (Richards et al. 2019, Table 1).

- kappa_order, gamma_order:

  ARIMA orders `c(p, d, q)` of the models of `kappa` and `gamma`
  (without a mean). `NULL` (default) chooses them as in Richards et al.
  (2019, Appendix E).

## Value

A model specification.

## Details

Following the paper:

- the deaths are Poisson with mean `m[x, y]` times the central exposure
  (eq. 4), and the parameters are estimated by maximum likelihood
  (Section 3);

- cohorts with four or fewer observations get `gamma = 0` and are not
  estimated (Appendix B);

- five identifiability constraints (Appendix C): `sum(kappa) = 0` (eq.
  17), `sum((y - ymin) kappa) = 0` (eq. 21), and `sum(w gamma) = 0`,
  `sum(w (c - cmin + 1) gamma) = 0` and
  `sum(w (c - cmin + 1)^2 gamma) = 0` (eq. 19, 20 and 22), with `w` the
  number of cells of cohort `c`. With the corner cohorts fixed at 0,
  only the constraints on `kappa` are needed; the constraints on `gamma`
  restrict the fit (the paper's over-constrained model, Section 6). The
  estimates maximise the likelihood subject to all five constraints,
  which hold exactly. For the paper's 55 ages and 45 years, the
  unsmoothed model then has its effective dimension of 241 (Table 2);

- with `smooth = TRUE` (the paper's APCI(S), Table 1), `alpha` and
  `beta` are P-splines over age (Eilers and Marx 1996) with penalties on
  the second-order differences of their coefficients (Appendix D;
  Appendix G calls the penalties quadratic), and the two smoothing
  parameters minimise the BIC, deviance plus `log(n)` times the
  effective dimension, with `n` the number of cells (eq. 24). `kappa`
  and `gamma` are not smoothed. With `smooth = FALSE`, every age has its
  own `alpha` and `beta`;

- `kappa` and `gamma` are projected by ARIMA models without a mean,
  fitted with [`stats::arima()`](https://rdrr.io/r/stats/arima.html) to
  the estimates as if they were observed (Section 5, Appendix E). By
  default the orders are chosen as in Appendix E: for `kappa`, the ARMA
  orders `p, q` in 0 to 3 with the smallest AICc for `d = 1` and for
  `d = 2`, then of these two the model with the smaller root mean
  squared error (Solo 1984); for `gamma`, the orders `p, q` in 0 to 3
  with the smallest AICc for `d = 1`. `gamma` is projected for the
  cohorts born after the last estimated one, which include the youngest
  cohorts with four or fewer observations (Appendix B). The two
  processes are independent;

- draws are simulated paths of both ARIMA processes, with their
  coefficients at the estimates. The paper adds parameter risk by a
  bootstrap of the mean of the process for `kappa` only, and the APCI
  processes have no mean (Appendix E). `alpha`, `beta` and the estimated
  `kappa` and `gamma` are fixed, so the deterministic part
  `beta[x] (y - ybar)` continues linearly. Predictive draws
  (`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
  given the rates.

Package-specific (reference unknown): the B-splines are cubic with
equally spaced knots about five years of age apart, on the range of the
ages widened by 1% at each end (as [`p_splines()`](p_splines.md); the
paper does not give the knots); the smoothing parameters are searched
between 1e-4 and 1e10 with the search of [`p_splines()`](p_splines.md);
the estimates come from penalised iteratively reweighted least squares
with the constraints as Lagrange conditions (the paper: the algorithm of
Currie 2013, with the constraints as rows of an augmented model matrix);
the AICc counts the ARMA coefficients and the innovation variance as
parameters, on the differenced series, as
[`forecast::Arima()`](https://pkg.robjhyndman.com/forecast/reference/Arima.html);
the root mean squared error that chooses `d` is taken over the one-step
residuals from the third year, which both orders of differencing have;
ARIMA paths are simulated from the state-space form of the fitted model,
from the filtered state at the last estimate. The paper's value-at-risk
framework and annuity factors are not implemented.

## References

Richards SJ, Currie ID, Kleinow T, Ritchie GP (2019). A stochastic
implementation of the APCI model for mortality projections. *British
Actuarial Journal* 24: e13.
[doi:10.1017/S1357321718000260](https://doi.org/10.1017/S1357321718000260)

Eilers PHC, Marx BD (1996). Flexible smoothing with B-splines and
penalties. *Statistical Science* 11(2): 89-121.
[doi:10.1214/ss/1038425655](https://doi.org/10.1214/ss/1038425655)

Currie ID (2013). Smoothing constrained generalized linear models with
an application to the Lee-Carter model. *Statistical Modelling* 13(1):
69-93.
[doi:10.1177/1471082X12471373](https://doi.org/10.1177/1471082X12471373)

Solo V (1984). The order of differencing in ARIMA models. *Journal of
the American Statistical Association* 79(388): 916-921.
[doi:10.1080/01621459.1984.10477111](https://doi.org/10.1080/01621459.1984.10477111)

## See also

[`stmomo()`](stmomo.md) for the age-period-cohort and Renshaw-Haberman
models.

## Examples

``` r
# \donttest{
if (requireNamespace("StMoMo", quietly = TRUE)) {
  ew <- mt_subset(as_mt_data(StMoMo::EWMaleData), ages = c(55, 90), periods = c(1971, 2011))
  fit <- mt_fit(apci(), ew)
  p <- mt_parameters(fit)
  p[parameter %in% c("order_kappa", "order_gamma")]
  fc <- mt_forecast(fit, h = 10, n_draws = 200, seed = 1)
  life_expectancy(fc, at = 65)
}
#>     stratum    sex  label period     mean   q0.025     q0.1     q0.5     q0.9
#>       <int> <char> <char>  <num>    <num>    <num>    <num>    <num>    <num>
#>  1:       1   male   male   2012 17.56500 17.37489 17.46074 17.55822 17.68231
#>  2:       1   male   male   2013 17.72639 17.55918 17.60802 17.73020 17.84358
#>  3:       1   male   male   2014 17.86312 17.66925 17.72444 17.85911 18.00434
#>  4:       1   male   male   2015 17.99612 17.73630 17.81373 17.99190 18.17277
#>  5:       1   male   male   2016 18.12505 17.81625 17.90695 18.13593 18.32173
#>  6:       1   male   male   2017 18.23230 17.82910 17.94403 18.24533 18.45270
#>  7:       1   male   male   2018 18.34464 17.90922 18.04602 18.36704 18.61837
#>  8:       1   male   male   2019 18.43827 17.94283 18.08382 18.44974 18.76477
#>  9:       1   male   male   2020 18.53830 17.96201 18.14273 18.56134 18.92158
#> 10:       1   male   male   2021 18.61627 17.91335 18.16544 18.64377 19.01162
#>       q0.975    at
#>        <num> <num>
#>  1: 17.75646    65
#>  2: 17.91052    65
#>  3: 18.06057    65
#>  4: 18.27048    65
#>  5: 18.40297    65
#>  6: 18.57092    65
#>  7: 18.75461    65
#>  8: 18.90759    65
#>  9: 19.06434    65
#> 10: 19.18684    65
# }
```
