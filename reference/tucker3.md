# Tucker3 tensor model across strata (three-way Lee-Carter model)

Joint model of the array of log rates `[age, period, stratum]`: the
three-way extension of the Lee-Carter model of Russolillo, Giordano and
Haberman (2011). As in Lee-Carter, the mean log rate over the years of
each age and stratum, `a[x, i]`, is removed (`centre = TRUE`), and the
centred array is decomposed by a Tucker3 decomposition
`X - a = G x1 A x2 B x3 C` with ranks `ranks = c(age, period, stratum)`:
age factors `A`, period factors `B` and stratum factors `C` linked by
the core array `G`. The period factors (columns of `B`) are forecast and
the rates rebuilt with `a`, `A`, `C` and `G` fixed.

## Usage

``` r
tucker3(
  ranks = c(3, 3, 2),
  transform = c("log", "coda"),
  ts = c("rwd", "arima", "ets"),
  residual_noise = FALSE,
  centre = TRUE,
  maxit = 200,
  tol = 1e-08
)
```

## Arguments

- ranks:

  Tucker ranks for age, period and stratum, or `"bic"` / `"aic"` to
  choose them by an information criterion on the training data
  (`n log(RSS / n)` plus the penalty times the number of parameters of
  the factor matrices and the core; ranks up to 6 for age and period and
  up to 4 for strata). For selection by forecast accuracy use
  [`tuned()`](mt_tune.md).

- transform:

  `"log"`, or the deprecated `"coda"`.

- ts:

  Forecast method for the period factors: `"rwd"`, `"arima"`, `"ets"`.

- residual_noise:

  Add the in-sample residual variation (per age and stratum, normal on
  the transformed scale) to the draws. The residuals contain the
  observation noise of the training years, which predictive draws add
  again as Poisson variation, so the noise is then counted twice.
  Default `FALSE` (earlier development versions: `TRUE`).

- centre:

  Remove the mean over the years of each age and stratum before the
  decomposition and add it back afterwards (Russolillo, Giordano and
  Haberman 2011). With `FALSE` the decomposition describes the log rates
  themselves.

- maxit, tol:

  Settings of the HOOI iterations.

## Value

A model specification.

## Details

The following parts are choices of this package without a published
source:

- estimation by higher-order orthogonal iteration (HOOI; De Lathauwer,
  De Moor and Vandewalle 2000), least squares on the log rates;

- the forecast of the period factors: a multivariate random walk with
  drift (bootstrap of joint year-to-year changes, which keeps their
  correlation), or ARIMA/ETS models per column (`ts = "arima"`, `"ets"`,
  package forecast);

- the ranks (default `c(3, 3, 2)`) and their choice by an information
  criterion (`ranks = "bic"`, `"aic"`);

- zero deaths are replaced by 0.5 before taking logarithms.

The model is experimental: compare it with established models such as
[`li_lee()`](li_lee.md) in a backtest before using it.

Draws carry the uncertainty of the forecast period factors; the other
factors are fixed. With fewer period factors than independent trends in
the data, the strata share their forecast uncertainty and intervals are
too narrow: on simulated Lee-Carter data with Poisson deaths and an own
random walk in each of 10 strata (the package's coverage benchmark), 95%
intervals covered about 75% with the default ranks and 95% with
`ranks = c(1, 10, 10)`. Predictive draws
(`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
given the rates and the exposures, and no other observation noise. If
intervals are too narrow on real data, that is a property of the model
(see [`recalibrated()`](recalibrated.md), [`ensemble()`](ensemble.md),
or overdispersed models such as `stanmomo(family = "nb")`).

**Deprecated:** `transform = "coda"` (centred log-ratio transformed
life-table deaths in place of log rates, with an auxiliary random walk
for the rate of the open age interval) is a combination of this package
without a published source; it still works and warns (class
`mt_warning_deprecated`). Use [`coda()`](coda.md) (`coherent = TRUE` for
several strata), the published compositional model.

Ranks are reduced to the array dimensions where needed, so the model
also runs with one stratum or one age group.

## References

Russolillo M, Giordano G, Haberman S (2011). Extending the Lee-Carter
model: a three-way decomposition. *Scandinavian Actuarial Journal*
2011(2): 96-117.
[doi:10.1080/03461231003611933](https://doi.org/10.1080/03461231003611933)

Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
*Journal of the American Statistical Association* 87(419): 659-671.

De Lathauwer L, De Moor B, Vandewalle J (2000). On the best rank-1 and
rank-(R1, R2, ..., RN) approximation of higher-order tensors. *SIAM
Journal on Matrix Analysis and Applications* 21(4): 1324-1342.
[doi:10.1137/S0895479898346995](https://doi.org/10.1137/S0895479898346995)

Kolda TG, Bader BW (2009). Tensor decompositions and applications. *SIAM
Review* 51(3): 455-500.
[doi:10.1137/07070111X](https://doi.org/10.1137/07070111X)

## Examples

``` r
x <- mt_example(deaths = "single")
fit <- mt_fit(tucker3(ranks = c(3, 3, 2)), x)
fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
```
