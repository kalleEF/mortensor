# Tucker3 tensor model across strata (three-way Lee-Carter model)

Joint model of the array of log rates `[age, period, stratum]`: the
three-way extension of the Lee-Carter model of Russolillo, Giordano and
Haberman (2011). Status: published, with the deviations listed below;
section and equation numbers refer to that paper.

## Usage

``` r
tucker3(
  ranks = c(2, 2, 2),
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

  Tucker ranks for age, period and stratum (default `c(2, 2, 2)`,
  Russolillo et al. 2011, Section 3.2), or `"bic"` / `"aic"` to choose
  them by an information criterion on the training data
  (package-specific: `n log(RSS / n)` plus the penalty times the number
  of parameters of the factor matrices and the core; ranks up to 6 for
  age and period and up to 4 for strata). For selection by forecast
  accuracy use [`tuned()`](mt_tune.md).

- transform:

  `"log"`, or the deprecated `"coda"` (use
  [`coda_multiway()`](coda_multiway.md)).

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

**Model (Section 2, eq. 7-8).** As in Lee-Carter, the mean log rate over
the years of each age and stratum, `a[x, i]`, is removed
(`centre = TRUE`), and the centred array is decomposed by a Tucker3
decomposition `X - a = G x1 A x2 B x3 C` with ranks
`ranks = c(age, period, stratum)`: column-wise orthonormal age factors
`A`, period factors `B` and stratum factors `C`, linked by the core
array `G`, estimated by least squares. Ranks `c(1, 1, 1)` give eq. (7),
`log m[x, t, i] - a[x, i] = b[x] k[t] c[i]` (one period index for all
strata); the default `c(2, 2, 2)` is the model of the application
(Section 3.2, ten European countries). The squared core elements over
the sum of squares of the centred array are the shares of variation of
each combination of components (Table 4); their sum is the fit
(`fit$fits$fit_pct`, `core_pct`).

**Forecasts.** The paper proposes a time-series model for the period
index, as in Lee-Carter, and leaves its form open (Section 3.2.2). Here
the period factors (columns of `B`) follow a multivariate random walk
with drift, as `k[t]` of Lee and Carter (1992); draws resample the joint
year-to-year changes, which keeps their correlation, and include the
uncertainty of the drift (`ts = "rwd"`). ARIMA or ETS models per column
(`ts = "arima"`, `"ets"`, package forecast) are alternatives.
Bergeron-Boucher, Simonacci, Oeppen and Gallo (2018) forecast this model
with one component per mode and a random walk with drift of the period
index, as `tucker3(ranks = c(1, 1, 1))` does; their draws resample the
errors of the index and keep the drift fixed. The rates are rebuilt with
`a`, `A`, `C` and `G` fixed. Draws carry the uncertainty of the forecast
period factors only. With fewer period factors than independent trends
in the data, the strata share their forecast uncertainty and intervals
are too narrow: in an example on simulated Lee-Carter data with Poisson
deaths and an own random walk in each of 10 strata, 95% intervals
covered about 74% with the default ranks and 95% with
`ranks = c(1, 10, 10)`. Predictive draws
(`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
given the rates and the exposures, and no other observation noise. If
intervals are too narrow on real data, that is a property of the model
(see [`recalibrated()`](recalibrated.md), [`ensemble()`](ensemble.md),
or overdispersed models such as `stanmomo(family = "nb")`).

**Deviations from the paper.**

- Estimation: the paper fits the model by alternating least squares
  (footnote 3; Kiers' Tucker3 program for Matlab, Appendix); here by
  higher-order orthogonal iteration (HOOI; De Lathauwer, De Moor and
  Vandewalle 2000), an alternating least-squares algorithm for the same
  criterion. The paper also rotates the solution to simple structure,
  which leaves the fitted values unchanged; signs are fixed so that the
  largest entry of each factor is positive.

- Centring: the text and eq. (7) centre over the years; the Matlab
  session in the Appendix also centres across the age groups. The
  package follows eq. (7), because the term removed by centring across
  ages varies over the years and the paper gives no forecast for it.

- Package-specific (reference unknown): the forecast model of the period
  factors (see above), the choice of the ranks by an information
  criterion (`ranks = "bic"`, `"aic"`), and the replacement of zero
  deaths by 0.5 before taking logarithms. The paper's data end with the
  closed age group 95-99; an open age interval is treated like the other
  age groups.

**Deprecated:** `transform = "coda"` (centred log-ratio transformed
life-table deaths in place of log rates) still works and warns (class
`mt_warning_deprecated`). Its forecast of all period factors by a
multivariate random walk, its auxiliary random walk for the rate of the
open age interval and its 0.5 deaths in empty cells are
package-specific. Use [`coda_multiway()`](coda_multiway.md), the
published multiway compositional model of Bergeron-Boucher, Simonacci,
Oeppen and Gallo (2018), or [`coda()`](coda.md) (`coherent = TRUE`) for
a common component with deviations per stratum.

Ranks are reduced to the array dimensions where needed (at most the
number of age groups minus 1, of years minus 1 and of strata, and at
least 1), so the model also runs with one stratum or one age group; a
warning (class `mt_warning_argument`) names the requested and the ranks
used, which [`mt_parameters()`](mt_parameters.md) reports (`rank`). The
model is experimental: compare it with established models such as
[`li_lee()`](li_lee.md) in a backtest before using it.

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

Bergeron-Boucher M-P, Simonacci V, Oeppen J, Gallo M (2018). Coherent
modeling and forecasting of mortality patterns for subpopulations using
multiway analysis of compositions: an application to Canadian provinces
and territories. *North American Actuarial Journal* 22(1): 92-118.
[doi:10.1080/10920277.2017.1377620](https://doi.org/10.1080/10920277.2017.1377620)

## Examples

``` r
x <- mt_example(deaths = "single")
fit <- mt_fit(tucker3(), x)
fit$fits$fit_pct
#> [1] 79.90136
fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
```
