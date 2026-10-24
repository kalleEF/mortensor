# Multiway compositional (CoDa) model of life-table deaths for subpopulations

`coda_multiway()` is the three-way compositional model of
Bergeron-Boucher, Simonacci, Oeppen and Gallo (2018): a Tucker3
decomposition of the centred log-ratios of the life-table deaths of
several populations, with time indexes shared by all of them. The strata
are the populations (for example the regions of a country). Status:
published, with the choices listed below where the paper leaves the
method open. Section and equation numbers refer to that paper. The model
is experimental: compare it with [`coda()`](coda.md) (`coherent = TRUE`)
and [`li_lee()`](li_lee.md) in a backtest before using it.

## Usage

``` r
coda_multiway(
  ranks = c(1, 1, 1),
  order = c(0, 1, 0),
  drift = TRUE,
  jump_off = c("fitted", "observed"),
  open = c("auto", "kannisto", "regression")
)
```

## Arguments

- ranks:

  Numbers of components `c(age, time, population)`, by default
  `c(1, 1, 1)`. They are reduced where the data allow fewer components
  (at most the number of parts of the composition minus 1, of years
  minus 1 and of populations, and no rank above the product of the other
  two), with a warning (class `mt_warning_argument`) that names the
  requested and the ranks used; the label of the fit and
  [`mt_parameters()`](mt_parameters.md) (`rank`) give the ranks used.

- order:

  ARIMA order `c(p, d, q)` of each time index; `c(0, 1, 0)` by default
  (random walk, Bergeron-Boucher et al. 2018, Section 5).

- drift:

  Include a drift when `d = 1`; ignored otherwise (`d = 0` includes a
  mean, `d = 2` no constant).

- jump_off:

  `"fitted"` (forecasts start from the fitted deaths) or `"observed"`
  (from the observed deaths of the last year).

- open:

  Treatment of the old ages: `"auto"` (the Kannisto extension where the
  age groups allow it, otherwise the regression of the open interval),
  `"kannisto"` or `"regression"`, as in [`coda()`](coda.md).

## Value

A model specification.

## Details

**Model (Section 2.2, eq. 5).** For population `i`:

1.  the death rates of each year give life-table deaths `d[x, t, i]`
    with radix 1 at the first age ([`life_table()`](life_table.md)
    conventions), as in [`coda()`](coda.md);

2.  `a[x, i]` is the geometric mean of `d[x, t, i]` over the years,
    closed to 1;

3.  the deaths are centred by perturbation (`d[x, t, i] / a[x, i]`,
    closed) and transformed by the centred log-ratio (clr), eq. (3);

4.  the array of centred log-ratios `[age, year, population]` follows
    the Tucker3 model
    `clr(d[x, t, i] / a[x, i]) = sum_pqr g[p, q, r] b[x, p] k[t, q] gamma[i, r]`,
    with age patterns `b`, time indexes `k` and population markers
    `gamma` (the paper's beta, kappa and gamma), each of unit length,
    and the core array `g`; `ranks = c(P, Q, R)` gives the numbers of
    components of age, time and population. The default `c(1, 1, 1)` is
    the model of the paper's forecasts (Sections 4 and 5);

5.  the time index follows a random walk with drift
    (`order = c(0, 1, 0)`, `drift = TRUE`), the ARIMA model that the
    paper selected by AIC for the index of Canadian women (Section 5);

6.  forecasts are transformed back (inverse clr, perturbation by
    `a[x, i]`, closure) to life-table deaths and then to death rates.

With one component per mode, population `i` changes by
`g gamma[i] b[x] k[t]`: all populations follow one index and one age
pattern, at a pace set by their markers. The forecasts are coherent in
this sense; populations with different markers keep changing at
different paces (Section 5.1). With one population, the fit and the
point forecasts are those of [`coda()`](coda.md) with the same index
model and jump-off; the draws of [`coda()`](coda.md) also carry the
uncertainty of the estimates.

**Estimation (Sections 2.2 and 4).** Least squares by alternating least
squares: higher-order orthogonal iteration (De Lathauwer, De Moor and
Vandewalle 2000), started from the higher-order singular value
decomposition. The paper uses the Tucker3 model for compositional data
of Gallo (2015), an alternating least-squares algorithm for the same
criterion. The clr of a composition sums to zero and the centring
removes the mean over the years, so `b` and `k` sum to zero. Signs: each
time index increases over the fitting period, each set of population
markers has a positive sum, and the largest core element of each age
pattern is positive; with one component, life-table deaths then move
from ages with negative `b` to ages with positive `b` as `k` rises.
`fit$fits$fit_pct` is the share of the sum of squares of the centred
array that the model explains (the explained variance of Section 4.1),
`fit$fits$core_pct` the share of each core element.

**Draws (Section 5).** The paper draws its prediction intervals from
simulations of the time index with resampled errors. Here each time
index is simulated from its ARIMA model with innovations resampled from
its residuals; `a`, `b`, `gamma`, `g` and the drift are fixed, and the
draws of the index are shared by all populations. Draws are expected
rates; predictive draws (`mt_forecast(predictive = TRUE)`) add Poisson
variation of the deaths given the rates and the exposures, and no other
observation noise.

**Choices where the paper is silent.**

- Jump-off: the paper states no adjustment, so forecasts start from the
  fitted deaths (`jump_off = "fitted"`). `jump_off = "observed"` starts
  each population from its observed deaths of the last year, perturbed
  by the change of the model, as in [`coda()`](coda.md)
  (Bergeron-Boucher, Canudas-Romo, Oeppen and Vaupel 2017, Appendix G).

- Several time indexes (`ranks[2] > 1`): the paper forecasts with one
  index. Each index then has its own ARIMA model, independent of the
  others, as for several components of a CoDa model in Shang and
  Haberman (2020).

- Other ARIMA orders (`order`, `drift`): the point forecast of an
  ARIMA(0,1,q) index with drift gets the first-step correction of
  [`coda()`](coda.md) (Bergeron-Boucher et al. 2017, eq. 20); its draws
  are shifted by the same amount.

- Old ages and zero deaths, as in [`coda()`](coda.md): with single-year
  age groups from 80 (`open = "auto"` or `"kannisto"`), the rates from
  age 80 come from a Kannisto model fitted per year and population and
  extended to age 120 (Bergeron-Boucher et al. 2017, Appendix C); with
  other age groups, the log death rate of the open interval of each
  population is regressed on the time indexes (package-specific). Zero
  deaths are replaced multiplicatively (Bergeron-Boucher et al. 2017,
  Appendix C). The paper uses five-year age groups to 100+ and five-year
  periods and combines the smallest populations to avoid zero counts;
  the models of this package need single calendar years.

- Strata of a hierarchy (see
  [`mt_aggregate_strata()`](mt_aggregate_strata.md)) enter the
  decomposition as populations of their own;
  [`reconciled()`](reconciled.md) makes their sums coherent.

The paper chooses the numbers of components by the DifFit criterion
(Timmerman and Kiers 2000), gives bootstrap confidence intervals of the
parameters (Kiers 2004) and compares models by the AIC of Poisson
deviances (Table 2); these are not implemented. It fits women and men
separately, as
`stratified("sex", female = coda_multiway(), male = coda_multiway())`
does.

## References

Bergeron-Boucher M-P, Simonacci V, Oeppen J, Gallo M (2018). Coherent
modeling and forecasting of mortality patterns for subpopulations using
multiway analysis of compositions: an application to Canadian provinces
and territories. *North American Actuarial Journal* 22(1): 92-118.
[doi:10.1080/10920277.2017.1377620](https://doi.org/10.1080/10920277.2017.1377620)

Bergeron-Boucher M-P, Canudas-Romo V, Oeppen J, Vaupel JW (2017).
Coherent forecasts of mortality with compositional data analysis.
*Demographic Research* 37(17): 527-566.
[doi:10.4054/DemRes.2017.37.17](https://doi.org/10.4054/DemRes.2017.37.17)

Gallo M (2015). Tucker3 model for compositional data. *Communications in
Statistics - Theory and Methods* 44(21): 4441-4453.
[doi:10.1080/03610926.2013.798664](https://doi.org/10.1080/03610926.2013.798664)

Shang HL, Haberman S (2020). Forecasting age distribution of death
counts: an application to annuity pricing. *Annals of Actuarial Science*
14(1): 150-169.
[doi:10.1017/S1748499519000101](https://doi.org/10.1017/S1748499519000101)

De Lathauwer L, De Moor B, Vandewalle J (2000). On the best rank-1 and
rank-(R1, R2, ..., RN) approximation of higher-order tensors. *SIAM
Journal on Matrix Analysis and Applications* 21(4): 1324-1342.
[doi:10.1137/S0895479898346995](https://doi.org/10.1137/S0895479898346995)

Timmerman ME, Kiers HAL (2000). Three-mode principal components
analysis: choosing the numbers of components and sensitivity to local
optima. *British Journal of Mathematical and Statistical Psychology*
53(1): 1-16.
[doi:10.1348/000711000159132](https://doi.org/10.1348/000711000159132)

Kiers HAL (2004). Bootstrap confidence intervals for three-way methods.
*Journal of Chemometrics* 18(1): 22-36.
[doi:10.1002/cem.841](https://doi.org/10.1002/cem.841)

## See also

[`coda()`](coda.md) for one population and for the coherent model with a
common component and deviations, [`tucker3()`](tucker3.md) for the
three-way model of log death rates.

## Examples

``` r
x <- mt_example_regions(years = 2000:2015, regions = 3)
fit <- mt_fit(coda_multiway(), x)
fit$fits$fit_pct
#> [1] 38.67032
fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
```
