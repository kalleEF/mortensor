# Compositional data (CoDa) model of life-table deaths

`coda()` is the Lee-Carter model in compositional form (Oeppen 2008), as
formulated by Bergeron-Boucher, Canudas-Romo, Oeppen and Vaupel (2017).
It models the life-table distribution of deaths by age (which sums to
1). With `coherent = TRUE` or `common_by`, it is their coherent model
for several populations. The code follows the reference implementation,
`model.Oeppen()` and `model.OeppenC()` of the R package
MortalityForecast 0.8.1 (Pascariu and Bergeron-Boucher), and reproduces
its fitted values and point forecasts.

## Usage

``` r
coda(jump_off = c("observed", "fitted"), coherent = FALSE, common_by = NULL)
```

## Arguments

- jump_off:

  `"observed"` (forecasts start from the observed deaths of the last
  year) or `"fitted"` (from the fitted deaths).

- coherent:

  Fit the coherent model with a common component per group.

- common_by:

  Stratum variable(s) defining groups with their own common component
  (implies `coherent = TRUE`); `NULL` with `coherent = TRUE` puts all
  strata in one group.

## Value

A model specification.

## Details

The model of one stratum has six steps:

1.  the death rates of each year give life-table deaths `d[x, t]` with
    radix 1 at the first age ([`life_table()`](life_table.md)
    conventions); for a closed age range, the survivors at its upper
    bound are an extra part, so that the parts sum to 1;

2.  `a[x]` is the geometric mean of `d[x, t]` over the years, closed to
    1;

3.  the deaths are centred by perturbation (`d[x, t] / a[x]`, closed)
    and transformed by the centred log-ratio (clr);

4.  the first component of a singular value decomposition gives
    `clr(d[x, t] / a[x]) = b[x] k[t]` (`b` sums to zero);

5.  `k[t]` follows a random walk with drift (ARIMA(0,1,0) with drift,
    the default of the reference implementation);

6.  forecasts are transformed back (inverse clr, perturbation by `a`,
    closure) to life-table deaths and then to death rates.

With `jump_off = "observed"` (the default, as `jumpchoice = "actual"` of
the reference implementation), forecasts are the observed deaths of the
last year `T` perturbed by the change of the model: `d[x, T + j]` is
proportional to `d_obs[x, T] exp(b[x] (k[T + j] - k[T]))`. With
`"fitted"`, forecasts start from the fitted deaths of the last year. The
forecast of year `T + j` uses the index `j` steps ahead; the reference
implementation uses `j + 1` steps, which makes no difference for the
observed jump-off with a random walk with drift.

**Coherent model.** With `coherent = TRUE` or `common_by` (for example
`"sex"`), the strata of each group (all strata, or each level of
`common_by`) share a common component (Bergeron-Boucher et al. 2017):
`clr(d[x, t, i] / a[x, i]) = B[x] K[t] + b[x, i] k[t, i]`. The common
component `B K` is the first component of the CoDa model of the average
distribution of the group, the closed geometric mean of the life-table
deaths of its strata (the centre of the compositions). The deviations of
each stratum are the first component of the remaining variation. `K`
follows a random walk with drift, and `k[t, i]` an AR(1) process with
mean (estimated by maximum likelihood), as in the reference
implementation (`model.OeppenC()`, defaults `order.B = c(0, 1, 0)` with
drift and `order.D = c(1, 0, 0)` without drift), so the deviations stay
bounded. If the data carry a hierarchy (see
[`mt_aggregate_strata()`](mt_aggregate_strata.md)), the average uses the
bottom strata only, as [`li_lee()`](li_lee.md) does.

**Open age interval.** Life-table deaths of an open interval are the
survivors at its lower bound; they do not determine its death rate, and
the reference implementation forecasts the deaths only. Here the log
death rate of the open interval gets loadings on the same indexes,
estimated by least squares: `log m[open, t] = a_o + b_o k[t]`, and in
the coherent model
`log m[open, t, i] = a_o[i] + B_o K[t] + b_o[i] k[t, i]`, with `B_o`
from the mean log rate of the group. The open interval thus follows the
rest of the age pattern (and the jump-off rule). An auxiliary random
walk on the observed log rate, as in [`tucker3()`](tucker3.md),
accumulates the noise of the few deaths at the highest ages: on French
data (ages 100+, forecasts of 10 years) a quarter to a half of its paths
exceeded a rate of 1.

**Zero deaths.** Cells with zero deaths have no logarithm. As
`replace.zeros()` of the reference implementation (method `"min"`), zero
rates are replaced by half the smallest positive rate of the stratum in
the training years; a warning is given when more than 10% of the cells
are affected.

Draws simulate the indexes: the random walk with drift includes the
uncertainty of the drift (normal, with the variance of the mean change),
the AR(1) processes their innovations; `a` and `b` are fixed. Draws are
expected rates; predictive draws (`mt_forecast(predictive = TRUE)`) add
Poisson variation of the deaths given the rates and the exposures, and
no other observation noise. If intervals are too narrow on real data,
that is a property of the model (see
[`recalibrated()`](recalibrated.md), [`ensemble()`](ensemble.md), or
overdispersed models such as `stanmomo(family = "nb")`). The model
weights all ages equally on the clr scale, so the noise of ages with few
deaths enters `k`; its intervals can be wider than those of
[`lee_carter()`](lee_carter.md).

## References

Oeppen J (2008). Coherent forecasting of multiple-decrement life tables:
a test using Japanese cause of death data. CODAWORK'08, Universitat de
Girona.

Bergeron-Boucher M-P, Canudas-Romo V, Oeppen J, Vaupel JW (2017).
Coherent forecasts of mortality with compositional data analysis.
*Demographic Research* 37(17): 527-566.
[doi:10.4054/DemRes.2017.37.17](https://doi.org/10.4054/DemRes.2017.37.17)

## See also

[`coda_causes()`](coda_causes.md) for causes of death,
[`lee_carter()`](lee_carter.md) for the model of death rates,
[`li_lee()`](li_lee.md) for its coherent version.

## Examples

``` r
x <- mt_example(deaths = "single", years = 1995:2012)
fc <- mt_forecast(mt_fit(coda(), x), h = 5, n_draws = 50, seed = 1)
fc2 <- mt_forecast(mt_fit(coda(coherent = TRUE), x), h = 5, n_draws = 50, seed = 1)
```
