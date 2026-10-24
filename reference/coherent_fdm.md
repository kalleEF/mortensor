# Coherent functional demographic model (product-ratio, demography)

The product-ratio method of Hyndman, Booth and Yasmeen (2013), with the
functional data models of the demography package: for a group of strata
(for example women and men, or the deprivation quintiles of one sex)

## Usage

``` r
coherent_fdm(
  order1 = 6,
  order2 = 6,
  common_by = NULL,
  smooth = TRUE,
  smooth_args = list(),
  K = 100,
  drange = c(0, 0.5),
  weighted = TRUE,
  beta = 0.05
)
```

## Arguments

- order1, order2:

  Number of principal components of the product and of the ratios.

- common_by:

  Stratum variable(s) defining groups with their own product.

- smooth:

  Smooth the rates of each stratum and year over age before the
  decomposition (see [`hyndman_ullah()`](hyndman_ullah.md)).

- smooth_args:

  List of arguments of
  [`demography::smooth.demogdata()`](https://pkg.robjhyndman.com/demography/reference/smooth.demogdata.html).

- K:

  Number of most recent years used to forecast the ratios.

- drange:

  Range of the fractional differencing parameter of the ARFIMA models of
  the ratios (`c(0, 0.5)`: stationary).

- weighted, beta:

  Weighted principal components of the product and of the ratios, with
  geometrically decaying weights `beta (1 - beta)^i` of the years
  (Hyndman, Booth and Yasmeen 2013: `beta = 0.05`). `weighted = FALSE`
  gives the unweighted decomposition, the default of
  [`demography::coherentfdm()`](https://pkg.robjhyndman.com/demography/reference/coherentfdm.html).

## Value

A model specification.

## Details

- the **product**, the geometric mean of the strata's rates, is modelled
  with a functional demographic model
  ([`hyndman_ullah()`](hyndman_ullah.md)) whose scores are forecast with
  ARIMA models, so it carries the common trend;

- the **ratio** of each stratum to the product is modelled with a
  functional demographic model whose scores are forecast with stationary
  ARFIMA models, so the differences between the strata stay bounded and
  revert towards their recent level.

Forecasts of the strata are products of the two. This is the method of
[`demography::coherentfdm()`](https://pkg.robjhyndman.com/demography/reference/coherentfdm.html)
and `forecast()` for its result; [`li_lee()`](li_lee.md) is a related
model with a common Lee-Carter factor.

Draws simulate the product and each ratio (uncertainty of the scores and
resampled model errors, as
[`simulate()`](https://rdrr.io/r/stats/simulate.html) for functional
forecasts in demography) and multiply them; the ratios of each path are
rescaled so that they multiply to one, as `forecast()` in demography
does for the point forecasts.

**Relation to the paper.** As in Hyndman, Booth and Yasmeen (2013,
Sections 2.1-2.4): the rates of each stratum are smoothed over age first
(weighted penalised regression splines, monotone above age 65, through
[`demography::smooth.demogdata()`](https://pkg.robjhyndman.com/demography/reference/smooth.demogdata.html));
the product is the geometric mean of the smoothed rates of the group
(their Eq. 3), and every stratum has a ratio model (their Eq. 4; they
note that J - 1 ratio models would suffice); the scores of the product
are forecast with ARIMA models chosen by the algorithm of Hyndman and
Khandakar (2008), those of the ratios with stationary ARFIMA models
(orders by the same algorithm, maximum likelihood); there is no jump-off
adjustment; `order1 = order2 = 6` are the numbers of components of their
applications. The principal components are weighted, with weight
`beta (1 - beta)^i` for the year `i` years before the last and
`beta = 0.05` as in their applications (Section 2.2 and Discussion);
[`demography::coherentfdm()`](https://pkg.robjhyndman.com/demography/reference/coherentfdm.html)
does not weight them by default (`weighted = FALSE`). The departures:

- the fractional differencing parameter of the ratio models lies in
  `drange = c(0, 0.5)` (the default of demography); the paper names
  `-0.5 < d < 0.5` as the stationary range;

- prediction intervals come from simulated paths with resampled model
  errors; the paper gives analytic variances for the rates (their Eq. 8)
  and simulates for life expectancy. Predictive draws add Poisson deaths
  in place of the observational variance of their Eq. 8;

- the ratios of each simulated path are rescaled to multiply to one, and
  with a hierarchy (see below) the aggregates get a ratio but do not
  enter the product. Both are package-specific (reference unknown);
  demography rescales the point forecasts in the same way, and the paper
  states the constraint for the observed ratios. Predictive draws
  (`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
  given the rates and the exposures, and no other observation noise. If
  intervals are too narrow on real data, that is a property of the model
  (see [`recalibrated()`](recalibrated.md), [`ensemble()`](ensemble.md),
  or overdispersed models such as `stanmomo(family = "nb")`).

With `common_by` (for example `"sex"`), each level of that stratum
variable forms its own group with its own product (for example the
quintiles of each sex). If the data carry a hierarchy (see
[`mt_aggregate_strata()`](mt_aggregate_strata.md)), the product is
computed from the bottom strata; aggregates get their own ratio to it.
Every group needs at least two strata
([`hyndman_ullah()`](hyndman_ullah.md) models a single stratum).

## References

Hyndman RJ, Booth H, Yasmeen F (2013). Coherent mortality forecasting:
the product-ratio method with functional time series models.
*Demography* 50(1): 261-283.
[doi:10.1007/s13524-012-0145-5](https://doi.org/10.1007/s13524-012-0145-5)

Hyndman RJ, Khandakar Y (2008). Automatic time series forecasting: the
forecast package for R. *Journal of Statistical Software* 27(3): 1-22.
[doi:10.18637/jss.v027.i03](https://doi.org/10.18637/jss.v027.i03)

## Examples

``` r
if (requireNamespace("demography", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 1990:2019)
  fit <- mt_fit(coherent_fdm(order1 = 2, order2 = 1), x)
  fc <- mt_forecast(fit, h = 10, n_draws = 50, seed = 1)
}
```
