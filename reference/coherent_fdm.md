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
  drange = c(0, 0.5)
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

Forecasts of the strata are products of the two. This is the original
method of
[`demography::coherentfdm()`](https://pkg.robjhyndman.com/demography/reference/coherentfdm.html);
[`product_ratio()`](li_lee.md) is a Lee-Carter version of it, and
[`li_lee()`](li_lee.md) a related model with a common Lee-Carter factor.

Draws simulate the product and each ratio (uncertainty of the scores and
resampled model errors, as
[`simulate()`](https://rdrr.io/r/stats/simulate.html) for functional
forecasts in demography) and multiply them; the ratios of each path are
rescaled so that they multiply to one, as in
[`demography::forecast.fdmpr()`](https://pkg.robjhyndman.com/demography/reference/forecast.fdmpr.html).
Predictive draws add the extra-Poisson and Poisson variation of the raw
deaths.

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

## Examples

``` r
if (requireNamespace("demography", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 1990:2019)
  fit <- mt_fit(coherent_fdm(order1 = 2, order2 = 1), x)
  fc <- mt_forecast(fit, h = 10, n_draws = 50, seed = 1)
}
```
