# Coherent multi-population models

Joint models for several strata (for example women and men, or regions)
whose forecasts do not diverge in the long run:

## Usage

``` r
li_lee(ratio_order = 1, common_by = NULL, jump_off = c("fitted", "observed"))

product_ratio(ratio_order = 1, jump_off = c("fitted", "observed"))
```

## Arguments

- ratio_order:

  Number of stratum-specific (ratio) factors.

- common_by:

  Stratum variable(s) defining groups with their own common factor
  (`li_lee()` only).

- jump_off:

  `"fitted"` (forecasts start from the fitted rates of the last year) or
  `"observed"` (forecasts are multiplied by the ratio of observed to
  fitted rates in the last year, per age and stratum; Lee and Miller
  2001). With few factors per stratum the fitted rates of the last year
  can be far from the observed ones; the observed jump-off then avoids a
  step between the last observed and the first forecast year.

## Value

A model specification.

## Details

- `li_lee()` (Li and Lee 2005):
  `log m[x,t,i] = a[x,i] + B[x] K[t] + b[x,i] k[t,i]`. The common factor
  `B K` is a Lee-Carter model of the pooled rates (sum of deaths over
  sum of exposures); the stratum factors `b k` are fitted to the
  residuals.

- `product_ratio()` (Hyndman, Booth and Yasmeen 2013): the geometric
  mean of the stratum rates (product) is modelled with a Lee-Carter-type
  factor, and the ratios of each stratum to the product with factors
  whose scores are stationary. This is a Lee-Carter-type version of the
  method; the original uses functional data models.

With `common_by` (for example `"sex"`), `li_lee()` fits one common
factor per level of that stratum variable (women and men each with their
own trend), and the other strata (for example regions and deprivation
quintiles) deviate from their group's common factor. If the data carry a
hierarchy (see [`mt_aggregate_strata()`](mt_aggregate_strata.md)), the
common factors are estimated from the bottom strata only, so aggregates
are not counted twice.

In both models the common period index follows a random walk with drift
(with drift uncertainty) and the stratum-specific scores follow AR(1)
processes without drift, so the differences between strata stay bounded.
Draws include the innovations of all indexes; predictive draws add
extra-Poisson and Poisson variation.

## References

Li N, Lee R (2005). Coherent mortality forecasts for a group of
populations: an extension of the Lee-Carter method. *Demography* 42(3):
575-594.
[doi:10.1353/dem.2005.0021](https://doi.org/10.1353/dem.2005.0021)

Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
method for forecasting mortality. *Demography* 38(4): 537-549.

Hyndman RJ, Booth H, Yasmeen F (2013). Coherent mortality forecasting:
the product-ratio method with functional time series models.
*Demography* 50(1): 261-283.

## Examples

``` r
x <- mt_example(deaths = "single")
fc <- mt_forecast(mt_fit(li_lee(), x), h = 10, n_draws = 100, seed = 1)
```
