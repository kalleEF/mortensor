# Lee-Carter version of the product-ratio method (deprecated)

**Deprecated.** `product_ratio()` replaces the functional data models of
the product-ratio method of Hyndman, Booth and Yasmeen (2013) by
Lee-Carter-type factors; this combination is a method of this package
without a published source. Use [`coherent_fdm()`](coherent_fdm.md), the
published product-ratio method. `product_ratio()` still works and warns
(class `mt_warning_deprecated`).

## Usage

``` r
product_ratio(ratio_order = 1, jump_off = c("fitted", "observed"))
```

## Arguments

- ratio_order:

  Number of stratum-specific (ratio) factors.

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

The geometric mean of the stratum rates (product) is modelled with a
Lee-Carter-type factor (SVD of its log rates, period index a random walk
with drift), and the ratio of each stratum to the product with factors
whose indexes follow AR(1) processes without constant, bounded to
`|phi| <= 0.99` (see [`li_lee()`](li_lee.md)). Predictive draws add
Poisson variation of the deaths given the rates and the exposures.

## References

Hyndman RJ, Booth H, Yasmeen F (2013). Coherent mortality forecasting:
the product-ratio method with functional time series models.
*Demography* 50(1): 261-283.

Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
method for forecasting mortality. *Demography* 38(4): 537-549.

## See also

[`coherent_fdm()`](coherent_fdm.md), [`li_lee()`](li_lee.md).

## Examples

``` r
sp <- suppressWarnings(product_ratio())
# the published method:
if (requireNamespace("demography", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 1990:2019)
  fit <- mt_fit(coherent_fdm(order1 = 2, order2 = 1), x)
}
```
