# Li-Lee coherent multi-population model

Joint model for several strata (for example women and men, or regions)
whose forecasts do not diverge in the long run (Li and Lee 2005):
`log m[x,t,i] = a[x,i] + B[x] K[t] + b[x,i] k[t,i]`. The common factor
`B K` is a Lee-Carter model of the pooled rates (sum of deaths over sum
of exposures, SVD with `K` refitted to the total deaths of each year; a
year where the search finds no such `K` keeps the SVD index, as with
`adjust = "dt"` in [`lee_carter()`](lee_carter.md)); the stratum factors
`b k` are the first `ratio_order` components of an SVD of the residuals
of each stratum.

## Usage

``` r
li_lee(ratio_order = 1, common_by = NULL, jump_off = c("fitted", "observed"))
```

## Arguments

- ratio_order:

  Number of stratum-specific factors. It is capped by the number of age
  groups and the number of years minus 2, with a warning (class
  `mt_warning_argument`) that names the requested and the fitted number;
  [`mt_parameters()`](mt_parameters.md) reports the fitted number
  (`ratio_order`).

- common_by:

  Stratum variable(s) defining groups with their own common factor.

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

With `common_by` (for example `"sex"`), `li_lee()` fits one common
factor per level of that stratum variable (women and men each with their
own trend), and the other strata (for example regions and deprivation
quintiles) deviate from their group's common factor. If the data carry a
hierarchy (see [`mt_aggregate_strata()`](mt_aggregate_strata.md)), the
common factors are estimated from the bottom strata only, so aggregates
are not counted twice.

The common period index follows a random walk with drift (with drift
uncertainty) and the stratum-specific indexes `k[t,i]` follow AR(1)
processes, so the differences between strata stay bounded. Two choices
of this package differ from Li and Lee (2005), who model `k[t,i]` as an
AR(1) process with a constant:

- the AR(1) processes have no constant (intercept): the coefficient is
  the least-squares estimate through the origin, and the stratum indexes
  revert to zero, their mean over the fitting period;

- the AR(1) coefficient is bounded to `|phi| <= 0.99`, a
  package-specific safeguard against unit roots and explosive deviations
  (reference unknown). A warning (class `mt_warning_li_lee`) names the
  strata where the bound binds; their deviations then revert very
  slowly, and the assumption of coherence may not suit them.

Draws include the innovations of all indexes and the drift uncertainty
of the common index; `a`, `B`, `b` and the AR(1) coefficients are fixed.
Predictive draws (`mt_forecast(predictive = TRUE)`) add Poisson
variation of the deaths given the rates and the exposures, and no other
observation noise. If intervals are too narrow on real data, that is a
property of the model (see [`recalibrated()`](recalibrated.md),
[`ensemble()`](ensemble.md), or overdispersed models such as
`stanmomo(family = "nb")`).

[`product_ratio()`](product_ratio.md), a Lee-Carter version of the
product-ratio method, is deprecated; use
[`coherent_fdm()`](coherent_fdm.md).

## References

Li N, Lee R (2005). Coherent mortality forecasts for a group of
populations: an extension of the Lee-Carter method. *Demography* 42(3):
575-594.
[doi:10.1353/dem.2005.0021](https://doi.org/10.1353/dem.2005.0021)

Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
*Journal of the American Statistical Association* 87(419): 659-671.

Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
method for forecasting mortality. *Demography* 38(4): 537-549.

## See also

[`coherent_fdm()`](coherent_fdm.md) for the product-ratio method.

## Examples

``` r
x <- mt_example(deaths = "single")
fc <- mt_forecast(mt_fit(li_lee(), x), h = 10, n_draws = 100, seed = 1)
```
