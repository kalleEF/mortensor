# Ensembles of models

A mixture (linear pool, Gneiting and Ranjan 2013) of the predictive
distributions of several models: each draw of the ensemble is a draw of
one component, chosen with the weights as probabilities. Mixing draws
keeps the paths coherent, so derived quantities such as life expectancy
can be computed from the ensemble draws.

## Usage

``` r
ensemble(
  ...,
  weights = "equal",
  inner_origins = 3,
  inner_h = 5,
  top_k = NULL,
  inner_n_draws = 100
)
```

## Arguments

- ...:

  Model specifications; names are used as labels.

- weights:

  `"equal"`, `"crps"`, `"stack"` or a numeric vector.

- inner_origins, inner_h:

  Inner backtest for `weights = "crps"` and `weights = "stack"`, and for
  `top_k`: the number of origins and the horizon of each inner forecast.

- top_k:

  Keep only the `top_k` components with the lowest CRPS in the inner
  backtest on the training years (the others get weight 0); `NULL` keeps
  all.

- inner_n_draws:

  Draws per inner forecast of the inner backtest.

## Value

A model specification.

## Details

Weights:

- `"equal"`;

- a numeric vector (normalised to sum to 1);

- `"crps"`: proportional to the inverse mean CRPS of each component in
  an inner rolling-origin backtest on the last years of the training
  window (`inner_origins`, `inner_h`), so the weights never use data
  after the forecast origin. This is a package heuristic (reference
  unknown), analogous to the inverse mean squared error weights of Bates
  and Granger (1969) with the mean CRPS in place of the mean squared
  error;

- `"stack"`: weights that minimise the squared error of the weighted
  mean of the components' median forecasts (log rates) in the same inner
  backtest, constrained to be non-negative and to sum to 1. This is
  stacking of point forecasts (Breiman 1996; for mortality, Kessy,
  Sherris, Villegas and Ziveyi 2022). Using these weights as the
  probabilities of the mixture is a deviation: they are optimal for the
  linear combination of the medians, not for the mixture distribution,
  whose median is not the weighted mean of the medians.

Each draw of the ensemble comes from one component for all strata, so
paths keep the dependence between strata and the coherence of reconciled
components ([`reconciled()`](reconciled.md)). A component that fails in
a stratum is replaced there by the other components in proportion to
their weights; this is recorded in the diagnostics. Any model
specification can be a component, including joint models, Bayesian
models, tuned, bootstrapped and reconciled models, and other ensembles.
The components must forecast the same strata: with
[`reconciled()`](reconciled.md) components, add the aggregates of the
hierarchy to the data first
([`mt_aggregate_strata()`](mt_aggregate_strata.md)). Weights from the
inner backtest (`"crps"`, `"stack"`, `top_k`) compare the components on
the cells that all of them forecast.

A mixture is wider than its components when they disagree, which often
corrects the too narrow intervals of single models; a linear pool of
different calibrated components is overdispersed (Gneiting and Ranjan
2013). Its CRPS is never larger than the weighted mean of the CRPS of
its components (the CRPS is convex in the forecast distribution). The
vignette "Choosing and combining models" compares ensembles and single
models in a backtest. In populations with few deaths per cell, the plain
random walk accumulates the Poisson noise of observed rates and its
long-horizon paths become implausible; the local level model with drift,
the random walk observed with noise, is
`arima_series(order = c(0, 1, 1), drift = TRUE)` (its reduced form,
Harvey 1989).

## References

Gneiting T, Ranjan R (2013). Combining predictive distributions.
*Electronic Journal of Statistics* 7: 1747-1782.

Bates JM, Granger CWJ (1969). The combination of forecasts. *Operational
Research Quarterly* 20(4): 451-468.

Breiman L (1996). Stacked regressions. *Machine Learning* 24(1): 49-64.

Kessy SR, Sherris M, Villegas AM, Ziveyi J (2022). Mortality forecasting
using stacked regression ensembles. *Scandinavian Actuarial Journal*
2022(7): 591-626.

Harvey AC (1989). *Forecasting, Structural Time Series Models and the
Kalman Filter*. Cambridge University Press.

## Examples

``` r
x <- mt_example(deaths = "single")
ens <- ensemble(lc = lee_carter(), rw = rwd())
fc <- mt_forecast(mt_fit(ens, x), h = 5, n_draws = 100, seed = 1)
```
