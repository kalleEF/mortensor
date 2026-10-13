# Ensembles of models

A mixture of the predictive distributions of several models: each draw
of the ensemble is a draw of one component, chosen in proportion to the
weights. Mixing draws (rather than averaging quantiles) keeps the paths
coherent, so derived quantities such as life expectancy can be computed
from the ensemble draws.

## Usage

``` r
ensemble(..., weights = "equal", inner_origins = 3, inner_h = 5, top_k = NULL)
```

## Arguments

- ...:

  Model specifications; names are used as labels.

- weights:

  `"equal"`, `"crps"`, `"stack"` or a numeric vector.

- inner_origins, inner_h:

  Inner backtest for `weights = "crps"` and `weights = "stack"`, and for
  `top_k`.

- top_k:

  Keep only the `top_k` components with the lowest CRPS in the inner
  backtest on the training years (the others get weight 0); `NULL` keeps
  all.

## Value

A model specification.

## Details

Weights:

- `"equal"`;

- a numeric vector (normalised to sum to 1);

- `"crps"`: inverse mean CRPS of an inner rolling-origin backtest on the
  last `inner_origins` years of the training window, so the weights
  never use data after the forecast origin;

- `"stack"`: weights on the simplex that minimise the squared error of
  the weighted median forecasts (log scale) in the same inner backtest
  (the non-negative stacking, fitted on training data only).

Each draw of the ensemble comes from one component for all strata, so
paths keep the dependence between strata and the coherence of reconciled
components ([`reconciled()`](reconciled.md)). A component that fails in
a stratum is replaced there by the other components in proportion to
their weights; this is recorded in the diagnostics. Any model
specification can be a component, including joint models, Bayesian
models, tuned, bootstrapped and reconciled models, and other ensembles.

A mixture is wider than its components when they disagree. In backtests
on HMD data (England & Wales, France) and on the ONS deprivation data,
Lee-Carter alone was too narrow (95% coverage 0.79 to 0.81), while an
equal-weight mixture of [`lee_carter()`](lee_carter.md) and
[`rwd()`](rwd.md) had 95% coverage of 0.92 to 0.98; its CRPS ranked
second (France), third (England & Wales) and fifth of eight (ONS), close
to the best model in each case (see the vignette "Comparison with
demography, bamp and StMoMo"). In populations with few deaths per cell,
use `rwd(observation_noise = TRUE)` as the random-walk component: the
plain random walk accumulates the Poisson noise of observed rates and
its long-horizon paths become implausible.

## Examples

``` r
x <- mt_example(deaths = "single")
ens <- ensemble(lc = lee_carter(), rw = rwd())
fc <- mt_forecast(mt_fit(ens, x), h = 5, n_draws = 100, seed = 1)
```
