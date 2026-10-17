# Models with reconciled forecasts

Wraps a model specification so that its forecasts are reconciled over a
hierarchy every time it is forecast, also inside
[`mt_backtest()`](mt_backtest.md). If the training data contain only the
bottom strata, the aggregate strata are added with
[`mt_aggregate_strata()`](mt_aggregate_strata.md) before fitting, so
that the base model forecasts every node; exposures for the forecast
(`exposure` of [`mt_forecast()`](mt_forecast.md)) can then be given for
the bottom strata only, and those of the aggregates are their sums.
Predictive draws (used by backtests) are coherent as well: the
observation noise of the base model (Poisson and extra-Poisson variation
for models of smooth rates, none for models of observed rates such as
[`rwd()`](rwd.md); for [`ensemble()`](ensemble.md) and
[`stratified()`](stratified.md) models that of the component or group
each draw comes from) is drawn for the bottom strata only, independently
per cell, and the aggregates are their sums.

## Usage

``` r
reconciled(
  spec,
  hierarchy,
  method = c("wls", "mint", "structural", "ols", "bottom_up", "proportional"),
  shrink = "auto",
  dependence = c("auto", "empirical", "none"),
  fixed = NULL
)
```

## Arguments

- spec:

  Base model specification (fitted to all strata).

- hierarchy:

  An [`mt_hierarchy()`](mt_hierarchy.md).

- method, shrink:

  See [`mt_reconcile()`](mt_reconcile.md).

- dependence:

  `"auto"` (default: `"none"` for models fitted jointly across strata,
  whose draws already carry their dependence, such as
  [`li_lee()`](li_lee.md) and [`tucker3()`](tucker3.md); `"empirical"`
  for models fitted per stratum), `"empirical"` or `"none"`. With
  `"empirical"` the paths of the strata are reordered before
  reconciliation (see [`mt_reconcile()`](mt_reconcile.md)), so that
  aggregates keep the observed co-movement of the strata; the
  observation noise of predictive draws is independent across cells, as
  for models without a hierarchy.

- fixed:

  Strata whose forecasts are kept (see
  [`mt_reconcile()`](mt_reconcile.md)).

## Value

A model specification.

## Details

Whether coherence costs accuracy, and which method and base model work
best, depends on the data: compare them in a backtest at every level of
the hierarchy (`mt_score_table(..., hierarchy = )`), as in the vignette
"Coherent forecasts: totals and their parts".

## Examples

``` r
x <- mt_example_regions(years = 2000:2019, regions = 3)
h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
fit <- mt_fit(reconciled(lee_carter(), h, method = "wls"), mt_aggregate_strata(x, h))
fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
```
