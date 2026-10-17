# Models with reconciled forecasts

Wraps a model specification so that its forecasts are reconciled over a
hierarchy every time it is forecast, also inside
[`mt_backtest()`](mt_backtest.md). If the training data contain only the
bottom strata, the aggregate strata are added with
[`mt_aggregate_strata()`](mt_aggregate_strata.md) before fitting, so
that the base model forecasts every node. Predictive draws (used by
backtests) are coherent as well: Poisson and extra-Poisson variation are
drawn for the bottom strata only and the aggregates are their sums.

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
  `"empirical"` the same training years order the paths before
  reconciliation and the predictive noise of the bottom strata, so
  aggregates of predictive draws keep the observed co-movement.

- fixed:

  Strata whose forecasts are kept (see
  [`mt_reconcile()`](mt_reconcile.md)).

## Value

A model specification.

## Details

On [`mt_example_regions()`](mt_example_regions.md) (4 regions, 90
strata, origins 2005 and 2008, h = 10), `li_lee(common_by = "sex")` with
`method = "bottom_up"` was the most accurate coherent forecast (CRPS
0.0425 bottom, 0.0269 aggregates; 95% coverage 0.95), and Lee-Carter per
stratum with `"wls"` reconciliation was calibrated (0.95-0.97) but less
accurate (0.0477, 0.0292). On the ONS and German regional data,
reconciliation changed the accuracy of the bottom strata little and
lowered the coverage of the totals (see the vignette "Coherent
forecasts: what reconciliation changes").

## Examples

``` r
x <- mt_example_regions(years = 2000:2019, regions = 3)
h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
fit <- mt_fit(reconciled(lee_carter(), h, method = "wls"), mt_aggregate_strata(x, h))
fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
```
