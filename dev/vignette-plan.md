# Vignette restructure (2026-10-17)

Owner's review of the documentation: too many comparisons and concepts at
once; Lee-Carter appears to be the core of the package although the aim is
any model, implemented or wrapped, and ensembles of any combination;
helpers defined in vignettes and hand-made result tables should be
package functions; leftover "plot of chunk" captions; scores such as CRPS
need clearer explanations.

## Principles

* One topic per concept vignette. A simple model (Lee-Carter, introduced
  as one example among many) keeps the focus on the topic.
* Visual comparisons first, tables second; every number in the text comes
  from the output.
* Concepts are explained where they are used (CRPS, coverage, PIT,
  jump-off, reconciliation), in plain words with a figure where possible.
* Applied vignettes show the whole pipeline once, with one decision per
  step and a reference to the concept vignette that explains it.
* Repeated vignette code becomes package functions.

## Structure (pkgdown Articles)

| Group | Vignette | Content | Source |
|---|---|---|---|
| Start | Getting started | data, one model, forecast, life expectancy, plot | `mortensor.Rmd`, shortened |
| Start | Your own data | reading tables, suppression, exposures, checks | `own-data.Rmd` |
| Concepts | Ungrouping age groups | methods compared against known single years; diagnostics (profile, errors with noise band, stability over time); grouped exposures; effect on a forecast | `ungrouping.Rmd` (was `granularity-hmd.Rmd`) |
| Concepts | Smoothing | methods over age, time and both; how much each smooths (figures); choice of the smoothing parameter; effect on one forecast; smoothing inside the fit | new |
| Concepts | Jump-off | fitted, observed and anchored start; the step in the first forecast year; effect on a forecast | new |
| Concepts | Backtesting and scores | rolling origins (figure), CRPS, interval score, coverage, PIT, MASE explained with figures; comparing models | new, from `forecasting-hmd.Rmd` |
| Concepts | Tuning and calibration | choosing settings by forecast accuracy (hold-out, rolling origins); calibrating the width of intervals; benchmark results | new |
| Concepts | Coherent forecasts | why totals and parts disagree; bottom-up, WLS, MinT, fixed totals, proportional; one figure per idea | `coherence.Rmd`, simplified |
| Concepts | Decomposition | age, cause and composition contributions; observed and forecast | `decomposition.Rmd`, shortened |
| Models | Choosing and combining models | the implemented and wrapped models side by side on one data set (forecasts and backtest), ensembles of any combination, weights | new |
| Models | Model reference | assumptions and uncertainty of every model | `models.Rmd` |
| Models | Bayesian APC with bamp | | `bayesian-apc.Rmd` |
| Models | Fidelity of the wrappers | demography, bamp, StMoMo against their native calls | `comparison.Rmd` |
| Applied | Inequality in England | ONS IMD quintiles: life expectancy gaps and cause-specific deaths, full pipeline | `deprivation.Rmd`, `causes.Rmd` |
| Applied | Differences between German Laender | Destatis: suppression, hierarchy, pooling, reconciliation, East-West gap | `regions.Rmd` |
| Help | FAQ, function map | | existing |

Removed after their content moved: `forecasting-hmd.Rmd`,
`deprivation.Rmd`, `causes.Rmd`, `regions.Rmd`.

## New functions (from vignette code)

* `mt_compare()`: backtests of a named list of specifications (each with
  its own preprocessing) on the same origins; returns one backtest with
  the names as models.
* `mt_score_table()`: scores of several targets side by side (rates and
  life expectancy), rounded for reading.
* `mt_compare_forecasts()`: one quantity (life expectancy, ASR, deaths) of
  several forecasts or approaches side by side, as "median (lower-upper)".
* `mt_plot_life_expectancy()` accepts a named list of forecasts.
* `mt_plot_backtest_design()`: the origins, training and test windows of a
  backtest.
* `mt_plot_jump_off()`: observed, fitted and forecast rates around the
  jump-off year.
* `mt_plot_tuning()`: scores of a tuning grid with standard errors.
* `mt_plot_contrast()`: gaps between strata over time, for one or several
  forecasts with the observed gaps.
* `mt_plot_ungroup()` (profile, time, error) and `mt_ungroup_check()`.
* `mt_backtest(contrast = )`: backtests of gaps; `mt_score_table(hierarchy = )`.
* `stratified()`: different models for groups of strata.

## Status (2026-10-17)

All vignettes of the structure above are written and knitted; the
applied ones choose every step by a backtest (model by the gap scores,
smoothing, jump-off, calibration only where they helped).
`mt_plot_jump_off()` was not needed: `mt_plot_compare(periods = )` shows
the jump-off.

## Open decisions for the owner

* StanMoMo (Bayesian mortality models in Stan, with stacking of model
  forecasts) as a further wrapped model, or StMoMo only (wrapped).
* Renaming `calibrated()` to `recalibrated()`.
