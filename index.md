# mortensor

Flexible probabilistic mortality forecasting from data observed at
different granularities: age groups and open age intervals, single or
grouped periods, population strata and causes.

> **Caution: experimental, vibe-coded software. Use at your own risk.**
> mortensor is written largely with AI coding assistants (“vibe-coded”).
> It has unit tests and a validation pipeline with thresholds (see
> `validation/`), but it has not been peer reviewed, its interfaces
> still change, and its results can be wrong in ways the tests do not
> catch. Check results independently before you use them in
> publications, reports or decisions.

**Status:** version 0.2 is a rewrite in development (see
`dev/progress.md`).

## Installation

``` r

# install.packages("remotes")
remotes::install_github("kalleEF/mortensor")
```

Requires R \>= 4.3.0. Hard dependencies are data.table and Matrix;
models and data adapters that build on other packages (forecast, bamp,
StMoMo, demography, ggplot2) use them when installed.

## What it does

| Layer | Functions |
|----|----|
| Data at any granularity | [`mt_data()`](reference/mt_data.md) (long or wide tables, symbols for suppressed cells, shared populations), [`mt_add()`](reference/mt_add.md), [`as_mt_data()`](reference/as_mt_data.md) (StMoMo, demography), [`mt_from_hmd()`](reference/mt_from_hmd.md), [`mt_read_hmd()`](reference/mt_from_hmd.md), [`mt_read_genesis()`](reference/mt_read_genesis.md), [`mt_read_gbe()`](reference/mt_read_gbe.md), [`mt_validate()`](reference/mt_validate.md), `ons_data`, `de_deaths` |
| Preparation | [`mt_as_exposure()`](reference/mt_as_exposure.md), [`mt_subset()`](reference/mt_subset.md), [`mt_group_levels()`](reference/mt_group_levels.md) (aggregation layers: cause groups, quintiles, regions), [`mt_relabel()`](reference/mt_relabel.md), [`mt_fill_suppressed()`](reference/mt_fill_suppressed.md) |
| Harmonisation | [`mt_regroup()`](reference/mt_regroup.md), [`mt_ungroup()`](reference/mt_ungroup.md) (penalized composite link model with exposure offset, 1D/2D; spline and GAM methods), [`mt_smooth()`](reference/mt_smooth.md) (over age, over time or both; P-splines, Whittaker, GAM), [`mt_smooth_check()`](reference/mt_smooth_check.md), [`mt_plot_smooth()`](reference/mt_smooth_check.md) |
| Preprocessing pipelines | [`mt_prep()`](reference/mt_prep.md) with [`step_ungroup()`](reference/mt_prep.md), [`step_smooth()`](reference/mt_prep.md), [`step_regroup()`](reference/mt_prep.md), [`step_grid()`](reference/mt_prep.md), [`step_shocks()`](reference/mt_prep.md), [`step_window()`](reference/mt_prep.md); models fitted to smoothed rates take their observation noise from the raw rates |
| Describing the past | [`mt_observed()`](reference/mt_observed.md), [`life_expectancy()`](reference/life_expectancy.md), [`asr()`](reference/derived.md), [`prob_dying()`](reference/derived.md) on data, [`mt_contrast()`](reference/mt_contrast.md) (gaps, ratios), [`mt_inequality()`](reference/mt_inequality.md) (SII, RII), [`standard_population()`](reference/standard_population.md) |
| Models | [`rwd()`](reference/rwd.md), [`lee_carter()`](reference/lee_carter.md) (Poisson or SVD with the adjustments and period choice of demography), [`tucker3()`](reference/tucker3.md), [`arima_series()`](reference/arima_series.md), [`ets_series()`](reference/arima_series.md), [`hyndman_ullah()`](reference/hyndman_ullah.md) (demography), [`stmomo()`](reference/stmomo.md) (APC, CBD, RH, M6, LC), [`li_lee()`](reference/li_lee.md), [`product_ratio()`](reference/li_lee.md), [`apc_bayes()`](reference/apc_bayes.md) and [`apc_select()`](reference/apc_select.md) (bamp), [`cause_shares()`](reference/cause_shares.md) |
| Wrappers for any model | [`ensemble()`](reference/ensemble.md) (weights from inner backtests, `top_k`), [`reconciled()`](reference/reconciled.md), [`calibrated()`](reference/calibrated.md) (spread fitted to the errors of the last 10 years), [`anchored()`](reference/anchored.md) (jump-off at the observed level), [`bootstrapped()`](reference/bootstrapped.md), [`tuned()`](reference/mt_tune.md) and [`mt_tune()`](reference/mt_tune.md) (full or two-stage grids) |
| Coherence | [`mt_hierarchy()`](reference/mt_hierarchy.md), [`mt_aggregate_strata()`](reference/mt_aggregate_strata.md), [`mt_reconcile()`](reference/mt_reconcile.md) and [`reconciled()`](reference/reconciled.md) (bottom-up, WLS, MinT, structural, OLS, proportional; draw by draw; non-negative; `fixed` nodes such as a trusted all-cause total), [`mt_coherence()`](reference/mt_coherence.md), [`mt_plot_coherence()`](reference/mt_coherence.md), [`mt_aggregate_draws()`](reference/mt_aggregate_draws.md), [`mt_order_strata()`](reference/mt_order_strata.md) |
| Forecasts | [`mt_fit()`](reference/mt_fit.md), [`summary()`](https://rdrr.io/r/base/summary.html), [`mt_parameters()`](reference/mt_parameters.md), [`residuals()`](https://rdrr.io/r/stats/residuals.html), [`mt_forecast()`](reference/mt_forecast.md) (simulated paths), [`mt_jump_off()`](reference/mt_jump_off.md), [`mt_summarise()`](reference/mt_draws.md), [`mt_deaths()`](reference/mt_deaths.md), [`mt_disaggregate()`](reference/mt_disaggregate.md) with [`mt_shape_rates()`](reference/mt_shape_rates.md) and [`mt_disaggregate_check()`](reference/mt_disaggregate_check.md) |
| Evaluation | [`mt_backtest()`](reference/mt_backtest.md) (rolling origin, no leakage; rates, life expectancy, ASR), [`mt_score()`](reference/mt_score.md) (CRPS, WIS, coverage, PIT, MASE) |
| Decomposition | [`mt_decompose()`](reference/mt_decompose.md) (differences by age, cause and subgroup: Arriaga, Horiuchi, stepwise replacement, Kitagawa, and Arriaga meets Kitagawa for subgroup composition; on observed data, forecast draws or both), [`mt_plot_decomposition()`](reference/mt_plot_decomposition.md) |
| Results and plots | [`life_table()`](reference/life_table.md), [`life_expectancy()`](reference/life_expectancy.md), [`prob_dying()`](reference/derived.md), [`asr()`](reference/derived.md), [`yll()`](reference/derived.md), [`mt_contrast()`](reference/mt_contrast.md), [`mt_inequality()`](reference/mt_inequality.md), [`mt_plot_forecast()`](reference/mt_plot.md), [`mt_plot_life_expectancy()`](reference/mt_plot.md), [`mt_plot_scores()`](reference/mt_plot_more.md), [`mt_plot_shares()`](reference/mt_plot_more.md), [`theme_mortensor()`](reference/mt_plot_style.md), [`mt_export()`](reference/mt_export.md), [`as_demogdata()`](reference/as_demogdata.md), [`as_stmomodata()`](reference/as_demogdata.md) |

## Function map

The objects of an analysis and the functions that connect them. The
[function map](vignettes/articles/function-map.Rmd) (on the website:
Articles, Function map) shows every function with its purpose, the
possible workflows and what the package does not cover.

``` mermaid
flowchart LR
  T["`**your tables**
  deaths, populations`"]
  D["`**mt_data**
  any resolution, strata`"]
  S["`**model specification**
  lee_carter(), li_lee(), rwd(),
  ensemble(), reconciled(), calibrated()`"]
  F["`**mt_fit**`"]
  P["`**mt_draws**
  simulated future rates`"]
  R["`**results**
  life expectancy, ASR, deaths,
  gaps, plots, CSV`"]
  B["`**backtest scores**`"]
  T -- "mt_data()" --> D
  D -- "mt_as_exposure(), mt_subset(), mt_group_levels(),
  mt_fill_suppressed(), mt_regroup(), mt_validate()" --> D
  S --> F
  D -- "mt_fit()" --> F
  F -- "mt_forecast()" --> P
  P -- "life_expectancy(), asr(), mt_contrast(),
  mt_decompose(), mt_plot_*(), mt_export()" --> R
  D -- "mt_backtest(), mt_score()" --> B
  D -- "mt_observed()" --> R
```

## Example

``` r

library(mortensor)

# deaths and populations by sex, age group and deprivation decile, England
# 2001-2019 (ONS; Open Government Licence)
x <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile"),
             age = "age_group", value = c(deaths = "all_causes", population = "population"))
x <- mt_as_exposure(x, type = "midyear")
# deciles to quintiles: an aggregation layer with its own name
x <- mt_group_levels(x, "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)),
                     name = "quintile")
mt_validate(x)

fit <- mt_fit(li_lee(common_by = "sex"), x)
fc <- mt_forecast(fit, h = 15, n_draws = 500, seed = 1)
life_expectancy(fc, at = 0)

bt <- mt_backtest(li_lee(common_by = "sex"), x, origins = 2008:2014, h = 5)
mt_score(bt)
```

Vignettes: *Getting started* (ONS deprivation data), *Your own data,
from CSV to report*, *Forecasting mortality by deprivation quintile*,
*Causes of death by deprivation*, *Coherent forecasts by sex, region and
cause* (German Laender), *Coherent forecasts: what reconciliation
changes*, *Decomposing differences in life expectancy*, *Model
reference*, *Data at different granularities* and *Forecasting and
backtesting* (HMD data), *Bayesian age-period-cohort models with bamp*,
*Comparison with demography, bamp and StMoMo*, *Frequently asked
questions*.

## Development

See `CLAUDE.md` for conventions and commands, `validation/README.md` for
the validation pipeline, and `dev/implementation-plan.md` for the plan.

## Licence

GPL-3.
