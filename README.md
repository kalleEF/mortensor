# mortensor

<!-- badges: start -->
[![R-CMD-check](https://github.com/kalleEF/mortensor/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kalleEF/mortensor/actions/workflows/R-CMD-check.yaml)
[![validation](https://github.com/kalleEF/mortensor/actions/workflows/validation.yaml/badge.svg)](https://github.com/kalleEF/mortensor/actions/workflows/validation.yaml)
<!-- badges: end -->

Flexible probabilistic mortality forecasting from data observed at different
granularities: age groups and open age intervals, single or grouped periods,
population strata and causes.

**Status:** version 0.2 is a rewrite in development (see `dev/progress.md`).

## Installation

```r
# install.packages("remotes")
remotes::install_github("kalleEF/mortensor")
```

Requires R >= 4.3.0. Hard dependencies are data.table and Matrix; models and
data adapters that build on other packages (forecast, bamp, StMoMo,
demography, ggplot2) use them when installed.

## What it does

| Layer | Functions |
|---|---|
| Data at any granularity | `mt_data()` (long or wide tables, symbols for suppressed cells, shared populations), `mt_add_deaths()`, `as_mt_data()` (StMoMo, demography), `mt_from_hmd()`, `mt_read_hmd()`, `mt_read_genesis()`, `mt_read_gbe()`, `mt_validate()`, `ons_data`, `de_deaths` |
| Preparation | `mt_subset()`, `mt_group_levels()` (aggregation layers: cause groups, quintiles, regions), `mt_fill_suppressed()`, `mt_as_exposure()` |
| Harmonisation | `mt_regroup()`, `mt_ungroup()` (penalized composite link model with exposure offset, 1D/2D; spline methods), `mt_smooth()`, `mt_as_exposure()`, `mt_disaggregate()` |
| Preprocessing pipelines | `mt_prep()` with `step_ungroup()`, `step_smooth()`, `step_regroup()`, `step_grid()`, `step_shocks()` |
| Models | `rwd()`, `lee_carter()`, `tucker3()`, `arima_series()`, `ets_series()`, `hyndman_ullah()` (demography), `stmomo()` (APC, CBD, RH, M6, LC), `li_lee()`, `product_ratio()`, `apc_bayes()` and `apc_select()` (bamp), `cause_shares()`, `mt_ensemble()`, `mt_tuned()`, `mt_bootstrap()` |
| Coherence | `mt_hierarchy()`, `mt_aggregate_strata()`, `mt_reconcile()` and `mt_reconciled()` (bottom-up, WLS, MinT, structural, OLS; draw by draw), `mt_aggregate_draws()`, `mt_cohere()` |
| Forecasts | `mt_fit()`, `mt_forecast()` (simulated paths), `mt_summarise()`, `mt_deaths()` |
| Evaluation | `mt_backtest()` (rolling origin, no leakage), `mt_score()` (CRPS, WIS, coverage, PIT, MASE) |
| Derived quantities and plots | `life_table()`, `life_expectancy()`, `prob_dying()`, `asr()`, `yll()`, `mt_plot_forecast()`, `mt_plot_scores()`, `mt_plot_shares()`, ... |

## Example

```r
library(mortensor)

# deaths and populations by sex, age group and deprivation decile, England
# 2001-2019 (ONS; Open Government Licence)
x <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile"),
             age = "age_group", value = c(deaths = "all_causes", population = "population"))
x <- mt_as_exposure(x)
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

Vignettes: *Getting started* (ONS deprivation data), *Causes of death by
deprivation*, *Coherent forecasts by sex, region and cause* (German
Laender), *Data at
different granularities* and *Forecasting and backtesting* (HMD data),
*Bayesian age-period-cohort models with bamp*, *Comparison with demography,
bamp and StMoMo*, *Frequently asked questions*.

## Development

See `CLAUDE.md` for conventions and commands, `validation/README.md` for the
validation pipeline, and `dev/implementation-plan.md` for the plan.

## Licence

GPL-3.
