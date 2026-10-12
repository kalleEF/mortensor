# mortensor

<!-- badges: start -->
[![R-CMD-check](https://github.com/kalleEF/mortensor/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kalleEF/mortensor/actions/workflows/R-CMD-check.yaml)
[![validation](https://github.com/kalleEF/mortensor/actions/workflows/validation.yaml/badge.svg)](https://github.com/kalleEF/mortensor/actions/workflows/validation.yaml)
<!-- badges: end -->

Flexible probabilistic mortality forecasting from data observed at different
granularities: age groups and open age intervals, single or grouped periods,
population strata and causes.

**Status:** version 0.2 is a rewrite in development (see `dev/progress.md`).
The 0.1.1 code is kept at commit 41058ba.

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
| Data at any granularity | `mt_data()`, `mt_add_deaths()`, `as_mt_data()` (StMoMo, demography), `mt_from_hmd()`, `mt_read_hmd()`, `mt_validate()`, `mt_example()` |
| Harmonisation | `mt_regroup()`, `mt_ungroup()` (native penalized composite link model, rates with exposure offset, 1D/2D), `mt_smooth()`, `mt_as_exposure()` |
| Preprocessing pipelines | `mt_prep()` with `step_ungroup()`, `step_smooth()`, `step_regroup()`, `step_grid()`, `step_shocks()` |
| Models | `rwd()`, `lee_carter()`, `tucker3()`, `arima_series()`, `ets_series()`, `hyndman_ullah()` (demography), `stmomo()` (APC, CBD, RH, M6, LC), `li_lee()`, `product_ratio()`, `apc_bayes()` (bamp), `cause_shares()` (multi-cause), `mt_ensemble()` (equal, CRPS or stacking weights), `mt_tuned()` |
| Constraints | `mt_reconcile()` (MinT, WLS, OLS; draw by draw), `mt_cohere()` |
| Forecasts | `mt_fit()`, `mt_forecast()` (simulated paths), `mt_summarise()`, `mt_deaths()` |
| Evaluation | `mt_backtest()` (rolling origin, no leakage), `mt_score()` (CRPS, WIS, coverage, PIT, MASE) |
| Derived quantities and plots | `life_table()`, `life_expectancy()`, `mt_plot_forecast()`, `mt_plot_life_expectancy()`, `mt_plot_coverage()`, ... |

## Example

```r
library(mortensor)

# simulated data: deaths in 5-year age groups and 90+, exposures by single year
x <- mt_example()
mt_validate(x)

# ungrouping is part of the model pipeline, so backtests refit it per origin
prep <- mt_prep(step_ungroup())
fit <- mt_fit(lee_carter(), x, prep = prep)
fc <- mt_forecast(fit, h = 15, n_draws = 500, seed = 1)
life_expectancy(fc, at = 65)

bt <- mt_backtest(lee_carter(), x, origins = 2004:2009, h = 10, prep = prep)
mt_score(bt)

# Human Mortality Database data for France (via demography)
fr <- as_mt_data(demography::fr.mort, years = 1950:2006)
```

Vignettes: *Getting started*, *Data at different granularities: ungrouping
HMD data*, *Forecasting and backtesting with HMD data*, *Bayesian
age-period-cohort models with bamp*.

## Development

See `CLAUDE.md` for conventions and commands, `validation/README.md` for the
validation pipeline, and `dev/implementation-plan.md` for the plan.

## Licence

GPL-3.
