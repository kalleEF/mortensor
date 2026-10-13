# mortensor 0.1.9000 (development of 0.2.0)

Clean break: the 0.1.x interface is replaced (see `dev/revamp-plan.md`).

* Data layer: `mt_axis` (`single_ages()`, `age_groups()`, `single_years()`,
  `parse_age_labels()`), `mt_data()` with tables at their own resolution,
  `mt_grid()`/`mt_set_grid()`, `mt_array()`, `mt_window()`,
  `mt_granularity()`, `mt_validate()` (including hierarchy sums),
  `mt_regroup()`, `mt_operator()`, `mt_as_exposure()`, `as_mt_data()` for
  data frames, `StMoMoData` and `demogdata`, and the simulated `mt_example()`.
* Native penalized composite link model and `mt_ungroup()`: deaths are
  ungrouped as rates with exposures as offset; the open interval is expanded
  to omega; 1D (default) and 2D; group totals preserved.
* `mt_smooth()` and preprocessing pipelines (`mt_prep()`, `step_*()`) that
  backtests refit on every training window.
* Interval-aware `life_table()` and `life_expectancy()` (also path by path
  on forecast draws).
* Model contract (`mt_fit()`, `mt_forecast()`, `mt_draws`) with models
  `rwd()`, `lee_carter()` (Poisson and SVD; agrees with StMoMo and
  demography), `tucker3()` (native HOOI; log and CoDa), `arima_series()`,
  `ets_series()`, and `mt_ensemble()`.
* `apc_bayes()`: seamless wrapper of the bamp package (Bayesian
  age-period-cohort models) with convergence checks and draws on the grid.
* Plots (`mt_plot_*()`, ggplot2), HMD readers (`mt_from_hmd()`,
  `mt_read_hmd()`), and four vignettes using HMD data.
* Ports of the 0.1 models and methods: `hyndman_ullah()`, `stmomo()` (APC, CBD, RH,
  M6, LC), `li_lee()`, `product_ratio()`, `cause_shares()` (Tucker4-ILR
  multi-cause forecasts), `mt_reconcile()` (MinT, WLS, OLS),
  `mt_cohere()`, stacking weights in `mt_ensemble()`, Whittaker and GAM
  smoothing,
  Kannisto closure in `life_table()`, monotone splitting of deaths, and the
  simulated multi-cause dataset `mt_example_causes()`.
* Grouped hierarchies for coherent forecasts (for example sex x region x
  deprivation quintile): `mt_hierarchy()`, `mt_aggregate_strata()`,
  multi-variable hierarchy checks in `mt_validate()`, `mt_reconcile()` on
  the death scale with WLS, MinT, structural, OLS or bottom-up
  reconciliation, an empirical copula for separately fitted strata and
  non-negativity; `mt_reconciled()` for reconciliation inside fits and
  backtests; `mt_aggregate_draws()`; `li_lee(common_by = "sex")`;
  simulated `mt_example_regions()` and a vignette.
* Rate-based spline ungrouping (`mt_ungroup(method = "hyman_kannisto")`,
  `"natural_gompertz"`, `"gam"`), `mt_disaggregate()` for forecast draws,
  `mt_bootstrap()` (parameter and preprocessing uncertainty by Poisson
  resampling and refitting) and information-criterion rank selection in
  `tucker3(ranks = "bic")`.
* `ons_data`: deaths by cause and populations by sex, age group and IMD
  decile, England 2001-2019 (ONS ad hoc release 3064, Open Government
  Licence), with the cause groups `ons_cause_groups`; built by
  `data-raw/ons_data.R`. `mt_read_ons()` reads it directly.
* `mt_read_ons()` (ONS Excel tables by sex, age and IMD decile, optional
  regions and cause groups), `apc_select()` (DIC-based choice of random-walk
  orders for `apc_bayes()`), and the plots `mt_plot_tucker()`,
  `mt_plot_shares()` and `mt_plot_scores()`.
* `mt_tune()` and `mt_tuned()`; `prob_dying()`, `asr()` (ESP 2013),
  `yll()`, `mt_export()`.
* `mt_backtest()` (rolling origin, scored at the observed resolution) and
  `mt_score()` (CRPS, WIS, coverage, PIT, MASE).
* Validation pipeline in `validation/` and GitHub workflows for R CMD check
  (R 4.3 to devel), coverage, lint and validation.
* Dependencies reduced to data.table and Matrix (plus base packages).

# mortensor 0.1.1

Maintenance release before the 0.2 rewrite (see `dev/implementation-plan.md`).

* `%||%` is imported from rlang, so the package works on R < 4.4 (review D1).
* `forecast()` and `accuracy()` are the generics from the generics package
  (re-exported); methods dispatch correctly after `library(forecast)` (D2).
* Licence changed to GPL-3.
* Undefined global variables declared; stray files removed from the repository.
* `dev/baseline/` records point forecasts of every model as a regression
  baseline for the rewrite.
