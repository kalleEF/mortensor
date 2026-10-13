# mortensor 0.1.9000 (development of 0.2.0)

Clean break: the 0.1.x interface is replaced (see `dev/revamp-plan.md`).

## Smoothing, jump-off and fixed totals (methods of the earlier UK analysis)

* Smoothing over time: `mt_smooth(along = "period", method = "whittaker")`
  (one smoothing parameter per stratum, chosen by GCV pooled over ages)
  and over the age-period surface with a quasi-Poisson GAM,
  `mt_smooth(dims = 2, method = "gam")` (`te(age, year)`, about one second
  per stratum; the 2D P-spline takes minutes). `mt_smooth_check()` reports
  how far the smoother moved the data (median absolute log difference,
  dispersion of the raw deaths around the smoothed ones, variance ratio of
  year-to-year changes, correlation); `mt_plot_smooth()` shows observed
  points against smoothed lines.
* Models fitted to smoothed data keep the raw deaths (`fit$raw_D`):
  observation noise of predictive draws, residuals, observed rates, the
  MASE scale of backtests and the empirical copula of reconciliation refer
  to the data before smoothing. Before, noise was estimated from the
  smoothed rates, which made intervals too narrow.
* Smoothing over time narrows the intervals of models with a period index;
  combined with `calibrated()` it gave the best calibrated forecasts for
  England & Wales men (Lee-Carter on the GAM surface, calibrated: CRPS
  0.0646, 95% coverage 0.91 against 0.0687 and 0.75). On the ONS quintiles,
  whose noise is close to Poisson, smoothing did not help.
* `anchored()` moves the start of any model's forecast towards the observed
  level of the last `years` (default 3) per age and stratum, with empirical
  Bayes shrinkage of sparse cells and a linear fade over `decay` years;
  ensembles are anchored per component and reconciled models stay
  coherent. `mt_jump_off()` reports the step from the last observed and the
  last fitted year into the first forecast year.
* Reconciliation with fixed nodes: `mt_reconcile(fixed = )` and
  `reconciled(fixed = )` keep trusted forecasts (for example the total, or
  all causes) and let the other strata absorb the gap. Bottom values that a
  projection would make negative are held at zero and the rest projected
  again (a simple form of non-negative reconciliation), so fixed nodes and
  aggregates are kept; before, zeroed values could move them.
  `mt_coherence()` and `mt_plot_coherence()` show how far forecasts are
  from adding up.
* Backtests score derived quantities: `mt_backtest(quantities =
  c("life_expectancy", "asr"))` adds rows with `target`, and
  `mt_score(target = )` scores them; `mt_plot_scores(by = "model")` draws
  a bar chart of models.
* `lee_carter(adjust = c("dt", "dxt", "e0", "none"))` re-estimates the
  period index as in `demography::lca()`, and `choose_period = TRUE` fits
  the longest recent period with a linear index (Booth, Maindonald and
  Smith 2002); `ensemble(top_k = )` keeps the best components of the inner
  backtest; `mt_tune(strategy = "two_stage")` searches a coarse grid first.
* Disaggregation: `tail` (`"kannisto"` or `"gompertz"`) and `tail_width`
  for the open interval in `mt_ungroup()` and `mt_disaggregate()`,
  year-specific shapes from `mt_shape_rates()`, `mt_disaggregate_check()`
  and `mt_plot_disaggregation()`.
* Vignettes: "Forecasting and backtesting with HMD data" (chosen period,
  anchoring, smoothing with calibration, scores of life expectancy),
  "Forecasting mortality by deprivation quintile" (smoothing check,
  jump-off, backtests of anchoring and smoothing), "Coherent forecasts"
  (coherence plot, a fixed total), "Causes of death by deprivation"
  (causes reconciled to a fixed all-cause forecast). New validation check
  `smoothing_jumpoff`.

## Calibration, deprivation and coherence vignettes

* `calibrated()` wraps any model specification (including `ensemble()`,
  `reconciled()` and `bootstrapped()`) and recalibrates the spread of its
  forecast draws, optionally also the drift, on the last `years` (default
  10) of the training data. The base model is refitted at several origins
  in that window. A horizon-dependent spread factor (and drift correction)
  is chosen by minimum CRPS, or by the interval score of one central
  interval (`target = "interval"`). Observation noise is modelled as in
  `mt_backtest()`, and draws of an ensemble get the noise of their
  component. Reconciled models are calibrated at the bottom level and
  stay coherent. `by` gives separate factors for groups of strata.
  `summary()` and `mt_parameters()` report the factors.
* Backtests (England & Wales men, origins 1996-2001, h = 10): calibration
  raised the 95% coverage of Lee-Carter from 0.75 to 0.88 and lowered its
  CRPS by 4%. An already calibrated ensemble was narrowed (0.95 to 0.90)
  at the same CRPS. Drift correction overfitted on the ONS data and is off
  by default. New validation check `calibration_ew`.
* Ensemble forecasts record the component of each draw (`$component`).
* Plot layout of the owner's earlier analysis as the package standard:
  `theme_mortensor()`, `scale_colour_mortensor()`, `scale_fill_mortensor()`,
  `mt_palette()` (ordered red-yellow-green scale for levels such as Q1-Q5,
  distinct colours otherwise). Forecast plots show observed points, fitted
  lines (stored in the draws as `$fitted`), dashed forecast medians with
  interval ribbons and a dotted line before the first forecast year;
  `mt_plot_forecast()` has panels of age by stratum and new arguments
  `colour` and `log_values`; `levels` defaults to 0.95. `fitted()` works
  for ensembles (weighted geometric mean of the components).
* `mt_inequality()` weights observed years by their own populations (it
  used the last year's for every year).
* Documentation review: corrected vignette titles and cross-references,
  Model reference (series models, calibration, sources of uncertainty of
  every model), FAQ, `mt_validate()` checks, help pages of `ensemble()`,
  `reconciled()`, `cause_shares()` and the derived quantities.
* Vignettes "Forecasting mortality by deprivation quintile" (ONS:
  inequality indicators, backtests by quintile, calibration, gap
  forecasts by model, fitting period) and "Coherent forecasts: what
  reconciliation changes" (ONS sex x quintile with totals: incoherence of
  separate forecasts, each reconciliation method's effect on forecasts,
  intervals and backtest accuracy per level, dependence between strata).
  The regional vignette adds forecasts for every Land and compares the
  reconciliation methods on them.

## Changes after two applied reviews and stress tests

* Renamed for clarity: model wrappers are nouns like the other model
  constructors, `ensemble()`, `reconciled()`, `tuned()`, `bootstrapped()`
  (were `mt_ensemble()`, `mt_reconciled()`, `mt_tuned()`, `mt_bootstrap()`);
  `mt_order_strata()` (was `mt_cohere()`); `mt_parameters()` replaces
  `lc_parameters()` and `apc_effects()`; `mt_window()` is internal
  (`mt_subset(periods = )`).
* Hierarchies: additivity per variable
  (`mt_hierarchy(region = "DE", cause = list(total = "all", additive =
  "rates"))`); exposures of crossed region-by-cause hierarchies are summed
  over regions and taken once over causes in aggregation, validation and
  reconciliation. Before, all-cause exposures were multiplied by the number
  of causes.
* `mt_fill_suppressed()` takes the bounds of suppressed counts (`lower`,
  `upper`), also without published totals; it stops when no total of the
  hierarchy is in the data and leaves uncovered cells missing (before, they
  kept the exposures as starting values).
* Inspection: `summary()` of fits (model, sources of uncertainty, fit),
  `mt_parameters()` for Lee-Carter, Li-Lee, product-ratio, random walks and
  Bayesian APC models (with R-hat), `fitted()`, `residuals()`,
  `mt_plot_residuals()`.
* Results: `mt_observed()` and observed indicators (`life_expectancy()`,
  `asr()`, `prob_dying()` accept data, with Poisson intervals),
  `mt_contrast()` (differences and ratios between strata, path by path),
  `mt_inequality()` (slope and relative index of inequality),
  `life_expectancy(paths = TRUE)`, `asr(per = )` with `ages`, `standard`
  and `per` columns, standards `"esp1976"` and `"who2000"`, `mt_export()`
  for any table.
* Future exposures (`mt_deaths()`, `mt_forecast()`, `mt_reconcile()`,
  `mt_aggregate_draws()`) are matched by age, year and stratum names and
  can be `mt_data` objects; unnamed arrays give a warning.
* Data input: broader age labels (`"1 to 4"`, `"85 years and over"`,
  `"90 & over"`, German labels), financial years (`"2001/02"`), total rows
  left out and unknown ages distributed with warnings, thousands
  separators and `decimal = ","`, trimmed labels, warnings for levels that
  differ between tables, detection of total levels among the parts,
  `mt_relabel()`, `mt_regroup(age = "common")`; `mt_as_exposure()` needs an
  explicit `type`; `mt_read_genesis()` makes no analysis decisions.
* Errors name the missing periods, strata or cells and the likely cause;
  duplicates name unused columns; overlapping age groups and rates above 1
  are located.
* Robustness: Lee-Carter handles age groups without deaths and sparse
  counts (ridge fallback); failed strata give one warning; forecasts with
  rates above 1 give a warning; `rwd(observation_noise = TRUE)` pools the
  level variance over ages where deaths are few; smoothing leaves age 0
  out and respects group widths; models with one age group work.
* Backtests: `n_origins` origins with complete horizons, seeded PIT,
  estimated cells not scored by default (`score_estimated`).
* `step_window()` (fitting period), `as_demogdata()`, `as_stmomodata()`,
  plots with the original data (`data`), stratum selection and facets for
  life expectancy, share plots that take their defaults from the
  hierarchy, a warning when both `h` and `to` are given.
* Documentation: reference index by analysis step, a function map with
  workflow diagrams, the "Model reference" vignette, a model-choice table
  and an aggregation table in the FAQ, no references to 0.1 or internal
  files in help pages.

* Data layer: `mt_axis` (`single_ages()`, `age_groups()`, `single_years()`,
  `parse_age_labels()`), `mt_data()` with tables at their own resolution,
  `mt_grid()`/`mt_set_grid()`, `mt_array()`,
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
  `ets_series()`, and `ensemble()`.
* `apc_bayes()`: seamless wrapper of the bamp package (Bayesian
  age-period-cohort models) with convergence checks and draws on the grid.
* Plots (`mt_plot_*()`, ggplot2), HMD readers (`mt_from_hmd()`,
  `mt_read_hmd()`), and four vignettes using HMD data.
* Ports of the 0.1 models and methods: `hyndman_ullah()`, `stmomo()` (APC, CBD, RH,
  M6, LC), `li_lee()`, `product_ratio()`, `cause_shares()` (Tucker4-ILR
  multi-cause forecasts), `mt_reconcile()` (MinT, WLS, OLS),
  `mt_order_strata()`, stacking weights in `ensemble()`, Whittaker and GAM
  smoothing,
  Kannisto closure in `life_table()`, monotone splitting of deaths, and the
  simulated multi-cause dataset `mt_example_causes()`.
* Grouped hierarchies for coherent forecasts (for example sex x region x
  deprivation quintile): `mt_hierarchy()`, `mt_aggregate_strata()`,
  multi-variable hierarchy checks in `mt_validate()`, `mt_reconcile()` on
  the death scale with WLS, MinT, structural, OLS or bottom-up
  reconciliation, an empirical copula for separately fitted strata and
  non-negativity; `reconciled()` for reconciliation inside fits and
  backtests; `mt_aggregate_draws()`; `li_lee(common_by = "sex")`;
  simulated `mt_example_regions()` and a vignette.
* Rate-based spline ungrouping (`mt_ungroup(method = "hyman_kannisto")`,
  `"natural_gompertz"`, `"gam"`), `mt_disaggregate()` for forecast draws,
  `bootstrapped()` (parameter and preprocessing uncertainty by Poisson
  resampling and refitting) and information-criterion rank selection in
  `tucker3(ranks = "bic")`.
* "Frequently asked questions" vignette; vdiffr snapshot tests of all plots;
  Tucker factor signs are fixed so that fits are identical across LAPACK
  builds.
* Vignettes "Getting started" and "Causes of death by deprivation" use
  `ons_data`; the 0.1-to-0.2 migration vignette is removed.
* `li_lee()` and `product_ratio()` gain `jump_off = "observed"`; Lee-Carter
  labels name their variant (SVD, observed jump-off, bootstrap).
* `asr()` and `prob_dying()` return the values of every path with
  `paths = TRUE`; `mt_plot_shares()` gains `groups` and `causes` filters.
* General preprocessing replaces the dataset-specific readers
  (`mt_read_ons()` is removed): `mt_data()` reads text values with symbols
  for suppressed or missing cells (`missing`, `zero`), wide tables with one
  value column per level of a stratum variable (one column per sex or
  cause), and tables without some stratum columns, which are shared over
  their levels (populations shared by causes of death).
  `mt_group_levels()` adds an aggregation layer to any stratum variable
  (cause groups, deciles to quintiles, districts to regions) with summed or
  shared exposures, a residual group (`other`, or the `total` minus the
  groups) and a new name for the layer. `mt_subset()` restricts strata,
  ages and periods. `mt_regroup(age = "deaths")` aggregates to the age
  groups of a table.
* The default modelling grid is the coarsest grid on which every table is
  exact (the age groups of the data when all tables share them), so models
  fit grouped data directly; `mt_ungroup()` ungroups to single years by
  default (new argument `age`). The error for deaths coarser than the grid
  names both remedies.
* Fixed: reconciled predictive draws of rate-additive hierarchies (causes)
  returned aggregated deaths instead of rates; the European Standard
  Population weights of 95+ were dropped for an open 90+ grid group;
  smoothing of 5-year age groups failed because of singular starting
  values.
* `rwd(observation_noise = TRUE)`: random walk with drift plus noise (local
  level model); sharper than the plain random walk. Interval calibration of
  random walks, Lee-Carter and ensembles on simulated, HMD and ONS data is
  documented in `dev/calibration.md`; the defaults stay.
* `ons_data`: deaths by cause and populations by sex, age group and IMD
  decile, England 2001-2019 (ONS ad hoc release 3064, Open Government
  Licence), with the cause groups `ons_cause_groups`; built by
  `data-raw/ons_data.R`.
* German regional data: `de_deaths` (deaths by Land, sex, age group and
  ICD-10 chapter, 1998-2024, GBE-Bund) and `de_population` (31 December
  populations, GENESIS-Online), Datenlizenz Deutschland by-2-0;
  `mt_read_gbe()` parses GBE-Bund exports, and `mt_fill_suppressed()`
  estimates suppressed death counts by raking them to the published totals. `mt_as_exposure(type = "dec31")` for populations on
  31 December. Missing cells are reported by `mt_validate()`.
* Ensembles draw one component per path for all strata (before, components
  were picked separately per stratum, which broke the dependence between
  strata and the coherence of reconciled components); every model type is
  tested as a component.
* `mt_read_genesis()` reads deaths and 31 December populations in the flat
  CSV format of GENESIS-Online (Destatis); `parse_age_labels()` understands
  the German age labels.
* `apc_select()` (DIC-based choice of random-walk
  orders for `apc_bayes()`), and the plots `mt_plot_tucker()`,
  `mt_plot_shares()` and `mt_plot_scores()`.
* `mt_tune()` and `tuned()`; `prob_dying()`, `asr()` (ESP 2013),
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
