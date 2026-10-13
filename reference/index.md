# Package index

## 

Read data

Deaths, exposures and populations in long or wide tables, at any age and
period resolution.

- [`mt_data()`](mt_data.md) : Mortality data observed at any granularity

- [`mt_add()`](mt_add.md) [`mt_add_deaths()`](mt_add.md)
  [`mt_add_exposure()`](mt_add.md) [`mt_add_population()`](mt_add.md) :

  Add observations to an `mt_data` object

- [`parse_age_labels()`](parse_age_labels.md) : Parse age-group labels

- [`as_mt_data()`](as_mt_data.md) : Convert other mortality data formats

- [`mt_from_hmd()`](mt_from_hmd.md) [`mt_read_hmd()`](mt_from_hmd.md) :
  Human Mortality Database tables

- [`mt_read_genesis()`](mt_read_genesis.md) : Read deaths and
  populations from GENESIS (Destatis) flat files

- [`mt_read_gbe()`](mt_read_gbe.md) : Read tables exported from GBE-Bund

## 

Prepare and check

Analysis decisions on exposures, ages, years, strata, aggregation layers
and suppressed cells.

- [`mt_as_exposure()`](mt_as_exposure.md) : Exposures from population
  counts

- [`mt_subset()`](mt_subset.md) : Restrict data to ages, periods and
  strata

- [`mt_group_levels()`](mt_group_levels.md) : Combine levels of a
  stratum variable

- [`mt_relabel()`](mt_relabel.md) : Rename levels of a stratum variable

- [`mt_fill_suppressed()`](mt_fill_suppressed.md) : Estimate suppressed
  death counts

- [`mt_validate()`](mt_validate.md) : Check mortality data for common
  problems

- [`mt_strata()`](mt_strata.md) [`mt_obs()`](mt_strata.md)
  [`mt_grid()`](mt_strata.md)
  [`as.data.table(`*`<mt_data>`*`)`](mt_strata.md) :

  Strata, grid and tables of an `mt_data` object

## 

Harmonise age groups and periods

Common resolution of all tables; ungrouping and smoothing; preprocessing
pipelines refitted in backtests.

- [`mt_regroup()`](mt_regroup.md) : Regroup data onto coarser age or
  period intervals
- [`mt_ungroup()`](mt_ungroup.md) : Estimate grid-resolution data from
  grouped observations
- [`mt_smooth()`](mt_smooth.md) : Smooth death rates on the grid
- [`mt_smooth_check()`](mt_smooth_check.md)
  [`mt_plot_smooth()`](mt_smooth_check.md) : Check and plot the result
  of smoothing
- [`mt_prep()`](mt_prep.md) [`step_ungroup()`](mt_prep.md)
  [`step_smooth()`](mt_prep.md) [`step_regroup()`](mt_prep.md)
  [`step_grid()`](mt_prep.md) [`step_window()`](mt_prep.md)
  [`step_shocks()`](mt_prep.md) [`mt_prep_apply()`](mt_prep.md) :
  Preprocessing pipelines
- [`mt_set_grid()`](mt_set_grid.md) : Set the modelling grid

## 

Describe the past

Observed indicators with Poisson intervals, gaps between strata,
inequality indices and decompositions (which also work on forecasts).

- [`mt_observed()`](mt_observed.md) : Observed rates as draws
- [`life_expectancy()`](life_expectancy.md) : Life expectancy from rates
  or forecast draws
- [`prob_dying()`](derived.md) [`asr()`](derived.md)
  [`yll()`](derived.md) : Derived quantities from forecast draws or
  observed data
- [`mt_contrast()`](mt_contrast.md) : Differences and ratios between
  strata, path by path
- [`mt_inequality()`](mt_inequality.md) : Slope and relative index of
  inequality, path by path
- [`mt_decompose()`](mt_decompose.md)
  [`summary(`*`<mt_decomposition>`*`)`](mt_decompose.md) : Decompose
  differences in mortality by age, cause and subgroup
- [`mt_plot_decomposition()`](mt_plot_decomposition.md) : Plot a
  decomposition
- [`standard_population()`](standard_population.md) : Standard
  populations
- [`life_table()`](life_table.md) : Life tables from death rates on any
  age axis

## 

Models for one population

- [`lee_carter()`](lee_carter.md) : Lee-Carter model
- [`rwd()`](rwd.md) : Random walk with drift on log death rates
- [`hyndman_ullah()`](hyndman_ullah.md) : Functional demographic model
  of Hyndman and Ullah (demography)
- [`stmomo()`](stmomo.md) : Generalised age-period-cohort models
  (StMoMo)
- [`arima_series()`](arima_series.md) [`ets_series()`](arima_series.md)
  : Per-age ARIMA and ETS models
- [`apc_bayes()`](apc_bayes.md) : Bayesian age-period-cohort model
  (bamp)
- [`apc_select()`](apc_select.md) : Choose the random-walk orders of a
  Bayesian APC model by DIC
- [`mcmc_control()`](mcmc_control.md) : MCMC settings for apc_bayes()

## 

Models for several populations and causes

- [`li_lee()`](li_lee.md) [`product_ratio()`](li_lee.md) : Coherent
  multi-population models
- [`tucker3()`](tucker3.md) : Tucker3 tensor model across strata
- [`cause_shares()`](cause_shares.md) : Coherent multi-cause forecasts
  with cause shares

## 

Combine models and make forecasts coherent

Model wrappers (ensembles, reconciliation, calibration, jump-off,
bootstrap, tuning) and hierarchies of strata.

- [`ensemble()`](ensemble.md) : Ensembles of models
- [`reconciled()`](reconciled.md) : Models with reconciled forecasts
- [`calibrated()`](calibrated.md) : Calibrate forecast distributions on
  recent years
- [`anchored()`](anchored.md) : Control the jump-off of any model
- [`bootstrapped()`](bootstrapped.md) : Bootstrap uncertainty of
  preprocessing and parameters
- [`mt_tune()`](mt_tune.md) [`tuned()`](mt_tune.md) : Tune model
  parameters on training data
- [`mt_hierarchy()`](mt_hierarchy.md) : Hierarchies of strata
- [`mt_aggregate_strata()`](mt_aggregate_strata.md) : Add aggregate
  strata to data
- [`mt_reconcile()`](mt_reconcile.md) : Reconcile forecast draws over a
  hierarchy
- [`mt_coherence()`](mt_coherence.md)
  [`mt_plot_coherence()`](mt_coherence.md) : Check the coherence of
  forecast draws
- [`mt_aggregate_draws()`](mt_aggregate_draws.md) : Aggregate forecast
  draws over strata or ages
- [`mt_order_strata()`](mt_order_strata.md) : Enforce an ordering
  between strata in forecast draws

## 

Fit, inspect and forecast

- [`mt_fit()`](mt_fit.md) : Fit a model to mortality data
- [`mt_parameters()`](mt_parameters.md) : Parameters of a fitted model
- [`fitted(`*`<mt_fit>`*`)`](fit-inspection.md)
  [`residuals(`*`<mt_fit>`*`)`](fit-inspection.md)
  [`summary(`*`<mt_fit>`*`)`](fit-inspection.md) : Fitted rates,
  residuals and summaries of fitted models
- [`mt_plot_residuals()`](mt_plot_residuals.md) : Residual heat map of a
  fitted model
- [`mt_forecast()`](mt_forecast.md) : Forecast draws from a fitted model
- [`mt_jump_off()`](mt_jump_off.md) : Step between the last observed
  year and the first forecast year
- [`mt_summarise()`](mt_draws.md)
  [`as.data.table(`*`<mt_draws>`*`)`](mt_draws.md) : Forecast draws
- [`mt_deaths()`](mt_deaths.md) : Expected deaths from draws and future
  exposures
- [`mt_disaggregate()`](mt_disaggregate.md) : Disaggregate forecast
  draws to a finer age grid
- [`mt_shape_rates()`](mt_shape_rates.md) : Age shapes within groups for
  disaggregating forecasts
- [`mt_disaggregate_check()`](mt_disaggregate_check.md)
  [`mt_plot_disaggregation()`](mt_disaggregate_check.md) : Check and
  plot disaggregated forecasts

## 

Evaluate

- [`mt_backtest()`](mt_backtest.md) : Rolling-origin backtest
- [`mt_score()`](mt_score.md) : Summarise backtest scores

## 

Report

- [`mt_plot_forecast()`](mt_plot.md) [`mt_plot_age()`](mt_plot.md)
  [`mt_plot_life_expectancy()`](mt_plot.md)
  [`mt_plot_rates()`](mt_plot.md) [`mt_plot_ungroup()`](mt_plot.md)
  [`mt_plot_coverage()`](mt_plot.md) : Plots
- [`mt_plot_tucker()`](mt_plot_more.md)
  [`mt_plot_shares()`](mt_plot_more.md)
  [`mt_plot_scores()`](mt_plot_more.md) : Further plots: Tucker factors,
  cause shares and backtest scores
- [`theme_mortensor()`](mt_plot_style.md)
  [`mt_palette()`](mt_plot_style.md)
  [`scale_colour_mortensor()`](mt_plot_style.md)
  [`scale_fill_mortensor()`](mt_plot_style.md) : Plot style of mortensor
- [`mt_export()`](mt_export.md) : Export forecasts and result tables

## Example data

- [`ons_data`](ons_data.md) [`ons_cause_groups`](ons_data.md) : Deaths
  and populations by deprivation decile, England 2001-2019
- [`de_deaths`](de_deaths.md) [`de_population`](de_deaths.md)
  [`de_cause_groups`](de_deaths.md) : Deaths and populations by Land,
  Germany
- [`mt_example()`](mt_example.md) : Simulated mortality data with known
  truth
- [`mt_example_causes()`](mt_example_causes.md) : Simulated deaths by
  cause
- [`mt_example_regions()`](mt_example_regions.md) : Simulated deaths by
  sex, region and deprivation quintile

## Exchange with other packages

- [`as_demogdata()`](as_demogdata.md)
  [`as_stmomodata()`](as_demogdata.md) : Convert data to demography and
  StMoMo objects

## Low level and extending

Building blocks for advanced use and for new models.

- [`mt_axis()`](mt_axis.md) [`single_ages()`](mt_axis.md)
  [`age_groups()`](mt_axis.md) [`single_years()`](mt_axis.md)
  [`period_groups()`](mt_axis.md) : Age and period axes
- [`axis_widths()`](axis_widths.md) [`axis_midpoints()`](axis_widths.md)
  [`is_open()`](axis_widths.md) : Interval widths and midpoints of an
  axis
- [`mt_array()`](mt_array.md) : Arrays on the modelling grid
- [`mt_granularity()`](mt_granularity.md) : Describe the granularity of
  each table
- [`mt_operator()`](mt_operator.md) : Aggregation operator from grid
  cells to observed cells
- [`new_mt_spec()`](new_mt_spec.md) : Model specifications
