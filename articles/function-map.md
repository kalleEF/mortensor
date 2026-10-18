# Function map

This page maps the functions of mortensor onto the steps of an analysis
and shows how their results feed into each other. Each box names a
function and its main purpose. Paths through the maps are workflows; the
section “What is not covered” lists steps that need other tools.

## The main path

Every analysis passes the same objects from top to bottom. The boxes are
objects; the labels on the arrows are the functions that connect them.

![Main path: your tables become an mt_data object with mt_data() and are
prepared and checked (map 1); mt_fit() fits a model specification to the
data, and mt_backtest() scores specifications on past years;
mt_forecast() gives mt_draws, simulated paths of future rates; tables
and plots of life expectancy, lifespan variation, rates, deaths, gaps,
decompositions and excess mortality come from the draws, from the
backtest scores (mt_score()) or from the observed data (map
3).](function-map-main.png)

## Map 1: from tables to analysis-ready data

The steps in their usual order. Each is a decision about the data; skip
those a data set does not need. Every step returns an `mt_data` object.

![Map 1, from tables to analysis-ready data: read the tables (mt_data(),
mt_read_genesis(), mt_read_gbe(), mt_from_hmd(), mt_read_hmd(),
as_mt_data()), then mt_as_exposure(), mt_subset(), mt_regroup(),
mt_fill_suppressed(), mt_group_levels() and mt_relabel(), mt_ungroup()
and mt_validate(), giving model-ready mt_data.](function-map-data.png)

Smoothing is not a step of this map.
[`mt_smooth()`](../reference/mt_smooth.md) with
[`mt_smooth_check()`](../reference/mt_smooth_check.md) and
[`mt_plot_smooth()`](../reference/mt_smooth_check.md) shows what a
smoother would do to the data; for forecasts, smooth inside the fit with
[`step_smooth()`](../reference/mt_prep.md) (map 2), so that observation
noise and backtests refer to the raw deaths.
[`mt_fit()`](../reference/mt_fit.md) warns about data smoothed
beforehand and [`mt_backtest()`](../reference/mt_backtest.md) refuses
them.

## Map 2: models, coherence and evaluation

A model specification is one model, optionally inside wrappers.
Preprocessing steps given as `prep` are repeated in every fit, so that
backtests do not see the test years.

![Map 2, models: models for one population per stratum (including
coda()) or for several strata jointly (including coda(coherent = TRUE)
and coda_causes()) form a model specification, optionally inside
wrappers (ensemble(), stratified(), reconciled(), recalibrated(),
anchored(), bootstrapped(), tuned()) that may use a hierarchy from
mt_hierarchy(); mt_fit() takes the specification, the model-ready data
and optional preprocessing from mt_prep(); the fit is inspected
(summary(), mt_parameters(), residuals()) or forecast with
mt_forecast().](function-map-models.png)

Choosing between specifications:

![Choosing between specifications: the data and the candidate
specifications, each with its preprocessing, go into mt_backtest();
mt_score(), the score plots (mt_plot_scores(), mt_plot_coverage()) and
mt_dm_test() (is the difference between two models larger than chance?)
lead to the chosen specification.](function-map-choice.png)

[`tuned()`](../reference/mt_tune.md) and
[`mt_tune()`](../reference/mt_tune.md) run the same comparison inside a
fit, over a grid of model and preprocessing settings (for example the
number of components and the smoothing of
[`hyndman_ullah()`](../reference/hyndman_ullah.md), or the fitting
window), with rolling origins or a single hold-out of the last years.
[`recalibrated()`](../reference/recalibrated.md) keeps the chosen model
and rescales the width of its intervals.

## Map 3: from draws to results

Indicators and comparisons are computed path by path, so their intervals
include the joint uncertainty of all strata. The same functions describe
observed data.

![Map 3, from draws to results: forecast draws (mt_draws) and observed
data (mt_data) give indicators (life_expectancy(), lifespan_variation(),
asr(), prob_dying(), yll(), mt_deaths()), comparisons (mt_contrast(),
mt_inequality(), mt_decompose(), and mt_excess() for observed data
against a baseline fit) and output (mt_summarise(), plots including
mt_plot_excess(), mt_export()); draws can be aggregated or refined
(mt_aggregate_draws(), mt_disaggregate(),
mt_coherence()).](function-map-results.png)

## Functions by purpose

| Step | Function | Purpose | Input | Output |
|----|----|----|----|----|
| Read | [`mt_data()`](../reference/mt_data.md) | read deaths, exposures, populations (long or wide, suppression symbols, shared populations) | data frames | `mt_data` |
| Read | [`mt_add()`](../reference/mt_add.md), [`mt_add_deaths()`](../reference/mt_add.md), [`mt_add_exposure()`](../reference/mt_add.md), [`mt_add_population()`](../reference/mt_add.md) | add or replace one table | `mt_data`, data frame | `mt_data` |
| Read | [`mt_read_genesis()`](../reference/mt_read_genesis.md), [`mt_read_gbe()`](../reference/mt_read_gbe.md) | read German official exports | files | `mt_data`, table |
| Read | [`mt_from_hmd()`](../reference/mt_from_hmd.md), [`mt_read_hmd()`](../reference/mt_from_hmd.md), [`as_mt_data()`](../reference/as_mt_data.md) | read or convert HMD, StMoMo, demography data | objects, files | `mt_data` |
| Read | [`parse_age_labels()`](../reference/parse_age_labels.md), [`age_groups()`](../reference/mt_axis.md), [`single_ages()`](../reference/mt_axis.md), [`single_years()`](../reference/mt_axis.md), [`period_groups()`](../reference/mt_axis.md) | age and period axes from labels or bounds | labels | `mt_axis` |
| Read | [`axis_widths()`](../reference/axis_widths.md), [`axis_midpoints()`](../reference/axis_widths.md), [`is_open()`](../reference/axis_widths.md) | interval widths, midpoints and the open interval of an axis | `mt_axis` | vector |
| Prepare | [`mt_as_exposure()`](../reference/mt_as_exposure.md) | exposures from populations (`type` required) | `mt_data` | `mt_data` |
| Prepare | [`mt_subset()`](../reference/mt_subset.md) | restrict ages, years, strata | `mt_data` | `mt_data` |
| Prepare | [`mt_strata()`](../reference/mt_strata.md), [`mt_obs()`](../reference/mt_strata.md), [`mt_grid()`](../reference/mt_strata.md), [`mt_set_grid()`](../reference/mt_set_grid.md), [`mt_granularity()`](../reference/mt_granularity.md) | look at strata (with their level in a hierarchy), tables, their resolution and the modelling grid | `mt_data` | table, grid |
| Prepare | [`mt_array()`](../reference/mt_array.md), [`mt_operator()`](../reference/mt_operator.md) | data as arrays on the grid; the matrix that sums grid cells to observed cells (for your own models) | `mt_data` | array, matrix |
| Prepare | [`mt_group_levels()`](../reference/mt_group_levels.md), [`mt_relabel()`](../reference/mt_relabel.md) | aggregation layers, recode labels | `mt_data` | `mt_data` |
| Prepare | [`mt_fill_suppressed()`](../reference/mt_fill_suppressed.md) | estimate suppressed counts | `mt_data`, hierarchy or bounds | `mt_data` |
| Check | [`mt_validate()`](../reference/mt_validate.md) | plausibility and hierarchy checks | `mt_data` | findings |
| Harmonise | [`mt_regroup()`](../reference/mt_regroup.md), [`mt_ungroup()`](../reference/mt_ungroup.md), [`mt_smooth()`](../reference/mt_smooth.md) | common resolution, single years (PCLM or monotone interpolation; the spline methods `"hyman_kannisto"`, `"natural_gompertz"` and `"gam"` are deprecated), smoothing over age, time or both | `mt_data` | `mt_data` |
| Harmonise | [`mt_plot_ungroup()`](../reference/ungroup-diagnostics.md), [`mt_ungroup_check()`](../reference/ungroup-diagnostics.md) | ungrouped against grouped rates, stability over time, errors against known single years; several methods side by side | `mt_data` and ungrouped data (or a list) | ggplot, table |
| Harmonise | [`mt_smooth_check()`](../reference/mt_smooth_check.md), [`mt_plot_smooth()`](../reference/mt_smooth_check.md) | how much smoothing removed (dispersion, trend ratio); several smoothers in one figure | `mt_data` pair or list | table, ggplot |
| Harmonise | [`mt_prep()`](../reference/mt_prep.md), [`mt_prep_apply()`](../reference/mt_prep.md) with [`step_ungroup()`](../reference/mt_prep.md), [`step_smooth()`](../reference/mt_prep.md), [`step_regroup()`](../reference/mt_prep.md), [`step_grid()`](../reference/mt_prep.md), [`step_shocks()`](../reference/mt_prep.md), [`step_window()`](../reference/mt_prep.md) | preprocessing inside fits and backtests | steps | pipeline |
| Model | [`lee_carter()`](../reference/lee_carter.md) (optionally rotated), [`rwd()`](../reference/rwd.md), [`li_lee()`](../reference/li_lee.md), [`hyndman_ullah()`](../reference/hyndman_ullah.md), [`coherent_fdm()`](../reference/coherent_fdm.md), [`coda()`](../reference/coda.md) (also coherent), [`stmomo()`](../reference/stmomo.md), [`stanmomo()`](../reference/stanmomo.md), [`arima_series()`](../reference/arima_series.md), [`ets_series()`](../reference/arima_series.md), [`apc_bayes()`](../reference/apc_bayes.md), [`tucker3()`](../reference/tucker3.md), [`coda_causes()`](../reference/coda_causes.md); deprecated: [`product_ratio()`](../reference/product_ratio.md) (use [`coherent_fdm()`](../reference/coherent_fdm.md)), [`cause_shares()`](../reference/cause_shares.md) (use [`reconciled()`](../reference/reconciled.md) or [`coda_causes()`](../reference/coda_causes.md)), `rwd(observation_noise = TRUE)` (use `arima_series(order = c(0, 1, 1), drift = TRUE)`), `tucker3(transform = "coda")` (use [`coda()`](../reference/coda.md)) | model specifications | settings | `mt_spec` |
| Model | [`apc_select()`](../reference/apc_select.md), [`mcmc_control()`](../reference/mcmc_control.md), [`mt_tune()`](../reference/mt_tune.md) | choose random-walk orders, MCMC settings; model and preprocessing settings by forecast accuracy on the last years (hold-out or rolling origins; rates, life expectancy or ASR) | data, grid | table, `mt_spec` |
| Combine | [`ensemble()`](../reference/ensemble.md) (optionally `top_k`), [`stratified()`](../reference/stratified.md), [`reconciled()`](../reference/reconciled.md) (also score-optimal, `method = "energy_score"`), [`bootstrapped()`](../reference/bootstrapped.md), [`tuned()`](../reference/mt_tune.md) | mixtures, different models for groups of strata, coherence, parameter uncertainty (and of smoothing, `type = "data"`), tuning | `mt_spec` | `mt_spec` |
| Combine | [`mt_order_strata()`](../reference/mt_order_strata.md) | keep an ordering between strata (for example women below men) path by path | `mt_draws` | `mt_draws` |
| Calibrate | [`recalibrated()`](../reference/recalibrated.md) | rescale the spread (optionally the drift) of any model’s draws to its forecast errors in the last years, overall or by age band, for rates, life expectancy or ASR | `mt_spec` | `mt_spec` |
| Jump-off | [`anchored()`](../reference/anchored.md), [`mt_jump_off()`](../reference/mt_jump_off.md) | start forecasts at the observed rates of the last year (optionally of several years, shrunk and fading out); report the first-year step | `mt_spec`, `mt_draws` | `mt_spec`, table |
| Coherence | [`mt_hierarchy()`](../reference/mt_hierarchy.md), [`mt_aggregate_strata()`](../reference/mt_aggregate_strata.md), [`mt_reconcile()`](../reference/mt_reconcile.md) (with `fixed` nodes; `"proportional"` scales the parts of fixed nodes), [`mt_aggregate_draws()`](../reference/mt_aggregate_draws.md) | define and apply totals of strata | hierarchy | `mt_data`, `mt_draws` |
| Coherence | [`mt_coherence()`](../reference/mt_coherence.md), [`mt_plot_coherence()`](../reference/mt_coherence.md) | sum of parts over aggregate, per year | `mt_draws`, hierarchy | table, ggplot |
| Fit | [`mt_fit()`](../reference/mt_fit.md) | estimate | `mt_spec`, `mt_data` | `mt_fit` |
| Inspect | [`summary()`](https://rdrr.io/r/base/summary.html), [`mt_parameters()`](../reference/mt_parameters.md), [`fitted()`](https://rdrr.io/r/stats/fitted.values.html), [`residuals()`](https://rdrr.io/r/stats/residuals.html), [`mt_plot_residuals()`](../reference/mt_plot_residuals.md) | model, assumptions, estimates, fit | `mt_fit` | text, table, plot |
| Forecast | [`mt_forecast()`](../reference/mt_forecast.md), [`mt_summarise()`](../reference/mt_draws.md) | simulated paths and their quantiles | `mt_fit` | `mt_draws`, table |
| Evaluate | [`mt_backtest()`](../reference/mt_backtest.md), [`mt_score()`](../reference/mt_score.md) | rolling-origin scores of rates, life expectancy or ASR, and of gaps between strata (`contrast`) | `mt_spec`, `mt_data` | table |
| Evaluate | [`mt_compare()`](../reference/mt_compare.md), [`mt_score_table()`](../reference/mt_compare.md) | backtests of several models on the same origins; scores of rates, life expectancy and gaps side by side, by level of a hierarchy | named list of `mt_spec`, `mt_data` | table |
| Evaluate | [`mt_dm_test()`](../reference/mt_dm_test.md) | Diebold-Mariano tests: is one model more accurate than another, per horizon, beyond chance? | backtest from [`mt_compare()`](../reference/mt_compare.md) | table |
| Evaluate | [`mt_plot_scores()`](../reference/mt_plot_more.md), [`mt_plot_coverage()`](../reference/mt_plot.md), [`mt_plot_pit()`](../reference/mt_plot_design.md), [`mt_plot_backtest_design()`](../reference/mt_plot_design.md), [`mt_plot_tuning()`](../reference/mt_plot_design.md) | scores by horizon or model, interval coverage, PIT histograms, the design of a backtest, tuning scores | backtest, tuning | ggplot |
| Describe | [`mt_observed()`](../reference/mt_observed.md) | observed rates as draws | `mt_data` | `mt_draws` |
| Report | [`life_expectancy()`](../reference/life_expectancy.md), [`asr()`](../reference/derived.md), [`prob_dying()`](../reference/derived.md) | indicators path by path, observed or forecast | `mt_draws` or `mt_data` | table |
| Report | [`lifespan_variation()`](../reference/lifespan_variation.md) | e-dagger, entropy, Gini, standard deviation, coefficient of variation and interquartile range of the age at death | `mt_draws`, `mt_data` or rates | table, array |
| Report | [`mt_excess()`](../reference/mt_excess.md), [`mt_plot_excess()`](../reference/mt_excess.md) | excess deaths, P-scores and loss of life expectancy against the expected values of a baseline model, by age, stratum and year | `mt_fit` or `mt_spec`, `mt_data` | table, ggplot |
| Report | [`yll()`](../reference/derived.md), [`mt_deaths()`](../reference/mt_deaths.md) | years of life lost, deaths with projected populations | `mt_draws` | table, array |
| Report | [`life_table()`](../reference/life_table.md), [`standard_population()`](../reference/standard_population.md) | life tables from rates; standard populations esp2013, esp1976, who2000 | rates, name | table |
| Report | [`mt_contrast()`](../reference/mt_contrast.md), [`mt_inequality()`](../reference/mt_inequality.md), [`mt_plot_contrast()`](../reference/mt_plot_contrast.md) | gaps, ratios, SII, RII; gaps over time | `mt_draws` or `mt_data` | table, ggplot |
| Report | [`mt_compare_forecasts()`](../reference/mt_compare.md), [`mt_plot_compare()`](../reference/mt_compare.md) | one quantity of several forecasts (models, jump-off choices, reconciliation methods) side by side | named list of `mt_draws` or `mt_data` | table, ggplot |
| Report | [`mt_decompose()`](../reference/mt_decompose.md), [`mt_plot_decomposition()`](../reference/mt_plot_decomposition.md) | differences in life expectancy (or other summaries) by age, cause and subgroup composition: Arriaga, Horiuchi, stepwise replacement, Kitagawa, Arriaga meets Kitagawa; jump-off share and changes of gaps | `mt_draws` or `mt_data` (or both; also census years with gaps) | `mt_decomposition`, table, ggplot |
| Report | [`mt_plot_forecast()`](../reference/mt_plot.md), [`mt_plot_age()`](../reference/mt_plot.md), [`mt_plot_life_expectancy()`](../reference/mt_plot.md), [`mt_plot_shares()`](../reference/mt_plot_more.md), [`mt_plot_rates()`](../reference/mt_plot.md), [`mt_plot_tucker()`](../reference/mt_plot_more.md), [`mt_export()`](../reference/mt_export.md) | figures and CSV files | draws, data, tables | ggplot, file |
| Data | `ons_data`, `ons_cause_groups`, `de_deaths`, `de_population`, `de_cause_groups`, [`mt_example()`](../reference/mt_example.md), [`mt_example_causes()`](../reference/mt_example_causes.md), [`mt_example_regions()`](../reference/mt_example_regions.md) | example data (published and simulated) |  | tables, `mt_data` |
| Report | [`theme_mortensor()`](../reference/mt_plot_style.md), [`scale_colour_mortensor()`](../reference/mt_plot_style.md), [`scale_fill_mortensor()`](../reference/mt_plot_style.md), [`mt_palette()`](../reference/mt_plot_style.md) | the package’s plot style for your own figures |  | theme, scales |
| Disaggregate | [`mt_disaggregate()`](../reference/mt_disaggregate.md), [`mt_shape_rates()`](../reference/mt_shape_rates.md), [`mt_disaggregate_check()`](../reference/mt_disaggregate_check.md), [`mt_plot_disaggregation()`](../reference/mt_disaggregate_check.md) | forecasts on finer ages with a fixed age shape within each group (the spline methods and `mt_shape_rates(method = "gam")` are deprecated) | `mt_draws` | `mt_draws`, table, ggplot |
| Exchange | [`as_demogdata()`](../reference/as_demogdata.md), [`as_stmomodata()`](../reference/as_demogdata.md) | cross-checks with demography and StMoMo | `mt_data` | `demogdata`, `StMoMoData` |
| Extend | [`new_mt_spec()`](../reference/new_mt_spec.md) | constructor for developers adding a model (methods for the internal generics `fit_model()` and `simulate_model()`) | model name, parameters | `mt_spec` |

## What is not covered

These steps of a mortality analysis need other tools or are open
development points:

- **Population projections.** mortensor forecasts death rates; future
  exposures (for expected deaths) come from a population projection, for
  example of the statistical office, passed as `exposure`.
- **Spatial models for small areas.** There is no model with spatial
  smoothing between neighbouring regions (such as BYM models), and no
  relational model for small populations (such as TOPALS); pooling works
  through common trends ([`li_lee()`](../reference/li_lee.md)); Bayesian
  APC models smooth over age, period and cohort within each stratum.
- **Models of life expectancy.** Models that forecast life expectancy
  itself (such as the Bayesian models of the UN or the double-gap model)
  are not implemented; mortensor forecasts rates and derives life
  expectancy from them.
- **Multiple imputation of suppressed cells.**
  [`mt_fill_suppressed()`](../reference/mt_fill_suppressed.md) returns
  one estimate; the uncertainty of suppressed counts is not carried into
  forecasts.
- **Causes of death with changing coding.** ICD revisions and coding
  changes (for example of dementia) are not bridged; cause groups that
  span the affected codes reduce the problem.
