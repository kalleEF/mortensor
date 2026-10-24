# mortensor

Flexible probabilistic mortality forecasting from data observed at
different granularities: age groups and open age intervals, single or
grouped periods, population strata and causes.

mortensor has no preferred model. It implements common mortality models
and wraps others from demography, StMoMo, StanMoMo, bamp and forecast
behind one interface: every model is fitted with
[`mt_fit()`](reference/mt_fit.md), forecast with
[`mt_forecast()`](reference/mt_forecast.md) and returns simulated paths
of the same form. Models can therefore be compared in the same backtest
and combined into ensembles of any composition, and the same steps
(smoothing, jump-off, recalibration, coherent totals) apply to all of
them.

> **Caution: experimental, vibe-coded software. Use at your own risk.**
> mortensor is written largely with AI coding assistants (“vibe-coded”).
> It has unit tests that check its methods against published results and
> reference implementations, but it has not been peer reviewed, its
> interfaces still change, and its results can be wrong in ways the
> tests do not catch. Check results independently before you use them in
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
StMoMo, StanMoMo, demography, ggplot2) use them when installed.

## What it does

| Layer | Functions |
|----|----|
| Data at any granularity | [`mt_data()`](reference/mt_data.md) (long or wide tables, symbols for suppressed cells, shared populations), [`mt_add()`](reference/mt_add.md), [`as_mt_data()`](reference/as_mt_data.md) (StMoMo, demography), [`mt_from_hmd()`](reference/mt_from_hmd.md), [`mt_read_hmd()`](reference/mt_from_hmd.md), [`mt_read_genesis()`](reference/mt_read_genesis.md), [`mt_read_gbe()`](reference/mt_read_gbe.md), [`mt_validate()`](reference/mt_validate.md), `ons_data`, `de_deaths` |
| Preparation | [`mt_as_exposure()`](reference/mt_as_exposure.md), [`mt_subset()`](reference/mt_subset.md), [`mt_group_levels()`](reference/mt_group_levels.md) (aggregation layers: cause groups, quintiles, regions), [`mt_relabel()`](reference/mt_relabel.md), [`mt_fill_suppressed()`](reference/mt_fill_suppressed.md) |
| Harmonisation | [`mt_regroup()`](reference/mt_regroup.md), [`mt_ungroup()`](reference/mt_ungroup.md) (penalized composite link model with exposure offset, 1D/2D; monotone interpolation; deprecated spline and GAM methods), [`mt_plot_ungroup()`](reference/ungroup-diagnostics.md) and [`mt_ungroup_check()`](reference/ungroup-diagnostics.md) (diagnostics), [`mt_smooth()`](reference/mt_smooth.md) (over age, over time or both; P-splines, Whittaker, GAM), [`mt_smooth_check()`](reference/mt_smooth_check.md), [`mt_plot_smooth()`](reference/mt_smooth_check.md) |
| Preprocessing pipelines | [`mt_prep()`](reference/mt_prep.md) with [`step_ungroup()`](reference/mt_prep.md), [`step_smooth()`](reference/mt_prep.md), [`step_regroup()`](reference/mt_prep.md), [`step_grid()`](reference/mt_prep.md), [`step_shocks()`](reference/mt_prep.md), [`step_window()`](reference/mt_prep.md); models fitted to smoothed rates take their observation noise from the raw rates |
| Describing the past | [`mt_observed()`](reference/mt_observed.md), [`life_expectancy()`](reference/life_expectancy.md), [`asr()`](reference/derived.md), [`prob_dying()`](reference/derived.md) on data, [`mt_contrast()`](reference/mt_contrast.md) (gaps, ratios), [`mt_inequality()`](reference/mt_inequality.md) (SII, RII), [`lifespan_variation()`](reference/lifespan_variation.md) (e-dagger, entropy, Gini, standard deviation, IQR), [`standard_population()`](reference/standard_population.md) |
| Models | [`rwd()`](reference/rwd.md), [`lee_carter()`](reference/lee_carter.md) (Poisson or SVD with the adjustments and period choice of demography; optional rotation of the age pattern, Li, Lee and Gerland 2013; negative binomial deaths, Delwarde, Denuit and Partrat 2007, with a dispersion by age, Li, Hardy and Tan 2009; averaged, filtered and constrained jump-off, Stoeldraijer et al. 2018, De Jong and Tickle 2006), [`coda()`](reference/coda.md) (compositional model of life-table deaths and its coherent version, Oeppen 2008; Bergeron-Boucher et al. 2017), [`p_splines()`](reference/p_splines.md) (two-dimensional P-splines, Currie, Durban and Eilers 2004; smooth dispersion by age, Djeundje and Currie 2011), [`cp_splines()`](reference/cp_splines.md) (constrained P-splines, Camarda 2019), [`bohk_rau()`](reference/bohk_rau.md) (age-specific rates of mortality improvement, Bohk-Ewald and Rau 2017), [`arima_series()`](reference/arima_series.md), [`ets_series()`](reference/arima_series.md), [`hyndman_ullah()`](reference/hyndman_ullah.md) (demography), [`stmomo()`](reference/stmomo.md) (APC, CBD, RH, M6, M7, M8, Plat, LC), [`apci()`](reference/apci.md) (age-period-cohort-improvement model, Richards et al. 2019), [`stanmomo()`](reference/stanmomo.md) (LC, RH, APC, CBD and M6, Bayesian, in Stan), [`li_lee()`](reference/li_lee.md), [`cae()`](reference/cae.md) (common age effect model, Kleinow 2015; non-diverging forecasts, Enchev, Kleinow and Cairns 2017), [`relative_model()`](reference/relative_model.md) (subpopulations relative to a reference population, Villegas and Haberman 2014), [`saint()`](reference/saint.md) (trend of a reference population and stationary spreads, Jarner and Kryger 2011; Jarner and Jallbjørn 2022), [`coherent_fdm()`](reference/coherent_fdm.md) (demography), [`apc_bayes()`](reference/apc_bayes.md) and [`apc_select()`](reference/apc_select.md) (bamp), [`nordpred_apc()`](reference/nordpred_apc.md) (Nordpred age-period-cohort model, Møller et al. 2003); experimental: [`tucker3()`](reference/tucker3.md) (three-way Lee-Carter), [`coda_multiway()`](reference/coda_multiway.md) (multiway compositional model, Bergeron-Boucher et al. 2018), [`coda_causes()`](reference/coda_causes.md) (also 2S-CoDA, Kjærgaard et al. 2019); deprecated: [`product_ratio()`](reference/product_ratio.md), [`cause_shares()`](reference/cause_shares.md) |
| Wrappers for any model | [`ensemble()`](reference/ensemble.md) (weights from inner backtests: stacking of predictive distributions with the CRPS, Yao et al. 2018; life expectancy bias, Kontis et al. 2017; inverse CRPS; `top_k`), [`stratified()`](reference/stratified.md) (different models for groups of strata), [`reconciled()`](reference/reconciled.md) (optionally with the Schaake shuffle of Clark et al. 2004, of every cell or of whole stratum paths, and MinT weights from past errors), [`coupled()`](reference/coupled.md) and [`mt_couple()`](reference/coupled.md) (separate forecasts given the dependence of a joint model by ensemble copula coupling, Schefzik et al. 2013), [`recalibrated()`](reference/recalibrated.md) (spread fitted to the errors of the last 10 years, optionally by age band or for life expectancy), [`anchored()`](reference/anchored.md) (jump-off at the observed level), [`bootstrapped()`](reference/bootstrapped.md), [`tuned()`](reference/mt_tune.md) and [`mt_tune()`](reference/mt_tune.md) (model and preprocessing settings chosen by forecast accuracy on the last years: hold-out or rolling origins, full or two-stage grids) |
| Coherence | [`mt_hierarchy()`](reference/mt_hierarchy.md), [`mt_aggregate_strata()`](reference/mt_aggregate_strata.md), [`mt_reconcile()`](reference/mt_reconcile.md) and [`reconciled()`](reference/reconciled.md) (bottom-up, WLS, MinT, structural, OLS, proportional; draw by draw; non-negative; `fixed` nodes such as a trusted all-cause total; score-optimal reconciliation learnt from past forecasts, Panagiotelis et al. 2023), [`mt_coherence()`](reference/mt_coherence.md), [`mt_plot_coherence()`](reference/mt_coherence.md), [`mt_aggregate_draws()`](reference/mt_aggregate_draws.md), [`mt_order_strata()`](reference/mt_order_strata.md) |
| Forecasts | [`mt_fit()`](reference/mt_fit.md), [`summary()`](https://rdrr.io/r/base/summary.html), [`mt_parameters()`](reference/mt_parameters.md), [`residuals()`](https://rdrr.io/r/stats/residuals.html), [`mt_forecast()`](reference/mt_forecast.md) (simulated paths), [`mt_jump_off()`](reference/mt_jump_off.md), [`mt_summarise()`](reference/mt_draws.md), [`mt_deaths()`](reference/mt_deaths.md), [`mt_disaggregate()`](reference/mt_disaggregate.md) with [`mt_shape_rates()`](reference/mt_shape_rates.md) and [`mt_disaggregate_check()`](reference/mt_disaggregate_check.md) |
| Evaluation and comparison | [`mt_backtest()`](reference/mt_backtest.md) (rolling origin, no leakage; rates, life expectancy, ASR), [`mt_score()`](reference/mt_score.md) (CRPS, WIS, coverage, PIT, MAE, RMSE, MASE; also of gaps between strata), [`mt_compare()`](reference/mt_compare.md) and [`mt_score_table()`](reference/mt_compare.md) (several models in one backtest, by level of a hierarchy), [`mt_dm_test()`](reference/mt_dm_test.md) (Diebold-Mariano tests of differences in accuracy), [`mt_compare_forecasts()`](reference/mt_compare.md), [`mt_plot_compare()`](reference/mt_compare.md), [`mt_plot_backtest_design()`](reference/mt_plot_design.md), [`mt_plot_tuning()`](reference/mt_plot_design.md), [`mt_plot_pit()`](reference/mt_plot_design.md) |
| Decomposition | [`mt_decompose()`](reference/mt_decompose.md) (differences by age, cause and subgroup: Arriaga, Horiuchi, stepwise replacement, Kitagawa, and Arriaga meets Kitagawa for subgroup composition; on observed data, forecast draws or both), [`mt_plot_decomposition()`](reference/mt_plot_decomposition.md), [`mt_pessimism()`](reference/mt_pessimism.md) (how much more pessimistic summed cause forecasts are, Wilmoth 1995) |
| Excess mortality | [`mt_excess()`](reference/mt_excess.md) (observed minus expected deaths or life expectancy, with intervals from a baseline model; by age, stratum and year, cumulative), [`mt_plot_excess()`](reference/mt_excess.md) |
| Results and plots | [`life_table()`](reference/life_table.md), [`life_expectancy()`](reference/life_expectancy.md), [`prob_dying()`](reference/derived.md), [`asr()`](reference/derived.md), [`yll()`](reference/derived.md), [`mt_contrast()`](reference/mt_contrast.md), [`mt_inequality()`](reference/mt_inequality.md), [`lifespan_variation()`](reference/lifespan_variation.md), [`mt_plot_forecast()`](reference/mt_plot.md), [`mt_plot_life_expectancy()`](reference/mt_plot.md), [`mt_plot_contrast()`](reference/mt_plot_contrast.md), [`mt_plot_scores()`](reference/mt_plot_more.md), [`mt_plot_shares()`](reference/mt_plot_more.md), [`theme_mortensor()`](reference/mt_plot_style.md), [`mt_export()`](reference/mt_export.md), [`as_demogdata()`](reference/as_demogdata.md), [`as_stmomodata()`](reference/as_demogdata.md) |

## Function map

The objects of an analysis and the functions that connect them. The
[function map](vignettes/articles/function-map.Rmd) (on the website:
Articles, Function map) shows every function with its purpose, the
possible workflows and what the package does not cover.

![Function map: your tables become an mt_data object with mt_data();
mt_fit() fits a model specification to it; mt_forecast() gives mt_draws,
simulated future rates; mt_backtest() gives backtest scores (mt_score(),
mt_dm_test()); results such as life_expectancy(), lifespan_variation(),
asr(), mt_contrast(), mt_decompose(), mt_excess(), plots and CSV files
come from the draws or from the observed
data.](reference/figures/README-function-map.png)

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

Vignettes:

- **Start**: *Getting started*, *Your own data, from CSV to report*.
- **Concepts, one at a time**: *Ungrouping age groups: methods and
  diagnostics*, *Smoothing death rates*, *Where a forecast starts: the
  jump-off*, *Backtesting and scores*, *Overdispersion: extra variation
  in death counts*, *Tuning and recalibration*, *Coherent forecasts:
  totals and their parts*, *Dependence between groups*, *Decomposing
  differences in life expectancy*, *Lifespan variation*.
- **Models**: *Choosing and combining models* (implemented and wrapped
  models side by side, ensembles of any combination), *Model reference*,
  *Bayesian age-period-cohort models with bamp*.
- **Benchmarks**: *Benchmarks: how well the pipeline works*, *Comparison
  with demography, bamp and StMoMo*.
- **Applications**: *Inequality in England: life expectancy and causes
  of death by deprivation* (ONS), *Causes of death by deprivation:
  coherent and combined forecasts* (ONS), *Differences between the
  German Laender* (Destatis), *Excess mortality*.
- **Help**: *Frequently asked questions*, *Function map*.

## Development

See `CLAUDE.md` for conventions and commands, and
`dev/implementation-plan.md` for the plan.

## Licence

GPL-3.
