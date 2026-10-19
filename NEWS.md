# mortensor 0.1.9000 (development of 0.2.0)

Clean break: the 0.1.x interface is replaced (see `dev/revamp-plan.md`).

## Testing

* The performance validation suite is removed: the `validation/`
  directory, its GitHub workflow, the benchmark helpers and the slow tests
  with thresholds on coverage, CRPS, accuracy and runtime. Forecast
  performance depends on the data; check it on your own data with
  `mt_backtest()`.
* The correctness tests stay: agreement with reference implementations
  and published values, exact identities, coherence and parameter
  recovery. Checks of the former suite that were not yet unit tests are
  now unit tests.

## Methods checked against the full papers; new published models

Work package E4. The full texts of papers that were earlier available only
as abstracts were read (sources and access in `dev/methods-audit.md`), and
the implementations follow them. Results change.

Corrected to the papers:

* `coda()` follows Bergeron-Boucher et al. (2017): the time index is
  ARIMA(0,1,1) with drift (the random walk is `order = c(0, 1, 0)`), the
  common component of the coherent model comes from the mean death rates,
  deviations are ARIMA(1,1,0), rates from age 80 are a Kannisto model to
  age 120 (`open = "auto"`; the regression of the open interval is
  `open = "regression"`), and intervals come from the paper's residual
  bootstrap. Fitted values and point forecasts reproduce the authors' code
  to 1e-9. 95% coverage on England & Wales men 0.888 (0.981 before); on
  French women the coherent model is more accurate (CRPS 0.0743 against
  0.0764; open interval 0.039 against 0.062).
* `coda_causes()` has `rank`, `order` and `drift` (Oeppen 2008 used rank 2
  and ARIMA(0,2,2)); `tucker3()` uses the ranks of Russolillo, Giordano
  and Haberman (2011), c(2, 2, 2), and reports their fit measure.
* Non-negative reconciliation is the algorithm of Wickramasuriya, Turlach
  and Hyndman (2020); the former heuristic, which could miss the optimum,
  is removed. With MinT and fixed nodes, the other nodes take covariances
  conditional on the fixed nodes. Results agree with a quadratic program
  to 4e-10.
* `life_table(open = "kannisto"/"gompertz")` fits the closure by maximum
  likelihood, as Thatcher, Kannisto and Vaupel (1998); their published
  parameters are reproduced within 0.11%, the former least-squares fit was
  off by up to 82%.
* `coherent_fdm()` uses weighted principal components (`weighted = TRUE,
  beta = 0.05`), as Hyndman, Booth and Yasmeen (2013): CRPS on France 0.0550
  against 0.0596, 95% coverage 0.802 against 0.760.
* `mt_ungroup(open_split = "pclm")` (default): the exposure of an open grid
  interval is split by the PCLM of the exposures, as Rizzi, Gampe and
  Eilers (2015); the stationary split is `open_split = "stationary"`.
* `mt_dm_test()`: the equal-weighted cosine variance is that of Lazarus et
  al. (2018); `varestimator = "daniell"` gives the test of Coroneo and
  Iacone (2020), oversized with overlapping test years.

New published models:

* `p_splines()`: two-dimensional Poisson P-spline forecasting (Currie,
  Durban and Eilers 2004), adapted from MortalitySmooth 2.3.4, whose fitted
  and forecast rates it reproduces to 1e-11. On England & Wales men CRPS
  0.035 against 0.071 for Lee-Carter; its life-expectancy intervals are
  narrow (95% coverage 0.40).
* `cp_splines()`: smooth constrained mortality forecasting (Camarda 2019),
  P-spline forecasts with asymmetric penalties on the rates of ageing and
  of change; reproduces the author's code to 2e-11. On French women CRPS
  0.066 against 0.085 for Lee-Carter; its intervals hold fit uncertainty
  only, as in the paper.
* `bohk_rau()`: varying age-specific survival improvements (Bohk-Ewald and
  Rau 2017), Bayesian, with the paper's forecast limits and quantile
  draws. On England & Wales men CRPS 0.045 (linear core 0.043) against
  0.071 for Lee-Carter, 95% coverage of life expectancy at 30 0.90.
* `stmomo("plat")`: the Plat (2009) model as specified in the StMoMo
  vignette (Villegas, Kaishev and Millossovich 2018); CRPS 3% below
  Lee-Carter at ages 55-89.
* `nordpred_apc()`: the Nordpred age-period-cohort model with a power link
  and cut trend (Møller et al. 2003), written for mortensor; projections
  agree with the nordpred package to 3e-8. Its draws hold coefficient
  uncertainty only and are far too narrow (95% coverage 0.11).

Fixes: quantiles of draws are interpolated as `a + w (b - a)`, so tied
draws no longer give a median above the upper quantile; errors of rejected
candidate models of the forecast package no longer appear in knitted
documents with `coherent_fdm()`.

Documentation: `dev/generalisation-plan.md` plans experimental
generalisations of the deprecated package methods, grounded in verified
literature.

## Published methods as defaults (audit of package-specific methods)

Every method was checked against its published source
(`dev/methods-audit.md`; article "Methods and sources"). Defaults are now
the published methods; package-specific options are labelled
"package-specific (reference unknown)" in their help pages; package
inventions with a published counterpart in the package are deprecated
(they work and warn with class `mt_warning_deprecated`). Results change.

Predictive draws:

* Deaths are Poisson given the model's rates (Brouhns, Denuit and
  Vermunt 2002). The extra-Poisson layer estimated from in-sample
  residuals is removed: it had no published source, counted observation
  noise twice for `tucker3()` and the functional data models, replaced
  the posterior predictive of the Bayesian models, and the validation
  generator added the same noise. `stanmomo(family = "nb")` uses its
  posterior negative-binomial predictive. Intervals of models of smooth
  rates are narrower: 95% coverage of `lee_carter()` on England & Wales
  men (ages 30-100, check `backtest_ew`) falls from 0.746 to 0.632; on
  simulated Lee-Carter data with Poisson deaths it is 0.944.

Models:

* `tucker3()` is centred as in the three-way Lee-Carter model of
  Russolillo, Giordano and Haberman (2011), without the residual noise
  layer; experimental. 95% coverage on simulated Lee-Carter data 0.748.
* Deprecated: `product_ratio()` (use `coherent_fdm()`), `cause_shares()`
  (use `reconciled()` over a cause hierarchy with `additive = "rates"`, or
  `coda_causes()`), `rwd(observation_noise = TRUE)` (use
  `arima_series(order = c(0, 1, 1), drift = TRUE)`, its published
  reduced form), `tucker3(transform = "coda")` (use `coda()`).
* The Li-Lee AR(1) without intercept clamped at 0.99 and the ridge and SVD
  fallbacks of `lee_carter()` for sparse counts are documented and warn.

Combining models and coherence:

* `reconciled()`, `mt_reconcile()`: `dependence = "none"` by default; the
  reordering of separately simulated strata (`"empirical"`) is opt-in and
  labelled package-specific. Aggregate intervals of separate models are
  narrower: 95% coverage of the England total (ONS, WLS) 0.434 against
  0.552. WLS and MinT weights from the base draws, the non-negativity
  heuristic, fixed nodes and proportional reconciliation are documented
  with their sources and deviations.
* `reconciled()` warns (class `mt_warning_hierarchy_exposure`) when the
  exposures of aggregate strata differ from those the hierarchy implies,
  for example causes of death declared without `additive = "rates"`.
* `recalibrated(scale = "constant")` by default (member-by-member
  post-processing, Van Schaeybroeck and Vannitsem 2015); the power law
  over horizons is an option. `quantity = "asr"` standardises as
  `mt_backtest()` does.
* `anchored()` defaults to the jump-off of Lee and Miller (2001):
  `years = 1, decay = Inf, shrink = FALSE`; averaging, linear decay and
  shrinkage (DerSimonian-Laird) are package-specific options. On the ONS
  quintiles the default is less accurate than the former settings (CRPS
  of Li-Lee 0.0708 against 0.0656).
* Attributions corrected: `ensemble(weights = "crps")` is package-specific;
  stacking weights optimise the medians; `bootstrapped(type = "grid")` is
  the residual bootstrap of Renshaw and Haberman (2008).

Harmonisation:

* `mt_ungroup(preserve_totals = FALSE)` by default (the published PCLM);
  the rescaling to group totals is opt-in. The open age interval inside
  the PCLM is iterated to convergence and documented as package-specific.
  Backtests keep exposure totals when they ungroup exposures.
* Deprecated: ungrouping methods `"hyman_kannisto"`, `"natural_gompertz"`
  and `"gam"`, `mt_shape_rates(method = "gam")`; `mt_disaggregate()`
  defaults to proportional allocation by a standard age schedule
  (`method = "shape"`).
* `mt_smooth()`: the weighted Whittaker smoother with smoothing parameters
  chosen by GCV per age, and a quasi-Poisson GAM in 1D as in 2D.
* `mt_fill_suppressed()`: cells without a total get the truncated-Poisson
  mean.

Evaluation and derived quantities:

* `mt_dm_test(varestimator = "ewc")` by default: the fixed-smoothing test
  with the equal-weighted cosine variance and t critical values (Lazarus
  et al. 2018). In simulations under the null
  (help page; check `dm_size`), the HLN-corrected test (`"acf"`) rejected
  in 10-19% of samples at nominal 5% with 7 overlapping origins and
  horizons 2-5, the EWC test in 4.8-7.4% with 7-16 origins. p-values are larger: random walk against Li-Lee at
  horizon 1, 0.049 against 0.00015.
* `asr()` standardises directly on the age groups of the standard
  (Eurostat 2013) with `asr_groups()` and `asr_columns()`; the former even
  split of standard weights over finer ages (`age_groups = "grid"`) gave
  ASRs 1.8-11% too high on French single ages.
* `yll()` uses life expectancy at the mean age at death in the interval,
  e(x + a_x) (Murray 1994): 8.5-9.9% lower with 5-year groups.
* `mt_decompose(composition = "pooled")` is labelled package-specific;
  the estimand of `mt_excess()` and its P-scores are specified.

Documentation:

* New article "Methods and sources": every method with its source and
  status (published, published with deviations, package-specific).
* All vignettes use the published defaults and were re-knitted; the
  package comparison was re-run. Results that changed: in "Choosing and
  combining models" the Diebold-Mariano test of life expectancy over 16
  origins separates only the random walk from M6 over all horizons
  (p = 0.029); the
  intervals of Lee-Carter, CBD and the functional data model are
  narrower (95% coverage 0.69, 0.71 and 0.46). In the England inequality
  analysis the Lee-Miller jump-off made rate forecasts worse and is not
  used. In the German Laender analysis Lee-Carter forecasts the East-West
  gap best. "Your own data" pairs Lee-Carter with the random walk in its
  ensemble.
* Warnings of joint models name the message once instead of an example
  stratum.

## Published methods added

Models:

* `coda()`: the compositional data model of life-table deaths (Oeppen
  2008; Bergeron-Boucher et al. 2017), adapted from MortalityForecast
  0.8.1, whose fitted values and point forecasts it reproduces. The
  coherent model of several populations (`coherent = TRUE` or
  `common_by`) has a common component and stratum deviations that follow
  AR(1) processes. The open age interval follows the indexes. Validation
  (checks `coda_ew`, `coda_coherent_fr`): CRPS on England & Wales men
  0.055 against 0.069 for `lee_carter()`; coherent model on France 0.077
  against 0.080 for `li_lee()`.
* `coda_causes()` (experimental): the multiple-decrement form of the CoDa
  model for causes of death, which add up to all causes in every draw.
  Less accurate than `cause_shares()` on simulated causes.
* `lee_carter(rotate = TRUE)`: rotation of the age pattern of decline as
  life expectancy rises (Li, Lee and Gerland 2013, as in MortCast); life
  expectancy at birth stays that of the unrotated forecast. On France,
  CRPS 1% lower than without rotation 10 years ahead and 2% lower 20
  years ahead (check `rotation_france`).
* `stmomo("m7")` and `stmomo("m8")`: the cohort extensions of the CBD
  model (Cairns et al. 2009), with `gc_order` and `xc`.
* `reconciled(method = "energy_score")`: score-optimal reconciliation
  (Panagiotelis et al. 2023), learnt from the base model's forecasts at
  inner origins; it reproduces `ProbReco::scoreopt()` to 4e-15. In
  backtests on simulated regions and on the ONS data it was less accurate
  than WLS (CRPS 11-32% higher), which the help page and the coherence
  vignette say.

Evaluation:

* `mt_dm_test()`: Diebold-Mariano tests (with the correction of Harvey,
  Leybourne and Newbold 1997) of the difference in accuracy between two
  models of the same backtest, per horizon, allowing for forecasts that
  overlap in their test years. Agrees with `forecast::dm.test()` to 3e-14.
  Horizons with too few origins for the autocovariances get no test.

Derived quantities:

* `lifespan_variation()`: life disparity (e-dagger), entropy, Gini
  coefficient, standard deviation, variance, coefficient of variation and
  interquartile range of the age at death, for rates, observed data and
  forecast draws, conditional on survival to an age `at`. Adapted from
  LifeIneq 00.04.09, which it matches to 1e-15 on tables with intervals of
  equal width. The help page quantifies how much the measures depend on
  the open age interval.
* Lifespan variation is a quantity of `mt_compare_forecasts()`,
  `mt_plot_compare()`, `mt_contrast()` and `mt_inequality()` (`quantity =
  "lifespan_variation"`, `measure`).
* `mt_excess()` and `mt_plot_excess()`: excess deaths, P-scores and the
  loss of life expectancy against the predictive draws of a baseline
  model, by age group, stratum and year, and summed over years. Cells with
  estimated deaths are left out and counted. France 2003: 15,100 excess
  deaths (95% interval 5,600 to 25,400). Intervals covered the true excess
  in 97-98% of simulated cases (check `excess_simulation`).

Vignettes and documentation:

* New vignettes "Lifespan variation" (France: observed, by age, the
  effect of the open age interval, forecasts of Lee-Carter, rotated
  Lee-Carter and CoDa) and "Excess mortality" (France 2003 and 1918,
  German Laender 2020-2024).
* "Backtesting and scores" tests the differences between models with
  `mt_dm_test()`; "Choosing and combining models" adds CoDa and M7 and a
  test from 16 origins; "Coherent forecasts" describes score-optimal
  reconciliation; "Model reference", FAQ, README and the function map
  (diagrams redrawn) cover the new functions.

## Fixes from a review of the whole package

Data:

* `mt_group_levels()` takes numeric `total`, `other` and `keep` levels as
  labels, keeps an unlisted level named like a group in that group (as
  `mt_relabel()` to an existing level does), and refuses a group named
  like the total.
* `mt_read_gbe()` reads values with several thousands separators
  (`"1.023.687"`) and footers in quotes.
* `mt_fill_suppressed()` no longer fills cells of strata published at
  another age resolution.
* `mt_data()` treats `NA` ages as unknown (distributed) and leaves out
  rows of unknown period with a warning; factor age bounds use their
  labels. `mt_from_hmd()` recognises open ages written as `"110+"`.

Harmonisation:

* `mt_ungroup()` keeps an observed first year of life out of the fit (its
  peak spilled into ages 1-6), works with the spline methods on grouped
  exposures, and warns when `dims = 2` is ignored
  (`mt_warning_ungroup_dims`). The spline and GAM methods keep the group
  totals next to a missing group and above a closed top group.
* `mt_smooth()`: cells with missing deaths no longer break the smooth;
  they stay missing (only the 2D GAM predicts them). Smoothing over time
  works with a single age group. `mt_smooth_check()` measures the smoothed
  cells only.
* `mt_disaggregate()`: `method = "shape"` with a shape of zeros gives a
  flat shape instead of `NaN`; missing group rates give missing fine
  rates.
* `step_shocks(method = "interpolate")` removes shock years at the start
  or end of the training data instead of failing.

Models:

* Models of observed rates (`rwd()`, `arima_series()`, `ets_series()`)
  fitted to smoothed rates (`step_smooth()`) now add the observation noise
  of the raw rates to predictive draws; their backtest intervals were too
  narrow.
* Each stratum is fitted with its own random number stream started from
  the current generator: fits that use random numbers follow `set.seed()`
  with any number of `cores` (bootstrap results differ from earlier
  versions).
* `rwd(drift_window = )` with `observation_noise = TRUE` no longer fails;
  `lee_carter(bootstrap = TRUE)` resamples the chosen fitting period;
  `lee_carter(adjust = "e0")` with age groups without deaths matches life
  expectancy again (it could silently leave the index unadjusted);
  `stmomo()` forecasts one year; `tucker3()`, `li_lee()` and
  `product_ratio()` run with a single age group; `coherent_fdm()` reports
  groups with fewer than two strata.
* `cause_shares()`: predictive draws are those of the total model, with
  its deaths split between the causes, so the causes add up to the total
  in every draw; the extra variation of the shares beyond multinomial
  sampling is estimated as such.

Wrappers and hierarchies:

* `reconciled(anchored(...))` and `reconciled(recalibrated(...))`
  reconcile again (the inner wrapper's fields replaced those of
  `reconciled()`); nested wrappers keep their own fields.
* Observation noise of reconciled predictive draws is independent across
  cells: the empirical copula reorders the paths only. The noise had been
  made comonotone across strata, which widened the intervals of the
  aggregates (validation `de_regions`: 95% coverage of the aggregates
  0.76, before 0.80).
* Exposures given for the bottom strata only are summed to the
  aggregates; a reconciled model backtested on bottom-level data scores
  the right strata (the draws of other strata had been recycled).
* Hierarchies of one variable are reconciled on the scale of deaths when
  the exposures of the parts add up to the total, on the scale of rates
  when each part has the total's exposure (causes of death).
* `ensemble()`: inner backtests are matched cell by cell, components that
  forecast different strata are refused, and the ensemble fit no longer
  carries the first component's calibration or jump-off. Ensembles and
  `tuned()` add observation noise according to their components.
* `stratified()` with a level named `other`, `recalibrated(seed = NULL)`,
  `bootstrapped(tuned(...))`, `mt_coherence()` with one year or one draw,
  and the component of each draw after reordering are fixed.
* `mt_hierarchy()` checks its groups; `mt_aggregate_strata()` warns about
  members missing from the data. `mt_tune()` refuses ungrouped data with
  the reason.

Evaluation and derived quantities:

* Backtests floor the predictive deaths at half a death after summing them
  to the observed cells, not per grid cell: with ungrouping inside the fit
  (`step_ungroup()`), the floor per single year had biased the scores of
  age groups with few deaths. Predictive draws may now hold zero deaths.
* Life expectancy and ASR in backtests are computed from the predictive
  deaths like the observed values, and only for years whose cells cover
  the whole age range.
* `mt_backtest()` refuses data whose deaths were ungrouped with
  `mt_ungroup()` (`mt_error_ungrouped_data`): use `step_ungroup()`.
  Test data are summed to the age groups of a fit with `step_regroup()`.
* `mt_compare()` chooses the default origins once for all models;
  `mt_score_table()` keeps a model that lacks one target and separates
  preprocessing pipelines for every target.
* `life_table()`: a rate of 0 in the open interval uses the last positive
  rate (instead of an infinite life expectancy).
* `mt_inequality()` leaves the total and groups of a stored hierarchy out
  of the default `order`.
* Plots: `mt_plot_compare()` with one year, `mt_plot_contrast()` with
  other probabilities, `mt_plot_age()` with unsorted `probs`.

Vignettes:

* The precomputed vignettes are re-knitted. With observation noise
  independent across cells, their backtests show that the intervals of
  life expectancy are too narrow for every model, most for national
  totals; the texts say so. In the German Laender, the ensemble of
  Lee-Carter and the coherent fdm now also forecasts the East-West gap
  best (with Lee-Carter).
* "Bayesian age-period-cohort models with bamp": the effects plot leaves
  out the R-hat rows of `mt_parameters()`, which get their own table; the
  MCMC fits are seeded.
* The README and the function map article show their diagrams as images
  (drawn from `dev/diagrams/*.mmd` by `dev/diagrams/render.js`), because
  the GitHub mobile app does not draw Mermaid diagrams; the website no
  longer loads Mermaid.

## Engineering

* `mt_parameters()` works for `ensemble()` (the parameters of each
  component, column `component`, and the weights), `stratified()` and
  `tuned(by = )` (column `group`), also inside other wrappers.
* `mt_backtest()` and `recalibrated()` run each origin with its own random
  number stream, as `mt_fit()` does for strata: results of models that use
  random numbers (bootstrap) no longer depend on `cores`.
* `coherent_fdm()` forecasts about three times faster: the ARFIMA models of
  the ratios are prepared once per coefficient instead of once per simulated
  path (`forecast::simulate()` refits them at every call). The draws are
  identical.

## StanMoMo, recalibration

* `stanmomo()` wraps the Bayesian Lee-Carter, Renshaw-Haberman, APC, CBD
  and M6 models of StanMoMo (Barigou et al. 2023), Poisson or negative
  binomial; draws are posterior draws of the expected rates, fits warn when
  the chains have not converged. Validation check `stanmomo_ew`.
* `calibrated()` is renamed `recalibrated()`, and the vignette "Tuning and
  calibration" is now "Tuning and recalibration": in demography,
  "calibrating" a model often means fitting it, while this wrapper adjusts
  the forecast distribution of a fitted model afterwards.
* `hyndman_ullah()` keeps smoothing the rates over age by default.

## Vignettes restructured, models and diagnostics added

* The vignettes are reorganised: one concept per vignette with a simple
  model, a model vignette, and two applied pipelines. New: "Choosing and
  combining models" (implemented and wrapped models side by side,
  ensembles of any combination), "Ungrouping age groups: methods and
  diagnostics" (replaces "Data at different granularities"), "Inequality
  in England" and "Differences between the German Laender" (replace the
  deprivation, causes, regions and HMD forecasting vignettes). Every step
  of the applied pipelines is chosen in a backtest of the quantity the
  question is about.
* `coherent_fdm()`: the product-ratio method with functional data models
  of demography (Hyndman, Booth and Yasmeen 2013), with simulated paths,
  groups (`common_by`) and hierarchies.
* `stratified()`: different models for groups of strata, for example
  Lee-Carter for all causes and Li-Lee for the causes of death, inside
  `reconciled()` with the causes scaled to the all-cause forecast.
* Backtests score gaps between strata path by path:
  `mt_backtest(contrast = c(var = , ref = ))` adds the targets
  `"life_expectancy_gap"` and `"asr_gap"`. `mt_score_table(hierarchy = )`
  scores the bottom series, aggregates and total separately.
* `mt_plot_ungroup()` compares several ungrouping results (a named list)
  with three views: age profile against the groups, stability over time,
  and errors against known single years with the Poisson noise band;
  `mt_ungroup_check()` gives the same as a table. `mt_regroup(measures = )`
  regroups only some tables (for example the deaths, to check ungrouping).
* Help pages no longer quote benchmark results; the vignettes show them.
* Fixes: StMoMo cohort models now have fitted rates for clipped cohorts
  (their fitted life expectancy dropped to zero in the last years); a
  missing rate gives a missing life expectancy instead of certain death;
  `mt_smooth(along = "period")` uses Whittaker smoothing by default;
  errors raised inside `do.call()` no longer print the whole function;
  `mt_inequality()` leaves out levels not in `order` (such as a hierarchy
  total); the legend order of `mt_plot_shares()` is fixed.

## Comparing models and forecasts

* `mt_compare()` backtests a named list of models (each with its own
  preprocessing) on the same origins; `mt_score_table()` shows the scores
  of rates, life expectancy and ASR side by side.
* `mt_compare_forecasts()` tabulates one quantity (life expectancy, ASR,
  deaths, rates) of several forecasts or of observed data as "median
  (lower-upper)"; `mt_plot_compare()` draws them in one figure, also
  around the jump-off year.
* `mt_plot_backtest_design()` shows the training and test years of a
  backtest; `mt_plot_tuning()` the scores of a tuning grid.
* Precomputed vignettes no longer show "plot of chunk" captions.
* `mt_plot_smooth()` compares several smoothers in one figure (a named
  list) and shows rates over age (`along = "age"`); `mt_plot_pit()` draws
  PIT histograms of a backtest.
* Smoothing over time (`mt_smooth(along = "period")`, `step_smooth()`)
  scales the weights within each age: one smoothing parameter now gives
  the same smoothness at every age. Before, the ages with many deaths were
  barely smoothed (effective degrees of freedom 15 of 19 at ages 80+ on
  the ONS data) and those with few deaths strongly.
* New concept vignettes, one topic each: "Smoothing death rates", "Where
  a forecast starts: the jump-off", "Backtesting and scores" (with the
  scores explained); the articles are grouped into getting started,
  concepts, models, applications and help.
* Vignette "Tuning and recalibration": the fitting window of Lee-Carter
  chosen by a hold-out, the width of intervals calibrated, both checked in
  a backtest.
* `hyndman_ullah()`: when the monotonic spline of `smooth.demogdata()`
  cannot be fitted (age ranges starting in middle age), the rates are
  smoothed without the constraint, noted in the fit diagnostics.
* `mt_tune()`: standard errors also for a hold-out of a single stratum
  (blocks of test years).
* `mt_plot_contrast()` draws gaps between strata over time (one or
  several forecasts, with observed gaps); `mt_plot_compare(periods = )`
  with one year shows medians and intervals per stratum side by side.
* "Getting started" is one short pass through the main steps;
  "Coherent forecasts: totals and their parts" (live, no longer
  precomputed) explains bottom-up, reconciliation and fixed totals with
  one figure each and a small backtest.

## Tuning and recalibration

* `hyndman_ullah()` smooths the rates of each year over age with
  `demography::smooth.demogdata()` before the decomposition (`smooth =
  TRUE`, settings in `smooth_args`), as in Hyndman and Ullah (2007).
  Before, it decomposed raw rates, contrary to its documentation, and the
  model errors that `simulate()` resamples contained the observation noise
  that predictive draws add a second time. `smooth = FALSE` keeps the old
  behaviour (label "unsmoothed").
* `mt_tune()` and `tuned()` choose settings by forecast accuracy on the
  last years of the training data:
  * `spec_fun` may return `list(spec = , prep = )`, so preprocessing
    (smoothing parameter, fitting window) is tuned together with the model
    (for example the number of components of `hyndman_ullah()`);
  * `inner_origins = 1` is a single hold-out of the last `inner_h` years;
    `target = "life_expectancy"` or `"asr"`; `metric = "rmse"` for point
    accuracy;
  * `select = "one_se"` takes the simplest setting within one standard
    error of the best; scores carry standard errors, scores by horizon
    (`by_horizon`) and the `chosen` setting;
  * `by` tunes and fits groups of strata separately (for example by sex).
* `mt_score()` reports the root mean squared error of the median (`rmse`).
* `recalibrated()`, `anchored()` and `reconciled()` around `tuned()` failed
  when forecasting ("No simulate_model() method for model 'tuned'"): the
  wrappers rebuilt the base model from its specification instead of the
  one the tuning chose. Fits now keep the specifications fitted underneath.
* `summary()` of tuned fits shows the chosen settings.
* `recalibrated()`:
  * `ages` estimates the spread factors (and drift) per age band, for
    models too narrow at some ages and too wide at others;
  * `quantity = "life_expectancy"` or `"asr"` calibrates the intervals of
    these quantities, which also depend on how the ages move together;
  * `n_boot` gives bootstrap intervals of the factors (forecasts of
    stratum x origin resampled);
  * the drift correction of `bias = TRUE` stays constant beyond the longest
    horizon of the calibration window instead of growing further;
  * a warning when the factors rest on the forecasts of one origin; the
    table reports the numbers of inner forecasts and scored values.
* Benchmark of tuning and calibration (`dev/tuning.md`, validation check
  `tuning_ew`): tuning the fitting window of Lee-Carter on a 10-year
  hold-out improved forecasts on England & Wales and France; tuning the
  components and smoothing of fdm the same way did not; calibration
  brought the 95% coverage of all models close to nominal.

## Function map

* The diagrams of the function map and the README fit the width of the
  website, and their labels are no longer clipped (Bootstrap's paragraph
  margin inside mermaid labels made the boxes too small). Map 1 lists the
  preparation steps in their usual order and leaves smoothing to
  `step_smooth()` inside fits; `mt_prep()` feeds `mt_fit()`; the function
  table lists every exported function.

## Changes after the third and fourth applied reviews

Two further reviews (`dev/review-applied-3.md`, `-4.md`; an epidemiologist
and a demographer working through the vignettes on their own questions).

* `mt_decompose()`:
  * Strata that are totals of others (`mt_group_levels(total = )`, a
    hierarchy) are no longer added to their parts: such a cause or group is
    refused with a message naming it. Causes passed as `groups` (levels
    sharing one population) give a warning.
  * `groups` needs `composition = "radix"` or `"pooled"`; the
    documentation states when each applies (radix only for membership fixed
    from the start age, such as sex or education).
  * `observed = "fitted"` uses the fitted rates for observed years, so that
    a forecast is compared with its own model; `"poisson"` adds the
    sampling noise of observed years. Decompositions of an observed against
    a forecast year report `jump_off`, the part due to the observed instead
    of the fitted rates.
  * `pairing = "draw"` pairs the draws of two objects (default: independent
    shuffling). `summary(change = )` gives the change of a gap between two
    years path by path, with intervals.
  * Non-finite draws are dropped with a warning; residuals above the
    tolerance give a warning. `print()` explains the components;
    `as.data.table()` and `mt_export()` work on decompositions.
  * `mt_plot_decomposition(comparisons = )` selects comparisons; titles and
    axes name the quantity, the unit and the method.
* Life tables: where `n m / (1 + (n - a) m)` exceeds 1 (extreme forecast
  rates), `q` is `1 - exp(-n m)` instead of 1.
* Reconciliation: `method = "proportional"` keeps fixed nodes and scales
  their parts in proportion to their forecasts, draw by draw (top-down,
  Athanasopoulos, Ahmed and Hyndman 2009). Non-negative projections warn
  when more than 5% of the draws hold strata at zero or fixed nodes move;
  `$reconciliation$share_draws_with_zeros`. `mt_plot_coherence()` shows the range of
  the paths.
* `anchored()`: the empirical Bayes shrinkage subtracted the sampling
  variance from the wrong term, so offsets of noisy cells were shrunk too
  little. ONS Li-Lee anchored: CRPS 0.0654 (before 0.0658) against 0.0666;
  life expectancy at birth 0.315 against 0.336. `summary()` of anchored fits
  shows the offsets.
* `mt_fit()` warns when data were smoothed before the fit (raw deaths are
  then missing and intervals too narrow); `mt_backtest()` refuses such data
  (the smoother has seen the test years). Smooth inside the fit with
  `prep = mt_prep(step_smooth(...))`.
* Backtests carry a `prep` column; `mt_score()` scores the same model with
  different preprocessing separately and reports the `unit` of the scores.
  `print()` of backtests with life expectancy or ASR names the other
  targets.
* Period axes may have gaps (census years 2001, 2011, 2021): observed
  indicators and decompositions work on such data; models ask for
  consecutive years.
* `mt_group_levels(total = )` records totals, which fits and forecast
  draws carry, so `mt_decompose()` recognises them.
* The plausibility warning of forecasts covers the open age interval.
* The categorical palette has 12 colours (more are interpolated).
* Vignettes: all-cause forecasts implied by the causes against all-cause
  models (Wilmoth 1995) in the causes and regional vignettes; proportional
  reconciliation; decomposition of forecasts (jump-off, the model's own
  change, change of the gap, comparison of two models); calibrated
  Lee-Carter without smoothing in "Forecasting and backtesting"; corrected
  interpretations and numbers in the deprivation vignette and the FAQ.

## Decomposition

* `mt_decompose()` splits the difference in a mortality summary between
  two sets of rates into contributions of ages, causes and population
  subgroups, as post-processing of observed data (optionally with Poisson
  draws), forecast draws (draw by draw, so contributions have intervals),
  or both (an observed year against a forecast year, or two objects).
  Comparisons: between levels of a stratum variable (`var`, `ref`), between
  years (`periods`), or between two objects (`y`). Methods: symmetric
  Arriaga for life expectancy with causes split by the implied sensitivity
  (exact), the line integral of Horiuchi, Wilmoth and Pletcher (2008) and
  stepwise replacement of Andreev, Shkolnikov and Begun (2002) for any
  quantity (life expectancy, probability of dying, age-standardised rate,
  or a function of the rates), and Kitagawa's rate and age-structure
  components of the crude death rate. With `groups`, the composition of a
  population of subgroups is separated from their mortality: as in Riffe,
  Tursun-zade and Trias Llimós (2024, "Arriaga meets Kitagawa"), with
  subgroups blended at the radix (`composition = "radix"`), or for the
  pooled rates (`"pooled"`, composition effects by age).
* `summary()` sums contributions path by path over ages, causes or
  subgroups; `mt_plot_decomposition()` draws them as bars, stacked by cause
  or subgroup, with 95% intervals for forecasts.
* Results agree with DemoDecomp (`horiuchi()`, `stepwise_replacement()`)
  and LEdecomp (`arriaga_sym()`) to machine precision on the same life
  table, and the Kitagawa step reproduces Table 1 of Riffe et al. (2024).
  New validation check `decomposition`; article "Decomposing differences in
  life expectancy" on the ONS data; plan in `dev/decomposition-plan.md`.

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
  combined with `recalibrated()` it gave the best calibrated forecasts for
  England & Wales men (Lee-Carter on the GAM surface, calibrated: CRPS
  0.0646, 95% coverage 0.91 against 0.0687 and 0.75). On the ONS quintiles
  smoothing over time made the forecasts worse.
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
  for the open interval in `mt_ungroup(method = "gam")`; in
  `mt_disaggregate()`,
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

* `recalibrated()` wraps any model specification (including `ensemble()`,
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
