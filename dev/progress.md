# Progress log

Status per work package of `dev/implementation-plan.md`. Branch for this
work: `claude/vigilant-albattani-2ssw1x` (plays the role of `dev-0.2`).

| WP | Status | Commit | Notes |
|---|---|---|---|
| P0.1, P0.2, P0.5 | done | 41058ba | 0.1.1: generics and `%\|\|%` fixed, GPL-3, baseline in `dev/baseline/` |
| P0.3 | done | this branch | Workflows: R-CMD-check (R 4.3, oldrel, release, devel; macOS, Windows), test-coverage, lint, validation, pkgdown |
| P0.4 | dropped | | 0.1.x vignettes not revived; the clean break replaces them |
| P0.6 | done | this branch | `CLAUDE.md`, this file |
| Q0.1 | done | this branch | `validation/` pipeline with quick and full suites |
| D1.1 | done | this branch | Skeleton, conditions, checks, name-clash test |
| D1.2 | done | this branch | `mt_axis`, `parse_age_labels()` |
| D1.3 | done | this branch | `mt_hierarchy()` (crossed variables, intermediate groups), `mt_aggregate_strata()`, multi-variable checks in `mt_validate()` |
| D1.4, D1.5 | done | this branch | `mt_data()`, `mt_add_*()`, grid, `mt_array()`, `mt_window()`, `mt_granularity()`, `mt_validate()` |
| D1.6 | done | this branch | `mt_regroup()`, `mt_operator()` |
| D1.7 | done | this branch | `mt_as_exposure()` (mid-year, January 1) |
| D1.8 | done | this branch | Native PCLM (1D, 2D), rate-based with stationary weights in the open interval |
| D1.9 | done | this branch | `mt_ungroup()`; E&W benchmark 0.015 / 0.024 (1D), 0.015 / 0.023 (2D) |
| D1.10 | done | this branch | `mt_smooth()` (Poisson P-splines, 1D/2D), `mt_prep()` with `step_ungroup()`, `step_smooth()`, `step_regroup()`, `step_grid()`, `step_shocks()` |
| D1.11 | done | this branch | `life_table()` (interval widths, reciprocal or Gompertz open interval), `life_expectancy()` |
| D1.12 | done | this branch | `as_mt_data()` for data frames, StMoMoData, demogdata; `mt_example()`; `mt_read_ons()` (Excel or data frames; deciles to quintiles; optional region and cause groups with a closing "rest"; population years filled with a warning); tested on a synthetic fixture |
| M2.1 | done | this branch | `new_mt_spec()`, `mt_fit()`, `mt_forecast()` (incl. predictive draws), per-stratum failures in diagnostics; contract test `expect_mt_model()` |
| M2.2 | done | this branch | `mt_draws`, `mt_summarise()`, `mt_deaths()` |
| M2.3 | done | this branch | `rwd()` (bootstrap or Gaussian), `arima_series()`, `ets_series()` |
| M2.4 | done | this branch | `lee_carter()` Poisson and SVD, jump-off, bootstrap; agrees with StMoMo (6e-5) and demography (1e-6 on fitted log rates) |
| M2.5 | done | this branch | `tucker3()` native HOOI, ranks fixed or chosen by BIC/AIC (`ranks = "bic"`, port of `tucker3_info_criterion`), log and CoDa (open interval via auxiliary series), RW/ARIMA/ETS factors, residual noise; regression vs 0.1.0 not reproduced (different preprocessing; see notes) |
| M2.6 | done | this branch | `mt_backtest()` scores at the observed resolution; leakage test |
| M2.7 | done | this branch | `mt_score()`: CRPS, WIS (Bracher et al. 2021), coverage, PIT, MASE |
| M2.8 | done | this branch | `mt_tune()` (inner rolling origin) and `tuned()` (tuning refitted inside each fit, so nested backtests do not leak) |
| C4.1 | done | this branch | `ensemble()` (equal, numeric, inner-backtest CRPS weights); interval calibration evaluated in `dev/calibration.md` (mixture of Lee-Carter and `rwd()` best calibrated at 95% on real data); `rwd(observation_noise = TRUE)` added |
| B3.1-B3.4 | done | this branch | `apc_bayes()` + `mcmc_control()`: fit at the death resolution (`periods_per_agegroup`), explicit convergence check with refits, draws mapped to the grid with shape factors, `apc_effects()`; `apc_select()` ranks random-walk orders by DIC (bamp's DIC = deviance at the posterior mean + 2 pD) |
| C4.2 | done | this branch | `life_expectancy()`, `prob_dying()`, `asr()` with `standard_population("esp2013")`, `yll()`, `mt_deaths()`; all path by path |
| C4.3 | done | this branch | `mt_plot_forecast()`, `mt_plot_age()`, `mt_plot_life_expectancy()`, `mt_plot_rates()`, `mt_plot_ungroup()`, `mt_plot_coverage()` (vdiffr snapshots of all plots in `tests/testthat/_snaps/plot-snapshots/`); `mt_plot_tucker()`, `mt_plot_shares()`, `mt_plot_scores()` |
| C4.4 | done | this branch | `mt_export()` (summary or long draws, CSV) |
| C4.5 | done | this branch | Vignettes: getting started (ONS deprivation data: Lee-Carter and Li-Lee by sex and quintile, life expectancy gap path by path, backtests), causes of death by deprivation (ONS, 12 cause groups at ages 50+: bottom-up Lee-Carter most accurate for all causes, Li-Lee by sex and cause for single causes, cause shares less accurate), earlier: granularity and ungrouping (HMD E&W, France), forecasting and backtesting (HMD France, E&W), Bayesian APC with bamp (HMD France; precomputed via `dev/precompute-vignettes.R`) |
| C4.8 | done | this branch | FAQ vignette (data, granularity, ungrouping, model choice, calibration, predictive draws, small counts, COVID-19 years, coherence, reproducibility, runtime) |
| C4.7 | partial | this branch | pkgdown reference index; the migration vignette was removed (0.1 is not used any more); tagging a release is the owner's decision |
| D1.13 | done | this branch | `ons_data` (ONS ad hoc 3064, part 1, 2001-2019, OGL; `data-raw/ons_data.R` downloads and checks the md5) and `ons_cause_groups` (12 groups per sex at ages 50+, chosen so that median cells by sex, quintile, age group and year have at least about 25 deaths) |
| D1.14 | done | this branch | `mt_read_genesis()` (Destatis flat CSV); German age labels in `parse_age_labels()`; `de_deaths`/`de_population` (GBE-Bund deaths by Land and ICD-10 chapter 1998-2024, GENESIS 12411-0013 populations; `data-raw/de_data.R`); `mt_read_gbe()`, `mt_fill_suppressed()` (raking of suppressed cells to Land and national totals: 8,954 cells at ages 50+ for 1998-2019, totals matched to 5e-7); `mt_as_exposure(type = "dec31")` |
| D1.15 | done | this branch | General preprocessing instead of dataset-specific readers (`mt_read_ons()`, `mt_read_destatis()` removed): symbols for suppressed cells (`missing`, `zero`), wide tables (one value column per level), shared tables (populations shared by causes) in `mt_data()`; `mt_group_levels()` (aggregation layer for any stratum variable; summed or shared exposures; residual group from `other` or `total`; new layer name; vectorised), `mt_subset()`; default grid = coarsest exact grid; `mt_ungroup(age = )`; `mt_regroup(age = "deaths")`. ONS pipeline reproduces the former reader exactly (deaths by cause identical) |
| D1.16 | done | this branch | Applied reviews (`dev/review-applied-1.md`, `-2.md`, two independent agents in the roles of an epidemiologist and a demographer without coding experience) implemented: renames (`ensemble()`, `reconciled()`, `tuned()`, `bootstrapped()`, `mt_order_strata()`, `mt_parameters()`), per-variable additivity of hierarchies (regions x causes), bounded `mt_fill_suppressed()`, inspection (`summary()`, `residuals()`, `mt_plot_residuals()`, R-hat), observed indicators (`mt_observed()`), `mt_contrast()`, `mt_inequality()`, exposures matched by names, `mt_relabel()`, `step_window()`, `as_demogdata()`, `as_stmomodata()`, backtest origins with complete horizons, Model reference and CSV-to-report vignettes, FAQ tables, reference index by analysis step |
| D1.17 | done | this branch | Stress tests with simulated defective data (`dev/stress-tests.R`, `dev/stress-tests.md`, `tests/testthat/test-robustness.R`): 31 scenarios; every scenario ends in a correct result or a message naming the cause and remedy; Lee-Carter robust to zero-death age groups (ridge fallback), plausibility warning for forecasts, pooled level variance in `rwd(observation_noise = TRUE)` |
| D1.18 | done | this branch | Function map (`vignettes/articles/function-map.Rmd`, Mermaid diagrams rendered by pkgdown via a header script and natively on GitHub; overview in the README) |
| D1.19 | done | this branch | `recalibrated()`: recalibration of the spread (and optionally the drift) of any model's draws on the last 10 training years by minimum CRPS or interval score, nested in backtests; works with ensembles (per-draw observation noise), reconciled (bottom level, coherent) and bootstrapped models; validation check `calibration_ew` |
| V4 | done | this branch | Vignettes on deprivation quintiles (ONS) and on the effect of reconciliation methods (ONS sex x quintile with totals), precomputed; regional vignette extended with forecasts for every Land by reconciliation method |
| U0 | done | this branch | Review of the owner's earlier UK analysis (`uk_analysis.R`, methodology document, plots) and integration plan `dev/integration-plan-uk-analysis.md`: plot layout as standard, jump-off control, reconciliation with a fixed total, diagnostics; list of elements not to port. Package labelled experimental and vibe-coded (README, DESCRIPTION, package help, startup message); documentation rule in `CLAUDE.md` |
| Q1 | done | this branch | Decomposition as post-processing (`dev/decomposition-plan.md`): `mt_decompose()` for observed data, forecast draws or both (comparisons between strata, years or objects); symmetric Arriaga with cause split, Horiuchi line integral, stepwise replacement, Kitagawa (crude rate), subgroup composition at the radix (Riffe, Tursun-zade and Trias Llimós 2024) or pooled; `summary()`, `mt_plot_decomposition()`; tests against DemoDecomp, LEdecomp and Table 1 of Riffe et al.; validation check `decomposition`; article on ONS data |
| Q2 | done | this branch | Applied reviews 3 and 4 (`dev/review-applied-3.md`, `-4.md`) implemented: totals among causes and groups refused in `mt_decompose()`, `composition` required for `groups`, `observed = "fitted"`/`"poisson"`, `jump_off`, `pairing`, `summary(change = )`, export; life table `q` capped by `1 - exp(-n m)`; `mt_reconcile(method = "proportional")` and warnings of non-negative projections; corrected shrinkage in `anchored()`; smoothed data warned in `mt_fit()` and refused in `mt_backtest()`; `prep` and `unit` in scores; period axes with gaps; vignettes on all-cause forecasts implied by causes, forecast decompositions and calibration without smoothing |
| Q3 | done | this branch | Tuning and recalibration review (`dev/tuning.md`): `hyndman_ullah(smooth = TRUE)` (smooth.demogdata); `mt_tune()`/`tuned()` with preprocessing in the grid, hold-out, life expectancy and ASR targets, RMSE, standard errors, one-SE rule, `by` groups (`mt_fit_split`); wrappers around `tuned()` fixed (spec stack); `recalibrated()` by age band, on life expectancy or ASR, bootstrap intervals, drift held beyond the window; validation check `tuning_ew`; CI lint fix for lintr 3.4 |
| Q4 | done | this branch | Vignette restructure (`dev/vignette-plan.md`): one concept per vignette (smoothing, jump-off, backtesting, tuning and calibration, coherence, ungrouping with diagnostics), "Choosing and combining models", applied pipelines for England (inequality, causes) and the German Laender with every step chosen in a backtest; `coherent_fdm()` (demography product-ratio) with validation check `coherent_fdm_france`; `stratified()`; gap scores in backtests (`contrast`); `mt_score_table(hierarchy = )`; `mt_plot_ungroup()` types and `mt_ungroup_check()`; `mt_regroup(measures = )`; comparison helpers (`mt_compare()`, `mt_plot_compare()`, `mt_plot_contrast()`, ...); benchmark numbers moved from help pages to vignettes; fixes in StMoMo fitted values, life tables with missing rates, condition calls |
| Q5 | done | this branch | Owner's decisions: `stanmomo()` (StanMoMo LC, RH, APC, CBD, M6; Poisson or NB; posterior draws of expected rates; convergence warnings) with validation check `stanmomo_ew`; `calibrated()` renamed `recalibrated()` (vignette "Tuning and recalibration"); `hyndman_ullah()` keeps smoothing by default |
| Q6 | done | this branch | Review of the whole package (six reviewers, about 55 confirmed findings), fixed by area with regression tests (`tests/testthat/test-fix-*.R`): data adapters and suppression, ungrouping (age 0 kept out of the fit, group totals of spline methods), smoothing with missing cells, backtest floor after aggregation, quantity scores from deaths, ungrouped data refused in backtests, wrapper field stack (`reconciled(anchored())`), independent observation noise in reconciled draws, bottom-level exposures, per-stratum random number streams, `cause_shares()` predictive counts, observation noise for observed-rate models fitted to smoothed rates; documentation corrected; precomputed vignettes re-knitted |
| E1 | done | this branch | Engineering after the review: `mt_parameters()` for ensembles and fits by groups; random number streams per origin in backtests and recalibration (bug fixed in `par_lapply(streams = TRUE)`, which passed positions instead of elements); `coherent_fdm()` forecasts 74 s to 26 s per German backtest origin with identical draws (ARFIMA simulation models built once per coefficient); README and function map diagrams as images (`dev/diagrams/`) for the GitHub app |
| E2 | done | this branch | Published methods: `coda()` (Oeppen 2008; Bergeron-Boucher et al. 2017; adapted from MortalityForecast 0.8.1) with the coherent model, `coda_causes()` (experimental), `lee_carter(rotate = TRUE)` (Li, Lee and Gerland 2013), `stmomo("m7"/"m8")`, `reconciled(method = "energy_score")` (Panagiotelis et al. 2023; less accurate than WLS in both backtests), `mt_dm_test()` (Diebold-Mariano with HLN correction), `lifespan_variation()` (adapted from LifeIneq 00.04.09) also as a quantity of comparisons, `mt_excess()`/`mt_plot_excess()`; vignettes "Lifespan variation" and "Excess mortality"; 12 validation checks added (`coda_*`, `rotation_france`, `dm_ew`, `stmomo_cbd_ew`, `lifespan_reference`, `excess_*`, `scoreopt_*`). Not implemented by decision: TOPALS, models of life expectancy, population projection |
| E3 | done | this branch | Audit of published methods against package-specific solutions (`dev/methods-audit.md`, article "Methods and sources"): defaults are the published methods (Poisson predictive draws without the extra-Poisson layer; `stanmomo(family = "nb")` posterior NB; centred `tucker3()`; `dependence = "none"`; `recalibrated(scale = "constant")`; Lee-Miller `anchored()`; `mt_ungroup(preserve_totals = FALSE)`; EWC Diebold-Mariano; ASR on the standard's age groups; YLL at e(x + a_x)); deprecated with warnings: `product_ratio()`, `cause_shares()`, `rwd(observation_noise = TRUE)`, `tucker3(transform = "coda")`, spline ungrouping methods; package-specific options labelled; exposure check of hierarchies in `reconciled()`; validation checks `dm_size`, `coverage_local_level`, `causes_reconciled` replace checks of deprecated methods; candidate published methods listed in the audit |
| E4 | done | this branch | Methods checked against full papers (literature verified against Crossref): `coda()` as Bergeron-Boucher et al. 2017 (ARIMA(0,1,1), Kannisto 80+, bootstrap), `coda_causes()` options of Oeppen 2008, `tucker3()` ranks of Russolillo et al. 2011, published non-negative reconciliation (Wickramasuriya et al. 2020), maximum-likelihood Kannisto closure (Thatcher et al. 1998), weighted `coherent_fdm()` (Hyndman et al. 2013), `open_split = "pclm"` (Rizzi et al. 2015), `varestimator = "daniell"` (Coroneo and Iacone 2020); new models `p_splines()`, `cp_splines()`, `bohk_rau()`, `stmomo("plat")`, `nordpred_apc()`; plan `dev/generalisation-plan.md`; models with abstract-only sources (NB Lee-Carter, APCI, CAE, SAINT) not implemented |
| T1 | done | this branch | Owner's decision: the package implements methods correctly; users check forecast performance with backtests. Performance validation suite removed (`validation/`, its workflow, `R/validation.R` benchmarks, slow tests with coverage, CRPS and accuracy thresholds); its correctness checks are unit tests (DM size grid, stepwise = Arriaga identity, ungrouping recovery at 90+, ProbReco path at 1e-8) |
| E5 | done | this branch | Owner's decisions of 2026-10-19 (`dev/generalisation-plan.md`, Section 0): correctness first, new methods opt-in. Implemented from the full papers: `lee_carter(family = "nb")` (Delwarde, Denuit and Partrat 2007; French estimates within 1% of the paper), `apci()` (Richards et al. 2019; equals a constrained `glm()`), `cae()` (Kleinow 2015; equals Lee-Carter for one population), `saint()` (Jarner and Kryger 2011 as revised by Jarner and Jallbjørn 2022), `coupled()`/`mt_couple()`/`reconciled(dependence = "ecc")` (ensemble copula coupling, Schefzik et al. 2013), `coda_multiway()` (Bergeron-Boucher et al. 2018; replacement named by the deprecated `tucker3(transform = "coda")`); CoDa review against the CoDa literature (no numerical departures; documentation); vignette "Dependence between groups"; performance study scripts removed from `dev/`; jump-off note in the Laender vignette |
| U1 | done | this branch | Methods of the earlier UK analysis (plan 4.2-4.8): `mt_smooth(along = "period")` (Whittaker over time, GCV pooled per stratum), `mt_smooth(dims = 2, method = "gam")` (quasi-Poisson `te(age, year)`), `mt_smooth_check()`, `mt_plot_smooth()`; models fitted to smoothed rates take residuals, observation noise, MASE scale and copula from the raw deaths (`fit$raw_D`); `anchored()` and `mt_jump_off()`; `fixed` nodes in `mt_reconcile()`/`reconciled()` with non-negative projection, `mt_coherence()`, `mt_plot_coherence()`; `mt_backtest(quantities = )` and `mt_score(target = )` for life expectancy and ASR, `mt_plot_scores(by = "model")`; `lee_carter(adjust = , choose_period = )` (demography's `lca` adjustments and BMS period); `ensemble(top_k = )`; `mt_tune(strategy = "two_stage")`; tail options and diagnostics of disaggregation (`mt_shape_rates()`, `mt_disaggregate_check()`, `mt_plot_disaggregation()`); validation check `smoothing_jumpoff`; vignettes updated (HMD, deprivation, coherence, causes) |
| V3 | done | this branch | Comparison vignette (fidelity and accuracy vs demography, bamp, StMoMo; stored results of `dev/compare-packages.R`, full run 2026-10-13); regional vignette on German Laender (precomputed); validation check `de_regions` |
| D1.12 (HMD) | done | this branch | `mt_from_hmd()`, `mt_read_hmd()` (HMDHFDplus, credentials from `HMD_USER`/`HMD_PASSWORD`) |

| P5.1 | done | this branch | `hyndman_ullah()` (demography; draws via `simulate()`); 95% coverage 0.90 on the benchmark, 0.73 on France |
| P5.2 | done | this branch | `stmomo()` APC, CBD, RH, M6, LC with integer ages, initial exposures for logit links, `simulate()` draws; RH fails to converge at some origins (recorded) |
| P5.3 | done | this branch | `cause_shares()` (Tucker4 on ILR cause shares, any all-cause model), `mt_example_causes()` |
| P5.4 | done | this branch | `mt_reconcile()` for grouped hierarchies on the death scale (WLS default, MinT, structural, OLS, bottom-up), empirical copula, non-negative; `reconciled()`; `mt_aggregate_draws()` |
| P5.5 | done | this branch | `mt_order_strata()`, `li_lee()`, `product_ratio()` (Lee-Carter-type version) |
| P5.6 | done | this branch | `ensemble(weights = "stack")` |
| P5.7 | done | this branch | Whittaker and Gaussian GAM (mgcv) smoothing in `mt_smooth()`; the GAM-Bayes uncertainty is generalised as `bootstrapped(type = "data")` (Poisson resampling of the observed deaths, preprocessing and model refitted) and `type = "grid"` (semi-parametric, Brouhns et al. 2005); Lee-Carter grid bootstrap 95% coverage 0.934 on the calibration benchmark |
| P5.8 | done | this branch | `mt_ungroup(method = "hyman_kannisto" / "natural_gompertz" / "gam")` (rate-based, exposure-weighted group rates preserved exactly) and `mt_disaggregate()` for forecast draws. Median absolute log error, E&W 30-89/90-100 and simulated 30-89/90+: PCLM 0.015/0.024, 0.012/0.006; Hyman-Kannisto 0.015/0.029, 0.013/0.028; natural-Gompertz 0.015/0.027, 0.013/0.030; GAM 0.015/0.031, 0.012/0.031. The 0.1 floor of the Kannisto tail at the open-group mean was dropped (simulated 90+ error 0.166 with it). PCLM stays the default |
| P5.9 | partial | this branch | Multi-cause vignette on `mt_example_causes()` (HMD has no causes); release 0.3.0 open |
| C4.6 | partial | this branch | 2D PCLM 87 s to 45 s (standard errors only for the final fit, Cholesky solves); no further work |

## Decisions taken during implementation

- 2026-10-12: Exported verbs use the `mt_` prefix; `mt_obs()` returns an
  observation table (the class name in the plan is not used as a class).
- 2026-10-12: The open age interval is expanded into single years up to
  omega inside the PCLM; for rates the exposure of the open interval is split
  with stationary-population weights from the current fit (three outer
  iterations). This replaces the artificial bin of `ungroup`.
- 2026-10-12: Default smoothing selection BIC on log10 lambda in [-2, 8],
  knots every 2 years, cubic B-splines, second-order differences
  (alternatives tested on E&W: third-order differences and wider knots gave
  the same accuracy within 0.003).
- 2026-10-12: 2D ungrouping takes about 90 s for 71 ages x 32 years; it is
  an option, 1D per year is the default.

- 2026-10-12: Backtests score predictive draws of observed log rates: smooth
  models (Lee-Carter, Tucker3) add extra-Poisson noise estimated from
  in-sample residuals and Poisson variation given the test-year exposures;
  models on observed log rates (random walk, ARIMA, ETS) are used as they
  are.
- 2026-10-12: Calibration on the simulation benchmark (6 replicates, h =
  1-10, 95% intervals): Lee-Carter 0.959 (min 0.924), Tucker3 0.947 (0.923),
  Tucker3-CoDa 0.943 (0.917), random walk 0.987 (over-covers: it treats iid
  noise as persistent). Tucker3 needs `residual_noise = TRUE` (default;
  0.87 without).
- 2026-10-12: England & Wales men 30-100, origins 1996-2001, h = 10:
  Lee-Carter CRPS 0.069 with 95% coverage 0.75; random walk 0.059 and 0.97.
  Lee-Carter fitted from 1981 with observed jump-off: CRPS 0.054, coverage
  0.74. Lee-Carter intervals are too narrow on these data (known limitation:
  only period-index uncertainty).
- 2026-10-12: The 0.1.0 regression baseline (`dev/baseline/`) used GAM
  smoothing and the old data structures; exact reproduction is not a goal.
  The Tucker3 core is tested by exact recovery of a low-rank array instead.

- 2026-10-12: bamp defaults changed to RW2 age, RW1 period and cohort: with
  RW2 period/cohort bamp's own forecasts exploded on French HMD data
  (ages 90-94, 15 years ahead, 90% interval 0.004-0.63 against 0.11-0.17
  with RW1). RW1 forecasts are biased when mortality keeps declining
  (backtest France 1996+10: CRPS 0.067 vs Lee-Carter 0.037, PIT mean 0.23).
- 2026-10-12: Splitting an open exposure interval (85+) into single years
  from its total alone gave errors of 27% at 85-99 and a factor of about 6
  at 100+ (France); `mt_ungroup()` now warns
  (`mt_warning_ungroup_open_counts`).
- 2026-10-12: lintr's `object_name_linter` is disabled (matrix notation such
  as `D`, `E`, `XtWX` is used throughout; lintr 3.1.1 failed with custom
  regexes); line length 120.

- 2026-10-12: `fdm()` renamed `hyndman_ullah()` (clash with `demography::fdm`).
- 2026-10-12: `cause_shares()` with noisy cause-specific rates over-covered
  (95% coverage 0.99) unless the rates are smoothed first
  (`step_smooth()`); with smoothing 0.94-0.97. Residual cause added only when
  the listed causes leave more than 0.1% unexplained (`residual = "auto"`).
- 2026-10-12: France 1986-1995 origins, h = 10 (CRPS, 95% coverage): Li-Lee
  0.058/0.95, product-ratio 0.059/0.91, Lee-Carter 0.061/0.95, StMoMo APC
  0.067/0.95, RH 0.067/0.97 (2 of 4 origins failed), FDM 0.069/0.73, CBD
  0.194/0.67 (CBD is designed for older ages; ages 30+ used here).

- 2026-10-12: First recorded full validation run (`validation/history.csv`,
  commit 803e6d8 code): all 32 threshold checks pass, 418 unit-test
  expectations. Open against the 0.2.0 criteria: the equal-weight ensemble of
  Lee-Carter, Tucker3 and the random walk over-covers on the simulation
  benchmark (95% coverage 0.994 against the target 0.90-0.98), because the
  random walk component is too wide on that generator.

- 2026-10-13: Coherence for sex x region x quintile (simulated, 90 strata,
  origins 2005/2008, h = 10; CRPS bottom/aggregate, 95% coverage
  bottom/aggregate): Lee-Carter per stratum 0.0482/0.0286, 0.97/0.98;
  reconciled WLS without copula 0.0466/0.0282 but aggregates 0.85; with
  the empirical copula 0.0477/0.0292, 0.97/0.95; MinT with the copula
  under-covers (0.88); `li_lee(common_by = "sex")` bottom-up 0.0425/0.0269,
  0.95/0.95. Defaults: WLS, copula for per-stratum models, none for joint
  models.

- 2026-10-13: `recalibrated()` in backtests (CRPS, 95% coverage). England &
  Wales men, origins 1996-2001, h = 10: Lee-Carter 0.0687/0.746, calibrated
  0.0657/0.877 (with drift 0.0653/0.877; interval target 0.0661/0.856);
  rwd() 0.0586/0.967, calibrated 0.0565/0.822; ensemble of both
  0.0593/0.945, calibrated 0.0599/0.898 (before per-draw noise of ensemble
  components: 0.0602/0.906). ONS quintiles, 6-year window, origins
  2012-2014, h = 5: Lee-Carter 0.0673/0.829, calibrated 0.0672/0.837,
  interval target 0.0665/0.897; Li-Lee 0.0687/0.836, interval target
  0.0685/0.889, with drift 0.0818/0.758. Calibration helps models that are
  persistently too narrow; it narrows mixtures that are already wide.
  Group-specific factors (`by`) widen contrasts between groups of joint
  models; the spread factor is held constant beyond the longest calibrated
  horizon.
- 2026-10-13: Reconciliation on ONS sex x quintile with totals (18 strata,
  origins 2008-2014, h = 5): bottom-level CRPS unchanged by any method
  (0.0676-0.0683; Li-Lee bottom-up 0.0665); aggregates best without
  reconciliation (England 0.0440 against 0.0449-0.0458) and 95% coverage of
  aggregates lower after reconciliation (England 0.73 against 0.57-0.64).
  The empirical copula raises the England coverage of bottom-up sums from
  0.46 to 0.58 (origins 2011-2014).

- 2026-10-13: Smoothing before forecasting (CRPS, 95% coverage). England &
  Wales men, origins 1996-2001, h = 10: Lee-Carter 0.0687/0.746; smoothed
  over time and calibrated 0.0658/0.873; the 2D GAM surface with
  calibration 0.0646/0.914; ensemble smoothed over time and calibrated
  0.0562/0.917. Smoothing over time alone narrows the intervals of models
  with a period index (its random walk is fitted to a smoother series), so
  smoothing is combined with `recalibrated()`. ONS quintiles (origins
  2008-2014, h = 5): Li-Lee 0.0666/0.847, smoothed over time 0.0669/0.774;
  life expectancy at birth 0.336/0.81 against 0.348/0.61 years. The ONS
  rates vary around the smoothed trend no more than Poisson noise
  (dispersion 0.4-1.0), so smoothing has little to remove there.
- 2026-10-13: Jump-off. England & Wales men: `anchored(lee_carter(), years
  = 3, decay = 10)` 0.0648 (first two years 0.033 against 0.040); observed
  jump-off 0.0583. ONS: anchored Li-Lee 0.0658 against 0.0666 (e0 0.313
  against 0.336, Q1 0.495 against 0.584); the observed jump-off of one year
  0.0709. Anchoring to three years is the safer default where the last year
  is noisy.
- 2026-10-13: Fixed nodes. With a total far below the sum of its parts
  (causes of ONS data reconciled to separately forecast all-cause rates,
  gap up to 84%), setting negative bottoms to zero moved the fixed totals by
  up to 51%. The projection now holds such bottoms at zero and projects the
  rest again; fixed totals are kept to 1.5e-7 (40 s for 130 strata x 135
  cells x 200 draws).

- 2026-10-13: Decomposition on ONS data. Gap Q5 - Q1 in life expectancy at
  birth, 2019: women 6.46 years (ages 40-84: 5.2), men 7.66 (6.3). Symmetric
  Arriaga and symmetric stepwise replacement are identical for life
  expectancy (difference 1.6e-14); the line integral (20 steps) differs by
  up to 0.007 years per age with a residual of 2e-4. England 2001-2019,
  quintile composition at the radix: women +2.72 years (composition
  -0.10), men +3.74 (-0.13); pooled rates: +2.98 (+0.10) and +3.95 (+0.03).
  The signs differ because the share of Q5 fell at birth and rose at ages
  75-84.

- 2026-10-17: Reviews 3 and 4. Corrected `anchored()` shrinkage (between-year
  variance minus the mean sampling variance): ONS anchored Li-Lee 0.0654
  against 0.0666 (e0 0.315 against 0.336, Q1 0.506 against 0.584); England
  & Wales men h 1-2 0.032 (before 0.040), overall 0.0639. Decomposition of
  forecasts, women Q5 - Q1 at birth: observed 2019 6.46 years, fitted 5.69,
  2035 4.58; the model narrows the gap by 1.11 years (0.03-2.04), 0.78 at
  ages 60-84; of the forecast gain of women in Q1 2019-2035 (3.5 years),
  the jump-off from observed to fitted rates is 0.4. Causes summed
  bottom-up give female Q1 e50 30.8 in 2034 against 31.7 observed in 2019;
  all-cause Lee-Carter 32.9 and Li-Lee 34.3 (cause forecasts more
  pessimistic, as in Wilmoth 1995).

- 2026-10-17: Function map. On the website the diagrams were scaled down
  to unreadable text (map 2 1812 px wide, README 1905 px, column 776 px)
  and their last lines were clipped (Bootstrap `p` margin inside mermaid
  HTML labels). Diagrams redrawn top-down within 776 px (map 2 853 px),
  CSS reset in `_pkgdown.yml`; map 1 no longer sends `mt_smooth()` output
  to models, `mt_prep()` points to `mt_fit()`. Checked by rendering the
  built pages in Chromium with mermaid 10.9.1.

- 2026-10-17: Tuning benchmark (`dev/tuning-benchmark.R`, `dev/tuning.md`).
  Lee-Carter with its fitting window tuned by a 10-year hold-out: E&W CRPS
  0.0573 against 0.0687, e30 0.29 against 0.47 years (window 15 at every
  origin). fdm tuned over components and knots by the same hold-out: no
  gain (E&W 0.0657 against 0.0641; France 0.104 against 0.096), choices
  alternate between origins. Smoothing fdm narrows its intervals (95%
  coverage 0.63 against 0.70). Tuned and calibrated fdm: 95% coverage
  0.90/0.91 (rates), 0.87/0.80 (life expectancy).

- 2026-10-18: Extensions (E2). Lifespan variation depends on the open age
  interval more than the first draft of the help page said: with 85+
  instead of 100+, e-dagger of French women in 2006 is 2.1 years (22%)
  higher, the Gini coefficient and the interquartile range up to 13%
  lower; from 95+ the differences are below 2%. Help page, FAQ and
  vignette state this. Excess deaths of a national total from separate
  models of its parts have too narrow intervals (Germany 2020: 11,000 to
  28,000 from 32 Land-by-sex Lee-Carter models, -1,000 to 55,000 from a
  model of the total), because independent trend errors cancel in the sum;
  the vignette recommends a model of the total or a joint model. With
  three backtest origins the Diebold-Mariano variance is not estimable;
  the model-choice vignette tests on 16 origins instead.

- 2026-10-18: Audit (E3). Removing the extra-Poisson layer narrows the
  intervals of models of smooth rates (Lee-Carter on England & Wales men,
  95% coverage 0.632 against 0.746); on simulated Lee-Carter data with
  Poisson deaths Lee-Carter covers 0.944, ARIMA(0,1,1) with drift 0.823
  and centred `tucker3()` 0.748. The Diebold-Mariano test with the HLN
  correction rejected a true null in 10-19% of simulated samples with 7
  overlapping origins (horizons 2-5), the EWC test in 4.8-7.4%. A cause hierarchy
  declared without `additive = "rates"` gave all-cause rates three times
  too low in `reconciled()`; it now warns.

- 2026-10-18: E4. The full papers changed several implementations: the
  CoDa model of Bergeron-Boucher et al. (2017) uses ARIMA(0,1,1), whose MA
  coefficient is near -1 on England & Wales men, so its 95% coverage fell
  to 0.888; the former non-negativity heuristic missed the optimum of the
  quadratic program in up to 20% of cells (MinT with a fixed total); the
  least-squares Kannisto closure was off by up to 82% from the published
  parameters. Among the new models, `p_splines()` and `bohk_rau()` were
  much more accurate than Lee-Carter on England & Wales men (CRPS 0.035
  and 0.045 against 0.071); intervals of `p_splines()`, `cp_splines()` and
  `nordpred_apc()` are narrow, as their published uncertainty covers part
  of the error only.

## Open questions

See `dev/implementation-plan.md` Section 17.
