# Integrating the methods of the earlier UK analysis into mortensor 0.2

Review of the owner's earlier AI-assisted analysis (`uk_analysis.R`, 5511
lines; `pipeline_methodology.Rmd`; 72 output files: plots and CSV tables),
received 2026-10-13. The analysis forecasts cause-specific mortality in
England by age group, sex and IMD quintile with Tucker3 models, an ensemble
of seven model families, WLS reconciliation of causes to all causes, and
age disaggregation to single years. This plan concerns its **methods**; no
data or results of the analysis enter the repository. The script was read,
not executed.

## 1. Summary

Most of the analysis has already been ported to mortensor 0.2, in a form
that is statistically sounder (Section 3): models, Tucker3 (log and CoDa),
Tucker4-ILR cause shares, age disaggregation, ensembles, tuning and the sex
ordering. What the 0.2 package lacks falls into four groups:

1. **the plot layout** the owner prefers (Section 4.1; implemented first);
2. **jump-off control** at the first forecast year for every model, with a
   diagnostic of the step (Section 4.2);
3. **reconciliation with a fixed total** (the all-cause forecast as a hard
   constraint, the default of the earlier analysis) and a coherence
   diagnostic plot (Section 4.3);
4. **diagnostics and conveniences**: smoothing checks and observed-versus-
   smoothed plots, a 2D GAM smoother over age and time, disaggregation
   diagnostics, trend-level accuracy, a model-comparison bar chart, Lee-Carter
   variants of demography, top-k ensembles and two-stage tuning (Sections
   4.4 to 4.9).

Several elements of the analysis should not be ported because they make
intervals or evaluations invalid (Section 5).

## 2. How the earlier analysis works

| Step | Method | Lines |
|---|---|---|
| Data | ONS deaths by IMD decile to quintiles; causes as sums of columns; residual cause = all causes minus named causes, clipped at 0; zeros and missing cells set to 1e-12 | 210-438 |
| Smoothing | Gaussian 2D P-spline GAM `te(age, year, bs = "ps")` per sex x quintile, REML, COVID years weight 0; fallback 1D Whittaker over time with LOO lambda and a trend-ratio floor (0.15); checks MAD <= 0.2, trend ratio, correlation | 445-700 |
| Models | Tucker3-CoDa and Tucker3-direct (HOOI, B columns forecast by ETS/ARIMA/RWD), Lee-Carter (demography, `adjust`, `chooseperiod`), FDM, StMoMo APC/CBD, RWD, ARIMA and ETS per series; all fitted to the smoothed surface | 866-1106 |
| Calibration | rolling or expanding origins on the smoothed series, grids per model (Tucker3: 315 rank and method combinations), MASE | 1108-1379 |
| Selection | holdout 2016-2019: combined MASE = 0.7 cell + 0.3 age-averaged trend; ensembles by inverse MASE (top k = 2-5) or NNLS stacking; ensemble used if within 5% of the best model | 1384-1650, 4828 |
| Intervals | default `gam_bayes`: posterior draws of the smoother, models refitted on each; alternative residual bootstraps; ensemble draws averaged across models; re-centred on the point forecast and capped at +-1.5 on the log scale | 2164-2472, 5010-5028 |
| Coherence | WLS distribution of the gap between the fixed all-cause forecast and the sum of 10 cause forecasts, weights from interval widths, cap of +-50% per cause, rest absorbed by the residual cause; optional Tucker4-ILR shares; optional sex constraint | 2100-2157, 3466-3919 |
| Disaggregation | GAM over age midpoints and year, parametric tail (Kannisto or Gompertz) blended in from age 95 with a logistic weight (tau = 2), iterative spline correction to group means; intervals by interpolated ratios | 1652-2093 |
| Output | CSV tables and 15 plot types per cause | 3925-4434 |

The jump-off is not handled explicitly: every model starts from its fitted
surface on smoothed data, and the plots join the last fitted value to the
first forecast value, which hides any step from the last observation.

## 3. Already in mortensor 0.2

| Earlier analysis | mortensor 0.2 | Difference |
|---|---|---|
| Tucker3-CoDa, Tucker3-direct | `tucker3(transform = "coda" / "log")` | open interval handled with an auxiliary series instead of a Gompertz step from 80-89 |
| Tucker4-ILR cause shares | `cause_shares()` | draw by draw, with total-cause uncertainty; CV without the oracle all-cause level |
| LC, FDM, StMoMo, RWD, ARIMA, ETS | `lee_carter()`, `hyndman_ullah()`, `stmomo()`, `rwd()`, `arima_series()`, `ets_series()` | simulated paths with period-index uncertainty instead of residual bootstraps of constant width |
| Ensembles (inverse MASE, NNLS) | `ensemble(weights = "crps" / "stack")` | mixture of draws; weights from an inner backtest on training years |
| Grid calibration | `mt_tune()`, `tuned()` | inner origins inside each training window (no leakage) |
| Holdout MASE | `mt_backtest()`, `mt_score()` | CRPS, WIS, coverage, PIT and MASE at the resolution of the observed deaths |
| WLS reconciliation | `mt_reconcile()`, `reconciled()` | draw by draw on expected deaths; methods bottom-up, OLS, structural, WLS, MinT |
| Sex constraint | `mt_order_strata()` | path by path |
| Hyman-Kannisto, natural-Gompertz, GAM disaggregation | `mt_ungroup(method = )`, `mt_disaggregate()` | exposure-weighted group means; forecasts disaggregated draw by draw |
| COVID years with weight 0 | `step_shocks(action = "exclude")` | inside the preprocessing pipeline |
| Residual cause | `mt_group_levels(other = "rest")` with an additive hierarchy | without clipping at zero |
| Interval recalibration | `calibrated()` | minimum-CRPS spread on the last 10 years |

## 4. Work packages

Each package follows the conventions of `CLAUDE.md`: tests, a validation
check with thresholds for statistical code, documentation and pkgdown
updates in the same change.

### 4.1 Plot layout as the package standard (P1, first)

The owner prefers the layout of the earlier plots
(`plot_fc_trend_byage_*`, `plot_fc_age_*`, `plot_obs_smooth_*`):

* `theme_minimal()`; base size 12 for single panels and 9 for grids of
  panels (strip text 8); legend at the bottom for grids, at the right
  otherwise;
* observed rates as small transparent points (size 0.5-0.8, alpha 0.35-0.4);
* fitted rates as solid lines (linewidth 0.5-0.7), with a faint ribbon
  (alpha 0.08) where an interval of the fit exists;
* forecast medians as dashed lines (linewidth 0.6-0.8) with a light ribbon
  (alpha 0.12-0.15) for the 95% interval;
* a dotted grey vertical line half a year before the first forecast year;
* `facet_grid(age ~ sex, scales = "free_y")` for trends by age;
* a title naming the quantity and a subtitle that explains the encoding
  ("Dots = observed, solid = fitted, dashed = forecast (95% interval)");
* for ordered strata such as deprivation quintiles, a diverging
  red-yellow-green palette (ColorBrewer RdYlGn).

Implementation:

* `theme_mortensor(base_size = 11)` and
  `scale_colour_mortensor()`/`scale_fill_mortensor()` for ordered levels,
  exported;
* `mt_plot_forecast()`, `mt_plot_age()`, `mt_plot_life_expectancy()`,
  `mt_plot_rates()`, `mt_plot_shares()`, `mt_plot_scores()`,
  `mt_plot_coverage()`, `mt_plot_tucker()` and `mt_plot_ungroup()` use the
  theme; forecast plots get the observed, fitted, forecast encoding
  (fitted rates are stored in the draws for this);
* y axis of rates: logarithmic, with the rates as labels (as now) or
  natural-log values as in the earlier plots (`log_values = TRUE`);
* vdiffr snapshots regenerated; precomputed vignettes re-knitted.

Corrections to the earlier layout:

* the middle colour of RdYlGn (`#FFFFBF`) is nearly invisible on white; it
  is replaced by a darker yellow;
* the earlier plots coloured Q1 (most deprived) green in the trend plots
  and red in the disaggregation plots; the package uses one direction
  throughout (first level red by default, `reverse = TRUE` for the
  opposite);
* the "Mean log(mx)" axis of the age-averaged trend showed the log of mean
  rates and averaged interval bounds over ages; aggregate trends are
  computed from the draws instead (`asr()`, `mt_aggregate_draws()`).

Decisions for the owner: the default direction of the palette and the
default y axis (rate labels or natural-log values).

### 4.2 Jump-off at the first forecast year (P2)

Problem: models that start from fitted rates produce a step between the
last observation and the first forecast year; the earlier analysis hid it
in the plots. In mortensor, `lee_carter()`, `li_lee()` and
`product_ratio()` offer `jump_off = c("fitted", "observed")`; the other
models start from fitted or smoothed values. On the ONS data, Li-Lee with
an observed jump-off forecast worse (CRPS 0.071 against 0.067), because the
last observed rates carry Poisson noise and year effects.

Plan:

* a wrapper `anchored(spec, to = c("observed", "smoothed"), years = 3,
  decay = 10)` for any model, including ensembles: per age and stratum,
  the ratio between the observed (mean of the last `years` years, or the
  smoothed) rates and the model's fitted rates in the last year scales the
  draws; the log ratio decays linearly to zero over `decay` years (`decay =
  Inf` keeps it, `decay = 0` gives the fitted jump-off). The ratio is
  shrunk towards 1 in cells with few deaths;
* `mt_jump_off(draws)`: a table of the step between the last observed,
  last fitted and first forecast year per stratum and age, on the log
  scale and for life expectancy, with a plot;
* tests: an exact jump-off with `decay = Inf`, the fitted forecast with
  `decay = 0`, coherent draws after `reconciled()`;
* validation: backtests of `fitted`, `observed` and blended jump-off for
  Lee-Carter, Li-Lee, Tucker3 and the ensemble on ONS, England & Wales and
  German Laender; thresholds on the CRPS of horizons 1 and 2;
* plots: the fitted line continues into the forecast, so a step stays
  visible.

### 4.3 Reconciliation with a fixed total and coherence diagnostics (P3)

The earlier analysis treats the all-cause forecast as a hard constraint
and distributes the gap to the causes by their forecast variance. In
`mt_reconcile()` every node moves, including the total.

Plan:

* `mt_reconcile(..., fixed = c("all"))` and `reconciled(..., fixed = )`:
  nodes whose forecasts are kept, by label or by level (`"total"`); the
  projection uses zero variance for these nodes (constrained least squares
  with the fixed nodes as additional constraints), draw by draw on the
  deaths scale; works with WLS, structural, OLS and MinT;
* non-negativity as now (bottom values set to zero and aggregates
  recomputed); a fixed total that cannot be met with non-negative parts
  is reported;
* `mt_plot_coherence(draws, hierarchy)`: ratio of the sum of the parts to
  the total, per year (box plots over cells and draws), as in
  `plot_consistency_coherent`;
* validation: ONS causes at ages 50+ and German Laender; coherence gap
  below 1e-8; CRPS and coverage of the fixed total equal to its own model;
  comparison with bottom-up and free WLS in the coherence vignette.

Not ported: the cap of +-50% per cause with absorption by the residual
cause, and shifting interval bounds by the point adjustment (Section 5).

### 4.4 Smoothing diagnostics and a 2D GAM smoother (P4)

* `mt_smooth(method = "gam", dims = 2)`: Gaussian or Poisson
  `mgcv::gam(log rate ~ te(age, year, bs = "ps"))` per stratum, REML, with
  weights from the deaths (the earlier analysis weighted all cells
  equally, which gives cells with few deaths too much influence);
* stored diagnostics: median absolute deviation between raw and smoothed
  log rates, ratio of the variance of year-to-year changes, correlation;
  reported by `summary()` of the smoothing step;
* `mt_plot_smooth(raw, smoothed)`: observed points, smoothed lines and the
  95% interval of the smooth by stratum, as in `plot_obs_smooth_*`;
* a Whittaker smoother over time (the package smooths over age).

### 4.5 Disaggregation options and diagnostics (P5)

* Gompertz tail as an alternative to Kannisto in the GAM method of
  `mt_ungroup()`, and the width of the blend (`tau`) as an argument;
* time-varying age shape for forecasts in `mt_disaggregate(shape = )`:
  the age pattern within groups taken from a GAM extrapolated in time;
* `mt_plot_disaggregation(type = )`: stages of the profile (GAM, tail,
  blend, corrected), exposure-weighted correction factors per group,
  heat map of second differences of log rates, single-age trends across
  the boundary between observed and forecast years.

### 4.6 Accuracy at the level of trends (P6)

* `mt_score(..., level = "trend")`: scores of age-standardised rates or
  life expectancy per stratum and year, computed from the draws; this
  replaces the earlier trend MASE on age-averaged rates (dominated by the
  oldest ages);
* `mt_plot_scores(type = "bar")`: the model comparison as horizontal bars
  with a reference line, as in `plot_model_mase_*`.

### 4.7 Lee-Carter variants of demography (P7)

* `lee_carter(adjust = c("none", "dt", "dxt", "e0"), choose_period =
  FALSE)`: re-estimation of k(t) to match total deaths (Lee and Carter
  1992), deaths by age (Booth, Maindonald and Smith 2002) or life
  expectancy (Lee and Miller 2001), and the choice of the fitting period by
  Booth, Maindonald and Smith (2002), as in `demography::lca()`; validated
  against `demography::lca()` like the existing wrappers.

### 4.8 Ensembles and tuning (P8)

* `ensemble(top_k = )`: the k components with the best inner-backtest
  scores, chosen on training years only;
* `mt_tune(strategy = "two_stage")`: a coarse grid followed by a fine grid
  around the best combination, useful for Tucker ranks.

### 4.9 Uncertainty of smoothing (P9, optional)

* `bootstrapped(type = "smooth")`: smoothed surfaces drawn from the
  posterior of the 2D GAM, the model refitted on each, and the forecast
  draws pooled. This adds the uncertainty of the smoothing to the forecast
  draws; it does not replace them, because draws from the smoother's
  posterior alone contain no future variation.

## 5. Not to be ported

| Element | Reason |
|---|---|
| `gam_bayes` intervals as prediction intervals | posterior draws of the smoother with refitted point forecasts contain no future innovations: they are confidence bands for the trend |
| Ensemble draws averaged across models | averaging independent draws shrinks the spread; mixtures keep it |
| Re-centring of intervals on the point forecast and the +-1.5 cap on the log scale | moves and truncates the predictive distribution; draws make the median the summary |
| Ensemble size, weights and stacking fitted on the holdout, 5% tolerance | the same holdout selects and evaluates; weights must come from training years |
| Calibration on the two-sided smoothed full series | smoothing over the test years leaks them into calibration |
| Cap of +-50% per cause in reconciliation, absorption by the residual cause, shifting interval bounds by the point adjustment | ad hoc; quantiles of a sum are not sums of quantiles; reconciliation works on draws |
| Cause rates capped at 0.95 times the all-cause forecast | redundant with coherent forecasts |
| Sex constraint as an automatic adjustment | wrong for some causes (breast cancer, dementia at old ages); kept as an option (`mt_order_strata()`), with a diagnostic |
| Missing cells set to 1e-12, failed forecasts set to 0 before weighting | indistinguishable from zero deaths; biases rates downward |
| Group-mean correction with unweighted means | biases single-year rates in the open interval downward; the package corrects with exposure weights |
| Interval disaggregation by interpolated ratios | assumes perfect correlation within groups; draws are disaggregated instead |
| Tucker4 cross-validation with the observed all-cause level | an oracle total makes the errors optimistic |
| Deleting plot files in a fixed directory | the package writes files only where the user asks |

## 6. Order and effort

| Order | Package | Effort | Depends on |
|---|---|---|---|
| 1 | 4.1 plot layout | 1 day; snapshots and vignettes re-knitted | none |
| 2 | 4.2 jump-off | 2 days with validation | fitted rates stored in draws (4.1) |
| 3 | 4.3 fixed total | 2 days with validation | none |
| 4 | 4.4 smoothing diagnostics | 1-2 days | none |
| 5 | 4.5 disaggregation diagnostics | 1-2 days | 4.1 |
| 6 | 4.6 trend scores | 1 day | none |
| 7 | 4.7 Lee-Carter variants | 1 day | none |
| 8 | 4.8 top-k, two-stage tuning | 1 day | none |
| 9 | 4.9 smoothing uncertainty | 1 day | 4.4 |

## References

Booth H, Maindonald J, Smith L (2002). Applying Lee-Carter under conditions
of variable mortality decline. *Population Studies* 56(3): 325-336.

Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
*Journal of the American Statistical Association* 87(419): 659-671.

Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter method
for forecasting mortality. *Demography* 38(4): 537-549.
doi:10.1353/dem.2001.0036

Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast
reconciliation for hierarchical and grouped time series through trace
minimization. *Journal of the American Statistical Association* 114(526):
804-819.
