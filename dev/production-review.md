# mortensor production readiness review

Review date: 2026-10-12. Package version reviewed: 0.1.0 (commit `09b80b6` on
`main`, plus the handoff notes in `dev/review-handoff.md`). This document
continues the first review pass recorded in the handoff notes and replaces its
"Still to do" list with results.

## 1. How this review was done

- Environment: Ubuntu 24.04, R 4.3.3 from apt, all Imports and Suggests
  installed (`rTensor` 1.5.0, `compositions` 2.0.9, `demography` 2.0.1,
  `StMoMo` 0.4.1 from CRAN; `forecast` 9.x from apt). No LaTeX, so the PDF
  manual check is not informative.
- The package was installed unmodified from the repository.
- Because `%||%` is missing on R < 4.4 (defect D1), every run after
  `R CMD check` defined `%||%` in the global environment (identical to the
  base R 4.4 definition, `if (is.null(x)) y else x`). This reproduces what a
  user on R >= 4.4 sees and exposes the defects hidden behind D1.
- Experiments (scripts kept outside the repository; numbers below are from
  single runs with fixed seeds):
  1. `R CMD build` and `R CMD check --as-cran`; `testthat`.
  2. Every chunk of every vignette, evaluated in order with `eval = TRUE`
     (the vignettes set `eval = FALSE` globally, see D3).
  3. The README quick start, timed step by step, plus all `plot()` types,
     `export()`, `accuracy()`, `disaggregate()` and `life_expectancy()`.
  4. Targeted checks of each code-reading finding from the handoff.
  5. A holdout leakage experiment and a prediction interval (PI) coverage
     backtest (example data and a simulated Lee-Carter-type process).

## 2. Summary

[[SUMMARY]]

## 3. Current state

### 3.1 R CMD check --as-cran (R 4.3.3, unmodified package)

Status: 3 ERRORs, 3 WARNINGs, 4 NOTEs. Build took 18 s, check 116 s.

| Result | Check | Cause | Package problem? |
|---|---|---|---|
| ERROR | examples | `print(example_mortality)` fails: `could not find function "%||%"` | Yes (D1) |
| ERROR | tests | 5 of 64 expectations fail, all from `%||%` | Yes (D1) |
| ERROR / WARNING | PDF manual | `pdflatex` not installed in the review container | No (environment) |
| WARNING | DESCRIPTION meta-information | `License: MIT LICENSE` is not a valid specification | Yes (D17) |
| WARNING | R files for syntax errors | locale `en_US.UTF-8` unavailable in container | No (environment) |
| NOTE | R code for possible problems | `%||%` undefined in 23 functions; undefined globals `group_label`, `line_id` | Yes (D1, D17) |
| NOTE | future file timestamps, HTML manual, `mortensor-manual.tex` | container limitations | No |

### 3.2 Tests

`testthat` runs 24 `test_that()` blocks with 64 expectations. On R 4.3.3, 5
fail (all D1). With `%||%` defined, all 64 pass. The suite never runs
`calibrate()` on more than toy settings, never checks numerical results of the
main workflow, and contains a malformed fixture (D15), so passing tests say
little about correctness.

### 3.3 Vignettes

[[VIGNETTES]]

### 3.4 README quick start and runtime

[[README]]

## 4. Verified defects

Severity: **Blocker** = install, check or the documented workflow fails;
**High** = default workflow silently returns wrong or misleading results;
**Medium** = wrong results in non-default paths or limited impact;
**Low** = hygiene and maintainability.

| ID | Severity | Area | Defect | Evidence (this review) |
|---|---|---|---|---|
| D1 | Blocker | Dependencies | `%||%` is used in 23 functions but never defined or imported; it exists in base R only from 4.4.0 while DESCRIPTION declares R >= 4.1.0 | R CMD check: examples ERROR, 5 test failures; `print()` of every main class fails on R 4.3.3 |
| D2 | Blocker | Generics | `forecast()` (`R/workflow.R:28`) and `accuracy()` (`R/metrics.R:40`) are new generics with the same names as the `forecast`/`generics` generics. After `library(forecast)` (forecast 9.x re-exports `generics::forecast`), `forecast(<mort_calibration>)` dispatches to the wrong generic | Reproduced: `non-numeric argument to binary operator` |
| D3 | Blocker | Documentation | All nine vignettes set `eval = FALSE` globally, so `R CMD build` succeeds without running any code. Evaluated, [[VIG_FAILCOUNT]] | Section 3.3 |
| D4 | High | Model selection | Holdout leakage: the holdout years are smoothed together with the training years, the last rolling-origin CV fold of every model evaluates exactly the holdout years, and the README path scores forecasts against smoothed (not observed) holdout values | Section 5.1, [[LEAK_SHORT]] |
| D5 | High | Ensemble | Ensemble size, inverse-MASE weights and NNLS stacking weights are chosen on the holdout and the same holdout score is reported and compared with single models (with a 5% preference margin for the ensemble). NNLS weights are fitted on log(mx) and applied to mx | Code: `R/comparison.R:151-237`, `R/workflow.R:250-262` |
| D6 | High | Ensemble | A component that fails (NA) contributes 0 to the weighted ensemble average without renormalising the weights per cell (`R/workflow.R:320`, `R/comparison.R:170`, `R/bootstrap.R:238,323`), and no warning is raised | Lee-Carter failing for one group: ensemble forecast for that group = 0.500 x Tucker3 forecast; all other groups 0.99-1.00 |
| D7 | High | PIs | `add_pi()` on an ensemble simulates only the first component (`extras$method` of the first model) | Ensemble of RWD and Tucker3: PI width equals RWD alone (0.180) or Tucker3 alone (0.184) depending only on component order |
| D8 | High | PIs | PIs do not widen with horizon for RWD/ARIMA/ETS/Lee-Carter/FDM (point forecast plus i.i.d. resampled in-sample residuals; Lee-Carter k_t uncertainty ignored), are re-centred on the point forecast and capped at +/-1.5 log units; `mx_median` is overwritten with the point forecast | 95% PIs covered 84-87% on the example data; on simulated random-walk data coverage fell from 80-85% at h = 1 to 37-68% at h = 10, with constant width for five of seven models; README run: width 0.021 log units at all horizons (Section 5.2) |
| D9 | Medium | PIs | `ensemble_bootstrap()` (used by `run_forecast()`) averages draws across models instead of sampling from the mixture, which shrinks the predictive spread | Two components with 10% spread: 95% log-width 0.27 (averaged draws) vs 0.48 (mixture) |
| D10 | High | Life tables | `life_expectancy()` and `life_table()` default to `width = 5`, also for single-year output of `disaggregate()` | Single-year Gompertz schedule: e30 = 201.6 with the default vs 59.0 with `width = 1`; README forecast disaggregated to single years: mean e30 = 207.3 with the default vs 62.7 with `width = 1` (63.4 from the grouped forecast) |
| D11 | High | CoDa | The open age group is lost in the CoDa transform: `compute_dx()` ignores the open-interval rate and `dx_to_mx()` replaces it with a Gompertz extrapolation of the two previous groups | Doubling the input 90+ rate changes the round-trip 90+ rate by a factor of 1.000; on `example_mortality` the round trip inflates m(90+) by 20% |
| D12 | Medium | Tucker3 | Data with one group fail: rank grid `2:min(4, n_g - 1)` becomes `2:0`; BIC selection returns `R_C = 2` with criterion `Inf` without a warning, and forecasting then errors (`ranks must be positive`) | Reproduced for Tucker3-Direct and Tucker3-CoDa |
| D13 | Medium | Calibration grids | ETS candidates `AAdN`, `MMN` and `ANN` with `damped = TRUE` always error and silently fall back to `auto.arima`; RWD drift windows `0` and `length(years)` give identical forecasts | 3 of 5 ETS specifications error; max abs difference between the two RWD windows = 0 |
| D14 | Medium | Pipeline | `run_forecast()` defaults (`holdout_train = holdout_test = NULL`) give a validation horizon of 0, so every CV score is `Inf` and all hyperparameters silently take their defaults | `calibrate_tucker3(h_val = 0)` returns ranks 3/3/2, `ets`, score `Inf` |
| D15 | Medium | Data input | `build_array()` sorts age labels as text, so `"10-14"` is placed before `"5-9"` and values are assigned to the wrong ages; a missing age row fills the whole age vector with 1e-12 (with only a message); the unit test fixture has duplicated age rows and only checks dimensions | Input rates for ages 0/5/10 = 0.001/0.0001/0.0002 come back as 0.001/0.0002/0.0001 |
| D16 | Medium | Scoring | `compute_wis()` divides by the sum of the weights instead of K + 1/2 as in Bracher et al. (2021), so values are not comparable with the literature or other software, and the scale changes with the quantile set | 7 quantiles: 4.00 x the published definition; 5 quantiles: 3.23 x |
| D17 | Low | Package metadata | `License: MIT LICENSE` is invalid (should be `MIT + file LICENSE` with the two-line CRAN LICENSE file); undefined globals `group_label`, `line_id`; copyright holder in LICENSE (`kalleEF`) differs from DESCRIPTION | R CMD check WARNING and NOTE |
| D18 | Medium | API | The `models` argument of `calibrate()` does not restrict the comparison: `compare_models()` always fits all models and `calibrate_all_models()` always tunes all non-Tucker3 models; user-facing names (`"lee_carter"`) differ from internal names (`"Lee-Carter"`) | [[MODELS_SHORT]] |
| D19 | Medium | Model selection | Models with NA cells are flagged `complete = FALSE` and excluded from the ensemble, but remain eligible as best single model; their MASE is computed on the subset of cells that did not fail | Lee-Carter and StMoMo with a failed group: finite combined MASE, `complete = FALSE` |
| D20 | Low | Life tables | `life_table()` documentation says Kannisto is the default open-age method; the code uses the reciprocal. Kannisto/Gompertz extrapolation uses group start ages and a sum that starts at exp(-mu_0) | e90 for a Gompertz hazard: exact 7.43, package Gompertz option 7.15; e30 differs by 2.4 years between reciprocal and Kannisto on `example_mortality` |
| D21 | Medium | Disaggregation | Single-year rates are constrained so that their unweighted mean equals the grouped rate; for the open interval the mean is taken over 90-100 only. The correct constraint is exposure weighted over the full interval | Stationary Gompertz population: single-year rates at ages 95 and 100 overestimated by 11% and 22%; e30 0.4 years too low |
| D22 | Medium | Life tables | e_x bounds are life tables of the pointwise upper and lower mx bounds (`R/life_table.R:320-327`), which is not a prediction interval for e_x | Code reading; requires simulated paths to fix |
| D23 | Medium | S3 workflow | `disaggregate(method = "gam")` stops with an error; `plot(<disaggregated forecast>, type = "disagg")` returns NULL with "No single-year disaggregated data available" because it reads `extras$mx_1yr_dt`, which only `run_forecast()` sets; single-year `export()` writes correct rows to a file named `*_forecast_grouped.csv`; `print()`/`summary()` report the holdout MASE of the selected model as "MASE"/"CV MASE" | README run, Section 3.4 |
| D24 | Medium | Multi-cause | `run_reconciled_forecast()` and `run_coherent_forecast()` call `run_forecast(<mort_config>)`, which only works with the ONS Excel input; `enforce_coherence()` hardcodes `sex == "Male"/"Female"` and `quintile` | [[MULTI_SHORT]] |
| D25 | Medium | Performance | Default `calibrate()` fits 315 Tucker3 configurations x 5 origins x 2 variants plus all other grids | [[RUNTIME_SHORT]] |
| D26 | Low | Hygiene | 170 `cat()` calls (not suppressible), 12 plot methods call `print(p)`, 17 `ggsave()` side effects in exported helpers, 40 `tryCatch(..., error = function(e) NULL)` that hide failures, 96 hard-coded `1e-12` floors; committed `.claude/settings.local.json` (contains a local Windows path) and `tests/testthat/Rplots.pdf`; no R CMD check CI workflow; 103 exported functions including internal helpers | `grep` counts on `R/` |

## 5. Methodological assessment

### 5.1 Model selection and evaluation

The holdout comparison in `calibrate()` is the package's central claim of
accuracy, and its design does not separate training from test data:

- `calibrate()` smooths all years, including the holdout, before splitting
  (`R/workflow.R:166-172`). `run_forecast()` does the same.
- Tucker3 rank and time-series selection (`cv_rolling_origin()`) and the
  grid search for all other models (`cv_expanding_model()`) receive the full
  series. With 20 years and the default 4-year holdout, the last CV fold of
  both routines trains on years 1-16 and evaluates years 17-20, which are the
  holdout years. Instrumenting `tucker3_direct_forecast()` during a default
  `calibrate()` call confirmed fits with training data ending at the last
  pre-holdout year (Section 5.1 table).
- When users follow the README and pass smoothed data to `calibrate()`,
  `data_raw` is the smoothed object, so the holdout "observations" are the
  GAM fit through the holdout years.
- Ensemble choice and weights are tuned on the holdout and then compared with
  single models on the same holdout (D5).

[[LEAK_TABLE]]

The MASE used throughout is computed on log rates and scaled by the in-sample
one-step naive error (Hyndman and Koehler 2006). That is reasonable as a
relative measure, but the package averages it with a "trend MASE" using fixed
weights 0.7/0.3 without justification, and reports no horizon-specific
accuracy or skill relative to a simple benchmark such as a random walk with
drift.

### 5.2 Prediction intervals

Two backtests ran each model through the package's own path
(`mort_calibration` -> `forecast()` -> `add_pi(n_sim = 200)`) with default
hyperparameters, and counted how often held-out observed log rates fell
inside the nominal 95% interval.

**(a) `example_mortality`, train 2001-2014, test 2015-2019.** The example data
are a deterministic log-linear trend plus i.i.d. noise (sd 0.04), the easiest
possible case for residual-bootstrap intervals.

| Model | Coverage (nominal 95%) | Log PI width, h = 1 | Log PI width, h = 5 |
|---|---|---|---|
| Tucker3-Direct | 0.85 | 0.134 | 0.134 |
| Tucker3-CoDa | 0.87 | 0.139 | 0.139 |
| Lee-Carter | 0.85 | 0.131 | 0.131 |
| FDM | 0.84 | 0.132 | 0.132 |
| RWD | 0.86 | 0.131 | 0.132 |
| ARIMA | 0.86 | 0.131 | 0.131 |
| ETS | 0.84 | 0.130 | 0.130 |

**(b) Simulated Lee-Carter-type data** (20 replicates; 13 age groups x 10
groups; log m = a_x + b_x k_t + group shift + noise, k_t a random walk with
drift, innovation sd 0.02, observation sd 0.03; train 20 years, forecast 10).
Here the true predictive spread grows with the horizon. Ignoring parameter
uncertainty, the correct 95% log width is about 0.14 at h = 1 and 0.27 at
h = 10.

| Model | Coverage h = 1 | Coverage h = 5 | Coverage h = 10 | Log width h = 1 | Log width h = 10 |
|---|---|---|---|---|---|
| Tucker3-Direct | 0.80 | 0.69 | 0.68 | 0.143 | 0.190 |
| Tucker3-CoDa | 0.81 | 0.70 | 0.67 | 0.154 | 0.194 |
| Lee-Carter | 0.84 | 0.66 | 0.53 | 0.110 | 0.110 |
| FDM | 0.80 | 0.61 | 0.46 | 0.111 | 0.111 |
| RWD | 0.85 | 0.69 | 0.55 | 0.114 | 0.114 |
| ARIMA | 0.82 | 0.60 | 0.48 | 0.108 | 0.108 |
| ETS | 0.83 | 0.56 | 0.37 | 0.109 | 0.109 |
| Correct (oracle, approx.) | 0.95 | 0.95 | 0.95 | 0.14 | 0.27 |

The cap at +/-1.5 log units never bound in these runs, and coverage computed
from the raw simulation quantiles (before re-centring and capping) was within
0.02 of the reported coverage. The under-coverage therefore comes from how
the simulations are generated, not from the post-processing.

**(c) README quick start.** The ensemble chosen by `calibrate()` was FDM,
ARIMA and ETS with weights 0.000, 0.977 and 0.023. `add_pi()` simulated from
FDM, the first listed component, which has weight 0 and the worst holdout
score (D7). Because `mx_raw` was the smoothed series, residuals were tiny, and
the 95% interval was 0.021 log units wide (about +/-1%) at every horizon from
1 to 10 years.

Mechanisms behind these results, all confirmed in the code paths exercised
above:

- RWD, ARIMA and ETS: `exp(point + resampled in-sample residual)`
  (`R/dispatch.R:256-298`). Residual variance does not grow with the horizon.
- Lee-Carter and FDM: the same construction (`R/bootstrap.R:159-162`), so the
  random-walk uncertainty of the period index k_t of the Lee-Carter model
  (Lee and Carter 1992), which is what makes its intervals widen with the
  horizon, is ignored.
- Tucker3: B-mode paths are simulated, but with the default ranks and ETS
  components the simulated spread is small relative to the residual term.
- Residuals are winsorised at +/-2 SD per series and come from fits to
  smoothed data, which narrows the tails further.
- `add_pi()` shifts all quantiles so that the median equals the point
  forecast and caps them at +/-1.5 log units.
- `gam_bayes` propagates only the smoothing uncertainty of the GAM fit, not
  the forecast innovations.

### 5.3 Demographic methods

- Ages and widths: several functions assume 5-year groups (`width = 5`
  defaults, `age_width` from the config instead of the data) and an open
  interval starting at 90 (`R/disaggregation.R:60,79`). Single-year output is
  then summarised with 5-year widths (D10).
- CoDa on d_x: the transformation is a sound idea for coherent age patterns,
  but the open interval is replaced by an extrapolation (D11), so the
  forecast of the oldest group is not learned from data.
- Smoothing is a Gaussian GAM on log(mx) with 0/1 weights. With deaths and
  exposures available, a Poisson P-spline model with a log-exposure offset is
  the standard approach (Currie et al. 2004) and gives correct weighting of
  small cells.
- StMoMo is fitted with an age index 1..n rather than ages
  (`R/models_stmomo.R:22,38`), logit-link models (CBD) are fitted with central
  exposures and return q rather than m, and fitted values always come from APC
  (`R/dispatch.R:92`).
- Disaggregation (D21): the penalized composite link model (Rizzi et al.
  2015) is a well-established alternative that respects the grouping
  constraint on the count scale.

### 5.4 Multi-cause forecasting

The reconciliation and coherence functions are tied to the ONS input path
(D24). WLS reconciliation closes the gap exactly only if a cause named
`"residual"` exists (`R/reconciliation.R:65`); `compute_cause_shares()` adds a
remainder part even when the cause list already contains `"residual"`,
creating a near-zero ILR component. `enforce_coherence()` imposes
female <= male mortality, which is a constraint on the output rather than a
coherent multi-population model in the sense of Li and Lee (2005).

### 5.5 Architecture

Two parallel implementations exist: `run_forecast()` and the S3 workflow
(`calibrate() -> forecast() -> add_pi() -> disaggregate()`). They differ in
ensemble PIs (D7 vs D9), coherence handling, GAM disaggregation, and which
objects carry single-year output, so fixing a defect in one path does not fix
the other. Many functions swallow errors, return NA or defaults, and print
progress with `cat()`, which makes failures invisible in batch use.

## 6. Corrections to the handoff notes

Running the code changed or sharpened several statements in
`dev/review-handoff.md`:

1. **`build_array` test.** The function does not return an all-NA array for
   the test fixture. The fixture itself is malformed (each year/group has two
   `30-34` rows and no `35-39` row), `build_array()` accepts the duplicates
   silently, and the test passes because it only checks dimensions. The real
   defects are the text sort of age labels and the silent floor fill (D15).
2. **Open-interval disaggregation bias.** The handoff expected old-age rates to
   be biased downward. In a stationary Gompertz population the bias is upward
   (+11% at 95, +22% at 100), because the 90+ rate also contains deaths and
   exposure above age 100 while the constraint averages only over 90-100. The
   sign depends on the age pattern and on how far the open interval extends;
   the constraint is wrong either way (D21).
3. **Kannisto off-by-one.** Present, but small (about 0.3 years in e90 in the
   test). The larger problems are the use of group start ages instead of
   midpoints and the mismatch between documentation and default (D20).
4. **Single-group data.** `n_g = 2` works (the rank is clamped internally);
   `n_g = 1` fails (D12).
5. **Generic masking.** `forecast` 9.x re-exports `generics::forecast` and
   `generics::accuracy`; the fix is to import and re-export those generics and
   register methods on them (D2).
6. **CRAN access.** In this container only `https://cran.r-project.org`
   is reachable; `cloud.r-project.org` and `packagemanager.posit.co` are
   blocked.

## 7. Prioritised roadmap

The order below follows dependencies: nothing downstream can be trusted until
the evaluation design is honest, and nothing can be evaluated until the
package runs and is tested. Effort estimates assume one developer familiar
with the code.

### Phase 0: make it install, check and run (1-2 weeks)

1. Import `%||%` from rlang (already in Imports) via
   `@importFrom rlang %||%` (D1).
2. Delete the package's own `forecast()` and `accuracy()` generics; import and
   re-export `generics::forecast` and `generics::accuracy` and register
   `forecast.mort_calibration`, `accuracy.mort_forecast` on them (D2).
3. Fix DESCRIPTION (`License: MIT + file LICENSE`, two-line LICENSE, full text
   in `LICENSE.md`), undefined globals, remove `.claude/settings.local.json`
   and `tests/testthat/Rplots.pdf` from git (D17, D26).
4. Add `r-lib/actions` R CMD check CI on release, devel, oldrel and the
   minimum declared R version, plus test coverage reporting.
5. Make vignettes executable (`eval = TRUE`) with fast settings
   (`rank_selection = "bic"`, few models, `n_sim = 50`); drop or ship data for
   chunks that need ONS files; precompute anything slower than about 1 minute
   (D3).
6. Add an end-to-end regression test on `example_mortality` that runs
   `calibrate() -> forecast() -> add_pi() -> disaggregate() ->
   life_expectancy()` with fast settings and checks numerical invariants
   (positivity, dimensions, e_x in a plausible range, PI ordering).

### Phase 1: honest evaluation and model selection (3-5 weeks)

1. Split before anything else: the holdout years must never enter smoothing,
   rank selection, hyperparameter tuning, ensemble weighting or stacking.
   CV origins must end at the last training year (D4).
2. Keep raw rates (and deaths and exposures) in `mort_data`; never overwrite
   `mx` with smoothed values; always score against observed data (D4).
3. Choose ensemble members and weights by inner CV on training years only;
   report the holdout score once, for the pre-chosen procedure; fit stacking
   weights on the scale on which they are applied (D5).
4. Renormalise ensemble weights per cell over available components and record
   component failures in a diagnostics table; never treat NA as 0; exclude
   incomplete models from "best model" or impute explicitly (D6, D19).
5. One model registry with user-facing names; honour `models` everywhere (D18).
6. Repair grids (valid ETS specifications, distinct RWD windows), screen
   Tucker3 ranks by BIC and cross-validate only the top few, add optional
   parallelism (`future.apply`) (D13, D25).
7. Fix single-group support and fail loudly when a model cannot be fitted (D12).

### Phase 2: probabilistic forecasts (3-4 weeks)

1. Replace "point forecast + i.i.d. in-sample residuals" with model-consistent
   simulation: simulate the time-series components forward (k_t random walk
   with drift for Lee-Carter, `simulate()` for ARIMA/ETS per series, B-mode
   paths for Tucker3), add observation noise, and include parameter
   uncertainty (D8).
2. Ensemble predictive distribution as a mixture of component distributions
   (sample components with probability equal to their weights) (D7, D9).
3. Remove re-centring and capping (or make them opt-in diagnostics); keep the
   simulated paths and compute every derived quantity (e_x, totals, ratios)
   pathwise, so that its interval is a real prediction interval (D8, D22).
4. Build calibration checks into the package: horizon-specific coverage, PIT
   histograms, CRPS and WIS (Bracher et al. 2021 normalisation, or delegate to
   `scoringutils`) on rolling-origin backtests (D16).

### Phase 3: demographic correctness (3-5 weeks)

1. Represent age intervals explicitly (start, width, open interval) and derive
   widths from data everywhere; single-year life tables must use width 1 (D10).
2. Keep the open age group inside the CoDa transform (d_x of the open interval
   from l_x at the open age; inverse via L_x of the open interval) (D11).
3. Life tables: documented default, midpoint-based extrapolation, correct
   numerical integration, and tests against published life tables (D20).
4. Disaggregation with an exposure-weighted constraint over the full open
   interval, or the penalized composite link model (Rizzi et al. 2015,
   `ungroup` package) using deaths and exposures (D21).
5. Smoothing as a Poisson model with a log-exposure offset when counts are
   available (Currie et al. 2004), instead of a Gaussian GAM on log rates.
6. StMoMo: real ages, correct exposure type for logit-link models,
   fitted values of the selected model.
7. Coherence across groups via established multi-population models (Li and
   Lee 2005; product-ratio method of Hyndman et al. 2013) instead of the
   female <= male constraint; group definitions from data, not hardcoded.

### Phase 4: architecture and API (in parallel with Phases 1-3)

1. One pipeline: `run_forecast()` becomes a thin wrapper around the S3
   workflow so both paths share code and behaviour (D23).
2. Model plug-in interface: each model supplies `fit()`, `forecast()`,
   `simulate()` and `fitted()` with declared data requirements (rates only,
   or deaths and exposures). This is also the integration point for external
   backends such as bamp (Section 8), StMoMo or BAPC.
3. Data class holding deaths, exposures, raw and smoothed rates, age
   intervals and arbitrary group variables; tidy long output of point
   forecasts and draws.
4. `cli`/`message()` with a `verbose` option instead of `cat()`; plot methods
   return ggplot objects; files are written only by `export()` (D26).
5. Fail loudly with informative errors; replace silent `tryCatch(NULL)`
   patterns with recorded failures (D26).
6. Reduce the exported API to the user-facing workflow; keep helpers internal.
7. Multi-cause: make reconciliation and coherence work with `mort_data`
   inputs; consider MinT reconciliation (Wickramasuriya et al. 2019) or
   compositional cause-of-death models (D24).

### Phase 5: validation and release

1. Benchmark against established implementations (`demography`, `StMoMo`,
   `BAPC`/`bamp`) on public data (e.g. Human Mortality Database via
   `HMDHFDplus`) with rolling-origin evaluation of point and interval
   accuracy, replicating the design of published comparisons such as Booth et
   al. (2006).
2. Document the validated methods with their limits; version 1.0.0 only after
   Phases 0-3 and the benchmark.

## 8. Integrating bamp

Question from the owner: can bamp (<https://github.com/volkerschmid/bamp>)
and its methods be integrated into a flexible mortality modelling package?

**Short answer: yes, as an optional model backend behind a plug-in interface,
not as copied code, and only after mortensor works on counts and on
predictive draws (Phases 1-2, Phase 4 item 2).** bamp supplies exactly what
mortensor lacks today: a likelihood-based model for deaths and population
with cohort effects and forecast intervals that widen with the horizon and
include parameter uncertainty.

### 8.1 What bamp is

bamp 3.0.0 (CRAN, released 2026-08-31; the GitHub repository is at the same
version) implements the Bayesian age-period-cohort (APC) model of Schmid and
Held (2007), building on Knorr-Held and Rainer (2001):

- Binomial likelihood with logit link: deaths ~ Bin(population, p);
  logit p = mu + age + period + cohort, with RW1 or RW2 smoothing priors,
  optional unstructured heterogeneity per effect and optional cell-level
  overdispersion.
- Unequal grids: `periods_per_agegroup = 5` handles 5-year age groups with
  annual periods through overlapping 1-year cohorts. This matches mortensor's
  typical data.
- Inference: a new Pólya-Gamma Gibbs engine (default, C implementation) and
  the legacy block Metropolis-Hastings engine (`method = "taylor"`); four
  chains, Gelman-Rubin checks (`checkConvergence()`), DIC, and greedy DIC
  model selection (`selectModel()`).
- Prediction (`predict_apc()`): period and cohort effects are extrapolated as
  random walks for every posterior draw, returning posterior draws of rates
  (and of cases when future population is supplied).

### 8.2 Test run on mortensor data

One group of `example_mortality` (male, Q1, 2001-2019; deaths =
round(mx x pop)), RW2 for age, period and cohort, 6,000 iterations x 4 chains
on 4 cores:

| Step | Result |
|---|---|
| `bamp()` | 141 s wall time; `checkConvergence()` TRUE (max R 1.008); DIC 129.7 |
| `predict_apc()`, 5 periods | 11 s; 4,000 posterior draws of rates for 13 ages x 24 periods |
| 2020 forecast vs observed (held out) | 30-34: 0.00040 (90% interval 0.00033-0.00049) vs 0.00047; 90+: 0.0656 (0.0619-0.0692) vs 0.0806. The 2020 shock in the example data is outside the intervals, as expected |
| Automatic MCMC length (40,000 iterations) | about 16 min per group, extrapolated from a calibration run |

### 8.3 Fit with mortensor and gaps to close

| Topic | bamp | Consequence for integration |
|---|---|---|
| Data | Counts: deaths and population matrices, periods x ages, one stratum per call; no missing values (one NA made all chains fail) | mortensor must keep deaths and exposures (Phase 4 item 3) and loop over groups; impute or reject missing cells explicitly |
| Groups | No joint multi-group model | 10 independent fits for sex x quintile; no borrowing of strength between groups (multivariate APC models exist in the literature but not in bamp) |
| Likelihood | Binomial with population as denominator | Close to Poisson for small rates; document the difference at the oldest ages |
| Output | Posterior draws of rates by age and period | Fits directly into draw-based PIs, ensembles as mixtures and pathwise e_x (Phase 2), once mortensor stores draws |
| Runtime | Minutes per group per fit | Unsuitable for the inner CV grid; evaluate with fixed settings on a few rolling origins, or select within the training data with DIC |
| Known issues in 3.0.0 | Prediction fails for `rw1+het`/`rw2+het` models (`predict_apc()` handles only `rw1`/`rw2`); `bamp()` runs the convergence check silently unless `verbose = TRUE`; the Pólya-Gamma variables are drawn from a moment-matched normal approximation, and whether the subsequent Metropolis-Hastings step restores exact invariance was not verified; limited test suite | Restrict the wrapper to RW1/RW2 with optional overdispersion, call `checkConvergence()` explicitly and fail on non-convergence, pin `bamp (>= 3.0.0)`, add mortensor-side integration tests, report the issues upstream |
| Licence | GPL-3 | Keep in `Suggests`, call only exported functions, guard with `requireNamespace()`; do not copy code into the MIT-licensed package (not legal advice) |

### 8.4 Proposed integration

1. Prerequisites: count-based `mort_data`, the model plug-in interface and
   draw-based outputs (Phases 1, 2 and 4).
2. Wrapper `model_bamp()` implementing the plug-in contract: transpose each
   group slice to periods x ages, fit with user-chosen RW orders and
   `overdisp`, check convergence, call `predict_apc()` with future exposures,
   and return draws as an array \[age x horizon x group x draw\].
3. Evaluation: include bamp in rolling-origin backtests with fixed MCMC
   settings and compare interval calibration with Lee-Carter and the
   ensemble. Riebler and Held (2017) found well-calibrated Bayesian APC
   projections when holding out the last 10 years of lung cancer mortality
   data; replicating that kind of check on mortensor's data is the relevant
   benchmark.
4. Effort: about 1-2 weeks for wrapper, tests and a vignette once the
   prerequisites exist; parallelise across groups.

## 9. References

All references below were checked against bibliographic records (Consensus
search index, CRAN package metadata, or the scoringutils source) during this
review.

- Booth H, et al. (2006). Lee-Carter mortality forecasting: a multi-country
  comparison of variants and extensions. *Demographic Research*.
- Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
  forecasts in an interval format. *PLoS Computational Biology* 17(2):
  e1008618. doi:10.1371/journal.pcbi.1008618
- Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting mortality
  rates. *Statistical Modelling* 4(4): 279-298.
- Hyndman RJ, Koehler AB (2006). Another look at measures of forecast
  accuracy. *International Journal of Forecasting*.
- Hyndman RJ, et al. (2013). Coherent mortality forecasting: the product-ratio
  method with functional time series models. *Demography*.
- Hyndman RJ, Ullah MS (2007). Robust forecasting of mortality and fertility
  rates: a functional data approach. *Computational Statistics & Data
  Analysis*.
- Knorr-Held L, Rainer E (2001). Projections of lung cancer mortality in West
  Germany: a case study in Bayesian prediction. *Biostatistics*.
- Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality. *Journal
  of the American Statistical Association*. doi:10.2307/2290201
- Li N, Lee R (2005). Coherent mortality forecasts for a group of populations:
  an extension of the Lee-Carter method. *Demography*.
- Riebler A, Held L (2017). Projecting the future burden of cancer: Bayesian
  age-period-cohort analysis with integrated nested Laplace approximations.
  *Biometrical Journal*.
- Rizzi S, Gampe J, Eilers PHC (2015). Efficient estimation of smooth
  distributions from coarsely grouped data. *American Journal of
  Epidemiology*. doi:10.1093/aje/kwv020
- Schmid VJ, Held L (2007). Bayesian age-period-cohort modeling and
  prediction - BAMP. *Journal of Statistical Software* 21(8).
  doi:10.18637/jss.v021.i08
- Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast
  reconciliation for hierarchical and grouped time series through trace
  minimization. *Journal of the American Statistical Association*.
  doi:10.1080/01621459.2018.1448825
