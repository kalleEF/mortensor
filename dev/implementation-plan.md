# mortensor implementation plan

Date: 2026-10-12. Companion to `dev/revamp-plan.md` (architecture and
decisions) and `dev/production-review.md` (defects D1-D26). This document turns
the revamp into work packages that can be implemented, reviewed and merged one
at a time.

## 1. How to use this plan

- Every **work package (WP)** is sized for one pull request (roughly 300-1,000
  lines including tests) and has an ID, dependencies, deliverables, tests and
  acceptance criteria. Larger WPs are split into numbered sub-PRs.
- WPs are grouped into phases that end in a release: **0.1.1** (stabilise the
  current code), **0.2.0** (first release of the new design), **0.3.0** (all
  models and methods of 0.1.0 ported), **1.0.0** (validated).
- Effort is given in developer-days for one person familiar with R package
  development. Totals: about 183 days, roughly 9 months for one full-time
  developer (Phase 0: 8, Phase 1: 47, Phase 2: 34, Phase 3: 12, Phase 4: 26,
  Phase 5: 34, Phase 6: 22).
- Progress is tracked in `dev/progress.md` (one line per WP: status, PR link,
  notes); each PR updates it.

### Definition of done (every PR)

1. `R CMD check --as-cran` with 0 errors, 0 warnings, 0 notes on the CI matrix.
2. New code covered by tests (target 90% line coverage for new files);
   invariant tests where the WP lists them.
3. Roxygen documentation with a runnable example (under 5 s) for every exported
   function; NEWS entry.
4. `lintr` clean with the project configuration; no `cat()`/`print()` for
   messages (use `cli`); no files written outside `tempdir()` in tests.
5. Defects listed under "Fixes" have a regression test named after the defect
   (e.g. `test_that("D10: life_expectancy uses interval widths", ...)`).
6. `dev/progress.md` updated.

## 2. Branching and release mechanics

| Branch | Purpose |
|---|---|
| `main` | 0.1.x maintenance only (Phase 0); tagged `v0.1.1` |
| `dev-0.2` | Integration branch for the rewrite; starts with an empty `R/` tree plus the Phase 0 infrastructure |
| `wp/<id>-<slug>` | One branch per WP, PR into `dev-0.2`, squash merge |

0.1.x code is kept for reference at tag `v0.1.1` and as the regression
baseline (WP P0.5). Code is ported by copying the relevant functions into the
new files and refactoring them under the new contracts; nothing is sourced
from the old tree. Release 0.2.0 merges `dev-0.2` into `main`.

## 3. Conventions

| Topic | Convention |
|---|---|
| Naming | Exported verbs carry the `mt_` prefix (`mt_data()`, `mt_regroup()`, `mt_ungroup()`, `mt_smooth()`, `mt_fit()`, `mt_forecast()`, `mt_backtest()`, `mt_score()`, `mt_tune()`, `mt_reconcile()`). Model constructors are plain nouns (`lee_carter()`, `tucker3()`, `apc_bayes()`, `rwd()`). Derived-quantity helpers are plain (`life_expectancy()`, `prob_dying()`). This avoids clashes found during planning: `smooth()` exists in `stats`, `ungroup()` in `dplyr` |
| Name-clash test | A test asserts that no export clashes with `base`, `stats`, `utils`, `dplyr`, `data.table`, `forecast`, `generics`, `demography`, `StMoMo`, unless it is a deliberate method on a re-exported generic |
| Classes | `mt_axis`, `mt_strata`, `mt_obs`, `mt_data`, `mt_spec` (+ `mt_spec_<model>`), `mt_fit`, `mt_draws`, `mt_prep`, `mt_backtest`, `mt_score`; each with a validator `validate_<class>()` (internal), `print()`, `format()`, and `as.data.table()` where tabular |
| Internals | `data.table` for tables, `Matrix` for sparse operators, base arrays for draws |
| Errors | `rlang::abort()` with classes `mt_error_<topic>` (e.g. `mt_error_grid`, `mt_error_model_fit`); messages built with `cli` |
| Messages | `cli` only, controlled by `options(mortensor.verbose = FALSE)`; no output by default |
| Randomness | Functions that simulate take `seed = NULL`; parallel code uses `future.apply` with `future.seed = TRUE` |
| Imports | `@importFrom` only; no `library()` in code; `requireNamespace()` guards for Suggests |
| Generics | `forecast()` and `accuracy()` are not defined; `mt_forecast()` is the package verb. If `forecast()` methods are offered for convenience, they are registered on `generics::forecast` |
| Plots | Return ggplot objects; never `print()` inside a function; no files written |
| Licence headers | Code adapted from MIT-licensed packages (`ungroup`) or GPL packages (`bamp`) carries a comment with source, version, licence and authors |

A `CLAUDE.md` file in the repository root (WP P0.6) records these conventions
for Claude Code sessions.

## 4. Data contracts

These contracts are fixed in Phase 1 and changed only through a dedicated PR.

### 4.1 `mt_axis`

```r
mt_axis(lower, upper, open_last = TRUE, unit = "year", labels = NULL)
# single_ages(from = 0, to = 100, open = TRUE) -> 0, 1, ..., 99, 100+
# age_groups(breaks, open = TRUE)             -> e.g. c(0, 1, 5, 10, ..., 90) with 90+
# single_years(from, to); period_groups(breaks)
# parse_age_labels(c("<1", "01-14", "15-19", "90+", "85 and over"))
```

- Fields: `lower` (numeric), `upper` (numeric, `Inf` for an open last
  interval), `width` (derived), `labels` (character), `open_last`.
- Invariants: strictly increasing, contiguous (`upper[i] == lower[i + 1]`),
  only the last interval may be open; `omega` (assumed maximum age) stored
  separately for computations that need a finite upper bound (default 110).
- `parse_age_labels()` covers the label styles seen so far: ONS (`"30-34"`,
  `"90+"`), the owner's German data (`"<1"`, `"01-14"`, `"90+"`), HMD
  (`"0"`, `"1-4"`, `"110+"`).

### 4.2 `mt_strata`

- A `data.table` of grouping variables (e.g. `sex`, `imd`, `region`,
  `cause`) with one row per stratum and a stratum id.
- Optional hierarchies: `hierarchy(cause = all_cause ~ chd + stroke + ...)`,
  `hierarchy(region = Germany ~ West + East + Berlin)`. Used for validation
  (sums) and reconciliation.

### 4.3 `mt_obs` (one observation set)

| Column | Type | Meaning |
|---|---|---|
| `age_lo`, `age_hi` | numeric | Age interval (`age_hi = Inf` for open) |
| `period_lo`, `period_hi` | numeric | Period interval in years (`[2010, 2011)`) |
| strata columns | factor | As declared in `mt_strata` |
| `value` | numeric | Count or rate |

Attributes: `measure` (`"deaths"`, `"exposure"`, `"population_midyear"`,
`"population_jan1"`, `"rate"`, `"qx"`), `source`, `notes`. Validation:
no overlapping intervals within a stratum and period, non-negative values,
integer-valued deaths (warning if not).

### 4.4 `mt_data`

```r
list(obs = list(<mt_obs>, ...),        # named by dataset id
     grid = list(age = <mt_axis>, period = <mt_axis>, strata = <mt_strata>),
     omega = 110,
     ops = <environment>)               # lazily built sparse operators, cache
```

`mt_data()` + `mt_add_deaths()`, `mt_add_exposure()`, `mt_add_population()`,
`mt_add_rates()`; `mt_set_grid()` (default single years of age and calendar
years); `mt_window(data, periods)` returns a copy restricted to a training
window (used by backtests).

### 4.5 Aggregation operator

`mt_operator(obs, grid)` returns a `dgCMatrix` with one row per observed cell
and one column per latent cell (age x period x stratum, column-major), entry 1
where the latent cell lies inside the observed interval. Built once per
observation set and cached in `ops`. Used by `mt_regroup()`,
`mt_aggregate_to()` and the PCLM engine.

### 4.6 `mt_spec`, `mt_fit`, `mt_draws`

```r
new_mt_spec(model, params, requirements, by = NULL)   # constructor helper
# requirements: list(input = "counts" | "rates", min_periods, joint_strata,
#                    allows_missing, needs_single_years)

mt_fit(spec, data, prep = NULL, ...)        # generic dispatch on spec class
# -> list(spec, prep_fitted, train_window, fits (per stratum or joint),
#         diagnostics = data.table(stratum, status, message, seconds))

mt_forecast(fit, h, n_draws = 500, exposure_future = NULL, seed = NULL)
# -> mt_draws
```

`mt_draws`: array `[age, period, stratum, draw]` of rates on the latent grid
with dimnames, plus `meta` (model, grid, seed, train window, horizon) and
optional `deaths` array when `exposure_future` is supplied. Memory guide: 71
ages x 15 years x 10 strata x 500 draws is about 43 MB; larger problems use
`mt_summarise()` early or chunk by stratum.

### 4.7 `mt_prep` (preprocessing inside the model pipeline)

Preprocessing steps are part of what is fitted, so backtests can refit them on
each training window (prevents review defect D4):

```r
prep <- mt_prep(
  step_ungroup(method = "pclm", dims = 1),
  step_smooth(method = "pspline"),
  step_shocks(periods = 2020:2022, action = "exclude"))
fit <- mt_fit(lee_carter(), data, prep = prep)
```

Each step has `prep_fit(step, data)` (learns on the training window) and
`prep_apply(step_fitted, data)`.

## 5. Phase 0: stabilise 0.1.x (release 0.1.1, on `main`)

| WP | Title | Depends | Deliverables | Tests and acceptance | Fixes | Days |
|---|---|---|---|---|---|---|
| P0.1 | Generics and `%\|\|%` | none | `@importFrom rlang %\|\|%`; remove own `forecast()`/`accuracy()` generics, import and re-export from `generics`, register methods | CI on R 4.1 and 4.3: `print(example_mortality)` works; `library(forecast); forecast(cal)` dispatches to `forecast.mort_calibration` | D1, D2 | 1 |
| P0.2 | Package metadata and licence | none | `License: GPL-3`; GPL text via `usethis::use_gpl3_license()`; remove MIT LICENSE; fix undefined globals; remove `.claude/settings.local.json` and `tests/testthat/Rplots.pdf`, add to `.gitignore`/`.Rbuildignore` | R CMD check without the DESCRIPTION warning and the globals note | D17, D26 (part) | 1 |
| P0.3 | CI | P0.1 | `r-lib/actions` check workflow (Ubuntu R devel/release/oldrel and R 4.1; macOS and Windows release), test-coverage workflow, lintr workflow | All green on `main` | D26 (part) | 1.5 |
| P0.4 | Executable vignettes | P0.1 | `eval = TRUE`; fast settings (`rank_selection = "bic"`, few models, `n_sim = 50`); remove chunks that need ONS data or rewrite them on `example_mortality`; replace hard-coded "example output" with real output | All vignettes build in CI in under 10 minutes; no chunk errors | D3 | 3 |
| P0.5 | Regression baseline | P0.1 | `dev/baseline/capture_baseline.R`: point forecasts of every 0.1.0 model with fixed hyperparameters on `example_mortality` (committed as `tests/testthat/fixtures/baseline_0.1.rds`) and, locally, on the owner's data (not committed) | Script reruns reproducibly (fixed seed); fixture documented | protects 0.1.0 behaviour | 1 |
| P0.6 | Repository hygiene for the rewrite | P0.2 | `CLAUDE.md` (conventions from Section 3), `dev/progress.md`, issue templates, `dev-0.2` branch created | Present and referenced from README | | 0.5 |

Phase 0 total: about 8 days.

## 6. Phase 1: data layer and harmonisation (on `dev-0.2`)

| WP | Title | Depends | Deliverables | Tests and acceptance | Fixes | Days |
|---|---|---|---|---|---|---|
| D1.1 | Skeleton | P0.6 | New DESCRIPTION (Imports: data.table, Matrix, rlang, cli, generics, ggplot2, forecast, rTensor, compositions; Suggests: bamp, demography, StMoMo, ungroup, MortalityLaws, scoringutils, posterior, future.apply, testthat, knitr, rmarkdown), `aaa-imports.R`, `utils-cli.R`, `utils-checks.R`, error classes, name-clash test | Package installs and checks with no functions | | 2 |
| D1.2 | `mt_axis` | D1.1 | Constructors, `parse_age_labels()`, validation, print, interval utilities (`covers()`, `overlaps()`, `widths()`, `midpoints()`) | Unit tests for every label style listed in 4.1; invariants (contiguity, single open interval) | basis for D10, D15 | 3 |
| D1.3 | `mt_strata` and hierarchies | D1.1 | Constructor from data, crossed and nested hierarchies, consistency checks | Hierarchy sums validated on the owner's German structure (regions and causes) using a synthetic fixture | basis for D24 | 3 |
| D1.4 | `mt_obs` and `mt_add_*()` | D1.2, D1.3 | Long-table ingestion with column mapping, measure types, exposure types | Round trip from long table to `mt_obs` and back; malformed inputs (duplicate cells, gaps, overlapping intervals) give classed errors | D15 | 4 |
| D1.5 | `mt_data`, `mt_set_grid()`, `mt_window()`, `mt_describe_granularity()`, `mt_validate()` | D1.4 | Container and checks: coverage, deaths > exposure, zero exposure, rate consistency, hierarchy sums | Validation report on fixtures; the owner's doubled all-cause counts are flagged by the hierarchy-sum check | | 4 |
| D1.6 | Aggregation operator, `mt_regroup()`, `mt_aggregate_to()` | D1.5 | Sparse operator builder and cache; exact regrouping of counts; rates recomputed from summed counts | Invariants: totals preserved; open interval mapped correctly; open interval changing over time; performance for 110 ages x 50 years x 50 strata under 1 s | D21 (constraint) | 4 |
| D1.7 | Exposures | D1.5 | `mt_as_exposure()` from mid-year or January-1 populations; `exposure_future` container for external projections | Known-value tests; declared approximation stored in metadata | | 2 |
| D1.8 | PCLM engine (`harmonise-pclm.R`) | D1.6 | Native 1D and 2D penalized composite link model adapted from `ungroup` 1.6.3 (MIT, attribution header) using the operator from D1.6; Poisson IRLS; B-spline basis; difference penalty; lambda by BIC/AIC on a log grid (default 1e-2 to 1e8) with a boundary warning; offset handled without the count-based artificial-bin reallocation; standard errors | (a) Counts without offset agree with `ungroup::pclm()` within 1%; (b) noise-free rate test recovers rates within 1% (the 6.1% shift of `ungroup` 1.6.3 is absent); (c) England & Wales benchmark (`StMoMo::EWMaleData`, grouped to 5-year + 90+): median absolute log error at most 0.02 (30-89) and 0.03 (90-100); (d) fitted vs observed group totals within 1% | finding of revamp plan 5.1 | 8 |
| D1.9 | `mt_ungroup()` user API | D1.8 | Methods: `"pclm"` (rates with exposure offset, default), `"monotone"` (cumulative Hyman spline, used for populations), dims 1 or 2; refuses 2D on death counts without exposures with an explanatory error; exposures ungrouped first when only grouped; optional bootstrap replicates; provenance recorded per cell | Reproduces the revamp-plan 5.1 results on the fixtures; refusal tested | D21 | 4 |
| D1.10 | `mt_smooth()` and `step_*()` framework (`mt_prep`) | D1.9 | Poisson P-spline smoothing with log-exposure offset (1D/2D); `step_ungroup()`, `step_smooth()`, `step_shocks()` with fit/apply semantics | Steps fitted on a window never read data outside it (test with a sentinel value in later years) | D4 | 5 |
| D1.11 | Interval-aware life tables | D1.2 | `life_table(rates, axis, ax = NULL, open = c("reciprocal", "kannisto", "gompertz"))` using interval widths and midpoints; documented default | Identities (l_x, d_x, L_x, T_x consistency); single-year vs abridged agreement on a smooth schedule within 0.1 years of e0; open-interval methods checked against an exact Gompertz integral | D10, D20 | 4 |
| D1.12 | Adapters and example data | D1.5 | `mt_read_long()` generic reader; `mt_read_ons()` (ported loaders); simulated multi-granularity dataset `mt_example` with known truth (generator in `data-raw/`); keep `example_mortality` | Adapters tested on small fixtures; `mt_example` documented | | 4 |

Phase 1 total: about 47 days. D1.2-D1.4 and D1.11 can run in parallel after
D1.1.

## 7. Phase 2: model contract, benchmarks, Lee-Carter, Tucker3, evaluation

| WP | Title | Depends | Deliverables | Tests and acceptance | Fixes | Days |
|---|---|---|---|---|---|---|
| M2.1 | Model contract and registry | D1.10 | `new_mt_spec()`, `mt_fit()`/`mt_forecast()` generics, requirement checks, per-stratum loop (`future.apply`), failures recorded in diagnostics (never NA or 0), `mt_register_model()` | A dummy model passes the contract test suite (`expect_mt_model()` helper reused by every model) | D6, D19, D26 | 4 |
| M2.2 | `mt_draws` and summaries | M2.1 | Class, `mt_summarise()` (mean, median, quantiles), `as.data.table()`, `mt_aggregate_to()` on draws, deaths from `exposure_future` with Poisson noise | Quantiles of known simulated draws; aggregation of draws preserves totals per draw | D22 (basis) | 3 |
| M2.3 | Benchmarks `rwd()`, `arima_series()`, `ets_series()` | M2.2 | Per-series fits on log rates; draws via `simulate()` of the fitted model (uncertainty grows with horizon); valid ETS specifications only; distinct drift windows | Coverage test on the simulated random-walk generator (from the review): 95% intervals cover 0.90-0.98 averaged over h = 1-10, no horizon below 0.85 | D8, D13 | 4 |
| M2.4 | `lee_carter()` | M2.2 | Native Poisson log-bilinear fit (Brouhns et al. 2002) and classic SVD fit (Lee and Carter 1992); k_t as random walk with drift; draws include k_t innovation and drift uncertainty, optional semi-parametric bootstrap of a_x, b_x (Brouhns et al. 2005) | Agreement with `StMoMo::lc()` (Poisson) and `demography::lca()` (SVD) on England & Wales within numerical tolerance; coverage test as in M2.3 | D8 | 6 |
| M2.5 | `tucker3()` port | M2.2 | Port `tucker3_coda_forecast()`, `tucker3_direct_forecast()`, `forecast_B_mode()`; CoDa variant keeps the open age group (auxiliary series for the open-interval rate); rank grids valid for 1 and 2 groups; draws from simulated B-mode paths plus residual noise | Regression vs 0.1.0 baseline for Tucker3-Direct (same ranks) within 1e-6; CoDa round trip reproduces all ages including the open interval; one-group data works; coverage test | D11, D12 | 6 |
| M2.6 | `mt_backtest()` | M2.1, D1.10 | Rolling-origin evaluation that refits `mt_prep` and the spec at each origin on `mt_window()` data; parallel over origins; caching of fits | Leakage test: changing data after the origin does not change the forecast; on `example_mortality` the honest MASE of the review (about 0.7) is reproduced | D4 | 5 |
| M2.7 | `mt_score()` | M2.6 | CRPS (sample-based), WIS (Bracher et al. 2021 normalisation), log score where available, coverage and PIT by horizon, MASE against a declared baseline; via `scoringutils` when installed, native fallback | WIS equals the published definition on hand-computed cases | D16 | 3 |
| M2.8 | `mt_tune()` | M2.6 | Grid search over spec parameters by inner rolling origin within the training window, or by information criterion where the model offers one | Nested use inside `mt_backtest()` touches only training data | D4, D5 | 3 |

Phase 2 total: about 34 days.

## 8. Phase 3: seamless bamp wrapper

| WP | Title | Depends | Deliverables | Tests and acceptance | Days |
|---|---|---|---|---|---|
| B3.1 | `apc_bayes()` spec and data conversion | M2.1 | Constructor (`age`, `period`, `cohort` in `"rw1"`, `"rw2"`; `overdispersion`; `mcmc = mcmc_control("fast" \| "final" \| custom)`); conversion per stratum to bamp matrices (periods x age groups), `periods_per_agegroup` from the axes; integer rounding recorded | Conversion tested on fixtures with 5-year and single-year ages | 3 |
| B3.2 | Fitting and convergence policy | B3.1 | Calls `bamp::bamp()` with captured output; explicit `bamp::checkConvergence()`; automatic refit with longer chains up to a limit; classed error with R-hat table if still not converged; core budget shared between strata and chains | Fixture fit converges; forced non-convergence (very short chains) produces the classed error | 3 |
| B3.3 | Forecast draws | B3.2 | `bamp::predict_apc()` samples to `mt_draws`; mapping from age groups to the single-year grid through the D1.9 ungrouping applied draw by draw; deaths with binomial noise from `exposure_future` | Draw dimensions and metadata; mean of draws matches bamp's own posterior summary | 3 |
| B3.4 | Effects, selection, documentation | B3.3 | `tidy()`/`plot()` of APC effects (bamp conventions stated); `mt_tune(apc_bayes(), method = "dic")` via `bamp::selectModel()`; heterogeneity variants blocked for forecasting with an explanatory error until bamp supports them | Tests skip without bamp; CI job with bamp installed runs a short fixture fit (under 60 s) | 3 |

Phase 3 total: about 12 days. Can start once M2.1 and M2.2 are merged, in
parallel with M2.3-M2.8.

## 9. Phase 4: composition, derived quantities, release 0.2.0

| WP | Title | Depends | Deliverables | Tests and acceptance | Fixes | Days |
|---|---|---|---|---|---|---|
| C4.1 | `mt_ensemble()` | M2.6 | Mixture of component draws; weights equal, inverse score, or stacking on inner-CV scores; components with failures are dropped per stratum with a diagnostic and weights renormalised | Mixture spread test (two components give the width of their mixture); weights never computed on the outer holdout | D5, D6, D7, D9 | 5 |
| C4.2 | Derived quantities | M2.2, D1.11 | `life_expectancy()`, `prob_dying()`, `deaths()`, `asr()` (standard populations shipped with source and licence noted), `yll()`; all pathwise over draws | Intervals for e_x come from draws; single-year and grouped results agree on a smooth schedule | D10, D22 | 5 |
| C4.3 | Plots | M2.2, M2.7 | Fan charts, Lexis-surface heatmaps, age profiles, backtest calibration (coverage by horizon, PIT), APC effects; all return ggplot objects | `vdiffr` snapshot tests | D26 | 4 |
| C4.4 | `mt_export()` | M2.2 | Tidy CSV/parquet export of summaries and, optionally, draws; file names reflect the grid | Round trip | D23 | 1 |
| C4.5 | Vignettes | all above | Getting started (end to end on `mt_example`); data and granularity; ungrouping (the evidence from revamp plan 5.1); models and backtesting; Bayesian APC with bamp; derived quantities. Slow chunks precomputed (`.Rmd.orig` pattern) | All vignettes build in CI in under 15 minutes | D3 | 6 |
| C4.6 | Performance pass | C4.1 | Profiling of fit, forecast and backtest; `matrixStats` for quantiles over draws; caching | Backtest of `lee_carter()`, `rwd()`, `tucker3()` on `mt_example` with 4 origins under 2 minutes on 4 cores | D25 | 3 |
| C4.7 | Release 0.2.0 | all | Migration note (0.1 calls to 0.2 calls), NEWS, pkgdown site, CRAN-ready check | Section 13 acceptance criteria met | | 2 |

Phase 4 total: about 26 days.

## 10. Phase 5: port everything that runs in 0.1.0 (release 0.3.0)

Each port reuses the M2.1 contract test suite and the 0.1.0 regression
baseline where behaviour should not change.

| WP | Title | Deliverables | Tests and acceptance | Fixes | Days |
|---|---|---|---|---|---|
| P5.1 | `fdm()` | Wrapper on `demography::fdm()` with draws from simulated coefficient series | Agreement with `demography` forecasts; coverage test | | 3 |
| P5.2 | `stmomo()` APC and CBD | Real ages, initial exposures for logit models, fitted values of the selected model, `simulate()`-based draws, gnm output captured | Agreement with direct StMoMo calls | (review 5.3) | 4 |
| P5.3 | Multi-cause: `tucker4()` and `mt_cause_shares()` | Port Tucker4/ILR coherent cause forecasting; no near-zero residual component when a residual cause exists | Regression vs baseline; cause shares sum to 1 per draw | (review 5.4) | 6 |
| P5.4 | `mt_reconcile()` | WLS (ported) and MinT (Wickramasuriya et al. 2019) reconciliation applied draw by draw over `mt_strata` hierarchies | Reconciled draws satisfy the hierarchy exactly; regression vs 0.1.0 WLS | D24 | 5 |
| P5.5 | `mt_cohere()`, `li_lee()`, `product_ratio()` | The 0.1.0 female <= male constraint as an explicit option on any stratum pair; Li-Lee (Li and Lee 2005) and product-ratio (Hyndman et al. 2013) models | Constraint holds per draw; models agree with reference implementations where available | D24 | 6 |
| P5.6 | NNLS stacking | Fixed stacking (fitted on the scale applied, on inner CV) as an `mt_ensemble()` weighting option | Weights computed only on training windows | D5 | 2 |
| P5.7 | Smoothing options | Gaussian GAM and Whittaker smoothers as `step_smooth()` methods; GAM-Bayes smoothing uncertainty as an option that adds to model draws | Regression vs 0.1.0 smoothers on fixtures | | 3 |
| P5.8 | Disaggregation options | Natural-Gompertz, Hyman-Kannisto and GAM disaggregation as `mt_ungroup()` methods with the exposure-weighted constraint | Benchmarked against PCLM on the D1.8 fixtures (documented as alternatives) | D21 | 3 |
| P5.9 | Release 0.3.0 | Vignette for multi-cause forecasting on a shipped synthetic multi-cause dataset | Section 13 criteria | | 2 |

Phase 5 total: about 34 days.

## 11. Phase 6: validation and release 1.0.0

| WP | Title | Deliverables | Days |
|---|---|---|---|
| V6.1 | Public benchmark | Rolling-origin benchmark of all models and ensembles on openly licensed data (and on the owner's data locally): point accuracy, interval calibration, runtime | 10 |
| V6.2 | Documentation | Reference pages reviewed, methods article draft, FAQ on granularity and ungrouping | 8 |
| V6.3 | CRAN submission and 1.0.0 | Reverse-dependency-free check, submission, fixes | 4 |

Phase 6 total: about 22 days.

## 12. Sequencing

```
Phase 0  P0.1 → P0.2 → P0.3 → P0.4 ─┐
                 P0.5, P0.6 ─────────┴→ dev-0.2
Phase 1  D1.1 → {D1.2, D1.3, D1.11} → D1.4 → D1.5 → D1.6 → D1.7
                                                     └→ D1.8 → D1.9 → D1.10
                                      D1.5 → D1.12
Phase 2  D1.10 → M2.1 → M2.2 → {M2.3, M2.4, M2.5} ; M2.1 → M2.6 → {M2.7, M2.8}
Phase 3  M2.2 → B3.1 → B3.2 → B3.3 → B3.4          (parallel to M2.3-M2.8)
Phase 4  {M2.x, B3.x} → C4.1 … C4.6 → C4.7 = release 0.2.0
Phase 5  P5.1 … P5.8 (largely independent) → P5.9 = release 0.3.0
Phase 6  V6.1 → V6.2 → V6.3 = release 1.0.0
```

Critical path: P0.1 → D1.1 → D1.2 → D1.4 → D1.5 → D1.6 → D1.8 → D1.9 →
D1.10 → M2.1 → M2.6 → C4.1 → C4.7. With two people, the data layer (D1.x) and
the model ports (M2.4, M2.5 prototyped on arrays) can overlap.

## 13. Release acceptance criteria

| Release | Criteria |
|---|---|
| 0.1.1 | R CMD check clean on R 4.1-devel; vignettes execute; review defects D1, D2, D3, D17 fixed |
| 0.2.0 | All Phase 1-4 WPs merged; R CMD check clean on the CI matrix; vignettes build under 15 minutes; on the simulated random-walk generator, 95% intervals of `rwd()`, `lee_carter()`, `tucker3()` and the ensemble cover 0.90-0.98 averaged over horizons 1-10 and at least 0.85 at every horizon; ungrouping meets the D1.8 benchmark thresholds; backtest of the honest example reproduced; regression baseline matched or deviations documented with the defect that caused them |
| 0.3.0 | Every model and method listed in revamp plan Section 13 available behind the contract; regression baseline for ported models matched or explained |
| 1.0.0 | Public benchmark published in the documentation; CRAN acceptance |

## 14. Testing strategy

| Level | Content | Runs |
|---|---|---|
| Unit | Each function's behaviour and errors | Every PR |
| Contract | `expect_mt_model(spec)` checks fit/forecast/draw shapes, failure handling, reproducibility with a seed | Every model PR |
| Invariant | Totals preserved by regrouping and aggregation; life-table identities; hierarchy sums after reconciliation; draws reproducible with a seed | Every PR |
| Benchmark | Ungrouping on England & Wales and simulation (thresholds in D1.8); interval coverage on the simulated generator (M2.3) | Nightly or on demand (`MORTENSOR_SLOW_TESTS=true`) |
| Regression | 0.1.0 baseline fixture for ported models | Every model PR |
| Snapshot | Plots (`vdiffr`) and printed summaries (`expect_snapshot()`) | Every PR |

Fixtures stay small (under 1 MB). Data from other packages (e.g.
`StMoMo::EWMaleData`) is loaded in tests via Suggests, never copied into the
repository.

## 15. Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Native PCLM 2D is slow for large grids (the review saw 6-12 minutes per series with `ungroup`) | Medium | Medium | 1D default; sparse operators; generalised linear array model (GLAM) structure for 2D as a later optimisation; cache fits |
| bamp runtime makes backtests impractical | High | Medium | `mcmc_control("fast")` presets, caching, fewer origins for APC; document runtime |
| bamp API changes (young 3.0.0 release) | Medium | Medium | Version pin, wrapper isolates bamp calls in one file, integration test in CI |
| Coverage targets not met for Tucker3 | Medium | High | Coverage test is part of acceptance; add parameter uncertainty (bootstrap) if needed; document if a model cannot meet it |
| Scope creep from 0.3.0 ports into 0.2.0 | Medium | Medium | Release criteria fixed in Section 13; ports only after 0.2.0 |
| Single developer availability | Medium | High | WPs small and self-contained; `dev/progress.md` and `CLAUDE.md` allow handover between people and Claude sessions |
| Owner's data issues (e.g. doubled all-cause counts) | Observed | High | Hierarchy-sum and plausibility checks in `mt_validate()` (D1.5) |

## 16. Working with Claude Code sessions

- One WP per session and PR; the session reads `CLAUDE.md`,
  `dev/progress.md` and the WP row in this plan first.
- Each session ends by updating `dev/progress.md` with status, decisions and
  open questions, so the next session can continue without the full history.
- Benchmarks that need R >= 4.4 (`ungroup` 1.6.3) or large data run in a
  separate session or locally; results are written to `dev/` as tables and
  figures.

## 17. Open questions for the owner

1. Minimum R version: keep R >= 4.1 (native pipe) or raise to 4.3?
2. CRAN as a target for 0.2.0, or GitHub-only until 1.0.0?
3. German single-year population (Destatis 12411-0005 or similar) available
   for the rate-based ungrouping and backtests on your data?
4. Lee-Carter: confirm a native implementation, validated against `StMoMo`
   and `demography`.

## 18. References

Verified during this review (see also `dev/production-review.md` and
`dev/revamp-plan.md`).

- Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
  forecasts in an interval format. *PLoS Computational Biology* 17(2):
  e1008618. doi:10.1371/journal.pcbi.1008618
- Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear regression
  approach to the construction of projected lifetables. *Insurance:
  Mathematics and Economics* 31: 373-393.
- Brouhns N, et al. (2005). Bootstrapping the Poisson log-bilinear model for
  mortality forecasting. *Scandinavian Actuarial Journal*.
- Hyndman RJ, et al. (2013). Coherent mortality forecasting: the product-ratio
  method with functional time series models. *Demography*.
- Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality. *Journal
  of the American Statistical Association*. doi:10.2307/2290201
- Li N, Lee R (2005). Coherent mortality forecasts for a group of populations:
  an extension of the Lee-Carter method. *Demography*.
- Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast
  reconciliation for hierarchical and grouped time series through trace
  minimization. *Journal of the American Statistical Association*.
  doi:10.1080/01621459.2018.1448825
