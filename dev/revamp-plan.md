# mortensor revamp plan

Date: 2026-10-12. Status: **agreed direction**; the owner's decisions are
recorded in Section 14. Builds on the findings in `dev/production-review.md`
(referred to as "the review"; defect IDs D1-D26).

## 1. Goal

Turn mortensor into a suite of composable functions that:

1. accepts mortality and population data at **different granularities**
   (age groupings, open intervals, period lengths, strata, causes) and keeps
   track of what was observed at which resolution;
2. builds **complex forecasts** from interchangeable parts: data
   harmonisation, one or more models, ensembles, coherence and reconciliation
   constraints, and derived demographic quantities;
3. includes the **functionality of BAMP** (Bayesian age-period-cohort
   modelling and prediction; Schmid and Held 2007) through a seamless wrapper
   around the `bamp` package, behaving like every other model in the suite;
4. produces **honest, calibrated probabilistic forecasts**: every forecast is
   a set of draws, every evaluation is out of sample, and every derived
   quantity carries its uncertainty.

## 2. Design principles

| Principle | Consequence |
|---|---|
| Counts are the primary data | Deaths and exposures are stored and modelled; rates are derived. Rate-only inputs are supported but flagged |
| Granularity is explicit | Every observation knows its age interval, period interval and strata; nothing assumes 5-year groups or a 90+ open interval |
| One latent grid, many observation grids | Models describe rates on a fine latent grid; each dataset is linked to it by an aggregation operator. Reporting can use any grid |
| Forecasts are draws | Point forecasts, intervals, e_x, deaths and totals are summaries of the same simulated paths |
| Strict temporal separation | Splitting happens before smoothing, tuning or weighting; evaluation code cannot see test data during fitting |
| Plug-in models | A small contract that every model implements; external engines (bamp, StMoMo, demography) are used through wrappers |
| Fail loudly | Errors and per-model failures are reported in a diagnostics table, never converted to NA, 0 or defaults |
| Quiet and pure by default | No `cat()`, no files written outside `export()`, plots return ggplot objects |

## 3. Architecture

```
  inputs (long tables, arrays, HMD/ONS/Eurostat adapters)
        │
        ▼
  ┌──────────────┐   validate, describe granularity, convert measures
  │ 1. DATA      │   mt_data(): observation sets + latent grid + operators
  └──────┬───────┘
         ▼
  ┌──────────────┐   regroup (exact), ungroup (PCLM), align, impute,
  │ 2. HARMONISE │   exposures from populations, shock handling
  └──────┬───────┘
         ▼
  ┌──────────────┐   model specs: APC-Bayes, Lee-Carter family, Tucker3/CoDa,
  │ 3. MODELS    │   P-spline, time series, multi-population; fit per stratum
  └──────┬───────┘   or jointly; tuning on training data only
         ▼
  ┌──────────────┐   forecast draws [age x period x strata x cause x draw]
  │ 4. FORECAST  │   with future exposures; ensembles as mixtures
  └──────┬───────┘
         ▼
  ┌──────────────┐   coherence across strata, reconciliation across
  │ 5. COMPOSE   │   hierarchies (causes, regions), scenarios
  └──────┬───────┘
         ▼
  ┌──────────────┐   life tables, e_x, q_x, deaths, ASR, YLL, summaries,
  │ 6. DERIVE    │   aggregation to any reporting grid, pathwise
  └──────┬───────┘
         ▼
  ┌──────────────┐   rolling-origin backtests, CRPS/WIS/log score,
  │ 7. EVALUATE  │   coverage by horizon, PIT, benchmark tables
  └──────────────┘
```

Each layer has its own S3 classes, validators, `print()`/`summary()`/`plot()`
methods and a tidy `as.data.table()` method.

## 4. Layer 1: data model and granularity

### 4.1 Objects

| Class | Content |
|---|---|
| `mt_axis` | Interval definition for one dimension: lower and upper bounds, open last interval, unit (years). Constructors: `age_axis(breaks, open = TRUE)`, `single_ages(0, 100, open = TRUE)`, `abridged_ages()`, `period_axis(breaks)`, `single_years(2001, 2019)` |
| `mt_strata` | Grouping variables with their levels and relations: crossed (`sex x imd`), nested (`region` within `nation`, `cause` within `all_cause`) |
| `mt_obs` | One observation set in long form: measure (`deaths`, `exposure`, `population_midyear`, `population_jan1`, `rate`, `qx`), interval bounds per cell, strata columns, value, source and flags |
| `mt_data` | A collection of `mt_obs` plus the latent grid (axes and strata) and, computed lazily, sparse aggregation operators that map latent cells to each observed cell |

### 4.2 The aggregation operator

For each observation set, a sparse matrix **C** (package `Matrix`) maps
latent cells to observed cells: an observed count is the sum of latent counts
over the age, period and strata it covers. Examples:

- deaths in 5-year age groups and a 90+ open interval, latent single years
  to 100+: each row sums 5 latent ages, the last row sums ages 90-100+;
- an open interval that changes over time (85+ before 2001, 90+ after):
  rows differ by period, no special code needed;
- cause-specific deaths available only by sex while all-cause deaths are
  by sex and deprivation quintile: rows sum over quintiles;
- 5-year periods in a cancer registry: rows sum over 5 latent years.

This is the composite link structure of Rizzi et al. (2015), applied to every
dimension. It gives two ways to use mixed-granularity data:

1. **Harmonise first:** ungroup coarse data to the latent grid with the
   penalized composite link model (1D, or 2D over age and period as in Rizzi
   et al. 2018), then use any model. Available for every model.
2. **Model through the operator (later):** count-likelihood models put **C**
   into the likelihood and fit the latent rates directly to the coarse data.
   This propagates the uncertainty of ungrouping into the forecast and
   generalises BAMP's unequal age/period grid (Knorr-Held and Rainer 2001).
   It requires native count models and is outside the current scope; the
   data layer is designed so it can be added without changing the API.

### 4.3 Functions

| Function | Purpose |
|---|---|
| `mt_data()`, `add_deaths()`, `add_exposure()`, `add_population()`, `add_rates()` | Build a data object from long tables or arrays; declare columns for age, period, strata, cause |
| `set_grid()` | Choose the latent grid (default: single years of age and calendar years; coarser grids on request) |
| `validate()` | Coverage, duplicates, overlapping intervals, gaps, negative or zero exposures, deaths > exposure, consistency between deaths, exposure and rates |
| `describe_granularity()` | Table of what is observed at which resolution, by dataset and period |
| `as_exposure()` | Person-years from mid-year or January-1 populations with a declared approximation |
| `read_hmd()`, `read_ons_*()`, `read_eurostat()` | Thin adapters into `mt_data`; source-specific code stays out of the core |

## 5. Layer 2: harmonisation

| Function | Method | Notes |
|---|---|---|
| `regroup()` | Exact summation of deaths and exposures through **C** | Rates are recomputed from summed counts (exposure-weighted by construction; fixes D21) |
| `ungroup()` | Penalized composite link model on **rates** (deaths with single-year exposures as offset), 1D per year by default, 2D on request (Rizzi et al. 2015, 2018); monotone spline of cumulative counts for exposures | Native implementation adapted from `ungroup` (MIT, compatible with GPL-3) with the fixes in Section 5.1; optional bootstrap replicates to propagate ungrouping uncertainty |
| `align()` | Bring several observation sets to the latent grid (regroup where finer, ungroup where coarser) | Records the operation per cell |
| `smooth()` | Poisson P-spline in age and period with log-exposure offset (Currie et al. 2004) | Applied to training data only; raw data are kept |
| `flag_shocks()` | Mark periods (e.g. 2020-2022) for exclusion, down-weighting or an indicator term | Explicit, model-aware, never silent |
| `impute_missing()` | Model-based imputation or exclusion | Models declare whether they accept missing cells |

### 5.1 Ungrouping: evidence and chosen approach

Single-year ages are the default grid (Section 14), so ungrouping 5-year
death data is on the critical path. The owner reported unstable output from
2D ungrouping of German data over time. Three checks during this review
explain why and what works (`ungroup` 1.6.3, R 4.5.3):

1. **Simulation with known single-year truth** (Kannisto-type hazards with
   age-specific improvement, cohort-size waves in the population, Poisson
   deaths grouped to 30-34 ... 85-89 and 90+; 5 replicates, large and small
   populations).
2. **England & Wales men 1980-2011** (`StMoMo::EWMaleData`, real single-year
   deaths and exposures), grouped the same way and ungrouped back, so every
   estimate can be checked against observed single-year rates.
3. **German deaths 1998-2024 by cause and region** (owner's data, deaths only):
   10-year groups were ungrouped and re-aggregated to the observed 5-year
   groups.

**What goes wrong:**

| Cause | Evidence |
|---|---|
| With an exposure offset, `pclm()`/`pclm2D()` return rates that are shifted by a common factor: the artificial bin appended after the open interval is removed with a count-based reallocation (`delete.artificial.bin()`) that is applied to rates | Noise-free test: every single-year rate 6.1% too low (factor 0.939) and every implied group total 6% short; in 2D the factor differs by year |
| `pclm2D()` defaults: `lambda = c(1, 1)` fixed (not estimated), knots every 7 ages (`kr = 7`), search range for lambda only 0.1-1000 | Defaults give 4.8-6.6% error in simulation; with estimation the lambda search hit the upper bound in simulation |
| 2D smoothing of **death counts** cannot represent cohort-size waves, which run diagonally through the age-period surface (Germany: deaths of women aged 85-89 fell 37% from 1999 to 2004 and recovered by 2009 as the 1915-1919 birth cohorts passed through) | England & Wales, deaths only, 2D with lambda estimated: median absolute log error 0.64 at 30-89 and 1.25 at 90-100 (rates off by factors of about 1.9 and 3.5), group totals off by a factor of about 2; German all-cause series: log errors of 0.37-0.89 (factors 1.4-2.4) and year-to-year roughness 3-4 times the data's |
| Deaths alone do not identify the split of the open 90+ group | England & Wales, deaths only, 1D per year: median absolute log error 0.69 at 90-100 (a factor of about 2), biased low |
| Point smoothers (`demography::smooth.demogdata()`) do not use the grouping constraint | England & Wales: 2.1% error at 30-89 but rates at 90-100 about 20% too low; group totals off by 1-2% |

**What works** (England & Wales, real single-year truth; errors are median absolute log errors, which for small values equal relative errors):

| Method | Error 30-89 | Error 90-100 | Year-to-year roughness vs grouped data | Group totals |
|---|---|---|---|---|
| PCLM on rates with single-year exposures, 1D per year, rescaled per year | 1.5% | 2.4% | 1.07 | exact |
| PCLM on rates with single-year exposures, 2D with lambda estimated (`kr = 2`, range 0.01-1e7), rescaled per year | 1.6% | 2.4% | 0.76 | within 0.6% |
| Same 2D without the rescaling | 6.6% | 6.1% | 0.13 | 6.5% off |

The simulation gives the same ranking (2D rate-based with rescaling: 0.8-1.3%
error at 30-89 and 0.9-1.5% at 90-100; 1D per year has 2-5 times the
year-to-year roughness of the truth).

**Decisions for `ungroup()`:**

1. Ungroup **rates**: deaths with single-year exposures as the
   offset. Cohort-size waves then live in the exposures and the rate surface
   is smooth in age and period.
2. If exposures are only grouped, ungroup the population first (monotone
   spline of cumulative population counts, which preserves totals exactly),
   then the deaths as rates. Refuse 2D ungrouping of death counts without
   exposures, with an explanatory error.
3. Correct the level shift: implement PCLM natively (adapted from `ungroup`,
   MIT-licensed, with attribution) so that the artificial bin is removed
   correctly for rates; until then, rescale each year so implied deaths
   match observed deaths.
4. Default 1D per year; 2D as an option for smoother time paths, with lambda
   estimated on a wide range, `kr = 2`, and a warning when lambda hits a
   bound.
5. Set the upper age (omega) explicitly and document the open-interval
   assumption; report fitted-versus-observed group totals as a diagnostic.
6. In backtests, ungroup within the training window only (2D smoothing over
   time would otherwise use future years; review D4).
7. Report the offset issue and the 2D defaults to the `ungroup` maintainers
   (owner's decision).

Figures: `dev/figures/ew_age_profile_2008.png`,
`dev/figures/ew_age95_over_time.png`,
`dev/figures/ungroup_age95_over_time_small.png`.

## 6. Layer 3: models

### 6.1 Contract

Every model is a specification object created by a constructor
(`apc_bayes()`, `lee_carter()`, ...) and implements:

```r
mt_fit(spec, data, ...)              # -> mt_fit object; uses training data only
mt_forecast(fit, h, n_draws, exposure_future = NULL)   # -> mt_draws
fitted(fit); residuals(fit); tidy(fit); diagnostics(fit)
requirements(spec)                   # counts or rates, min periods, joint strata,
                                     # missing cells allowed, supports operator C
```

Models run independently per stratum (`by = c("sex", "imd")`) or jointly
(`joint = TRUE`) when they support it. A registry (`register_model()`) lets
other packages add models. Hyperparameters are tuned with `tune()` on the
training data only (grid, information criteria or inner rolling-origin CV).

### 6.2 Model families

| Family | Constructor | Implementation | Joint strata | Uses **C** |
|---|---|---|---|---|
| Bayesian APC (BAMP functionality, Section 7) | `apc_bayes()` | `bamp` wrapper | no (one fit per stratum) | no (bamp's own unequal-grid handling) |
| Lee-Carter and variants | `lee_carter()` | Poisson log-bilinear fit natively; `demography` wrapper for comparison | via `li_lee()` | later |
| Functional demographic model | `fdm()` | `demography` wrapper | via `product_ratio()` | no |
| Multi-population coherent | `li_lee()`, `product_ratio()` | Li and Lee (2005); Hyndman et al. (2013) | yes | no |
| Generalised APC family | `stmomo()` | `StMoMo` wrapper with real ages and correct exposures | no | no |
| 2D P-spline smoothing and forecasting | `pspline2d()` | Currie et al. (2004), Poisson | no | yes |
| Tensor models (current core idea) | `tucker3(transform = c("coda", "log"))`, `tucker4()` | Port from 0.1.0 after fixing D11, D12 | yes | no |
| Benchmarks | `rwd()`, `arima_series()`, `ets_series()` | `forecast` package, simulated with `simulate()` | no | no |

## 7. Bayesian APC module: a seamless bamp wrapper

Decision (Section 14): BAMP functionality comes from wrapping the `bamp`
package; no native engine is built now. Because mortensor moves to GPL-3,
adapting bamp code with attribution is legally possible where the wrapper
needs it, but the default is to call bamp's exported functions and to report
problems upstream.

### 7.1 What "seamless" means

`apc_bayes()` must behave like every other model: the same constructor style,
the same `mt_fit()` / `mt_forecast()` / `backtest()` / `mt_ensemble()` calls,
the same draws as output. The wrapper does all bamp-specific work:

| Requirement | Wrapper behaviour |
|---|---|
| Interface | `apc_bayes(age = "rw2", period = "rw2", cohort = "rw1", overdispersion = TRUE, mcmc = mcmc_control())` |
| Data conversion | Deaths and exposures per stratum from `mt_data` to bamp's matrices (rows = periods, columns = age groups); exposures rounded to integers with the rounding recorded |
| Age and period resolution | Fit bamp at the resolution the deaths were observed: 5-year age groups with annual periods use `periods_per_agegroup = 5` (Knorr-Held and Rainer 2001); single-year deaths use single-year ages. Forecast draws are mapped to the single-year latent grid by the harmonisation layer, draw by draw |
| Strata | One bamp fit per stratum; strata run in parallel under one core budget shared with bamp's own chains, so the machine is not oversubscribed |
| Convergence | `checkConvergence()` called explicitly after every fit; if it fails, refit with longer chains up to a limit, then stop with the R-hat table. Never return unconverged results silently |
| Forecast draws | `predict_apc()` posterior samples converted to `mt_draws` \[age x period x stratum x draw\]; deaths drawn with binomial noise from external future exposures |
| Missing cells | bamp cannot fit with NA cells: the wrapper stops with a clear message and points to `impute_missing()` or `regroup()` |
| Heterogeneity models | Prediction for `rw1+het`/`rw2+het` fails in bamp 3.0.0: blocked with an informative error until fixed upstream, or served by an adapted prediction routine (possible under GPL-3; decide in Phase 3) |
| Likelihood | bamp is binomial with population as denominator; documented, and the fitted p is treated as the central rate only where rates are small (checked at old ages) |
| Model selection | `tune(apc_bayes(), method = "dic")` via `bamp::selectModel()`, or out-of-sample selection through `backtest()` |
| Effects | `tidy()` and `plot()` of age, period and cohort effects via `bamp::effects()` with the convention stated |
| Reproducibility | Seeds passed through; parallel RNG streams (`L'Ecuyer-CMRG`) for chains and strata |
| Messages | bamp output captured into the diagnostics table; nothing printed by default |
| Dependency | `bamp (>= 3.0.0)` in Suggests; tests skip if not installed; an integration test on a small fixture runs in CI |

### 7.2 Runtime budget

The review's test fit took 141 s for one stratum with short chains
(6,000 iterations, 4 chains on 4 cores), and about 16 minutes with bamp's
automatic chain length. Because the chains already use all cores, 10 strata
on the same machine take roughly 25 minutes (short chains) to about 2.7 hours
(automatic length) per fit. Backtests multiply this by the number of
origins, so the wrapper needs explicit `mcmc_control()` presets ("fast" for
backtests with an enforced convergence check, "final" for the reported
forecast) and caching of fits.

### 7.3 Not in the current scope (possible later)

A native engine (e.g. TMB with optional full MCMC via `tmbstan`; Kristensen
et al. 2016, Monnahan et al. 2018) would add Poisson likelihood, joint
multi-stratum APC models (Riebler et al. 2009; Riebler et al. 2011), missing
cells and fitting coarse counts directly through the aggregation operator.
These are recorded as options for later.

### 7.4 Known modelling risks

- APC extrapolation depends on how the linear trend is split between period
  and cohort; forecasts of RW2 models widen quickly. Riebler and Held (2017)
  discuss criticism that Bayesian APC intervals are too wide and report
  well-calibrated intervals in their application. Calibration must therefore
  be checked in mortensor's own backtests.
- bamp 3.0.0 is a recent release with a limited test suite; pin the version
  and keep mortensor-side integration tests.

## 8. Layers 4-6: forecasting, composition and derived quantities

### 8.1 Draws as the common currency

`mt_draws` stores simulated paths as a long `data.table` or array
\[age x period x strata x cause x draw\] with metadata (model, grid, seed).
Conventions: 1,000 draws by default, reproducible seeds, `posterior`-package
compatible summaries.

### 8.2 Composition functions

| Function | Purpose |
|---|---|
| `mt_ensemble(..., weights)` | Mixture of model draws. Weights: equal, inverse score, or stacking on out-of-sample scores from inner CV; never fitted on the evaluation holdout. Precedent for probabilistic model ensembles in mortality: Kontis et al. (2017) |
| `cohere()` | Coherence across strata (e.g. sexes, deprivation quintiles) via joint models (Li-Lee, product-ratio, Tucker3 across groups) or post-hoc constraints, documented as such |
| `reconcile()` | Make draws add up over a hierarchy (causes to all-cause, regions to nation) with bottom-up, top-down or MinT reconciliation (Wickramasuriya et al. 2019), applied draw by draw |
| `cause_shares()` | Compositional forecasts of cause shares given an all-cause forecast (successor of the Tucker4/ILR code) |
| `scenario()` | Explicit user adjustments (e.g. shock periods, target improvement rates), recorded in the forecast metadata |

### 8.3 Derived quantities (always pathwise over draws)

| Function | Output |
|---|---|
| `life_table()` | Interval-aware (single-year, abridged, arbitrary), open-age methods documented and tested |
| `life_expectancy(age = c(0, 30, 65))` | e_x with intervals from draws (fixes D10, D22) |
| `prob_dying(from, to)` | e.g. probability of dying between 30 and 70 |
| `deaths()` | Expected and simulated deaths given future exposures |
| `asr(standard)` | Age-standardised rates for a chosen standard population |
| `yll(reference)` | Years of life lost against a reference life table |
| `aggregate_to(grid)` | Any reporting grid via **C**, e.g. 5-year groups for publication |

## 9. Layer 7: evaluation

| Function | Purpose |
|---|---|
| `backtest(spec, data, origins, h)` | Rolling-origin evaluation that refits the whole pipeline (smoothing, tuning, weighting) at every origin |
| `score(bt, metrics)` | CRPS, WIS (Bracher et al. 2021 definition), log score, coverage and PIT by horizon, MASE against a declared baseline; via `scoringutils` |
| `compare()` | Benchmark tables with differences and their uncertainty |
| `calibration_plot()` | Coverage by horizon and PIT histograms |

## 10. Example of the target API (proposal, not implemented)

```r
library(mortensor)

d <- mt_data() |>
  add_deaths(ons_deaths, age = "age_group", period = "year",
             strata = c("sex", "imd"), cause = "cause") |>   # 5-year groups, 90+
  add_exposure(ons_pop, age = "age", period = "year",
               strata = c("sex", "imd"), type = "midyear") |> # single years
  set_grid(age = single_ages(30, 100, open = TRUE), period = single_years(2001, 2019)) |>
  flag_shocks(periods = 2020:2022, action = "exclude")

spec <- mt_ensemble(
  apc_bayes(age = "rw2", period = "rw2", cohort = "rw1", overdispersion = TRUE,
            by = c("sex", "imd"), mcmc = mcmc_control("fast")),   # bamp wrapper
  lee_carter(by = c("sex", "imd")),
  tucker3(transform = "coda"),
  weights = stacking(score = "crps", cv = rolling_origin(h = 5, n = 4))
)

bt <- backtest(spec, d, origins = 2010:2014, h = 5)
score(bt, metrics = c("crps", "coverage"))           # is the ensemble calibrated?

fit <- mt_fit(spec, d)
fc  <- mt_forecast(fit, h = 15, n_draws = 1000, exposure_future = ons_pop_projection)
fc  <- reconcile(fc, hierarchy = cause ~ all_cause, method = "mint")

summarise(fc, life_expectancy(age = 30), deaths(), asr(standard = "esp2013"),
          grid = age_axis(seq(30, 90, 5), open = TRUE))
```

## 11. What happens to the current code

| 0.1.0 component | Fate |
|---|---|
| `mort_data`, `mort_forecast`, `mort_calibration` | Replaced by `mt_data`, `mt_draws`, `mt_fit` (clean break, no converters; a migration note maps old calls to new ones) |
| Tucker3/CoDa and Tucker4/ILR cores | Ported behind the model contract after fixing the open age group (D11) and single-group handling (D12) |
| `calibrate()`, `compare_models()`, `calibrate_all_models()` | Replaced by `tune()`, `mt_ensemble()` and `backtest()` |
| Bootstrap PI functions | Replaced by model-specific `simulate()` paths |
| Life tables, disaggregation | Rewritten interval-aware; PCLM replaces spline disaggregation |
| `run_forecast()`, `run_reconciled_forecast()`, `run_coherent_forecast()` | Removed; examples rebuilt from composable steps |
| ONS loaders | Moved to `read_ons_*()` adapters |
| Plotting | Kept as ideas; reimplemented on tidy outputs, returning ggplot objects |

## 12. Engineering

- **Licence:** GPL-3 (decision 1). DESCRIPTION `License: GPL-3`; the MIT
  LICENSE file is removed. Only the copyright holder(s) of the existing code
  can relicense it; confirm that no external contributions under MIT need
  separate permission.
- **Dependencies (core):** `data.table`, `Matrix`, `rlang`, `cli`,
  `generics`, `ggplot2`, `forecast`, `rTensor`, `compositions`.
  **Suggests:** `bamp (>= 3.0.0)`, `demography`, `StMoMo`, `ungroup`,
  `MortalityLaws`, `scoringutils`, `posterior`, `future.apply`.
  (`ungroup` is used only as a reference in tests; the PCLM code is adapted
  into mortensor, see Section 5.1.)
- **One package** (decision 4), organised by layer in `R/` (`data_*.R`,
  `harmonise_*.R`, `model_*.R`, `draws_*.R`, `compose_*.R`, `derive_*.R`,
  `evaluate_*.R`, `adapters_*.R`).
- **Regression baseline from 0.1.0:** before porting, record the point
  forecasts of every 0.1.0 model on `example_mortality` and on the owner's
  real data with fixed hyperparameters. Ported models must reproduce them
  (within tolerance) except where a documented defect fix changes them (e.g.
  D11). This protects what worked well in 0.1.0 while the structure changes.
- **Testing:**
  - Unit tests per function.
  - Invariant tests, e.g. regrouping preserves totals, life-table identities
    hold, and draws aggregate consistently.
  - Posterior predictive and convergence checks for the bamp wrapper on a
    simulated fixture with known rates.
  - Regression tests against `demography`, `StMoMo` and `bamp` on fixtures.
  - A small public benchmark run in CI. Ship only data whose licence allows
    redistribution (check the terms of use of sources such as the HMD before
    including any extract); otherwise use simulated data.
- **CI:** R CMD check on release, devel, oldrel and the minimum R version;
  coverage; pkgdown site built from executable vignettes, with precomputed
  results for slow models.
- **Reproducible environment:** a `Dockerfile` based on a pinned
  `rocker/r-ver` image plus an `renv` lockfile for development and
  benchmarks.
- **API stability:** `lifecycle` badges; 0.2.0 is a deliberate breaking
  release; 1.0.0 after the benchmark (Phase 6).

## 13. Phased roadmap

Effort assumes one developer familiar with R package development; ranges
reflect uncertainty.

| Phase | Content | Deliverable | Effort |
|---|---|---|---|
| 0 | Stabilise 0.1.x (review Phase 0: D1, D2, D17, CI, executable vignettes); relicense to GPL-3; record the 0.1.0 regression baseline | 0.1.1 that installs, checks and runs | 1-2 weeks |
| 1 | Data layer and granularity: axes, strata, observation sets, operator **C**, validation, `regroup()`, `ungroup()` with the methods chosen in Section 5.1, exposures, interval-aware life tables; single-year default grid | `mt_data` API with tests | 4-6 weeks |
| 2 | Model contract, draws, `backtest()`/`score()`, benchmarks (RWD, ARIMA, ETS), Lee-Carter, Tucker3-CoDa and Tucker3-Direct ported | Honest benchmark harness; first comparison of Tucker3 against standard models | 4-6 weeks |
| 3 | Seamless bamp wrapper `apc_bayes()` (Section 7) | BAMP functionality available | 2-3 weeks |
| 4 | Composition and derived quantities for the release-1 models: ensembles (mixture, inverse score, stacking on inner CV), pathwise life tables, e_x, deaths with external population projections, ASR | **Release 0.2.0** (first release) | 4-6 weeks |
| 5 | Port every remaining model and method that runs in 0.1.0 (committed, see list below) | **Release 0.3.0** | 5-8 weeks |
| 6 | Validation: public benchmark on openly licensed data, executable vignettes, methods paper draft | **Release 1.0.0** | 4-6 weeks |

Total: roughly 5-8 months of focused work.

**Committed for 0.3.0 (all models and methods running in 0.1.0):**
FDM; StMoMo APC and CBD (with real ages and correct exposures); Tucker4/ILR
coherent multi-cause forecasting and cause shares; WLS reconciliation
(alongside MinT) in `reconcile()`; sex/group coherence in `cohere()` (as a
documented constraint, plus Li-Lee and product-ratio models); NNLS stacking
(fitted on the scale it is applied to, on inner CV); GAM (Gaussian) and
Whittaker smoothing as options next to the Poisson P-spline; GAM-Bayes
smoothing-uncertainty propagation as an option on top of model simulation;
the natural-Gompertz, Hyman-Kannisto and GAM disaggregation methods as
alternative `ungroup()` methods with the exposure-weighted constraint;
Kannisto, Gompertz and reciprocal open-age methods in `life_table()`.

## 14. Decisions (agreed with the owner, 2026-10-12)

| # | Topic | Decision |
|---|---|---|
| 1 | Licence | Change to GPL-3 |
| 2 | Bayesian APC | Wrapper around `bamp` only for now, but seamless (Section 7); native engine not planned |
| 3 | Population scope | External population (exposures and external projections); no cohort-component projection |
| 4 | Package structure | One package |
| 5 | Default latent grid | Single years of age |
| 6 | Backward compatibility | Clean break at 0.2.0 |
| 7 | Scope | Release 0.2.0: APC (bamp) + Lee-Carter + Tucker3 + benchmarks; release 0.3.0 definitively adds all models and methods running in 0.1.0 (Section 13) |

## 15. References

Checked against bibliographic records (Consensus search index, CRAN package
metadata) during this session; see also the reference list of the review.

- Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
  forecasts in an interval format. *PLoS Computational Biology* 17(2):
  e1008618. doi:10.1371/journal.pcbi.1008618
- Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting mortality
  rates. *Statistical Modelling* 4(4): 279-298.
- Hyndman RJ, et al. (2013). Coherent mortality forecasting: the product-ratio
  method with functional time series models. *Demography*.
- Knorr-Held L, Rainer E (2001). Projections of lung cancer mortality in West
  Germany: a case study in Bayesian prediction. *Biostatistics*.
- Kontis V, et al. (2017). Future life expectancy in 35 industrialised
  countries: projections with a Bayesian model ensemble. *The Lancet* 389:
  1323-1335.
- Kristensen K, et al. (2016). TMB: Automatic differentiation and Laplace
  approximation. *Journal of Statistical Software* 70(5).
  doi:10.18637/jss.v070.i05
- Li N, Lee R (2005). Coherent mortality forecasts for a group of populations:
  an extension of the Lee-Carter method. *Demography*.
- Monnahan CC, et al. (2018). No-U-turn sampling for fast Bayesian inference
  in ADMB and TMB: introducing the adnuts and tmbstan R packages. *PLoS ONE*
  13. doi:10.1371/journal.pone.0197954
- Riebler A, et al. (2009). The analysis of heterogeneous time trends in
  multivariate age-period-cohort models. *Biostatistics*.
- Riebler A, et al. (2011). Estimation and extrapolation of time trends in
  registry data: borrowing strength from related populations. *The Annals of
  Applied Statistics*.
- Riebler A, Held L (2017). Projecting the future burden of cancer: Bayesian
  age-period-cohort analysis with integrated nested Laplace approximations.
  *Biometrical Journal*.
- Rizzi S, Gampe J, Eilers PHC (2015). Efficient estimation of smooth
  distributions from coarsely grouped data. *American Journal of
  Epidemiology*. doi:10.1093/aje/kwv020
- Rizzi S, et al. (2018). How to estimate mortality trends from grouped vital
  statistics. *International Journal of Epidemiology*.
- Schmid VJ, Held L (2007). Bayesian age-period-cohort modeling and
  prediction - BAMP. *Journal of Statistical Software* 21(8).
  doi:10.18637/jss.v021.i08
- Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast
  reconciliation for hierarchical and grouped time series through trace
  minimization. *Journal of the American Statistical Association*.
  doi:10.1080/01621459.2018.1448825
