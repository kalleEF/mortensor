# mortensor revamp plan

Date: 2026-10-12. Status: proposal for discussion. Builds on the findings in
`dev/production-review.md` (referred to as "the review"; defect IDs D1-D26).

## 1. Goal

Turn mortensor into a suite of composable functions that:

1. accepts mortality and population data at **different granularities**
   (age groupings, open intervals, period lengths, strata, causes) and keeps
   track of what was observed at which resolution;
2. builds **complex forecasts** from interchangeable parts: data
   harmonisation, one or more models, ensembles, coherence and reconciliation
   constraints, and derived demographic quantities;
3. includes the **functionality of BAMP** (Bayesian age-period-cohort
   modelling and prediction; Schmid and Held 2007) as a first-class model
   family, and extends it to Poisson likelihoods, several strata and mixed
   granularity;
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
| Plug-in models | A small contract that every model implements; external engines (bamp, StMoMo, demography) are wrapped, not copied |
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
2. **Model through the operator:** count-likelihood models (Bayesian APC,
   Poisson P-spline, Poisson Lee-Carter) put **C** into the likelihood and fit
   the latent rates directly to the coarse data. This propagates the
   uncertainty of ungrouping into the forecast and is the generalisation of
   BAMP's unequal age/period grid (Knorr-Held and Rainer 2001).

### 4.3 Functions

| Function | Purpose |
|---|---|
| `mt_data()`, `add_deaths()`, `add_exposure()`, `add_population()`, `add_rates()` | Build a data object from long tables or arrays; declare columns for age, period, strata, cause |
| `set_grid()` | Choose the latent grid (default: finest resolution present, single years where any input has single years) |
| `validate()` | Coverage, duplicates, overlapping intervals, gaps, negative or zero exposures, deaths > exposure, consistency between deaths, exposure and rates |
| `describe_granularity()` | Table of what is observed at which resolution, by dataset and period |
| `as_exposure()` | Person-years from mid-year or January-1 populations with a declared approximation |
| `read_hmd()`, `read_ons_*()`, `read_eurostat()` | Thin adapters into `mt_data`; source-specific code stays out of the core |

## 5. Layer 2: harmonisation

| Function | Method | Notes |
|---|---|---|
| `regroup()` | Exact summation of deaths and exposures through **C** | Rates are recomputed from summed counts (exposure-weighted by construction; fixes D21) |
| `ungroup()` | Penalized composite link model, 1D or 2D (Rizzi et al. 2015, 2018) | Wrap the `ungroup` package first; native implementation later if needed. Returns draws or bootstrap replicates of the ungrouped counts so the uncertainty can be propagated |
| `align()` | Bring several observation sets to the latent grid (regroup where finer, ungroup where coarser) | Records the operation per cell |
| `smooth()` | Poisson P-spline in age and period with log-exposure offset (Currie et al. 2004) | Applied to training data only; raw data are kept |
| `flag_shocks()` | Mark periods (e.g. 2020-2022) for exclusion, down-weighting or an indicator term | Explicit, model-aware, never silent |
| `impute_missing()` | Model-based imputation or exclusion | Models declare whether they accept missing cells |

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
| Bayesian APC (BAMP functionality, Section 7) | `apc_bayes()` | `bamp` wrapper, then native engine | native engine | yes (native) |
| Lee-Carter and variants | `lee_carter()` | Poisson log-bilinear fit natively; `demography` wrapper for comparison | via `li_lee()` | later |
| Functional demographic model | `fdm()` | `demography` wrapper | via `product_ratio()` | no |
| Multi-population coherent | `li_lee()`, `product_ratio()` | Li and Lee (2005); Hyndman et al. (2013) | yes | no |
| Generalised APC family | `stmomo()` | `StMoMo` wrapper with real ages and correct exposures | no | no |
| 2D P-spline smoothing and forecasting | `pspline2d()` | Currie et al. (2004), Poisson | no | yes |
| Tensor models (current core idea) | `tucker3(transform = c("coda", "log"))`, `tucker4()` | Port from 0.1.0 after fixing D11, D12 | yes | no |
| Benchmarks | `rwd()`, `arima_series()`, `ets_series()` | `forecast` package, simulated with `simulate()` | no | no |

## 7. Bayesian APC module (BAMP functionality)

### 7.1 Feature parity with bamp 3.0.0 and extensions

| bamp feature | Plan in mortensor |
|---|---|
| Binomial likelihood, logit link | `family = binomial()`; plus `family = poisson()` with log-exposure offset as the default for central death rates; optional negative binomial |
| RW1/RW2 smoothing priors for age, period, cohort | Same, as sparse precision matrices; optional sum-to-zero and drift constraints; Sørbye-Rue scaling option |
| Unstructured heterogeneity per effect (`rw2+het`) | Same, including prediction (broken in bamp 3.0.0) |
| Cell-level overdispersion | Same (Gaussian term in the linear predictor) |
| Unequal age and period grids (`periods_per_agegroup`) | Generalised through the latent grid and **C**: any combination of age and period resolutions, varying open intervals |
| Period and cohort covariates | Same, plus age-varying covariate effects |
| MCMC, Gelman-Rubin checks, DIC | Laplace approximation by default, full MCMC on request; R-hat and effective sample size reported and enforced; DIC and WAIC for parity, rolling-origin CV recommended for forecasting |
| `selectModel()` greedy DIC search | `tune()` over RW orders and options, by CV or DIC |
| `predict_apc()` random-walk extrapolation | `mt_forecast()` returning draws of rates and, with future exposures, of deaths including Poisson/binomial noise |
| `simulate_apc()` | `simulate()` method, used for simulation-based calibration tests |
| `effects()`, `plot()` with identifiability conventions | `tidy()` and `plot()` of age, period and cohort effects with documented conventions |
| (not in bamp) several strata | Multivariate APC with shared and stratum-specific effects (Riebler et al. 2009; Riebler et al. 2011), the natural model for sex x deprivation data |
| (not in bamp) missing cells | Supported: missing cells drop out of the likelihood |
| (not in bamp) mixed granularity | Through **C** (Section 4.2) |

### 7.2 Engine options

| Option | Pros | Cons | Use |
|---|---|---|---|
| A. Wrap `bamp` (Suggests) | Parity within weeks; reference implementation | One stratum per call, binomial only, no missing data, GPL-3 boundary, young release | Stage 1 and as validation reference |
| B. Native engine on TMB (Kristensen et al. 2016) | Fast Laplace approximation for sparse latent Gaussian models; Poisson/binomial; **C** and multiple strata straightforward; CRAN-compatible; MIT-compatible when written from the papers | C++ templates to maintain; Laplace approximation treats hyperparameters by empirical Bayes unless combined with MCMC | Stage 2 default engine |
| C. Full MCMC via `tmbstan` (Monnahan et al. 2018) on the same TMB template | Full Bayes on the identical model; checks the Laplace approximation | Slower | Optional `engine = "nuts"` |
| D. INLA (as in BAPC; Riebler and Held 2017) | Fast and accurate for this model class | Not on CRAN | Not as a hard dependency; optional comparison only |

**Recommendation:** A first, then B with C as an option, validated against A
on shared examples (posterior predictive agreement of rates) and by
simulation-based calibration. The native engine must be written from the
published model descriptions, not from bamp's GPL-3 source (see decision 1).

### 7.3 Known modelling risks

- APC extrapolation depends on how the linear trend is split between period
  and cohort; forecasts of RW2 models widen quickly. Riebler and Held (2017)
  discuss criticism that Bayesian APC intervals are too wide and report
  well-calibrated intervals in their application. Calibration must therefore
  be checked in mortensor's own backtests, not assumed.
- The bamp test run in the review took 141 s for one group with short chains
  (about 16 minutes with automatic chain length). Joint models over 10 strata
  and single-year ages need the sparse TMB engine to be practical.

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
| `cohere()` | Coherence across strata (e.g. sexes, deprivation quintiles) via joint models (Li-Lee, product-ratio, multivariate APC) or post-hoc constraints, documented as such |
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
  apc_bayes(age = rw2(), period = rw2(), cohort = rw1(), overdispersion = TRUE,
            family = poisson(), joint = TRUE, engine = "tmb"),
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
| `mort_data`, `mort_forecast`, `mort_calibration` | Replaced by `mt_data`, `mt_draws`, `mt_fit`; converters provided |
| Tucker3/CoDa and Tucker4/ILR cores | Ported behind the model contract after fixing the open age group (D11) and single-group handling (D12) |
| `calibrate()`, `compare_models()`, `calibrate_all_models()` | Replaced by `tune()`, `mt_ensemble()` and `backtest()` |
| Bootstrap PI functions | Replaced by model-specific `simulate()` paths |
| Life tables, disaggregation | Rewritten interval-aware; PCLM replaces spline disaggregation |
| `run_forecast()`, `run_reconciled_forecast()`, `run_coherent_forecast()` | Removed; examples rebuilt from composable steps |
| ONS loaders | Moved to `read_ons_*()` adapters |
| Plotting | Kept as ideas; reimplemented on tidy outputs, returning ggplot objects |

## 12. Engineering

- **Dependencies (core):** `data.table`, `Matrix`, `rlang`, `cli`,
  `generics`, `ggplot2`. **Suggests:** `bamp`, `demography`, `StMoMo`,
  `ungroup`, `scoringutils`, `posterior`, `tmbstan`, `future.apply`.
  **LinkingTo:** `TMB`, `RcppEigen` (Stage 2 engine).
- **Package split (decision 4):** keep one package until the TMB engine
  lands; split into `mortensor` (core) and `mortensor.apc` (compiled engine)
  if compile times or CRAN checks become a burden.
- **Testing:**
  - Unit tests per function.
  - Invariant tests, e.g. regrouping preserves totals, life-table identities
    hold, and draws aggregate consistently.
  - Simulation-based calibration for the Bayesian engine.
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
reflect uncertainty, especially for the native engine.

| Phase | Content | Deliverable | Effort |
|---|---|---|---|
| 0 | Stabilise 0.1.x (review Phase 0: D1, D2, D17, CI, executable vignettes) | 0.1.1 that installs, checks and runs | 1-2 weeks |
| 1 | Data layer and granularity: axes, strata, observation sets, operator **C**, validation, `regroup()`, `ungroup()` (wrapper), exposures, interval-aware life tables | `mt_data` API with tests | 4-6 weeks |
| 2 | Model contract, draws, `backtest()`/`score()`, benchmarks (RWD, ARIMA, ETS), Lee-Carter (Poisson), Tucker3 ported | Honest benchmark harness; first comparison of Tucker3 against standard models | 4-6 weeks |
| 3 | Bayesian APC Stage 1: `apc_bayes(engine = "bamp")` wrapper with convergence enforcement, draws, future exposures | BAMP functionality available | 2 weeks |
| 4 | Bayesian APC Stage 2: native TMB engine (Poisson/binomial, RW1/RW2, heterogeneity, overdispersion, **C**, missing cells), optional NUTS; validation against bamp and by simulation-based calibration | Engine handling mixed granularity | 6-10 weeks |
| 5 | Composition and derived quantities: ensembles (mixture, stacking), Li-Lee/product-ratio, multivariate APC, reconciliation on draws, cause shares, pathwise summaries | Full forecasting workflow | 5-8 weeks |
| 6 | Validation and release: public benchmark on openly licensed data, executable vignettes, methods paper draft, 1.0.0 | Release | 4-6 weeks |

Total: roughly 6-10 months of focused work. Phases 1-2 can start in
parallel with Phase 3 once the model contract is fixed.

## 14. Decisions needed from the owner

| # | Decision | Options | Recommendation |
|---|---|---|---|
| 1 | Licence | Stay MIT (wrap bamp via Suggests, write the native engine from the papers) or move to GPL-3 (allows reusing bamp code with attribution) | Stay MIT unless reusing bamp's C code saves substantial time; this is not legal advice |
| 2 | Native Bayesian engine | TMB (+ tmbstan), Stan only, INLA, or bamp only | TMB with optional tmbstan |
| 3 | Population scope | Exposures and external projections only, or a full cohort-component population projection (fertility, migration) | Exposures and external projections; cohort-component later if needed |
| 4 | Package structure | One package or core + compiled engine package | One package until Phase 4 |
| 5 | Default latent grid | Single years (needs ungrouping) or finest observed resolution | Finest observed resolution; single years when any input has them |
| 6 | Backward compatibility | Converters and deprecation period, or clean break | Clean break at 0.2.0 with a migration note, since 0.1.0 has no stable release |
| 7 | Scope cut for the first release | All model families, or APC + Lee-Carter + Tucker3 + benchmarks | The reduced set; add FDM, StMoMo and multi-cause after the benchmark |

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
