# Applied review of mortensor (reviewer B)

Reviewer profile: demographer and biostatistician; user of demography, StMoMo
and bamp through their documented interfaces; no package development
background. Goal: regional mortality forecasts (Laender by sex, with national
totals, plus causes of death at older ages) for a peer-reviewed paper, with
model choice justified by a backtest and checked against demography/StMoMo.

Material used: README, the vignettes (`mortensor`, `causes`, `faq`,
`forecasting-hmd`, `granularity-hmd`, `bayesian-apc.Rmd.orig`,
`regions.Rmd.orig`), the help pages and `_pkgdown.yml`. All code was run with
`devtools::load_all()` on 13 October 2026, with at most 200 draws and two
origins. The R sources changed while I worked (a warning about forecast rates
above 1 and a ridge fallback in `lee_carter()` appeared between my runs);
results below are from the state at the time of each run, and the P1 findings
were re-run on the later state.

## 1. Summary verdict

mortensor has a clear and attractive design: data at their own resolution,
preprocessing refitted inside leakage-free backtests, forecasts as joint paths,
and derived quantities computed path by path. The Lee-Carter core reproduces
StMoMo and demography. For my paper, however, the package is not yet usable
without reading source code, and two problems would produce wrong or absurd
published numbers without an error: (a) a crossed region by cause hierarchy
gives national and regional all-cause rates that are off by a factor equal to
the number of causes, and (b) the random walk with drift, and the ensembles
that the vignettes recommend, generate life expectancy paths far below any
plausible value in small Laender (for women in Bremen in 2035 under `rwd()`: a
median of 61 years and a 2.5% quantile of 0.5 years). Methods reporting (equations, defaults, sources of uncertainty,
parameters, residuals, MCMC diagnostics) also needs a dedicated reference.

## 2. What worked well

* Reading `de_deaths` and `de_population` with `mt_data()` worked first time,
  including German labels such as "1 Jahr bis unter 15 Jahre", deaths in
  mixed widths (0, 1-14, 5-year groups, 90+) and single-year populations
  shared by all causes.
* `mt_as_exposure(type = "dec31")` states its approximation and records it.
* `mt_fill_suppressed()` matched all published totals (`max_gap` 8e-7, 33,225
  cells filled) and marks the estimates with provenance `"estimated"`.
* `mt_backtest()` refits preprocessing and model on each training window and
  scores at the resolution of the observed deaths. Five specifications for 57
  strata and two origins ran in 5-9 seconds each.
* Paths are joint across ages, years and strata, so e0 gaps and rate ratios
  carry their own uncertainty; `mt_ensemble()` mixes whole paths.
* Cross-check: `lee_carter()` on `StMoMo::EWMaleData` (ages 50-89, 1971-2011)
  agreed with `StMoMo::fit(lc())` to 8e-6 in a_x, 2e-7 in b_x and 8e-5 in
  k_t; forecast medians in 2031 were within 0.3% of StMoMo's central
  forecast; `lee_carter(method = "svd")` matched `demography::lca(adjust =
  "dt")` exactly. The 95% interval at age 70 was slightly wider in mortensor
  (0.00927-0.01573 against 0.00968-0.01502), consistent with the documented
  drift uncertainty.
* The vignettes are candid about what a backtest cannot settle (gaps under
  Li-Lee versus separate models, the dementia trend, the bias of RW1 priors
  in `apc_bayes()`); most references are correct and complete.

## 3. Walkthrough log

**Step 1, data.** Followed `regions.Rmd.orig`: `mt_data()`,
`mt_as_exposure(type = "dec31")`, `mt_subset(periods = c(1998, 2019))`,
`mt_regroup(age = "deaths")`, `mt_fill_suppressed(x, mt_hierarchy(region =
"Deutschland", cause = "all"))` (55 s). The vignette then says "All totals are
matched", but the very next call `mt_validate(x, hierarchy = h_fill)` returns
"34 finding(s): 32 error": every Land's all-cause exposure is flagged with
"total / sum(parts): median 0.056" (1/18), because causes share one exposure
while the hierarchy is additive on deaths.

**Step 2, suppressed cells.** The published counts in `de_deaths` contain no
1s or 2s (values 0, 3, 4, 5, ...), so suppressed cells are 1 or 2. The raked
estimates range from 0.60 to 2.75: the method ignores these known bounds.

**Step 3, hierarchy.** `mt_hierarchy(region = list(total = "Deutschland",
groups = list(West = ..., East = ...)), sex = "both")` and
`mt_aggregate_strata()` gave 57 strata and a clean `mt_validate()`. To label
bottom and aggregate strata in backtest output I had to copy the vignette's
`mortensor:::hierarchy_structure()`, an internal function.

**Step 4, backtest** (origins 2011 and 2014, h = 5, 200 draws), `mt_score(bt,
by = "spec")`:

| Specification | CRPS | cov95 | mean PIT |
|---|---|---|---|
| `lee_carter()` separate | 0.0857 | 0.77 | 0.67 |
| `mt_reconciled(lee_carter(), h, "wls")` | 0.0861 | 0.74 | 0.67 |
| `mt_reconciled(li_lee(common_by = "sex"), h, "bottom_up")` | 0.0859 | 0.78 | 0.69 |
| `mt_reconciled(rwd(), h, "bottom_up")` | 0.0966 | 0.88 | 0.63 |
| ensemble Li-Lee + `rwd()`, bottom-up | 0.0838 | 0.91 | 0.66 |

On these scores the ensemble wins, as in the vignette. The three single
models differ in the third decimal, and the package gives no standard error or
paired test, so I could not say whether the differences are real.

**Step 5, forecast and e0.** `mt_forecast(mt_fit(ensemble), to = 2035,
n_draws = 200, seed = 1)` and `life_expectancy(fc, at = 0)`:
Deutschland.female e0 in 2035 had mean 82.0, median 84.7 and 2.5% quantile
43.3. One path had a death rate of 1.47 at ages 20-24 for Germany, driven by
Bremen women aged 20-24 (151 deaths per person-year). With `rwd()` alone,
Bremen.female e0 2035 had median 61.2 and 2.5% quantile 0.5; Germany female
2.5% quantile 17.3; the largest simulated rate was 162,061. With
`rwd(observation_noise = TRUE)` the Bremen 2.5% quantile was still 21.5.
`li_lee(common_by = "sex")` gave sensible results (Bremen 85.5, 82.3-87.9).
Even at h = 5 (origin 2014), the ensemble's 95% interval for Bremen female e0
in 2019 started at 57.9. In the later code state, `mt_forecast()` warns
"forecast death rates above 1 in 34.5% of the draws" for `rwd()` and suggests
"an ensemble" as a remedy, the very construct that propagated the problem.
The CRPS of log rates averaged over cells did not reveal any of this.

**Step 6, national e0 calibration.** By hand (the backtest cannot score e0),
from origin 2014 the observed Deutschland.both e0 in 2019 (81.20) lay below
the 95% intervals of `rwd()` (81.35-82.29) and of the ensemble (81.51-82.48),
and at the edge of Li-Lee's (81.22-82.94). Mean PIT of 0.63-0.79 in all
backtests reflects the post-2011 slowdown; the package has no fitting-period
choice (as in Booth, Maindonald and Smith) to test against it.

**Step 7, parameters and diagnostics.** `lc_parameters()` works for
`lee_carter()`; for `li_lee()` it stops with "Not a Lee-Carter fit." There is
no accessor for B_x, K_t, b_x,i, k_t,i or the AR(1) coefficients; I found
them in `fit$fits$commons` and `fit$fits$strata` by reading
`R/model-multipop.R`. `fitted()`, `residuals()` and `coef()` return `NULL`,
`logLik()` and `predict()` have no method, `summary(fit)` prints the list
structure. For `apc_bayes()` the diagnostics show only "ok"; R-hat values and
the bamp hyperpriors are not reported.

**Step 8, causes at 65+ with national totals.** After `mt_group_levels(xc,
"cause", de_cause_groups, other = "rest", total = "all", additive = "rates")`
and `mt_hierarchy(region = "Deutschland", cause = "all")`:

* `mt_aggregate_strata()` gives Deutschland.female.all an exposure of
  14,993,856 for ages 80-84 in 2019, eight times the correct 1,874,232
  (exposures summed over the eight cause groups), so its rate is 0.0059
  instead of 0.047.
* Bottom-up forecasts with the default `additive = "deaths"` give
  Bayern.female.all 0.0061 at 80-84 in 2024 (observed 2019: 0.046).
* With `additive = "rates"`, regional rates are summed: Deutschland.female
  neoplasms 0.168 (observed 0.011) and all causes 0.81.
* `mt_hierarchy(..., additive = c("deaths", "rates"))` is silently accepted
  and uses `"deaths"`.

The regional vignette avoids this by never crossing regions with causes, so my
target design (national cause-specific totals coherent with the Laender)
cannot be built. A smaller trap: `de_population` has no national rows, so the
Deutschland strata in `de_deaths` get no exposures; `mt_data()` accepts this
silently and `mt_fit()` later fails with "Deaths and exposures must be
complete on the grid ... can be estimated with mt_fill_suppressed()", which
points the wrong way.

**Step 9, causes.** Li-Lee with `common_by = c("sex", "cause")` and
`cause_shares()` ran quickly at 65+ (origin 2014); all-cause 95% coverage was
0.69-0.76.

**Step 10, reproducibility.** `mt_forecast(seed = 1)` is reproducible. Two
identical `mt_backtest()` calls (default `seed = 1`) differ in the `pit`
column (mean relative difference 0.03) because the randomised PIT uses
`runif()` outside the seeded block. With `seed = NULL` the seed is not
recorded in the draws. The default `n_origins = 5` produced 13 origins with
h = 10 (1999-2011) on `mt_example(years = 1990:2012)`, with 1,846 scored
cells at horizon 1 and 568 at horizon 10.

## 4. Prioritised recommendations

### P1: blocking for publication-quality use

**P1.1 Mixed additivity in crossed hierarchies.** Problem and evidence: Step
8 and Step 1. Change: let `mt_hierarchy()` take additivity per variable, for
example `mt_hierarchy(region = "Deutschland", cause = list(total = "all",
additive = "rates"))`, so that deaths sum over both, exposures sum over
regions only, and reconciliation projects deaths over regions and rates over
causes. Until then, `mt_aggregate_strata()`, `mt_reconcile()` and
`mt_validate()` should stop with an error when a hierarchy variable shares
exposures. Add a test with region by cause, and a section of the regional
vignette that builds national cause-specific forecasts.

**P1.2 Explosive random-walk paths in sparse cells.** Problem and evidence:
Step 5. Change: (a) treat `rwd()` on log rates as unsuitable for cells with
few deaths: estimate the drift on pooled or smoothed rates, or make
`observation_noise = TRUE` the default and shrink its drift; (b) let the
plausibility check cover derived quantities (for example e0 below the
observed minimum minus 10 years) and stop recommending "an ensemble" in that
warning when a component is the cause; (c) report per-stratum e0 intervals
for the smallest Laender in the regional vignette before recommending
`mt_ensemble(ll = li_lee(), rw = rwd())`.

**P1.3 A methods reference.** Problem: what I must write in a methods section
is scattered over help pages or only in the code. Examples: Li-Lee common
factor is the SVD Lee-Carter of pooled rates with second-stage k_t; stratum
factors come from an SVD of log residuals with zero deaths set to 0.5; AR(1)
without intercept, phi bounded to [-0.99, 0.99], no uncertainty in phi, a_x
or b_x; stratum innovations independent across strata; extra-Poisson noise
estimated from in-sample residual variance minus the Poisson part;
reconciliation and aggregate rates use exposures of the last training year;
bamp hyperpriors are bamp's defaults. Change: a vignette "Model reference"
with, per constructor, the equation, the estimator, the time-series model of
each index, the defaults, a table of included uncertainty sources (index
innovations, drift, parameters, Poisson and extra-Poisson noise, ungrouping)
and references. Add `summary.mt_fit()` that prints these settings for the
actual fit.

**P1.4 Parameters, residuals and convergence.** Problem: Step 7. Change: a
generic `mt_parameters(fit)` for every model (replacing `lc_parameters()` and
`apc_effects()`), `fitted()`, `residuals()` (deviance or Pearson) and
`coef()` methods, `mt_plot_residuals()` (age-period heat map, as StMoMo's
`plot(residuals())`), and R-hat and effective sample sizes in the
`apc_bayes()` diagnostics, with the hyperpriors exposed as arguments.

### P2: important

**P2.1 Backtest targets and score uncertainty.** Evidence: Steps 4-6. Change:
`mt_backtest(target = c("log_rate", "e0", "asr", "deaths"))`, weights (by
deaths or exposure) and an age range in `mt_score()`, standard errors over
origins and strata, and `mt_compare(bt, reference = )` with paired
differences. Score cells with provenance `"estimated"` only on request.

**P2.2 Backtest defaults and reproducibility.** Make `n_origins` mean the
number of origins and default to origins with complete horizons; draw the PIT
uniforms inside the seeded block; store the RNG state when `seed = NULL`.

**P2.3 Bounded estimation of suppressed counts.** Add `lower` and `upper` to
`mt_fill_suppressed()` (here 1 and 2) and an option to return several
imputations, so that backtests and forecasts can carry this uncertainty.

**P2.4 Data checks.** `mt_data()` should warn when strata have deaths and no
exposure source, and the error in `mt_fit()` should name those strata.
`?de_deaths` should state that `de_population` has no national total.

**P2.5 Path-level helpers.** `life_expectancy()` lacks `paths = TRUE` (which
`asr()` and `prob_dying()` have); `as.data.table(fc)` returns the summary,
not draws; contrasts need array indexing (`e[, st$label == ..., ]`). Add
`paths = TRUE` everywhere, `as.data.table(fc, draws = TRUE)` and
`mt_contrast(x, numerator, denominator, type = c("difference", "ratio"))`.

**P2.6 Export hierarchy levels.** Add a `level` column (bottom, intermediate,
total) to `mt_strata()` and to backtest output, and drop `:::` from the
vignette.

**P2.7 Remedies for trend change.** Add a fitting-period step
(`step_window(from = )` or `lee_carter(fit_from = )`), tunable with
`mt_tuned()`, and correlated stratum innovations in `li_lee()` (or document
that `dependence = "auto"` gives `"none"` although regional deviations are
independent).

**P2.8 Model-choice guidance.** Replace the four bullets of the FAQ with a
decision table: data size per cell, number of strata, need for coherence,
cohort effects, horizon, and the model and preprocessing to try first.

### P3: nice to have

* Remove user-facing references to "mortensor 0.1" (in `?apc_select`,
  `?cause_shares`, `?mt_bootstrap`, `?mt_cohere`, `?mt_disaggregate`,
  `?mt_ensemble`, `?mt_smooth`, `?tucker3`, `?mt_ungroup`) and to
  `dev/calibration.md` and "the validation pipeline", which users cannot see;
  give the evidence in a vignette.
* README lists a vignette "Comparison with demography, bamp and StMoMo" that
  does not exist; "Where next" in the getting-started vignette uses old
  titles; the shipped `regions.Rmd` is still the simulated version.
* FAQ on exposures omits `type = "dec31"`; `?mt_ensemble` omits `"stack"`
  from the `weights` argument; `?stmomo` gives the authors as "Villegas,
  Millossovich and Kaishev" (the paper lists Villegas, Kaishev and
  Millossovich).
* `cause_shares()` ("Tucker4-ILR") and `tucker3()` have no published
  application to mortality; label them as experimental.
* Export to `demogdata` and `StMoMoData` for cross-checks, and document
  where the underlying `fitStMoMo`, `fdm` and bamp objects sit.
* Life table: offer an a0 rule (Coale-Demeny or Andreev-Kingkade); align
  `omega = 110` with the Gompertz integration to age 130; add the WHO world
  standard population.

## 5. Proposed structure of functions, objects and documentation

The object chain `mt_data` -> `mt_prep` -> `mt_spec` -> `mt_fit` ->
`mt_draws` -> summaries, with `mt_backtest` -> `mt_score`, is sound and easy
to explain. The friction lies in naming and in where the boundaries sit.

**Naming.**

* `mt_reconcile()` (draws) and `mt_reconciled()` (specification) differ by one
  letter. The convention says model constructors are nouns, yet
  `mt_reconciled()`, `mt_bootstrap()` and `mt_tuned()` are wrappers with
  verb-like `mt_` names. Use `reconciled()`, `bootstrapped()`, `tuned()` and
  `ensemble()` for wrappers, and keep `mt_` verbs for operations.
* "Coherent" has three meanings: non-divergent (Li-Lee, product-ratio),
  additive over a hierarchy, and ordered (`mt_cohere()`). Rename
  `mt_cohere()` to `mt_order_strata()` and reserve "coherent" for
  non-divergence; use "consistent" or "additive" for hierarchies.
* Four aggregation verbs do different things: `mt_regroup()` (ages and
  periods), `mt_group_levels()` (levels of one variable), `mt_aggregate_strata()`
  (adds hierarchy nodes) and `mt_aggregate_draws()`. Keep them, but add a
  table in `?mt_hierarchy` and the FAQ saying which to use when.
* `mt_window()` duplicates `mt_subset(periods = )`; keep one.

**Tailoring.** Most functions are general. Points tied to a specific source
or analysis: `mt_read_genesis()` hard-codes the 31 December exposure and
creates a `"rest"` cause, both analysis decisions that belong to
`mt_as_exposure()` and `mt_group_levels()`; the `mt_validate()` message
"total equals twice the sum of parts" targets one data quirk;
`de_cause_groups` and `ons_cause_groups` are vignette choices exported as data
(fine, if documented as examples). Move the GENESIS, GBE, HMD, StMoMo and
demography readers into one "Adapters" group.

**Reference index (`_pkgdown.yml`).** Suggested groups: Read data (`mt_data`,
`mt_add*`); Adapters; Prepare and check (`mt_subset`, `mt_group_levels`,
`mt_as_exposure`, `mt_fill_suppressed`, `mt_validate`); Harmonise
(`mt_regroup`, `mt_ungroup`, `mt_smooth`, `mt_prep` and steps); Models,
single population; Models, several populations; Bayesian models; Model
wrappers (ensemble, reconciled, bootstrap, tuned); Hierarchies; Forecast and
draws; Inspect fits (new `mt_parameters`, `residuals`, `summary`); Evaluate;
Derived quantities; Plots; Low-level and extending (`new_mt_spec`,
`mt_operator`, axis helpers), which currently sit among user functions.

**Vignettes.** Keep the current set and add (1) "Model reference" (P1.3),
(2) "Choosing a model" (P2.8) and (3) a finished regional vignette that crosses
regions with causes, shows e0 intervals for Bremen and Saarland, and backtests
e0 as well as log rates.
