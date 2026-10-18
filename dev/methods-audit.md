# Audit: published methods and package-specific solutions (2026-10-18)

The owner's principle: mortensor implements published methods; it does not
develop new ones. Four read-only reviews (data and harmonisation, models,
wrappers and coherence, evaluation and derived quantities) classified every
methodological component as published (P), published with deviations (PD),
package-specific engineering (T-eng: input checks, numerical floors) or
package-specific method that changes results (T-method). This file records
the main findings and the decisions taken (work package E3).

## Main findings

| Area | Component | Class | Default? | Finding |
|---|---|---|---|---|
| Models | Extra-Poisson layer of predictive draws (`extra_noise_sd()`) | T-method | yes, every model of smooth rates | Moment estimator without source; double counts noise for `tucker3()` and the fdm models; replaces the posterior predictive of the Bayesian models; the validation generator added the same noise, so coverage checks confirmed the package's own assumption |
| Models | `tucker3(residual_noise = TRUE)` | T-method | yes | Counts observation noise twice; the default was chosen because coverage on the package's generator was 0.87 against a threshold of 0.90 |
| Models | `tucker3()` | T-method | opt-in | Closest published model (three-way Lee-Carter of Russolillo, Giordano and Haberman 2011) not cited or followed (no centring) |
| Models | `cause_shares()`, `product_ratio()`, `rwd(observation_noise = TRUE)`, `tucker3(transform = "coda")` | T-method | opt-in | No published source; published counterparts exist in the package (`reconciled()`, `coda_causes()`, `coherent_fdm()`, ARIMA(0,1,1) with drift) |
| Models | Li-Lee AR(1) without intercept clamped to 0.99; Lee-Carter ridge and SVD fallbacks | PD / T-eng | yes | Silent; not documented in the help pages |
| Coherence | Dependence by an "empirical copula" | T-method | yes | Presented as the Schaake shuffle, but one template value per stratum and year; needed to pass `regions_coherence`, whose generator has exactly that dependence |
| Coherence | WLS and MinT weights from the base draws | PD | yes | Published versions use in-sample residual (co)variances; not stated |
| Coherence | Non-negativity by holding negatives at zero and projecting again | T-method | yes | Described as a form of Wickramasuriya et al. (2020), which solves a quadratic program |
| Wrappers | `recalibrated(scale = "horizon")` power law, drift correction, age bands, targets | T-method | power law default | Closest published method (member-by-member post-processing) not cited |
| Wrappers | `anchored()`: 3-year offsets, 10-year linear fade, own between-stratum variance estimator | T-method | yes (of the wrapper) | Only `years = 1, decay = Inf, shrink = FALSE` is the published Lee-Miller jump-off |
| Wrappers | `ensemble(weights = "crps"/"stack")`, `bootstrapped(type = "grid")` | T-method / PD | opt-in / default | Inverse-CRPS weights unsourced; point stacking used as mixture weights; bootstrap attributed to the wrong paper |
| Harmonisation | `mt_ungroup(preserve_totals = TRUE)` | T-method | yes | Rescaling after the PCLM; steps at group boundaries; in 2D re-imposes noise |
| Harmonisation | Open interval inside the PCLM (stationary-population weights, 3 fixed iterations) | T-method | yes | Validated only on data generated with the same assumption |
| Harmonisation | `hyman_kannisto`, `natural_gompertz`, `gam` ungrouping; `mt_shape_rates(method = "gam")` | T-method | `hyman_kannisto` default of `mt_disaggregate()` | No source; undocumented constants and a hidden switch |
| Harmonisation | Whittaker over time with weights rescaled per age and pooled GCV; Gaussian 1D GAM | T-method / PD | opt-in | Deviate from the weighted Whittaker smoother and the Poisson GAM |
| Harmonisation | Suppressed cells: clipping inside IPF; clipped pooled means | T-method | opt-in | Truncated-Poisson mean is the standard alternative |
| Evaluation | Diebold-Mariano test with few overlapping origins | PD | yes | Size 10-27% at nominal 5% (7 origins, h 2-5); lag-0 fallback rejects a true null about half the time |
| Derived | ASR with standard weights split evenly over finer ages | T-method | yes | 2-7% higher than ESP 2013 standardisation |
| Derived | `yll()` with e(x) at the start of the interval | PD | yes | 10-12% too high with 5-year groups |
| Derived | `mt_decompose(composition = "pooled")` | T-method | required choice | Presented with a reference that only supplies the generic line integral |
| Derived | Excess: estimand, P-score NA with a zero draw, `p_positive` | PD | yes | Under-specified |

## Decisions

Defaults become the published methods. Package-specific options stay where
they are useful, labelled "package-specific (reference unknown)" in the help
pages and in the article "Methods and sources". Package inventions with a
published counterpart in the package are deprecated: they keep working,
warn (class `mt_warning_deprecated`) and name the replacement.

1. Predictive draws: Poisson deaths given the model's rates (Brouhns, Denuit
   and Vermunt 2002); the extra-Poisson layer is removed;
   `stanmomo(family = "nb")` uses its posterior negative-binomial
   predictive. Validation generator: Lee-Carter with Poisson deaths only.
2. `tucker3()`: centring as in Russolillo et al. (2011), `residual_noise =
   FALSE`; `transform = "coda"` deprecated.
3. Deprecated: `product_ratio()`, `cause_shares()`, `rwd(observation_noise
   = TRUE)`, ungrouping methods `hyman_kannisto`, `natural_gompertz`,
   `gam`, `mt_shape_rates(method = "gam")`.
4. Reconciliation: `dependence = "none"` by default; the copula is opt-in
   and labelled; WLS/MinT, non-negativity, fixed nodes and proportional
   documented with their published sources and deviations.
5. `recalibrated(scale = "constant")`; `anchored(years = 1, decay = Inf,
   shrink = FALSE)` (Lee and Miller 2001) with DerSimonian-Laird shrinkage
   when chosen; ensemble weights and bootstrap attributions corrected.
6. `mt_ungroup(preserve_totals = FALSE)`; open-interval treatment documented
   and iterated to convergence; `mt_disaggregate()` default changed;
   standard weighted Whittaker smoother and quasi-Poisson GAM; truncated
   Poisson mean for suppressed cells.
7. `mt_dm_test(varestimator = "ewc")`: fixed-smoothing test with the
   equal-weighted cosine variance (Lazarus et al. 2018; Coroneo and Iacone
   2020); simulated size 0.048-0.066 at nominal 5% for 7-16 origins and
   h up to 5.
8. ASR by direct standardisation on the standard's age groups (Eurostat
   2013); `yll()` with e(x + a_x) (Murray 1994); "pooled" composition and
   excess outputs labelled and specified.

## Candidates to add (published methods)

Verified references; not implemented in E3.

| Method | Reference | Fit for mortensor |
|---|---|---|
| Nordpred: APC with a power link and damped drift, used for cancer incidence and mortality projections | Møller et al. (2003) Statistics in Medicine 22(17): 2751-2766 | Public-health projections; simple Poisson GLM |
| 2D P-spline forecasting | Currie, Durban and Eilers (2004) Statistical Modelling 4(4): 279-298 | Implemented in E4 as `p_splines()`, adapted from MortalitySmooth 2.3.4 |
| Smooth constrained forecasting | Camarda (2019) Demographic Research 41(38): 1091-1130 | P-spline forecasts constrained by asymmetric penalties on the rates of ageing and of change (corrected after reading the full text; implemented in E4 as `cp_splines()`) |
| Plat model | Plat (2009) Insurance: Mathematics and Economics 45(3): 393-404 | A StMoMo specification; small addition to `stmomo()` |
| APCI model | Richards, Currie, Kleinow and Ritchie (2019) British Actuarial Journal 24: e13 | CMI model used in UK practice |
| Negative-binomial Lee-Carter | Delwarde, Denuit and Partrat (2007) ASMBI 23(5): 385-401 | Published overdispersion in place of the removed layer |
| Common age effect model for several populations | Kleinow (2015) Insurance: Mathematics and Economics 63: 147-152 | Applied to socio-economic groups in England |
| SAINT model for small populations | Jarner and Kryger (2011) ASTIN Bulletin 41(2): 377-418 | Small populations against a reference population |
| Varying age-specific improvements | Bohk-Ewald and Rau (2017) Genus 73: 1 | Bayesian; heavier |
