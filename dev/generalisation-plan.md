# From package-specific solutions to experimental methods: development plan (draft)

Draft of 2026-10-18 for the owner. It builds on the audit of work package E3
(`dev/methods-audit.md`). Nothing in it is implemented. Numbers are quoted
from the repository files named next to them.

## 1. Purpose and principles

The audit replaced package-specific solutions by published methods. Several
of the set-aside solutions answered problems that remain in mortensor's own
results. This plan turns each problem into experimental functionality built
on published methods.

1. **Published method first.** Implement the published method as published
   and record deviations in "Methods and sources". Add package-specific
   parts only where the literature leaves a gap, and name them.
2. **Experimental means opt-in and labelled.** An experimental option is
   never a default. Help page and methods article give the status
   "experimental", and the model label ends in "(experimental)", so score
   tables show it. A registry `mt_methods()` (function, option, source,
   status) keeps help pages, article and tests consistent.
3. **Known truth without the method's assumption.** Earlier package methods
   were checked on generators that contained their own assumption: the
   extra-Poisson noise; the regional dependence that the empirical
   reordering copies (`regions_coherence`: aggregate coverage 0.955 with it,
   0.843 without); the stationary open interval of the PCLM; the log-linear
   cause shares of `mt_example_causes()` (`R/model-causes.R`), on which
   `cause_shares()` beat `coda_causes()` (`NEWS.md`). Every new check uses
   three generators: **G0**, the method's own model (tests the code); **G1**,
   a different plausible model, where the method must not lose to the
   current default; **G2**, a stress case (trend break, small counts, short
   series).
4. **Real data decide defaults.** A default changes only after G1 and
   backtests on at least three of England & Wales men (StMoMo), France
   (demography), `ons_data` and `de_deaths`, with 16 or more origins where
   possible, tested with `mt_dm_test()` and the EWC variance (Lazarus et al.
   2018). Scores: CRPS, coverage by horizon, and the
   quantity of the question (e0, gaps, totals).
5. **Evidence status.** "(abstract)" marks works read only as abstracts; the
   plan reports only what these state, and implementation waits for the
   full text.

## 2. Problems

### a. Intervals too narrow: overdispersion and unmodelled variation

**Problem.** Draws of models of smooth rates carry period-index uncertainty
plus Poisson deaths. Real deaths vary more around a rigid fitted surface.

**Evidence.** Without the extra-Poisson layer, 95% coverage of
`lee_carter()` on E&W men (ages 30-100, origins 1996-2001, h = 10) fell from
0.746 to 0.632, and CRPS rose from 0.0687 (`dev/progress.md`, 2026-10-13) to
0.0709 (`validation/thresholds.csv`, `backtest_ew`). The Bayesian Lee-Carter
covers 0.729 at equal accuracy (`stanmomo_ew`). ONS rates vary around a
smoothed trend no more than Poisson noise (dispersion 0.4-1.0;
`dev/progress.md`). Overdispersion depends on data and model.

**Set aside.** `extra_noise_sd()` added moment-estimated lognormal noise per
cell to every model of smooth rates. It had no source, counted noise twice
for `tucker3()` and the functional data models, and replaced Bayesian
posterior predictives.

**Literature.** A negative binomial Lee-Carter accounts for overdispersion
(Delwarde, Denuit and Partrat 2007, abstract). Lee-Carter intervals are
excessively narrow mainly because of the rigid structure; heterogeneity in
each age-period cell improves the fit and widens intervals (Li, Hardy and
Tan 2009, abstract). Variance proportional to the mean, with smooth
dispersion estimated by quasi-likelihood, reduces the volatility of a
forecast table (Djeundje and Currie 2011, abstract). A Bayesian Lee-Carter
with all sources of variability gives wider, appropriate intervals (Pedroza
2006, abstract).

**Proposal.**
- `lee_carter(family = c("poisson", "nb"), dispersion = c("global", "age"))`
  by maximum likelihood (Delwarde et al. 2007). The family belongs to the
  fit: `mt_forecast(predictive = TRUE)` and `mt_backtest()` draw from it, as
  for `stanmomo(family = "nb")`. Models of observed rates get no extra
  noise. The age-specific form follows Li et al. (2009) once read in full.
- `mt_dispersion(fit)`: Pearson dispersion by age and stratum (McCullagh and
  Nelder 1989), to show whether the problem applies.
- Phase 3: smooth dispersion (Djeundje and Currie 2011) for `p_splines()`
  (Currie, Durban and Eilers 2004, implemented in E4 with MortalitySmooth's
  constant overdispersion as an option).

**Validation.** Check `overdispersion_sim`. G0: negative binomial
Lee-Carter. G1: Poisson-lognormal noise with age-varying variance; a rank-2
Lee-Carter truth fitted with rank 1. G2: populations of 2e4. Metrics:
coverage (80%, 95%) per horizon, CRPS, PIT. Thresholds: G0 cov95 0.90-0.98;
G1 cov95 at h 1-3 at least 0.90 with CRPS no worse than Poisson Lee-Carter;
on Poisson truth (`simulate_lc_dgp()`) CRPS within 1%. Real (`nb_lc_ew`):
cov95 at least 0.70, CRPS at most 0.0709; no loss on ONS.

**Risks.** The dispersion absorbs structural misfit and widens all horizons
alike; the trend error that grows with the horizon remains (b). Cell noise
averages out over ages, so e0 intervals barely change.


### b. Calibration that worsens with the horizon

**Problem.** Forecast errors grow faster with the horizon than model
spreads.

**Evidence.** One factor for all horizons raised the backtest 95% coverage
of Lee-Carter rates on E&W from 0.63 to 0.85; coverage still falls with the
horizon (`vignettes/tuning-recalibration.Rmd`). `calibration_ew`: 0.847 with
one factor, 0.877 with the former power law and noise layer together.
Lee-Carter e0 intervals were too broad up to 40 years and too narrow after
50 (Lee and Miller 2001, abstract).

**Set aside.** `recalibrated(scale = "horizon")`: a spread factor s1 h^γ
without a source.

**Literature.** Member-by-member post-processing works per lead time,
location and variable and keeps rank correlations (Van Schaeybroeck and
Vannitsem 2015, abstract). Distinct training periods may suit distinct
horizons (Gneiting et al. 2005). Extrapolating errors of past forecasts is
an established method, limited by short error series and by errors that
shrink over forecast rounds (Keilman, Pham and Hetland 2002). The scaled
model for error sums increments S(j,t)(η_j + δ(j,t)): a persistent trend
error η with variance k and independent noise δ with variance 1 - k,
correlated over ages; the error variance grows like k t² + (1 - k) t (Alho
and Spencer 2005: 280-283, as given by Li et al. 2007). Shang and Hyndman
(2017) build mortality intervals from bootstrapped in-sample h-step errors,
scaled to nominal coverage. Conformal intervals for time series need
short-term independent or strongly mixing errors; without feedback their
multi-step coverage was poor after a shift (Xu and Xie 2023).

**Proposal.**
- `recalibrated(scale = "per_horizon")`: one factor per horizon, the
  published form. Horizons with fewer than five inner forecasts take the
  nearest estimable factor (package-specific). Needs `years = 20` or more.
- `recalibrated(scale = "scaled_error")`: cell factor = scaled-model error,
  σ (k h² + (1 - k) h)^½ per age band, over the model's own spread; σ and k
  by minimum CRPS. Published parts, package-specific combination.
- `recalibrated(method = "errors")`: median paths plus resampled past h-step
  error vectors of the model (Keilman et al. 2002; Shang and Hyndman 2017).
  A vector holds all ages and strata of one inner origin and horizon, which
  keeps their dependence (joint bootstrap, Panagiotelis et al. 2023). No
  conformal guarantee is claimed.
- `bias` and `ages` stay opt-in, with their published analogues named (mean
  correction, post-processing per location). The drift correction worsened
  Li-Lee on ONS (CRPS 0.0818 against 0.0687; `dev/progress.md`).
- Coverage by horizon, with a slope, in `mt_score_table()`.

**Validation.** Check `horizon_calibration_sim`. G1: drift that changes at a
random year (variable decline, Booth, Maindonald and Smith 2002, abstract);
AR(1) deviations around a trend; a random walk with uncertain drift. G2:
15-year series, where factors must fall back to one constant. Threshold on
G1: |cov95(h) - 0.95| at most 0.05 inside the window, CRPS no worse than
`"constant"`. Real: E&W with 20-year windows (origins 1981-2001), France
(1960-1996).

**Risks.** With the defaults, horizon 10 rests on one origin. One stall
(England after 2011) can dominate past errors.


### c. Dependence between strata: totals, gaps, inequality

**Problem.** Sums of independently simulated strata cancel their errors.
Joint models that tie strata together pull gaps back to past levels.

**Evidence.** England total from separate models by sex and quintile,
reconciled by WLS: 95% coverage 0.434, 0.552 with the reordering
(`NEWS.md`). England e0 coverage
0.24; Li-Lee and the coherent fdm covered the female gap in 4% and 18% of
cases, Lee-Carter per stratum in all (`vignettes/inequality-england.Rmd`).
Lee-Carter e0 intervals for Germany covered none of 15 values
(`vignettes/laender-germany.Rmd`). German excess deaths in 2020: 11,000 to
28,000 from 32 Land models, -1,000 to 55,000 from a model of the total
(`dev/progress.md`).

**Set aside.** `dependence = "empirical"`: one template value per stratum and
year, whole paths ranked by mean level, years resampled with replacement.

**Literature.** The Schaake shuffle matches, per margin, the ranks of
ensemble members to those of observations from as many distinct dates;
consecutive dates keep persistence. It assumes a stationary correlation
structure, and its authors note that sums of separately calibrated parts
underestimate extremes (Clark et al. 2004). Ensemble copula coupling takes
the ranks of a raw ensemble and assumes its dependence correct (Schefzik,
Thorarinsdottir and Gneiting 2013). Reordering to the copula of residuals
needs as many draws as training observations and did not improve on base
forecasts in simulations; MinT with shrinkage did best for jointly generated
base forecasts (Panagiotelis et al. 2023). Reconciled intervals need the
covariance of base forecast errors (Hyndman et al. 2011); Raftery et al.
(2013) aggregate with perturbations drawn with residual correlations. The
common age effect (CAE) model (Kleinow 2015, abstract) fitted better than
variants of Li-Lee, and gave non-diverging joint scenarios (Enchev, Kleinow
and Cairns 2017, abstract); a special case fitted English deprivation deciles
best, with group period effects more important than group age effects (Wen,
Cairns and Kleinow 2021, abstract). Villegas et al. (2017, abstract) compare
two-population models and give guidance on data and calibration. A
one-factor copula couples filtered populations (Chen, MacMinn and Sun 2015,
abstract).

**Proposal.**
- `reconciled(dependence = "schaake")`: template per margin (age, stratum,
  horizon) of past h-step errors or year-to-year changes from distinct
  years; draws in blocks of the template size (labelled deviation).
- `reconciled(dependence = "ecc", template = li_lee(common_by = "sex"))`:
  margins of per-stratum models, ranks of a joint model's draws.
- `reconciled(method = "mint", covariance = "errors")`: W from inner h-step
  errors, the published input (Wickramasuriya, Athanasopoulos and Hyndman
  2019).
- `cae()`: stratum period indexes as a multivariate random walk with drift
  and correlated innovations; `index = "coherent"` for stationary
  deviations.

**Validation.** Check `dependence_sim` with dependence no method assumes:
common shocks whose correlation varies by age, spatially correlated regional
innovations, a break shared by half the strata; correlation 0 to 0.9.
Metrics: coverage and CRPS of totals, gaps (`contrast`) and SII/RII
(`mt_inequality()`), energy and variogram scores; margins unchanged.
Thresholds: totals cov95 0.85-0.98 at correlations of 0.5 and above; no loss
against `"none"` at 0. Real: `dependence_ons`, German Laender; DM tests.

**Risks.** Templates are short (ONS: 19 years). Dependence can change. ECC
inherits the template model's pull on gaps. CAE needs identifiability
constraints (Enchev et al. 2017).


### d. Jump-off bias and the first forecast years

**Problem.** Fitted rates of the last year can miss the observed level;
observed rates of one year carry their noise into every forecast year.

**Evidence.** The Lee-Miller start improves E&W Lee-Carter in years 1-2
(CRPS 0.0305 against 0.0418) and worsens Li-Lee on ONS (0.0714 against
0.0671); the former fading three-year offsets gave a ratio of 0.982
(`smoothing_jumpoff`). On ONS the observed start improves Li-Lee e0 (0.316
against 0.338 years) and worsens its rates (0.071 against 0.067) and those of
Lee-Carter (0.079 against 0.068; `vignettes/jump-off.Rmd`).

**Set aside.** `anchored(years = 3, decay = 10, shrink = TRUE)`: three-year
mean offset, linear fade, DerSimonian-Laird shrinkage.

**Literature.** Stoeldraijer et al. (2013) use observed jump-off rates after
Lee and Miller (2001), describe a conflict between matching the first year
and robustness to data updates, and find the choice moves Dutch male e0 in
2050 by one year. Their Table 1 lists a Danish short-term correction of the
jump-off error and a Dutch method with an 8-year trend for the near future
and a 20-year trend for the eventual level. Janssen (2018) summarises a
later study of six jump-off options: the last observed year was most
accurate for most countries, an average over more years the most robust.
Projecting half the trend after 10 years improved long-term cancer
predictions (Møller et al. 2003, abstract). SAINT models deviations from a
trend by a time series model (Jarner and Kryger 2011, abstract).

**Proposal.**
- Relabel `anchored(years = k)` as an averaged observed jump-off, after
  reading the study that Janssen (2018) summarises; it is not in the
  collections.
- `anchored(decay = "ar1")`: the offset decays with the AR(1) coefficient of
  the in-sample deviations (observed minus fitted log rates), the deviation
  process of Li-Lee and SAINT; it generalises the fixed fade.
- `anchored(to = "filtered")`: offset from a local level filter of the
  deviations (Harvey 1989; state-space Lee-Carter, De Jong and Tickle 2006,
  abstract).
- Backtest metric `revision`: change of the forecast for a target year
  between consecutive origins (robustness as in Stoeldraijer et al. 2013).

**Validation.** Check `jumpoff_sim`. G1: a level shift one to three years
before the origin; independent year effects; AR(1) deviations with φ 0.3,
0.7 and 0.95. Metrics: CRPS at h 1-3 and 4-10, e0 error, revision.
Thresholds: `"ar1"` within 2% of the better of fitted and observed start on
every generator; on ONS no worse than the fitted start. Real: E&W, France,
ONS, Germany.

**Risks.** Short deviation series bias AR(1) estimates towards fast decay.
Accuracy and robustness can conflict.


### e. Combining models

**Problem.** Single models are too narrow or wrong in different ways;
mixtures need weights.

**Evidence.** On ONS the ensemble covers 0.926 against 0.810 for Li-Lee, at
CRPS 0.0699 against 0.0671 (`ons_deprivation`). On simulated data it
over-covers (0.974; `coverage_ensemble`).

**Set aside.** Weights inverse to the mean inner CRPS; stacking weights of
medians used as mixture weights.

**Literature.** Stacking of predictive distributions maximises a proper
score (the CRPS qualifies) of held-out predictions, does better than
Bayesian model averaging when no candidate is true, and is unstable in small
samples (Yao et al. 2018). Barigou et al. (2023), cited in the package,
average mortality models with leave-future-out validation. Kontis et al.
(2017) weighted 21 models by exp(-|bias|) of e0 over 13 withheld years and
pooled draws in proportion; 90% intervals covered 90.2% (women) and 83.9%
(men). Weights inverse to the past mean absolute errors of the two best
models gave the lowest one-step errors, equal weights the highest (Shang
2012, Table 18). Janssen (2018) summarises a study that gives equal weights
to the models of a model confidence set. Linear pools of calibrated
components are overdispersed (Gneiting and Ranjan 2013, as cited in the
package). No text of Hall and Mitchell (2007) was available. The CRPS of a
point forecast is its absolute error, so inverse-CRPS weights extend
Shang's weights to distributions.

**Proposal.**
- `ensemble(weights = "crps_stack")`: simplex weights that minimise the mean
  CRPS of the mixture on inner rolling-origin forecasts (Yao et al. 2018,
  leave-future-out); optional shrinkage towards equal weights
  (package-specific).
- `ensemble(weights = "bias")`: the weights of Kontis et al. (2017).
- Relabel `"crps"` (Bates and Granger 1969; Shang 2012); deprecate `"stack"`
  once `"crps_stack"` passes; bootstrap intervals of the weights.

**Validation.** Check `stacking_sim`: M-closed (truth among the components)
and M-open (Lee-Carter with a drift break and a cohort effect; components
Lee-Carter, `rwd()`, `coda()`, ARIMA(0,1,1)). Thresholds: weight of the true
model at least 0.6 with 10 inner origins; M-open CRPS no worse than equal
weights. Real: E&W (16 origins), France, ONS, Germany; DM test against equal
weights.

**Risks.** Three inner origins by default give unstable weights; weights are
learnt in one regime; recalibrating a wide pool can narrow it.


### f. Coherent forecasts of causes of death

**Problem.** Cause forecasts must add up to all causes, and their sum tends
to be pessimistic.

**Evidence.** Summed causes give female Q1 e50 of 30.8 in 2034, below the
31.7 observed in 2019, against 32.9 (Lee-Carter, all causes) and 34.3
(Li-Lee) (`dev/progress.md`, 2026-10-17)
and the worst e50 gap (CRPS 0.499 against 0.153 when scaled to the
Lee-Carter total; `vignettes/inequality-england.Rmd`). `coda_causes()`: CRPS
0.0816 against 0.0631 for reconciled Lee-Carter (`coda_causes_coherence`).

**Set aside.** `cause_shares()` (Tucker model of log-ratio cause shares);
`tucker3(transform = "coda")`.

**Literature.** Wilmoth (1995, abstract) measures the pessimism of
cause-based forecasts, approximates it from cause distributions in two
periods, and attributes it mainly to cancer and heart disease, with a caveat
on diagnostic practice. Oeppen (2008) forecasts multiple-decrement life
tables as compositions; for Japanese women the cause model was not more
pessimistic than the all-cause model, and more testing is needed. Kjærgaard
et al. (2019, abstract) weight ages, years and causes and separate joint
from cause-specific variation, which improved cancer forecasts in France and
the Netherlands. Bergeron-Boucher et al. (2018, abstract) adapt Tucker3 to
compositions of life-table deaths by time, age and population: the structure
of the deprecated `tucker3(transform = "coda")`, whose help says it has no
published source.

**Proposal.**
- `mt_pessimism()`: per path and year, the e0 difference between an
  all-cause forecast and summed causes; Wilmoth's measure after the full
  text.
- Review `tucker3(transform = "coda")` against Bergeron-Boucher et al.
  (2018); restore it as published with deviations if it matches.
- `coda_causes()` weights and joint components (Kjærgaard et al. 2019) after
  the full text. `cause_shares()` stays deprecated.

**Validation.** Check `causes_triad`. G0: current log-linear shares. G1a:
independent Lee-Carter per cause. G1b: causes sharing a period index.
Methods: reconciled Lee-Carter (WLS; proportional to a fixed total),
`coda_causes()`, CoDa Tucker3. Metrics: CRPS of cause rates, all-cause e50,
`mt_pessimism()`. Real: ONS causes at 50+, `de_deaths` chapters. A default
changes only for a method best on G1a, G1b and ONS.

**Risks.** Coding changes; zeros in compositions; small causes.


### g. Small populations and pooling

**Problem.** Small strata give noisy rates and trends; pooling must keep
real differences.

**Evidence.** ARIMA(0,1,1) with drift per age covers 0.823 on simulated
30-year series (`coverage_local_level`), `rwd()` 0.974. For the German
Laender, an ensemble of Li-Lee and that model covers 0.865 (`de_regions`). Pooling the sexes helps in France (coherent fdm CRPS
ratio 0.889, `coherent_fdm_france`). On ONS, Li-Lee bounds some AR(1)
coefficients at 0.99 and misses the gap.

**Set aside.** `product_ratio()` with Lee-Carter components;
`rwd(observation_noise = TRUE)`, a moment estimator of a local level model.

**Literature.** The product-ratio method keeps ratio functions stationary for
coherence; coherent forecasts were at least as accurate overall (Hyndman,
Booth and Yasmeen 2013). Bottom-up grouped functional forecasts were most
accurate for Japanese prefectures (Shang and Hyndman 2017). A relative model
of English deprivation groups against a reference population fitted and
forecast best (Villegas and Haberman 2014, abstract). SAINT models a
reference trend and deviations from it (Jarner and Kryger 2011, abstract;
Jarner and Jallbjørn 2022, abstract); Cairns et al. (2011, abstract) model a
small population against a large one with a mean-reverting spread.
Alexander, Zagheni and Barbieri (2017) pool across areas and smooth over
time; their intervals were calibrated in simulation, and a Brass benchmark
did well because their data came from a Brass model.

**Proposal.**
- `relative_model(reference = lee_carter(), deviation = c("ar1", "rw",
  "rwd"))`: a model of the reference stratum (England, Germany), and each
  stratum's log ratio to it as an age pattern times an index (Villegas and
  Haberman 2014; SAINT). `"ar1"` is coherent; `"rw"` and `"rwd"` let gaps
  persist or trend. It generalises `product_ratio()`.
- Test `bootstrapped(arima_series(order = c(0, 1, 1), drift = TRUE))` for the
  missing parameter uncertainty; if coverage stays low, a state-space
  `local_level()` with signal-to-noise ratios pooled over ages (Harvey 1989;
  pooling package-specific).

**Validation.** Check `small_pop_sim`. G1: a reference with AR(1)
deviations; a reference with drifting deviations; populations 1e3 to 1e5.
Real: `relative_ons`, Laender, regions. Metrics: CRPS, coverage, gap scores,
implausible draws. Thresholds: ONS gap coverage at least 0.8 with `"rw"`;
rate CRPS within 3% of Li-Lee.

**Risks.** The reference's own error dominates; drifting deviations widen
gaps without limit.


### h. Old ages, open intervals and ungrouping

**Problem.** Open age intervals and grouped old ages bias rates, life
expectancy and lifespan variation.

**Evidence.** PCLM error at ages 90-100 on E&W: 0.024 against a noise floor
of 0.021 (`ungroup_ew_1d`); the deprecated spline methods 0.027-0.031
(`dev/progress.md`, P5.8). An open 85+ exposure split from its total gave
errors of 27% at 85-99 and a factor of about 6 at 100+ (France). E-dagger of
French women in 2006 is 22% higher with 85+ than with 100+
(`dev/progress.md`).

**Set aside.** Spline ungrouping with Kannisto or Gompertz tails;
`preserve_totals = TRUE`; stationary-population weights for the open
interval inside the PCLM (default, package-specific).

**Literature.** Rizzi, Gampe and Eilers (2015) close an open interval with a
maximum age and a zero-count bin, and ungroup grouped exposures first and
insert them in the composition matrix; dividing separately ungrouped deaths
and exposures was less accurate, most in the open class. Their standard
errors understate the variability. Thatcher, Kannisto and Vaupel (1998) fit
six models at ages 80-98, test extrapolations at 99 and above, and find the
logistic and Kannisto models best. On 360 cohorts aged 80+, models with
deceleration fit better, the log-quadratic predicts well most consistently,
and no model fits everywhere (Feehan 2018). The collected literature has no
method that forces exact group totals after a PCLM fit; `preserve_totals`
stays labelled.

**Proposal.**
- Done in E4: `mt_ungroup(open_split = "pclm")`, the two-step method of
  Rizzi et al. (2015), is the default; stationary weights are an option.
  `life_table(open = "kannisto"/"gompertz")` is fitted by maximum
  likelihood (Thatcher et al. 1998).
- `life_table(open = c(..., "logquad", "beard"))`: Feehan's further models,
  and `open = "aic"`.
- `mt_open_check()`: fit below a cut age, score observed ages above it (the
  design of Thatcher et al. 1998).
- Ungrouping uncertainty through `bootstrapped(type = "data")`.

**Validation.** Check `open_interval_holdout`: France (ages to 110) and E&W
(to 100), cuts at 85 and 90; errors of log rates at 90-100, e85 and
e-dagger. G1 truths: log-quadratic and gamma-Gompertz. Threshold: the
default closure within 5% of the best on both truths and both countries.

**Risks.** Feehan studies cohorts, mortensor period data. Data quality above
age 100.


### i. Life expectancy as the target

**Problem.** Users ask about life expectancy and its gaps; rate models give
miscalibrated intervals for both.

**Evidence.** e0 intervals are too narrow on ONS (England 0.24) and for the
German Laender (0.67-0.78); fdm covers 0.30 for E&W men aged 55-90
(`vignettes/model-choice.Rmd`). On E&W, e0 intervals cover well while rate
intervals are too narrow, because errors of single ages partly cancel
(`vignettes/tuning-recalibration.Rmd`); calibrating rates does not
calibrate e0.

**Scope.** Work package E2 excluded models of life expectancy
(`dev/progress.md`). The options are for the owner.

**Literature.** Raftery et al. (2013) model e0 gains by a double-logistic
function in a Bayesian hierarchy; out of sample, 80% intervals held 82% and
95% intervals 92% of outcomes; for two sexes they suggest a gap model.
Best-practice e0 can be forecast by ARIMA and the gap to it by geometric
Brownian motion or mean reversion (Torri and Vaupel 2012, abstract); female
e0 from the gap to record e0, male e0 from the sex gap (Pascariu,
Canudas-Romo and Vaupel 2018, abstract). Camarda (2019) constrains P-spline
rate forecasts to observed ranges of age profiles and time trends; it was
best in 39 of 72 out-of-sample comparisons on e0, e-dagger and rates.
Bohk-Ewald and Rau (2017) forecast age-specific improvement rates with other
countries' trends.

**Options.**
1. Evaluation only (in scope): SAPE (Raftery et al. 2013) and e0
   coverage by horizon in `mt_score()`.
2. External e0 drives rates (in scope):
   `lee_carter(rotate = TRUE, e0 = <draws>)`, as MortCast takes e0 from a
   separate model (Li, Lee and Gerland 2013).
3. e0 benchmarks (scope change): double-gap and best-practice gap
   specifications that return e0 draws for backtests, with international
   series loaded from other packages.
4. Constrained P-spline rates (Camarda 2019): implemented in E4 as
   `cp_splines()`, which constrains the rates of ageing and of change by
   asymmetric penalties. A version constrained by external forecasts of e0
   or lifespan disparity has no source in the collected literature.

## 3. Roadmap

| Phase | Content | Validation checks added |
|---|---|---|
| 0 Foundation | G1/G2 generators in `R/validation.R`; coverage by horizon; variogram score; `revision`; SAPE; `mt_dispersion()`; `mt_methods()` | none with thresholds; current defaults scored on G1 and recorded in `history.csv` |
| 1 Published, light | `scale = "per_horizon"`; `weights = "crps_stack"`, `"bias"`; `dependence = "schaake"`, `"ecc"`; MinT from errors; further closures; relabels | `horizon_calibration_*`, `stacking_*`, `dependence_*`, `open_interval_holdout` |
| 2 Models | `lee_carter(family = "nb")`; `cae()`; `relative_model()`; `anchored(decay = "ar1", to = "filtered")`; bootstrapped local level | `overdispersion_sim`, `nb_lc_ew`, `cae_ons`, `relative_ons`, `jumpoff_sim` |
| 3 After full texts or decisions | `scaled_error`; `method = "errors"`; causes; Tucker3 CoDa review; P-splines with dispersion; e0 options; `local_level()` | `causes_triad`, further checks as built |

Phase 0 comes first. The inner-error machinery of b serves c (MinT from
errors, Schaake templates) and e. `cae()` and `relative_model()` share the
joint simulation of period indexes. Recalibration and stacking checks are
re-baselined after `family = "nb"`, because base spreads change. Each item
is one work package with NEWS, progress log, help, vignettes and pkgdown
(`CLAUDE.md`).

**Decisions for the owner.**
1. Adopt Section 1: three generators per check, the label "experimental".
2. May `lee_carter(family = "nb")` and `recalibrated(scale = "per_horizon")`
   become defaults if they pass? Both change results again.
3. Which dependence method first: Schaake shuffle, ECC or `cae()`?
4. Offer reverting (coherent) deviations for deprivation groups, although
   they forecast the ONS and German gaps worst? Proposed: yes, with the
   deviation process as an explicit argument.
5. Review the deprecation of `tucker3(transform = "coda")`; relabel
   `anchored(years)` and `ensemble(weights = "crps")`.
6. Life expectancy: options 1-4 of i.
7. Obtain full texts first: Delwarde et al. 2007; Li, Hardy and Tan 2009;
   Djeundje and Currie 2011; Van Schaeybroeck and Vannitsem 2015; Alho and
   Spencer 2005; Kleinow 2015; Enchev et al. 2017; Wen et al. 2021; Villegas
   and Haberman 2014; Jarner and Kryger 2011; Kjærgaard et al. 2019;
   Bergeron-Boucher et al. 2018; Wilmoth 1995; Hall and Mitchell 2007; and
   the two studies on jump-off rates and model confidence sets that Janssen
   (2018) summarises (not yet verified).
8. (Done in E4: the audit's description of Camarda (2019) is corrected.)
9. Run `cae()` and negative binomial Lee-Carter backtests in the full suite.

## 4. References

"(abstract)": only the abstract was read. Citations follow the manifests of
the literature collections or the package's roxygen references.

- Alexander M, Zagheni E, Barbieri M (2017). A flexible Bayesian model for estimating subnational mortality. Demography 54(6): 2025-2041. doi:10.1007/s13524-017-0618-7
- Alho JM, Spencer BD (2005). Statistical Demography and Forecasting. Springer Series in Statistics. Springer. doi:10.1007/0-387-28392-7 (pages 280-283 as cited by Li et al. 2007)
- Barigou K, Goffard P-O, Loisel S, Salhi Y (2023). Bayesian model averaging for mortality forecasting using leave-future-out validation. International Journal of Forecasting 39(2): 674-690. doi:10.1016/j.ijforecast.2022.01.011
- Bates JM, Granger CWJ (1969). The combination of forecasts. Operational Research Quarterly 20(4): 451-468.
- Bergeron-Boucher M-P, Simonacci V, Oeppen J, Gallo M (2018). Coherent modeling and forecasting of mortality patterns for subpopulations using multiway analysis of compositions: an application to Canadian provinces and territories. North American Actuarial Journal 22(1): 92-118. doi:10.1080/10920277.2017.1377620 (abstract)
- Bohk-Ewald C, Rau R (2017). Probabilistic mortality forecasting with varying age-specific survival improvements. Genus 73(1): 1. doi:10.1186/s41118-016-0017-8
- Booth H, Maindonald J, Smith L (2002). Applying Lee-Carter under conditions of variable mortality decline. Population Studies 56(3): 325-336. doi:10.1080/00324720215935 (abstract)
- Cairns AJG, Blake D, Dowd K, Coughlan GD, Khalaf-Allah M (2011). Bayesian stochastic mortality modelling for two populations. ASTIN Bulletin 41(1): 29-59. (abstract; DOI not confirmed)
- Camarda CG (2019). Smooth constrained mortality forecasting. Demographic Research 41(38): 1091-1130. doi:10.4054/DemRes.2019.41.38
- Chen H, MacMinn R, Sun T (2015). Multi-population mortality models: a factor copula approach. Insurance: Mathematics and Economics 63: 135-146. doi:10.1016/j.insmatheco.2015.03.022 (abstract)
- Clark M, Gangopadhyay S, Hay L, Rajagopalan B, Wilby R (2004). The Schaake shuffle: a method for reconstructing space-time variability in forecasted precipitation and temperature fields. Journal of Hydrometeorology 5(1): 243-262. doi:10.1175/1525-7541(2004)005<0243:TSSAMF>2.0.CO;2
- Coroneo L, Iacone F (2020). Comparing predictive accuracy in small samples using fixed-smoothing asymptotics. Journal of Applied Econometrics 35(4): 391-409. doi:10.1002/jae.2756
- Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting mortality rates. Statistical Modelling 4(4): 279-298. doi:10.1191/1471082X04st080oa (abstract)
- De Jong P, Tickle L (2006). Extending Lee-Carter mortality forecasting. Mathematical Population Studies 13(1): 1-18. doi:10.1080/08898480500452109 (abstract)
- Delwarde A, Denuit M, Partrat C (2007). Negative binomial version of the Lee-Carter model for mortality forecasting. Applied Stochastic Models in Business and Industry 23(5): 385-401. doi:10.1002/asmb.679 (abstract)
- Djeundje VAB, Currie ID (2011). Smoothing dispersed counts with applications to mortality data. Annals of Actuarial Science 5(1): 33-52. doi:10.1017/S1748499510000047 (abstract)
- Enchev V, Kleinow T, Cairns AJG (2017). Multi-population mortality models: fitting, forecasting and comparisons. Scandinavian Actuarial Journal 2017(4): 319-342. doi:10.1080/03461238.2015.1133450 (abstract)
- Feehan DM (2018). Separating the signal from the noise: evidence for deceleration in old-age death rates. Demography 55(6): 2025-2044. doi:10.1007/s13524-018-0728-x
- Gneiting T, Raftery AE, Westveld AH, Goldman T (2005). Calibrated probabilistic forecasting using ensemble model output statistics and minimum CRPS estimation. Monthly Weather Review 133(5): 1098-1118. doi:10.1175/MWR2904.1
- Gneiting T, Ranjan R (2013). Combining predictive distributions. Electronic Journal of Statistics 7: 1747-1782.
- Hall SG, Mitchell J (2007). Combining density forecasts. International Journal of Forecasting 23(1): 1-13. doi:10.1016/j.ijforecast.2006.08.001 (title only)
- Harvey AC (1989). Forecasting, Structural Time Series Models and the Kalman Filter. Cambridge University Press.
- Hyndman RJ, Ahmed RA, Athanasopoulos G, Shang HL (2011). Optimal combination forecasts for hierarchical time series. Computational Statistics & Data Analysis 55(9): 2579-2589. doi:10.1016/j.csda.2011.03.006
- Hyndman RJ, Booth H, Yasmeen F (2013). Coherent mortality forecasting: the product-ratio method with functional time series models. Demography 50(1): 261-283. doi:10.1007/s13524-012-0145-5
- Janssen F (2018). Advances in mortality forecasting: introduction. Genus 74(1): 21. doi:10.1186/s41118-018-0045-7
- Jarner SF, Jallbjørn S (2022). The SAINT model: a decade later. ASTIN Bulletin 52(2): 483-517. doi:10.1017/asb.2021.37 (abstract)
- Jarner SF, Kryger EM (2011). Modelling adult mortality in small populations: the SAINT model. ASTIN Bulletin 41(2): 377-418. doi:10.2143/AST.41.2.2136982 (not registered with Crossref; abstract)
- Keilman N, Pham DQ, Hetland A (2002). Why population forecasts should be probabilistic - illustrated by the case of Norway. Demographic Research 6(15): 409-454. doi:10.4054/DemRes.2002.6.15
- Kjærgaard S, Ergemen YE, Kallestrup-Lamb M, Oeppen J, Lindahl-Jacobsen R (2019). Forecasting causes of death by using compositional data analysis: the case of cancer deaths. Journal of the Royal Statistical Society Series C 68(5): 1351-1370. doi:10.1111/rssc.12357 (abstract)
- Kleinow T (2015). A common age effect model for the mortality of multiple populations. Insurance: Mathematics and Economics 63: 147-152. doi:10.1016/j.insmatheco.2015.03.023 (abstract)
- Kontis V, Bennett JE, Mathers CD, Li G, Foreman K, Ezzati M (2017). Future life expectancy in 35 industrialised countries: projections with a Bayesian model ensemble. The Lancet 389(10076): 1323-1335. doi:10.1016/S0140-6736(16)32381-9
- Lazarus E, Lewis DJ, Stock JH, Watson MW (2018). HAR inference: recommendations for practice. Journal of Business & Economic Statistics 36(4): 541-559. doi:10.1080/07350015.2018.1506926 (abstract)
- Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter method for forecasting mortality. Demography 38(4): 537-549. doi:10.1353/dem.2001.0036 (abstract)
- Li JS-H, Hardy MR, Tan KS (2009). Uncertainty in mortality forecasting: an extension to the classical Lee-Carter approach. ASTIN Bulletin 39(1): 137-164. doi:10.2143/AST.39.1.2038060 (abstract)
- Li N, Lee R, Gerland P (2013). Extending the Lee-Carter method to model the rotation of age patterns of mortality decline for long-term projections. Demography 50(6): 2037-2051. doi:10.1007/s13524-013-0232-2
- Li Q, Reuser M, Kraus C, Alho J (2007). Aging of a giant: a stochastic population forecast for China, 2001-2050. MPIDR Working Paper WP-2007-032. doi:10.4054/MPIDR-WP-2007-032
- McCullagh P, Nelder JA (1989). Generalized Linear Models (2nd ed.). Chapman and Hall.
- Møller B, Fekjær H, Hakulinen T, Sigvaldason H, Storm HH, Talbäck M, Haldorsen T (2003). Prediction of cancer incidence in the Nordic countries: empirical comparison of different approaches. Statistics in Medicine 22(17): 2751-2766. doi:10.1002/sim.1481 (abstract)
- Oeppen J (2008). Coherent forecasting of multiple-decrement life tables: a test using Japanese cause of death data. CODAWORK'08, Universitat de Girona. hdl:10256/742
- Panagiotelis A, Gamakumara P, Athanasopoulos G, Hyndman RJ (2023). Probabilistic forecast reconciliation: properties, evaluation and score optimisation. European Journal of Operational Research 306(2): 693-706. doi:10.1016/j.ejor.2022.07.040
- Pascariu MD, Canudas-Romo V, Vaupel JW (2018). The double-gap life expectancy forecasting model. Insurance: Mathematics and Economics 78: 339-350. doi:10.1016/j.insmatheco.2017.09.011 (abstract)
- Pedroza C (2006). A Bayesian forecasting model: predicting U.S. male mortality. Biostatistics 7(4): 530-550. doi:10.1093/biostatistics/kxj024 (abstract)
- Raftery AE, Chunn JL, Gerland P, Ševčíková H (2013). Bayesian probabilistic projections of life expectancy for all countries. Demography 50(3): 777-801. doi:10.1007/s13524-012-0193-x
- Rizzi S, Gampe J, Eilers PHC (2015). Efficient estimation of smooth distributions from coarsely grouped data. American Journal of Epidemiology 182(2): 138-147. doi:10.1093/aje/kwv020
- Schefzik R, Thorarinsdottir TL, Gneiting T (2013). Uncertainty quantification in complex simulation models using ensemble copula coupling. Statistical Science 28(4): 616-640. doi:10.1214/13-STS443
- Shang HL (2012). Point and interval forecasts of age-specific life expectancies: a model averaging approach. Demographic Research 27(21): 593-644. doi:10.4054/DemRes.2012.27.21
- Shang HL, Hyndman RJ (2017). Grouped functional time series forecasting: an application to age-specific mortality rates. Journal of Computational and Graphical Statistics 26(2): 330-343. doi:10.1080/10618600.2016.1237877
- Stoeldraijer L, van Duin C, van Wissen L, Janssen F (2013). Impact of different mortality forecasting methods and explicit assumptions on projected future life expectancy: the case of the Netherlands. Demographic Research 29(13): 323-354. doi:10.4054/DemRes.2013.29.13
- Thatcher AR, Kannisto V, Vaupel JW (1998). The Force of Mortality at Ages 80 to 120. Monographs on Population Aging 5. Odense University Press.
- Torri T, Vaupel JW (2012). Forecasting life expectancy in an international context. International Journal of Forecasting 28(2): 519-531. doi:10.1016/j.ijforecast.2011.01.009 (abstract)
- Van Schaeybroeck B, Vannitsem S (2015). Ensemble post-processing using member-by-member approaches: theoretical aspects. Quarterly Journal of the Royal Meteorological Society 141(688): 807-818. doi:10.1002/qj.2397 (abstract)
- Villegas AM, Haberman S (2014). On the modeling and forecasting of socioeconomic mortality differentials: an application to deprivation and mortality in England. North American Actuarial Journal 18(1): 168-193. doi:10.1080/10920277.2013.866034 (abstract)
- Villegas AM, Haberman S, Kaishev VK, Millossovich P (2017). A comparative study of two-population models for the assessment of basis risk in longevity hedges. ASTIN Bulletin 47(3): 631-679. doi:10.1017/asb.2017.18 (abstract)
- Wen J, Cairns AJG, Kleinow T (2021). Fitting multi-population mortality models to socio-economic groups. Annals of Actuarial Science 15(1): 144-172. doi:10.1017/S1748499520000184 (abstract)
- Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast reconciliation for hierarchical and grouped time series through trace minimization. Journal of the American Statistical Association 114(526): 804-819. doi:10.1080/01621459.2018.1448825
- Wilmoth JR (1995). Are mortality projections always more pessimistic when disaggregated by cause of death? Mathematical Population Studies 5(4): 293-319. doi:10.1080/08898489509525409 (abstract)
- Xu C, Xie Y (2023). Conformal prediction for time series. IEEE Transactions on Pattern Analysis and Machine Intelligence 45(10): 11575-11587. doi:10.1109/TPAMI.2023.3272339
- Yao Y, Vehtari A, Simpson D, Gelman A (2018). Using stacking to average Bayesian predictive distributions (with discussion). Bayesian Analysis 13(3): 917-1007. doi:10.1214/17-BA1091
