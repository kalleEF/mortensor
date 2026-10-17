# Applied review of mortensor, round 4 (reviewer B)

* Persona: demographer and biostatistician; user of demography, StMoMo,
  bamp, DemoDecomp and LEdecomp through their documented interfaces; no
  package development background. Goal: forecasts for the German Laender
  by sex with national totals, causes of death at older ages, model choice
  by backtest, decomposition of forecast changes and gaps by age and cause,
  cross-checks against established packages, for a peer-reviewed paper.
* Version reviewed: mortensor 0.1.9000, branch
  `claude/vigilant-albattani-2ssw1x`, commit `c531ece`, loaded with
  `devtools::load_all()`. R 4.3.3, LEdecomp 1.0.14, DemoDecomp 1.14.1,
  demography 2.0.1, StMoMo 0.4.1, mgcv 1.9.1. Date: 17 October 2026.
* Material: README, vignettes `decomposition`, `regions`, `coherence`,
  `forecasting-hmd`, `models`, `faq`, `causes`; help pages of
  `mt_decompose()`, `anchored()`, `calibrated()`, `mt_reconcile()`,
  `mt_smooth()`, `lee_carter()`, `life_table()`; source of
  `R/derive-decompose.R`, `R/model-anchor.R`, `R/hierarchy.R`,
  `R/derive-lifetable.R`. Data: `ons_data`, `de_deaths`/`de_population`,
  `demography::fr.mort`, `StMoMo::EWMaleData`. At most 300 draws.
* Scripts and outputs: this directory (`s01_*.R` to `s23_*.R`, `*.out`).
  The previous round is `dev/review-applied-2.md`; its P1 items (crossed
  hierarchies, exploding random walks, model reference, parameters and
  residuals) are addressed in this version and are not repeated.

## 1. Summary verdict

The decomposition module is correct where it matters most: symmetric
Arriaga, Horiuchi and stepwise replacement reproduce LEdecomp and
DemoDecomp to machine precision on abridged and single ages, at `at > 0`
and with causes; the Kitagawa step and the radix and pooled composition
formulas are implemented as defined; draw-by-draw decomposition with
path-wise sums is the right design, and all numbers in the decomposition
vignette reproduce. Lee-Carter `adjust` and `choose_period` agree with
`demography::lca()`. Calibration claims reproduce.

For my paper, four problems would still produce wrong or indefensible
published numbers without an error:

1. the documented workflow for national cause-specific forecasts (regional
   vignette) implies that all-cause mortality of German women at 50+ rises
   by 2035, against an all-cause model that forecasts a 17% fall, and the
   vignette does not say so (P1.1);
2. `mt_decompose()` silently counts total levels (cause `"all"`, aggregate
   regions) as parts when `hierarchy` is omitted, although the data and the
   draws carry their hierarchy (P1.2);
3. decompositions of change from an observed year to a forecast year mix
   the jump-off step and observation noise into the age contributions; in
   the decomposition vignette 41% of the forecast narrowing of the
   deprivation gap is the jump-off, not the model's dynamics (P2.1);
4. `anchored()` shrinks every age to one offset per stratum, even in
   Bayern, because its empirical Bayes variance estimate is biased to zero
   (P2.2).

The radix composition of Riffe et al. is used in the vignette outside the
validity conditions stated by its authors (area-based groups, radix at
birth), and decompositions break on paths with death probabilities capped
at 1 (P2.3, P2.4). With these fixed, the package would support the
decomposition part of the paper; the cause-specific national forecasts
need the fixed-node reconciliation that the package already offers.

## 2. Findings

### P1: wrong results or blocking

**P1.1 National cause-specific forecasts imply rising all-cause mortality,
undiscussed.**
Location: `vignettes/regions.Rmd`, "Forecasts by cause" (lines 631-720),
table at line 671 and text from line 708.

Evidence (`s20_asr.R`, `s10_de_causes.R`): age-standardised death rate,
women 50+, Germany, per 100,000 (ESP 2013):

| | ASR 50+ | e50 |
|---|---|---|
| observed 1998 / 2010 / 2019 | 2643 / 2214 / 2037 | 2019: 34.65 |
| 2035, Li-Lee by sex and cause, bottom-up (median) | 2219 (1783-3123) | 34.65 |
| 2035, vignette ensemble by cause, bottom-up (median) | 2424 (2054-3375) | |
| 2035, Li-Lee on all causes, bottom-up (median) | 1692 (1339-2255) | 36.38 |

The cause-specific bottom-up forecasts, including the specification the
vignette recommends ("most accurate and best calibrated"), put the
national median all-cause rate in 2035 9-19% above 2019 and give no gain
in female e50 in 16 years, while the all-cause model forecasts a 17% fall
and +1.7 years. The driver is the mental and nervous group (ASR 223 in
2019, median 661 in 2035 under Li-Lee, 793 in the vignette), whose rise is
largely a coding shift. The text mentions wide intervals but not that the
median all-cause trajectory reverses. This is the known result of Wilmoth
(1995) that projections by cause are more pessimistic than all-cause
projections. The same Li-Lee forecasts make the East-West e50 gap of women grow
from 0.08 to 0.53 years by 2035; the mental and nervous group contributes
0.49 years of the 2035 gap (0.15 in 2019) (`mt_decompose(..., cause = "cause")`). A referee would
reject these numbers. The 5-year backtest cannot detect the problem.

Fix: in the vignette, reconcile the causes to the all-cause forecasts of
each Land (`reconciled(..., method = "wls", fixed = <labels of the
Land.sex.all strata>)`, already available), or use `cause_shares()`, and
show the all-cause comparison table above next to the cause-specific one;
cite Wilmoth (1995); state which assumption the national all-cause figure
rests on. A helper that reports the bottom-up all-cause forecast against a
separate all-cause model (and warns when their medians diverge by more
than, say, 5%) would catch this in other analyses.

**P1.2 `mt_decompose()` double counts totals when `hierarchy` is omitted.**
Location: `R/derive-decompose.R`, `mt_decompose()` (argument `hierarchy =
NULL`, line 151) and `decomp_populations()` (line 323).

Evidence (`s18_nohier.R`): ONS causes at 50+, Q5 against Q1, 2019, women:

```r
mt_decompose(xc, var = "quintile", ref = "Q1", cause = "cause")
#>  from_value 25.36  to_value 31.71  difference 6.35
mt_decompose(xc, var = "quintile", ref = "Q1", cause = "cause",
             hierarchy = mt_hierarchy(cause = list(total = "all", additive = "rates")))
#>  from_value 31.69  to_value 37.26  difference 5.57
```

Without the hierarchy, `"all"` is added as a thirteenth cause, rates are
doubled and e50 is wrong by 6 years, with no warning. The same happens
with `groups = "quintile"` on `mt_aggregate_strata()` output, where
"England" becomes a subgroup (radix e0 80.53 instead of 80.42). Yet the
aggregated data hold `$hierarchy` and reconciled draws hold
`fc$hierarchy` (`s19_coh.R`: `TRUE`), so the information is there.

Fix: default `hierarchy` to `x$hierarchy %||% y$hierarchy`; abort when a
level of `cause` or `groups` equals the sum of the other levels (the check
that `mt_validate()` already has for total levels among parts); add a
test.

### P2: important

**P2.1 Changes from an observed year to a forecast year include the
jump-off and observation noise.**
Location: `R/derive-decompose.R`, `decomp_period()` (observed years taken
from `draws$observed`); `vignettes/decomposition.Rmd` lines 190-215; help
section "Observed years inside forecast draws are the same in every draw".

Observed years are the raw rates, fixed in every draw; forecast years are
the model's expected (smooth) rates. The difference therefore contains the
step from observed to fitted rates at the jump-off and the sampling noise
of the observed year, both attributed to ages as if they were forecast
change.

Evidence (`s16_jumpoff_gap.R`, `s15_vig_fc.R`, ONS, Li-Lee by sex):

| Female Q5 - Q1 gap in e0 | total | ages 0-59 | 60-84 | 85+ |
|---|---|---|---|---|
| observed 2019 | 6.46 | 1.90 | 3.86 | 0.70 |
| fitted 2019 (Li-Lee) | 5.69 | 1.82 | 3.41 | 0.46 |
| forecast 2035 | 4.57 | 1.48 | 2.63 | 0.47 |

The vignette says "Li-Lee narrows the female gap from 6.5 to 4.6 years
mainly at ages 60 to 84, where its deviations of Q1 from the common trend
return towards their average". Of the 1.89-year narrowing, 0.77 years
(41%) is the step between observed and fitted 2019 rates; at 85+ the
whole narrowing is that step. With `li_lee(jump_off = "observed")` the
2035 gap is 5.34 instead of 4.57. `mt_jump_off()` shows a step of 1.07
years in e0 from observed 2019 to forecast 2020 for women in Q1.

German Laender (`s06_de_decomp.R`): of the forecast gain 2019-2035 of
Bremen women (2.32 years), 0.31 years is observed-to-fitted; for Bremen
men the 2019-2020 "change" of 0.46 years has absolute age contributions
summing to 0.93 years, which is mostly noise of the 2019 rates. The
intervals of these changes ignore that the observed e0 of Bremen women in
2019 has a Poisson 95% interval of 82.62-83.56 (width 0.94 years;
`s23_bremen_se.R`).

Fix: add `observed = c("raw", "fitted", "poisson")` (the fitted rates are
already stored in `draws$fitted`; Poisson draws exist in
`mt_observed()`), report the jump-off step in `summary(by = NULL)` and
message when it exceeds a threshold, and rewrite the vignette paragraph to
separate jump-off from forecast dynamics (or decompose from the first
forecast year, or from fitted 2019).

**P2.2 `anchored()` shrinkage removes all age variation.**
Location: `R/model-anchor.R`, `anchor_offsets()`, line 99:

```r
tau2 <- max(sum(w * (r[, s] - rbar)^2) / sum(w) - mean(v), 0)
```

The between-age variance is the deaths-weighted variance of the offsets
minus the unweighted mean of the sampling variances `v = 1/D`. With
weights `w = D` the expected weighted variance is `tau2 + n / sum(D)`, not
`tau2 + mean(1/D)`; `mean(1/D)` is dominated by the few deaths of young
ages and is 10-20 times too large. Result (`s13_anchor.R`, Lee-Carter,
last 3 years):

| Stratum | weighted var | `mean(1/D)` (code) | `n / sum(D)` (consistent) | tau2 code | tau2 consistent |
|---|---|---|---|---|---|
| Bayern female | 0.00085 | 0.00166 | 0.00009 | 0 | 0.00077 |
| Bremen female | 0.00335 | 0.03197 | 0.00151 | 0 | 0.00184 |
| Saarland female | 0.00151 | 0.02275 | 0.00086 | 0 | 0.00064 |

tau2 is 0 in all six strata tested, so every age gets the same offset: in
Bayern women the raw offsets range from -0.075 to +0.112 (70-74: +0.086,
based on thousands of deaths), all replaced by one number. The help says
the offset is "per age and stratum".

Fix: subtract `sum(w * v) / sum(w)` (equal to `n / sum(w)`), or use a
DerSimonian-Laird or REML estimate; add a test that offsets vary by age in
a large stratum; re-run the backtests quoted in `?anchored`. Also state
in the help that with `decay = 10` forecasts more than 10 years ahead are
unchanged (the 2035 gap above is identical with and without anchoring) and
that the fade itself adds an apparent trend in the fade years.

**P2.3 Radix composition used outside its validity conditions.**
Location: `?mt_decompose` section "Subgroups"; `vignettes/decomposition.Rmd`
lines 130-180.

The journal version of Riffe, Tursun-zade and Trias Llimós states that the
method is valid when subgroup weights are part of the initial conditions
of the combined life table and later group prevalence is determined only
by mortality. The vignette applies it at birth to IMD quintiles, which are
areas: people move between them and areas are reclassified with each IMD
version (the vignette mentions the latter only). The same call with
`groups = "region"` for Laender assumes no internal migration. On ONS data
the two variants give composition effects of opposite sign (radix -0.10,
pooled +0.10 years for women, 2001-2019; `s01_vignette_decomp.R`), so the
choice is not cosmetic. The vignette nonetheless calls the radix version
"the one of Riffe et al. and exact" and presents its numbers first.

Fix: state the condition in the help and the vignette; recommend
`composition = "pooled"` for area-based groups and regions, and radix for
groups fixed after the radix age (education from 30 or 35, country of
birth); use an example that meets the condition, or label the IMD example
as illustrative.

**P2.4 Paths with death probabilities capped at 1 break the decomposition.**
Location: `R/derive-lifetable.R` line 57 (`Q[Q > 1 | !is.finite(Q)] <- 1`)
with `arriaga_dir()` (line 537) and `cause_split()` (line 556).

When a path has `n m / (1 + n m / 2) > 1` at some age, `l_x` becomes 0
above it and Arriaga's `L_x / l_x` is `0/0`.

* Without causes (`s12_lx0.R`, one rate set to 2.5 at age 85 in one of 20
  draws): `sum_contributions` and `max_residual` are `NaN`, and
  `summary(d, by = "age")` fails with "non-numeric argument to binary
  operator" from `row_quantiles()`.
* With causes (`s10_de_causes.R`, separate Lee-Carter per cause, bottom-up,
  which `mt_forecast()` flags with rates above 1 in 51% of draws): the
  `NaN` sensitivities are silently replaced by the numerical derivative, so
  results look valid but are not additive: `max_residual` 5.6 years for
  Bremen women and 17.2 years for Bremen men (e50), and
  `summary(d, by = "cause")` returns plausible-looking contributions.

Fix: check `l_x > 0` and finite contributions per path; abort (or drop the
paths with a warning that counts them) and point to the implausible draws;
apply the numerical-derivative fallback only where `dM` is near zero, not
where `dl` is not finite; consider `q = 1 - exp(-n m)` when the
linear formula exceeds 1 (as LEdecomp does) so that `l_x` stays positive.

**P2.5 Inconsistent fixed nodes are moved without warning.**
Location: `R/hierarchy.R`, `fixed_nodes()` (line 371) and the help of
`mt_reconcile(fixed = )` (line 284: "they move only when every bottom
stratum below them would have to be negative").

Evidence (`s11_fixed.R`, ONS sex by quintile, Lee-Carter per node):
`fixed = c("both.England", "female.England", "male.England")`, whose base
forecasts do not add up (median gap 2.1% in 2035), returns without warning
and `max_relative_change_fixed = 0.335`: "fixed" nodes move by up to a
third in some draws. A user fixing the national and the Land all-cause
forecasts (the natural set-up for P1.1) would fall into this.

Fix: check that the fixed set is coherent in each draw (or contains no
node together with all its children), and warn when
`max_relative_change_fixed` exceeds a small tolerance; correct the help.

**P2.6 Non-consecutive years cannot be used.**
Location: period axis construction (`R/data-axis.R`), `mt_data()`,
`as_mt_data.demogdata()`.

Evidence (`s09_edges.R`, `s02_crosscheck.R`): `mt_data()` on ONS rows for
2001, 2011 and 2019 builds the period grid `2001, 2002-2010, 2011,
2012-2018, 2019`; `as_mt_data(extract.years(fr.mort, c(1950, 2000)))`
builds `1950, 1951-1999, 2000`. Every later call fails with "No deaths in
2002-2010 ... restrict the years with mt_subset(periods = )", a remedy
that cannot work. Comparing two published years (two census-based life
tables, or a 1990 and a 2019 table) is a common use of decomposition and
is impossible.

Fix: let the period axis have gaps (one interval per observed year) or
add an argument for it; at least give a correct error.

### P3: minor

* **a0 and life-table conventions.** `a_x = n/2` at age 0
  (`lt_columns()`), and `mt_decompose()` has no `ax` or `open` arguments
  although `life_expectancy()` accepts them. e0 of French women 1950:
  69.156 (mortensor) against 69.188 (`demography::lifetable`) and 69.183
  (LEdecomp, Andreev-Kingkade a0). The NEWS claim of agreement with
  LEdecomp "to machine precision" holds only without an age-0 row (the
  test uses `at = 30`): with age 0 the age-0 contribution differs by 4e-4
  and other ages by up to 4e-6. Add an a0 rule (Andreev-Kingkade, as in
  LEdecomp) and pass `ax`/`open` through. `lee_carter(adjust = "e0")`
  differs from `demography::lca(adjust = "e0")` by up to 2.8% in fitted
  rates for the same reason (`s21_lca_fit.R`); `"dt"` and `"dxt"` agree to
  0.09% (rounding of reconstructed deaths). Parameters differ by a
  constant shift of `k` (0.32) because mortensor re-centres `k` after the
  second stage and demography does not; say so in `?lee_carter`.
* **Stepwise order with causes.** Mortensor replaces age by age (all
  causes of an age, then the next age); LEdecomp and DemoDecomp replace
  cause by cause. By-cause totals differ by up to 0.008 years
  (`s03_causes.R`). Document the order and offer an argument.
* **Forecast composition is zero by default.** With `groups` and forecast
  years, composition comes from the last observed exposures, so the
  composition row is exactly 0 (`s06_de_decomp.R`). Print a note
  suggesting `exposure` (a population projection); with a projection the
  effect is estimated as expected (`s18_nohier.R`: -0.20 years).
* **Runtime of pooled composition with causes.** 16 Laender, 8 causes, 9
  ages, 200 draws: 2.0 minutes at `n_steps = 5` (`max_residual` 0.024);
  at the default 20 steps and 1000 draws about 40 minutes. Mention in the
  help, or use the exact stepwise route by default for this case.
* **References.** `?mt_decompose` and the vignette cite only the working
  paper of Riffe et al.; the journal version exists (Geographical Analysis
  58(2), 2026, doi:10.1111/gean.70033, listed in
  `dev/decomposition-plan.md`). Say that equation numbers refer to the
  working paper. "Vaupel 2002 calls these 'current conditions' and
  'current rates'" (`decomposition.Rmd` line 153) stretches the paper,
  whose distinction concerns tempo effects of changing mortality (and, at
  most, heterogeneity), not radix versus pooled weighting of subgroup life
  tables; cite it as an analogy or replace it (for example by Vaupel and
  Yashin 1985 on heterogeneity). I could verify the bibliographic details
  of all decomposition references but not the equation numbers (the
  working paper PDF was not reachable).
* **`forecasting-hmd.Rmd`.** Line 154 announces "Four modifications" and
  lists three. Line 192, "The coverage of life expectancy is at or above
  95% for all four", is false for plain Lee-Carter (0.917 in my run with
  the vignette's settings, `s14b_ew.R`). Booth, Maindonald and Smith
  (2002) is cited in the text but missing from the references.
* **Smoothing plus calibration: the gain is mostly calibration.** Ablation
  on England & Wales men (`s14b_ew.R`, `s22_gam.R`; origins 1996-2001,
  h = 10, 200 draws):

  | Lee-Carter variant | CRPS | cov95 | CRPS e30 |
  |---|---|---|---|
  | plain | 0.0687 | 0.746 | 0.473 |
  | Whittaker over time only | 0.070 | 0.694 | 0.501 |
  | 2D GAM only | 0.0742 | 0.586 | |
  | calibrated only | 0.0657 | 0.877 | 0.442 |
  | Whittaker over time + calibrated | 0.066 | 0.873 | 0.442 |
  | 2D GAM + calibrated | 0.0646 | 0.914 | |

  The documented numbers reproduce, but calibration alone gives most of
  the gain, Whittaker smoothing over time adds nothing, and the GAM adds
  0.0011 in CRPS on one population with six overlapping origins. The FAQ,
  NEWS and `?mt_smooth` credit "smoothing with calibration"; show the
  ablation and present smoothing over time as optional.
* **Backtest labels ignore preprocessing.** `mt_backtest()` labels a model
  by its specification only, so Lee-Carter with and without `prep` are both
  "Lee-Carter" and are pooled by `mt_score(by = "model")`. Add the
  preprocessing to the label or a `prep` column.
* **HMD data to 110+.** `mt_decompose()` and `mt_fit()` stop with
  "Exposures must be positive on the grid ... male, age 109, 1990" on
  `as_mt_data(demography::fr.mort)`; suggest `mt_regroup(x, age = 0:100)`
  in the message.
* **Previous P2.1 only partly done.** `mt_score()` still has no standard
  errors, paired comparisons, weights or age ranges; the regional vignette
  still ranks models on CRPS differences in the third decimal
  (0.088 against 0.089).

## 3. What works well

* **Exactness and agreement with reference implementations.** Symmetric
  Arriaga, Horiuchi (20 steps) and symmetric stepwise replacement agree
  with LEdecomp `arriaga_sym()`, DemoDecomp `horiuchi()` and
  `stepwise_replacement()` to 1e-14 or better on abridged ages from age 1,
  at `at = 65`, and on HMD single ages (1e-6, caused only by a0). With
  causes, Arriaga and Horiuchi agree with `LEdecomp::LEdecomp(...,
  n_causes = )` to 2.5e-16 and 5e-14. The identity "symmetric Arriaga =
  symmetric stepwise from the youngest age" holds exactly, as the help
  explains.
* **Cause split.** The proportional (implied-sensitivity) split is
  numerically stable on real data: in 59 of 180 age cells of ONS
  2018-2019 changes the absolute cause terms exceed three times the age
  total (where LEdecomp warns), yet the cause totals agree with Horiuchi
  to 1.4e-5 years.
* **Kitagawa, radix and pooled.** The Kitagawa rate and composition terms
  match a hand calculation (female Q1: -138.1 and -56.4 per 100,000). The
  radix version reproduces the weighted mean of subgroup e0; the pooled
  version reproduces e0 of the summed deaths and exposures; summing
  composition terms over subgroups by age is the identified quantity.
* **Draw-by-draw design.** Contributions and their sums have intervals
  formed path by path; the East-West gap by age and cause with intervals,
  observed against forecast years, two models against each other, and a
  forecast against hold-out data (`mt_decompose(fc, observed, periods =
  2006)`) all work; Arriaga with causes for 342 strata and 200 draws takes
  under a second.
* **Documentation reproduces.** All numbers in the decomposition vignette
  that I checked reproduce (gap 6.5 and 7.7 years; 5.2 and 6.3 at 40-84;
  England 2.7 and 3.7; pooled 0.2-0.5 higher; Q1 women +3.5 years with 2.5
  at 50-84; gap 6.5 to 4.6; Lee-Carter 1.9 below Li-Lee; Kitagawa 140-360
  and the sign statements).
* **Errors.** Missing (suppressed) cells, wrong `at`, wrong `ref` and
  method-quantity mismatches give clear errors with remedies.
* **Lee-Carter options.** `choose_period = TRUE` selects the same first
  year as `demography::lca(chooseperiod = TRUE, breakmethod = "bms")`
  (1968 for French women, 1979 for men).
* **Calibration, coherence, reconciliation.** `calibrated()` claims
  reproduce exactly (0.0646, 0.914); backtests score e0 from predictive
  draws, which is consistent; a consistent fixed total is kept to 1e-9
  with WLS, MinT and OLS; non-negativity works (share set to zero 1e-4);
  `mt_coherence()` reports exact coherence for bottom-up draws and the
  expected spread for separate models.
* **The previous round's P1 items are fixed**: crossed region by cause
  hierarchies aggregate correctly, `rwd(observation_noise = TRUE)` and
  `li_lee()` give plausible e0 in Bremen, and the model reference and
  parameter accessors exist.

## 4. Cross-checks performed

| # | Check | Script | Result |
|---|---|---|---|
| 1 | Decomposition vignette numbers, observed part | `s01` | all reproduce |
| 2 | Arriaga sym. vs LEdecomp, abridged from age 0 | `s02` | age 0 differs by 4e-4 (a0), other ages by up to 4e-6 |
| 3 | Same from age 1, and at 65 | `s02` | 1.7e-15, 1.6e-15 |
| 4 | Horiuchi and stepwise vs DemoDecomp (age 1+) | `s02` | 5.7e-14, 7.1e-15 |
| 5 | e0 France 1950/2000 vs demography and LEdecomp | `s02` | 1950: 0.03 years lower (a0); 2000: within 0.007 |
| 6 | Causes (12 groups, 50+): Arriaga, Horiuchi, stepwise vs LEdecomp | `s03` | 2.5e-16, 5e-14, stepwise 0.008 (order) |
| 7 | Stability of cause split vs Horiuchi | `s04` | max difference 1.4e-5 years |
| 8 | German Laender: gap, change, radix/pooled, jump-off | `s05`, `s06` | runs; composition 0 in forecasts; jump-off 0.3 years in Bremen |
| 9 | `lee_carter(adjust, choose_period)` vs `demography::lca` | `s07`, `s21` | period identical; dt/dxt fitted 0.09%; e0 2.8% |
| 10 | Edge cases: gaps in years, suppressed cells | `s09` | gaps fail (P2.6); suppressed: clear error |
| 11 | Explosive draws, causes by Land | `s10`, `s12` | NaN or non-additive (P2.4) |
| 12 | Fixed nodes, consistent and inconsistent | `s11` | 1e-9; 0.335 silently (P2.5) |
| 13 | `anchored()` offsets and shrinkage | `s13` | tau2 = 0 everywhere (P2.2) |
| 14 | E&W backtest variants, ablation | `s14b`, `s22` | numbers reproduce; ablation in P3 |
| 15 | Forecast part of decomposition vignette, Kitagawa by hand | `s15` | reproduce |
| 16 | Share of gap narrowing due to jump-off | `s16` | 0.77 of 1.89 years (P2.1) |
| 17 | HMD single ages, hold-out vs forecast | `s17` | 1e-6 vs LEdecomp; hold-out works |
| 18 | Omitted hierarchy; projected exposures | `s18` | silent double counting (P1.2); projection works |
| 19 | `mt_coherence()`; help example | `s19` | as documented |
| 20 | All-cause vs cause bottom-up for Germany | `s20` | P1.1 |
| 21 | Poisson interval of observed e0, small Laender | `s23` | Bremen women width 0.94 years |
| 22 | References (web search) | | Riffe et al. WP and journal version exist; Vaupel 2002 bibliographic details correct |

## 5. What a paper reviewer would ask that the package cannot yet provide

1. Uncertainty of score differences: standard errors over origins and
   strata, paired tests (Diebold-Mariano type) and deaths-weighted scores,
   before claiming that one model is "most accurate" on third-decimal
   CRPS differences.
2. A strategy for the pandemic years: `de_deaths` runs to 2024 but all
   vignettes stop in 2019; a 2026 paper must forecast from 2024 and needs
   year exclusion or down-weighting and a jump-off after 2022.
3. Decomposition of forecast uncertainty by source (common trend, Land
   deviations, Poisson noise, ensemble component), and of differences
   between models by parameter (jump-off, drift, age pattern), not only by
   age.
4. Lifespan variation (e-dagger, Gini) as built-in quantities with their
   decomposition, and continuous-change decomposition over the forecast
   path (Vaupel and Canudas-Romo 2003).
5. Sensitivity of e0 and decompositions to life-table conventions (a0,
   open-age closure), reported in one call.
6. Uncertainty from the imputation of suppressed cells (several
   imputations propagated into forecasts), which matters for rare causes
   in small Laender.
7. Handling of cause-coding breaks (dementia between F and G chapters):
   bridge coding or comparability ratios before extrapolation, not only a
   warning in the text.
8. An adapter for official population projections (Destatis variants) as
   `exposure`, so that composition effects and deaths in forecast years
   are not frozen at 2019.

## References added by this review

Riffe T, Tursun-zade R, Trias Llimós S (2026). Arriaga meets Kitagawa:
life expectancy decomposition with population subgroups. *Geographical
Analysis* 58(2). doi:10.1111/gean.70033.

Vaupel JW, Yashin AI (1985). Heterogeneity's ruses: some surprising effects
of selection on population dynamics. *The American Statistician* 39(3):
176-185.

Wilmoth JR (1995). Are mortality projections always more pessimistic when
disaggregated by cause of death? *Mathematical Population Studies* 5(4):
293-319.
