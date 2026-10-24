# Review of mortensor 0.2.0

Independent review of the whole package at the state of pull request #15
(commit 2eec7aa, merged as 0.2.0), 2026-10-24. The reviewer read the code
against the sources cited in "Methods and sources", ran targeted scripts on
the package data and the HMD extracts of demography and StMoMo, and ran
`R CMD check`. The scripts and their output stayed outside the repository.
The plan that follows from the review, with the verification of each
finding, is `dev/review-0.2-plan.md`.

## Summary

17 findings: 0 critical, 1 high, 4 medium, 12 low.

The methods checked against their sources are implemented as published or
with documented deviations. The high finding is a wrong table from
`mt_parameters()` for `lee_carter(choose_period = TRUE)`. The medium
findings are a crash and silent reductions of requested orders, a silent
fallback in the SVD second stage, a wrong error from `mt_coherence()`, and
two-dimensional ungrouping and smoothing that take minutes on the package
example. The low findings are documentation gaps, dead code, deprecated
code that can go, a test gap, the scoring of floored cells and two
`R CMD check` details.

Read only, not checked equation by equation: `stmomo()`, `stanmomo()`,
`apc_bayes()`, `hyndman_ullah()`, the SAINT estimators beyond the
documented misprints, the P-spline forecast extension, the anchored
shrinkage and the Lee-Miller offset.

## Findings

### F01. `mt_parameters()` recycles and mislabels `k` for a chosen period (high, confirmed)

`R/model-inspect.R` (`yrs <- as.character(fit$periods)` and
`row(lab[s], "k", yrs, f$k)`). With `choose_period = TRUE` the fit holds
`k` for the chosen years only, while `fit$periods` holds all years.
data.table recycles the shorter vector with a warning. On StMoMo's
`EWMaleData` (ages 0-88, 1961-2011) the table reports 51 values of `k`
labelled 1961 to 2011, where the 21 values of 1991 to 2011 repeat.

Fix: index the years from the stored first year of the fit; a length check
in `row()` turns any future recycling into an error; a test (F13).

### F02. `lc_svd(adjust = "dt")` keeps the SVD index silently when no root is found (medium, confirmed from code)

`R/model-lee_carter.R`, second stage of Lee and Carter (1992): on a
`uniroot()` error the unadjusted `k[t]` is kept. A code comment states the
choice; no warning or diagnostics entry tells the user that the fitted
total deaths of that year differ from the observed total. The Poisson
fallbacks of the same file warn.

Fix: count the years without a root, add them to the diagnostics, warn
with class `mt_warning_model_fit`, and document it under `adjust`.

### F03. `mt_coherence()` on draws without aggregate strata fails with an internal error (medium, confirmed)

`R/hierarchy.R`: with no aggregate strata the rows are empty and
`merge(..., by = "stratum")` fails with a data.table error.

Fix: `mt_abort()` with class `mt_error_hierarchy` naming the cause
(aggregate the data with `mt_aggregate_strata()` or fit `reconciled()`).

### F04. Two-dimensional ungrouping and smoothing take minutes on the package example (medium, confirmed)

`R/harmonise-ungroup.R` (dense tensor basis `kronecker(bp$B, fb$bas$B)`
and dense penalties, one PCLM per stratum over all years),
`R/harmonise-pclm.R` (`pclm_fit()`: two rounds of two golden-section
searches, each with about 15 IRLS fits of the dense system) and
`R/harmonise-smooth.R`. On `mt_example()` (2 strata, 13 age groups, 30
years, 71 single ages): `mt_ungroup()` 6.8 s, `mt_ungroup(dims = 2)`
166 s, `mt_smooth()` 4.7 s, `mt_smooth(dims = 2)` 655 s.
`step_ungroup()` and `step_smooth()` with `dims = 2` run once per
backtest origin.

Fix: the array arithmetic of generalised linear array models (as in
`p_splines()`) for the 2D fits, which avoids the Kronecker products; state
the run time in the help.

### F05. `li_lee(ratio_order)` crashes above the number of ages and is reduced silently otherwise; Tucker ranks likewise (medium, confirmed)

`R/model-multipop.R` (`r <- min(p$ratio_order, ncol(D) - 2)`),
`R/model-tucker3.R` and `R/model-coda-multiway.R` (ranks). `li_lee(ratio_order
= 10)` on 6 age groups and 10 years stops with "subscript out of bounds";
on 61 ages the order is cut to 8 without a warning and the label stays
"Li-Lee". `tucker3(ranks = c(5, 5, 5))` on 6 x 10 x 20 is reduced without
a warning.

Fix: cap by the number of ages as well, warn with class
`mt_warning_argument` whenever the cap binds, and report the effective
order or ranks in `mt_parameters()` and the label.

### F06. The bootstrap pool of the CoDa index includes the diffuse residuals of `stats::arima()` (low, confirmed)

`R/model-coda-core.R`: `innov` drops the first `d` residuals of a model
with `d` differences, `pool` keeps them. They come from the diffuse
initialisation and are near zero, so the resampled spread is slightly too
small (ARIMA(0,1,1) with drift, 30 years: sd 0.914 with them, 0.930
without).

Fix: `pool = res[(d + 1):n] - mean(res[(d + 1):n])`.

### F07. The normal approximation in the calibration sample is undocumented (low, confirmed from code)

`R/model-recalibrate.R`: cells with 50 or more expected deaths use the
normal approximation of the Poisson quantiles; the help does not say so.

Fix: one sentence in the help, labelled package-specific.

### F08. The open-interval PCLM reuses the smoothing settings of the death ungrouping (low, confirmed from code)

`R/harmonise-ungroup.R` (`open_weights_pclm()`): the split of the open
interval passes the user's `lambda` for the deaths; with `dims = 2` its
first element is taken. The help does not mention the coupling.

Fix: choose the smoothing of the open-interval split by BIC on its own, or
document the coupling.

### F09. `summary()` and `print()` of a chosen-period fit do not report the chosen period (low, confirmed)

`R/model-inspect.R`: the period shown is the range of all periods.

Fix: report the fitted period when it differs, and a row `first_year` in
`mt_parameters()`.

### F10. Deprecated functions and options can be removed (low, confirmed)

`product_ratio()`, `cause_shares()`, `ensemble(weights = "stack")`,
`recalibrated(scale = "horizon")`, the spline methods of `mt_ungroup()`
(`"hyman_kannisto"`, `"natural_gompertz"`, `"gam"`), and the deprecated
methods of `mt_disaggregate()` and `mt_shape_rates()` keep their code,
tests and help. The warnings already name the published replacements.

Fix: remove them in a later version with a NEWS entry.

### F11. Internal functions that nothing calls (low, suspected)

A static scan names `axis_from_bounds()`, `oldage_internal()`,
`Matrix_diag()`, `coda_multiway_point()`, `nordpred_projection()`,
`relative_kernel()`, `saint_trend_rates()` and `mt_inform()`;
`simulate_lc_dgp()` and `pkg_data()` are used in tests only.

Fix: delete the unused ones or move test helpers to
`tests/testthat/helper-*.R`.

### F12. `lee_carter(bootstrap = TRUE)` has no `seed`, while `bootstrapped()` has one (low, confirmed)

Two fits of `lee_carter(bootstrap = TRUE, n_boot = 5)` followed by
`mt_forecast(seed = 1)` give different draws. `CLAUDE.md` asks for
`seed = NULL` with `with_seed()` in functions that simulate.

Fix: `seed = NULL` in `lee_carter()` for the bootstrap, or document the
difference in both help pages.

### F13. No test of the period choice against `demography::lca(breakmethod = "bms")`, none of `mt_parameters()` with a chosen period (low, confirmed)

The help claims the same period as demography. The reviewer found the
mean-deviance ratios equal to four decimals and the same choice (1991) on
`EWMaleData`, 1961-2011, ages 0-89. A parameter test on that fit would
have caught F01.

Fix: add both checks to the tests.

### F14. Zero-death cells dominate the mean scores of a backtest (low, confirmed)

`R/eval-backtest.R`: observed and predicted log rates are floored at 0.5
deaths, as documented. With age 30 set to zero deaths from 2015, the
floored cells had a CRPS of 4.73 against 0.03 to 0.05 elsewhere; one such
cell outweighs about 100 ordinary cells in the mean.

Fix: report the number of floored observations in `mt_score()`, or score
sparse cells on deaths.

### F15. `R CMD check` details (low)

A NOTE on the hidden file `.git` (in a git worktree `.git` is a file that
`R CMD build` keeps), and a WARNING "checking for unstated dependencies in
'tests'", attributed to `MortalitySmooth::Mort2Dsmooth()` in a test.

### F16. The default `mt_example()` cannot be fitted as is (low, confirmed)

`mt_fit(lee_carter(), mt_example())` stops with a clear message, because
the default deaths are grouped and the exposures single-year. The help
examples pass `deaths = "single"`.

Fix: one sentence in the README's first example.

### F17. The CoDa zero replacement caps `delta` beyond the paper without saying so (low, confirmed from code)

`R/model-coda-core.R`: Bergeron-Boucher et al. (2017, Appendix C) define
`delta_t = min(D) / 2 / sum_x D_tx`; the package adds the cap `0.5 / z`.
The help describes only the paper's rule.

Fix: document the cap as a package safeguard.

## Checks that passed

Area A (methods): Poisson Lee-Carter equal to StMoMo (parameters within
1e-6, French women 1980-2005); SVD Lee-Carter equal to `demography::lca()`
up to the handling of zero rates; `sum(b) = 1`, `sum(k) = 0`; the averaged
jump-off of Stoeldraijer et al. (2018); the constrained fit and bootstrap
of Li, Hardy and Tan (2009, Appendix 2); the negative binomial score and
weights of Delwarde et al. (2007); the state-space jump-off of De Jong and
Tickle (2006); the period choice of Booth, Maindonald and Smith (ratios
equal to demography to four decimals). CoDa: the bootstrap of Appendix B,
the zero replacement of Appendix C, the common factor, the drift
correction (eq. 20) and the centring of the authors' code
(Bergeron-Boucher et al. 2017); `coda_causes()` sums to the all-cause draws
(relative difference 2e-16); the 2S-CoDA weights of Kjærgaard et al.
(2019). `cae(index = "coherent")` is Model 3 of Enchev et al. (2017,
Table 10). APCI constraints and cohort weights of Richards et al. (2019).
CP-spline penalties of Camarda (2019). SAINT recursion and EM counts with
the documented misprint corrections. Evaluation: CRPS, WIS with the
weights of Bracher et al., randomised PIT, the HLN correction, the EWC and
Daniell variance estimators of Coroneo and Iacone (2020); backtests refit
on the training window only. Derived quantities: life tables (e65 equal to
`demography::lifetable()` within 1e-5), all ten closures of the open
interval, Wilmoth's m-ratio equal to a hand computation of eq. (12).
Hierarchies: every reconciliation method coherent to machine precision,
fixed nodes moved by at most 2e-12, predictive draws of reconciled models
coherent. Harmonisation: `mt_ungroup()` recovers single-year deaths on
`mt_example()` within 4% at ages 30-89 and 2% at 90+. Randomness: the same
seed reproduces forecasts, backtests (also with `cores = 2`) and
reconciled draws.

Area B (edge cases): short series, one age group, missing cells, age
groups without deaths, `h = 1`, `n_draws = 1`, predictive draws,
single-stratum joint models and nested wrappers either run or stop with a
clear, classed error.

Area D (documentation): `pkgdown::check_pkgdown()` finds no problems;
every export is in the function map except `mt_axis()` and `mt_couple()`;
no em-dashes in `R/`, the vignettes, README or NEWS.
