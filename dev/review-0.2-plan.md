# Plan: findings of the review of 0.2.0

Source: `dev/review-0.2.md` (independent review of the whole package,
2026-10-24). Each finding was checked again before it entered this plan.
The owner's principles apply (`dev/generalisation-plan.md`, Section 0):
correctness first, new methods opt-in, no options for problems that may
not occur, tests against the published method or exact identities, and
every work package updates NEWS, the progress log, the help and the
vignettes that use the changed code (`CLAUDE.md`).

## 1. Verification

| Finding | Severity | Checked by | Result | Decision |
|---|---|---|---|---|
| F01 `mt_parameters()` with a chosen period | high | reproduced on `EWMaleData`: 51 rows of `k`, recycling warning, 21 values in the fit | confirmed | fix (R1) |
| F02 SVD second stage without a root | medium | code: the unadjusted `k[t]` is kept without a warning | confirmed | warn and document (R1) |
| F03 `mt_coherence()` without aggregates | medium | reproduced: data.table error "Elements listed in `by`" | confirmed | classed error (R1) |
| F04 2D ungrouping and smoothing | medium | reviewer's timings on `mt_example()`: 166 s and 655 s | confirmed | faster algorithm (R3) |
| F05 `li_lee(ratio_order)`, Tucker ranks | medium | reproduced: "subscript out of bounds" on 6 age groups | confirmed | cap, warn, report (R1) |
| F06 CoDa bootstrap pool | low | code: `pool` keeps the first `d` residuals that `innov` drops | confirmed | fix (R2) |
| F07 normal approximation in recalibration | low | code | confirmed | document (R2) |
| F08 smoothing of the open-interval split | low | code | confirmed | own BIC choice or document (R2) |
| F09 summary of a chosen period | low | reviewer's output | confirmed | report the fitted period (R1) |
| F10 deprecated code | low | code | confirmed | owner decision (R4) |
| F11 unused internal functions | low | grep: `axis_from_bounds()`, `Matrix_diag()` and `relative_kernel()` are unused; `oldage_internal()`, `coda_multiway_point()`, `nordpred_projection()` and `saint_trend_rates()` are used by tests; `mt_inform()` is the messaging function `CLAUDE.md` prescribes | partly confirmed | delete the three unused ones (R2) |
| F12 seed of the Lee-Carter bootstrap | low | reproduced: two fits with `mt_forecast(seed = 1)` differ | confirmed | `seed = NULL` (R2) |
| F13 tests of the period choice | low | test files | confirmed | add tests (R1) |
| F14 floored cells in mean scores | low | reviewer's backtest | confirmed | owner decision (R4) |
| F15 `R CMD check` details | low | `tools:::.check_packages_used_in_tests()` finds no unstated dependency; the WARNING comes from the unreachable Bioconductor index of the container; `.git` is a directory in the primary checkout and is excluded | not reproduced | no change |
| F16 default `mt_example()` | low | reviewer's output | confirmed | README sentence (R2) |
| F17 cap of the CoDa zero replacement | low | code | confirmed | document (R2) |
| Function map without `mt_axis()` and `mt_couple()` | low | reviewer's scan | confirmed | add (R2) |
| `relative_model()` did not converge in one trial backtest (an aggregate stratum under `reconciled()`) | low | E8 vignette agent; the package warns with class `mt_warning_relative_model` | reported, not yet reproduced | reproduce and fix or document (R1) |

## 2. Work packages

### R1. Correct output and errors

- F01: `mt_parameters()` labels `k` with the years of the fitted period;
  `row()` stops on vectors of different lengths.
- F09: `summary()` and `print()` report the fitted period when
  `choose_period = TRUE` changes it; `mt_parameters()` gains `first_year`.
- F13: tests on `StMoMo::EWMaleData` (1961-2011, ages 0-89): the chosen
  first year equals `demography::lca(breakmethod = "bms")`, and the
  parameter table of the fit has one row of `k` per fitted year with the
  right labels.
- F03: `mt_coherence()` stops with class `mt_error_hierarchy` and names the
  cause when the draws hold no aggregate strata; a test checks the class.
- F05: `li_lee()` caps `ratio_order` by the number of ages and of years,
  `tucker3()` and `coda_multiway()` cap their ranks as now; all three warn
  with class `mt_warning_argument` when the cap binds and report the
  effective order or ranks in `mt_parameters()`. Tests: the warning class,
  the crash case of the review, and the reported effective values.
- F02: `lc_svd(adjust = "dt")` counts the years without a root of the
  total-deaths equation, records them in the diagnostics, warns with class
  `mt_warning_model_fit`, and the help of `adjust` says what happens.
  Test: data with a year whose equation has no root in reach.
- `relative_model()`: reproduce the non-converging fit (aggregate stratum
  under `reconciled()`, ages 50 and over); fix the iterations if the cause
  is in the code, otherwise document the case in the help.

### R2. Small corrections and labels

- F06: the bootstrap pool of the CoDa index drops the first `d` residuals,
  as `innov` does. Test: pool and innovations have the same length and the
  pool has mean zero. The vignettes that use the CoDa bootstrap are
  re-knitted if their numbers change.
- F12: `lee_carter()` gains `seed = NULL`, used with `bootstrap = TRUE`
  through `with_seed()`, as `bootstrapped()` does. Test: two fits with the
  same seed give identical bootstrap parameters.
- F07, F17: help sentences, labelled package-specific: the normal
  approximation of `recalibrated()` for cells with 50 or more expected
  deaths, and the cap `0.5 / z` of the CoDa zero replacement.
- F08: the split of the open interval chooses its smoothing by BIC on its
  own (it smooths exposures, not deaths); if this changes ungrouped
  results, the change goes into NEWS and the ungrouping vignette is
  re-knitted.
- F16: the README's first example says why `deaths = "single"` is passed.
- Function map: `mt_axis()` and `mt_couple()`.

### R3. Two-dimensional ungrouping and smoothing

- F04: the 2D PCLM of `mt_ungroup(dims = 2)` and the 2D P-spline of
  `mt_smooth(dims = 2)` use the array arithmetic of generalised linear
  array models (Currie, Durban and Eilers 2006), as `p_splines()` already
  does, with the composition matrix applied by year. The search of the
  smoothing parameters stays as it is.
- Tests: the new code gives the results of the current dense code within
  1e-8 on `mt_example()` and on a grid with an open interval, for fixed
  and for chosen smoothing parameters.
- Help: the run time per stratum on the package example.
- Reference (verified on Crossref): Currie ID, Durban M, Eilers PHC
  (2006). Generalized linear array models with applications to
  multidimensional smoothing. Journal of the Royal Statistical Society
  Series B 68(2): 259-280. doi:10.1111/j.1467-9868.2006.00543.x

### R4. Decisions of the owner

The owner decided on 2026-10-24:

- F10: remove the deprecated functions and options now, in the
  development version after 0.2.0 (next release 0.3.0): `product_ratio()`,
  `cause_shares()`, `ensemble(weights = "stack")`,
  `recalibrated(scale = "horizon")`, the spline methods of `mt_ungroup()`,
  and the deprecated methods of `mt_disaggregate()` and
  `mt_shape_rates()`, with a NEWS entry listing each replacement.
  `rwd(observation_noise = TRUE)` and `tucker3(transform = "coda")` were
  not part of this decision; on 2026-10-25 the owner decided to remove
  them as well, since nothing else uses them and published replacements
  exist (`arima_series(order = c(0, 1, 1), drift = TRUE)`,
  `coda_multiway()`).
- F11 (moved into R2): delete `axis_from_bounds()`, `Matrix_diag()` and
  `relative_kernel()`. Keep the test helpers and `mt_inform()`.
- F14: report the number of floored observations in `mt_score()` (column
  `n_floored`), with a help sentence on the derived quantities and coarser
  age groups.
- From E8, not a review finding: add the robustness measure of
  Stoeldraijer et al. (2018), the standard deviation of the changes of the
  forecast for one target year between successive fitting periods, as an
  opt-in measure of `mt_score()`, tested against a hand computation. The
  mean absolute revision stays, labelled package-specific.

## 3. Order

R1, then R3, then R2, then R4 (after the owner's decisions). Each work
package is one pull request with its tests, lint, the help, NEWS and the
progress log; R3 also re-knits the ungrouping and smoothing vignettes.
