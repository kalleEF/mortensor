# Plan: ONS deprivation data in the vignettes, and comparison with fdm and bamp

Status: draft for discussion (2026-10-13). Nothing below is implemented yet
except the reader `mt_read_ons()`, which already reads the files (see 2.3
for the changes they need).

## 1. The data

Two workbooks of an ONS user-requested tabulation, "Death Registrations and
Populations by Index of Multiple Deprivation, England":

| | Part 1 | Part 2 |
|---|---|---|
| Deaths (sheet 1) | 2001-2024 | 2021-2024 |
| Populations (sheet 2) | 2001-2020 | 2021-2024 |
| IMD version | 2004, 2007, 2010, 2015, 2019 (changes in 2006, 2009, 2013, 2018) | 2025 |
| Layout | 2 sexes x 20 age groups (`<1`, `01-04`, 5-year groups, `90+`) x 10 IMD deciles x 30 cause columns plus all causes | same |

Checks already run:

* Populations match the published England totals (49.45 million in 2001,
  56.55 million in 2020 and 2021, 58.62 million in 2024). Deaths per year
  range from 452,863 (2011) to 569,700 (2020).
* No all-cause cell at ages 40+ has zero deaths (minimum 40); younger ages
  and single causes have small counts.
* Deaths 2021-2024 appear in both parts with different IMD versions. Per
  decile the counts differ by up to 4% (decile 3: +2.5% to +4.0% with IMD
  2025; decile 1: -0.6% to -1.8%). Deaths and populations are consistent
  within a year only when both come from the same part.
* The cause columns overlap: heart failure includes I11.0, I13.0 and I13.2,
  which are also in hypertension (I10-I15); alcohol-specific deaths include
  K70 (alcohol-related liver disease) and X65 (also in suicide). The listed
  columns add up to 70% of all deaths.
* The COVID-19 column is headed "Deaths involving COVID-19", while the
  title of Table 1 says "underlying cause" and the ICD-10 lookup lists the
  codes U07.1, U07.2, U09.9 and U10.9. The column has 69,299 deaths in 2020.
  Whether it counts underlying cause or any mention on the certificate is
  to be confirmed with the colleagues who requested the table.
* Neonatal deaths have no cause code and appear only in "All causes".
* Licence: the cover sheet states Crown copyright 2025, re-use under the
  Open Government Licence, with the accreditation "Source: Office for
  National Statistics licensed under the Open Government Licence".

## 2. Data handling in the package

### 2.1 Series construction across IMD versions (decision needed)

* **A (recommended):** each year from the part that has both deaths and
  populations for it: 2001-2020 from part 1, 2021-2024 from part 2. Rates
  are consistent within every year; the switch to IMD 2025 in 2021 is a
  break of the same kind as the four earlier version changes.
* B: part 1 only, 2001-2020. Consistent, but loses the four most recent
  years and the post-pandemic period.
* C: part 1 deaths for 2021-2024 with part 2 populations. Mixes IMD 2019
  numerators with IMD 2025 denominators; I would avoid it.

The 2021-2024 overlap can serve as a short sensitivity analysis in a
vignette (how much do decile-specific rates move when only the IMD version
changes).

### 2.2 Storage (decision needed)

* **Recommended:** a compact processed dataset in the package
  (`data/ons_imd.rda`: year, sex, age group, decile, deaths for the cause
  groups of 2.4, population; 9,600 rows, roughly 100-200 KB), documented
  with source, licence and accreditation. A script `data-raw/ons_imd.R`
  rebuilds it from the two workbooks with `mt_read_ons()`. The workbooks
  themselves (2.2 MB) go into `data-raw/` (excluded from the built package)
  or stay outside the repository with a download link, if ONS publishes the
  tabulation at a stable URL.
* The vignettes then load `ons_imd` without network access or Excel files,
  which keeps CI and CRAN checks simple.

### 2.3 Reader changes

* Accept several deaths files; for every year take deaths and populations
  from the same file (rule A), with an explicit `years` mapping as an
  override.
* Record the IMD version per year (sheet "IMD") as metadata of the result.
* Keep all-cause deaths for every age; restrict cause groups to the ages
  used for cause models (no cause for neonatal deaths).

### 2.4 Cause groups (decision needed)

Disjoint groups covering the main causes, with shares of all deaths
2001-2019:

| Group | Columns | Share |
|---|---|---|
| cancer | breast, colorectal, lung, prostate, all other cancers (C00-C97, D05) | 27.2% |
| circulatory | chronic ischaemic heart disease, stroke incl. TIA | 20.7% |
| dementia | dementia | 7.3% |
| respiratory | COPD, asthma | 5.5% |
| rest | all causes minus the above (includes COVID-19) | 39.3% |

Heart failure, hypertension, alcohol-specific deaths and suicide stay
outside the groups because they overlap; they can be shown as separate
single-cause analyses. COVID-19 stays in "rest": a group with zero deaths
before 2020 has no log-ratio coordinate.

### 2.5 The pandemic years (decision needed)

* Forecasts with jump-off 2024 and 2020-2022 marked as shock years
  (`step_shocks()`), compared with forecasts trained through 2019.
* Backtests with origins 2008-2014 and horizon 5 (scored before 2020), plus
  one stress test from origin 2019 scored on 2020-2024.

## 3. Vignettes

ONS data become the main worked example. The HMD vignettes on ungrouping
stay, because the accuracy of ungrouping can only be measured where true
single-year rates exist (England & Wales, France).

1. **Getting started with ONS deprivation data** (end to end): read and
   validate, granularity, ungrouping to single ages (PCLM and the spline
   alternatives), Lee-Carter and Li-Lee by sex and quintile, forecasts,
   life expectancy at birth and at 65 by quintile, the gap between the most
   and least deprived quintile with intervals, probability of dying between
   40 and 75, age-standardised rates (ESP 2013), export.
2. **Coherent forecasts by sex and deprivation**: hierarchy sex x decile,
   deciles grouped into quintiles as an intermediate level, plus totals
   (`mt_hierarchy(sex = "both", decile = list(total = "all", groups =
   quintiles))`). Base, reconciled and hierarchical Li-Lee forecasts with
   backtests at every level. The simulated regional vignette stays as a
   short example for regions, since the ONS data have none.
3. **Causes of death by deprivation**: `cause_shares()` on the groups of
   2.4, coherent with all-cause forecasts, `mt_plot_shares()`.
4. **Model choice and evaluation**: backtests of `rwd()`, `lee_carter()`,
   `li_lee()`, `tucker3(ranks = "bic")`, `hyndman_ullah()`, `apc_bayes()`
   (with `apc_select()`), `mt_bootstrap()` and ensembles; scores by horizon
   and age (`mt_plot_scores()`), coverage, and the pandemic stress test.
5. **Comparison with the original packages** (section 4).

Expensive parts (bamp chains, many backtests) are precomputed by
`dev/precompute-vignettes.R` as for the current bamp vignette; the vignettes
read the stored results and state how to rerun them.

## 4. Comparison with fdm (demography) and bamp

Two questions, kept apart because they need different designs.

### 4.1 Fidelity: do the wrappers reproduce the original packages?

Same input, same settings, original call next to the mortensor call:

* **fdm:** `demography::fdm()` and `forecast.fdm()` on a `demogdata` object
  against `hyndman_ullah()`. Common input: single-year rates for one sex and
  quintile after PCLM ungrouping (fdm smooths over single ages). Compare
  point forecasts of log rates (maximum absolute difference) and the 80% and
  95% bounds. demography computes intervals analytically, mortensor from
  simulated paths; differences are reported together with their Monte Carlo
  error.
* **bamp:** `bamp::bamp()` and `bamp::predict_apc()` on the 5-year groups
  against `apc_bayes()` with the same priors and MCMC settings. Posterior
  medians and interval bounds should agree within MCMC error; differences
  are reported relative to the posterior standard deviation.
* **Lee-Carter:** `StMoMo::lc()` and `demography::lca()` against
  `lee_carter()` (the unit tests already show agreement to 6e-5 and 1e-6);
  shown as a table for completeness.

### 4.2 Added value: how do the workflows compare on forecast accuracy?

Rolling-origin backtests in which each package is used as its documentation
recommends (fdm with its own smoothing on single ages; bamp on 5-year groups
with its own prediction), next to mortensor pipelines (ungrouping, predictive
draws, reconciliation, ensembles). All forecasts are scored on the observed
5-year-group deaths with CRPS, WIS and coverage, so the comparison does not
depend on any package's grid. Points to show explicitly:

* intervals for the smooth rate (fdm, bamp without overdispersion) against
  predictive intervals for observed rates (Poisson and extra-Poisson
  variation), which explains most coverage differences;
* runtime;
* cases where an original package forecasts better are reported as such.

Optionally the stored mortensor 0.1 baseline (`dev/baseline/`) as a third
column, if the colleagues' earlier results used 0.1.

## 5. Steps after agreement

1. Reader changes (2.3), `data-raw/ons_imd.R`, dataset and documentation.
2. Vignettes 1 to 3 (ONS), then 4 and 5 with precomputed results.
3. Validation checks on the ONS data (backtest scores with thresholds,
   coherence, fidelity tolerances for fdm and bamp) in `validation/`.
4. Update README, pkgdown articles, NEWS and `dev/progress.md`.

## 6. Questions for the owner

1. Series construction: rule A (recommended), B or C?
2. Storage: processed dataset in the package with the workbooks in
   `data-raw/`, or workbooks outside the repository? Is there a public ONS
   URL for the tabulation to cite?
3. Cause groups: the five groups of 2.4, or the groupings the colleagues
   used (for example IHD and stroke separately)? Is the COVID-19 column
   underlying cause or "involving"?
4. Pandemic years: jump-off 2024 with shock years, plus training to 2019?
5. Comparison: fdm and bamp plus Lee-Carter, fidelity and added value, and
   whether to include mortensor 0.1 results.
6. Should the ONS vignettes replace the simulated examples in "Getting
   started", with the regional vignette kept on simulated data?
