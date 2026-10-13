# Applied usability review 1: a regional public health user

Reviewer persona: Dr. A, epidemiologist at a regional public health institute
(rates, person-years, life tables, age standardisation, Poisson regression;
RStudio user who adapts vignette code; no data.table or S3 experience).
Goal: forecasts to 2035 of all-cause and cause-specific mortality by sex and
deprivation quintile (and district groups) from office CSV tables, with
report-ready life expectancy, age-standardised rates (ASR), expected deaths
and inequality gaps with intervals.

Version tested: commit `5421c28` with the uncommitted working tree as it
changed between 09:55 and 10:45 UTC on 2026-10-13. The developer edited the
import code during the review; items that changed while I worked are marked
"(changed during review)". Scripts and simulated office files are in the
session scratchpad (`review1/`), outside the repository.

## 1. Summary verdict

mortensor has the right architecture for my task: one data container that
accepts my age groups, shared populations for causes, 31 December
populations, COVID years as shocks, coherent cause forecasts, and life
expectancy and ASR computed path by path with intervals. The getting-started
and regional vignettes explain modelling assumptions better than most
packages. I could not reach a trustworthy report from my own CSV files
without reading source code. Three problems produced wrong numbers without an
error (suppressed cells filled with population counts, exposures matched to
strata by position, preprocessed values plotted as observed). Sparse young
ages made the default models produce death rates far above 1. And the
quantities a report needs most (observed life expectancy and ASR,
inequality gaps, expected deaths under a population projection, exported
tables) require array code copied from vignettes.

## 2. What worked well

- **Reading long tables.** `mt_data()` read my German file with causes in a
  column and a population table without causes; the populations were shared
  by all causes as the Details section promises. Labels such as
  `"unter 1 Jahr"` and `"1 bis unter 5 Jahre"` parsed.
- **Helpful import warnings.** An undeclared symbol produced:
  `Values of deaths that are neither numbers nor declared markers are treated
  as missing: "<3". Declare them with 'missing' or 'zero'.` This is exactly
  what I needed to read.
- **Exposure choices are explicit.** `mt_as_exposure(type = "dec31")` and the
  FAQ entry "I have populations on 31 December" answered my question at once.
- **COVID years.** The FAQ snippet
  `mt_prep(step_shocks(2020:2022, action = "interpolate"))` worked first time.
- **Suppressed cells with published totals.** For causes with a published
  all-cause total, `mt_fill_suppressed(x, mt_hierarchy(Todesursache =
  "Insgesamt", additive = "rates"))` filled 990 cells with values between 1
  and 2, and `mt_validate()` confirmed the sums.
- **Derived quantities.** `life_expectancy(fc, at = 0)`, `asr(fc)` and
  `prob_dying(fc, 30, 70)` return tidy tables with my stratum columns and
  quantiles. The idea of computing them path by path is well explained.
- **Coherent causes.** `mt_reconciled(li_lee(common_by = c("Geschlecht",
  "Todesursache")), hc, method = "bottom_up")` on 40 strata fitted and
  forecast in 2.5 seconds; `mt_plot_shares()` gave a clear figure.
- **Speed.** Li-Lee for 10 strata with 200 draws ran in under 1 second.
- **Decision table.** The "Summary of the decisions" table at the end of the
  regional vignette is the most useful page for an applied user; every step
  is named with its function.

## 3. Walkthrough log of my own task

I simulated two office files (about 1.7 million inhabitants, 2008-2023,
5-year groups with `<1`, `1-4` and an open `90+`): a German file
(`Jahr;Geschlecht;Quintil;Todesursache;Altersgruppe;Gestorbene`, sex as
`"männlich"`/`"weiblich"`, `"."` for 1-2 deaths, `"-"` for zero, rows
`"Insgesamt"` for all ages and all causes, populations on 31 December from
2007) and an English file (sex `1`/`2`, `"<3"`, `"All ages"` rows, mid-year
populations).

**Step 1. Reading the English file.**

```r
d <- read.csv("deaths_en.csv"); p <- read.csv("pop_midyear_en.csv")
x <- mt_data(deaths = d, population = p, by = c("Sex", "IMD_quintile"), age = "Age",
             period = "Year", value = c(deaths = "Deaths", population = "MidYearPop"))
#> Error: Cannot parse age label(s): "All ages".
```

I removed the rows by hand. (Changed during review: `mt_data()` now drops
total rows with a warning.) The strata were labelled `"1.1"`, `"2.5"`. I
found the way to recode sex only by reading in `?mt_group_levels` that "a
group with one member renames the level":
`mt_group_levels(x, "Sex", list(male = "1", female = "2"))`.

**Step 2. Reading the German file.** `"90 und älter"` failed to parse at
first; I recoded it to `"90+"`. (Changed during review: it now parses.)

**Step 3. Suppressed `"<3"` cells without stratum totals.** `mt_fit()`
stopped with `Deaths and exposures must be complete on the grid (no missing
cells); suppressed cells can be estimated from published totals with
mt_fill_suppressed().` My file has only age totals, so I guessed:

```r
r <- mt_fill_suppressed(x, mt_hierarchy(IMD_quintile = "all"))
attr(r, "fill")   # filled 432, inconsistent 0, iterations 0, max_gap 0
```

No warning. The 432 filled death counts ranged from 1,760 to 11,306: they
are the population counts of those cells (rates of 1.0 at ages 5-9).
`mt_validate(r)` reported "no findings", and `mt_fit(lee_carter(), r)` then
failed with `The model failed in every stratum: Lee-Carter fit diverged.`
I replaced `"<3"` by 1.5 by hand to go on. (Verified in source: cells that no
total covers keep their starting value, the exposure; the warning "could not
be estimated" never fires because those cells are no longer `NA`.)

**Step 4. 31 December populations.** With populations for 2007-2023,
`mt_as_exposure(x, type = "dec31")` worked, but `mt_fit()` gave the same
"suppressed cells ... mt_fill_suppressed()" error: the grid still contained
2007, which has populations and no deaths. The regional vignette calls
`mt_subset(periods = ...)` at this point without saying why. I found the
cause only by printing the object and comparing period ranges.

**Step 5. Fitting with the COVID years interpolated.**

```r
prep <- mt_prep(step_shocks(2020:2022, action = "interpolate"))
fc <- mt_forecast(mt_fit(li_lee(common_by = "Geschlecht"), xa, prep = prep),
                  to = 2035, n_draws = 200, seed = 1)
life_expectancy(fc, at = 0)
#> weiblich, Q1, 2035: q0.025 47.6, q0.5 76.2, q0.975 76.8
```

With `rwd()` the 2.5% quantile of e0 was 0.5 years; with `step_smooth()`
(P-splines) added, single draws reached a death rate of 292 per
person-year. The cause is sparse counts at ages 1-14 (0-5 deaths per cell).
(Changed during review: `mt_forecast()` now warns "forecast death rates
above 1 in 7.0% of the draws ... consider wider age groups, pooling,
smoothing or an ensemble". This is a good message.)
`step_smooth(method = "whittaker")` gave plausible e0 intervals, but it
smooths across age 0 and across groups of width 1, 4 and 5: the 2019 infant
rate of men in Q1 fell from 0.00396 to 0.00113 and the rate at 1-4 rose
from 0.00014 to 0.00102. The forecast infant rate for 2024 (median 0.0014)
is then about 40% of the observed 2008-2023 mean (0.0036).

**Step 6. Observed indicators for the report.** `life_expectancy(x)` on my
data object gave `` `age` must be an object of class <mt_axis> `` and
`asr(x)` gave `` `draws` must be an object of class <mt_draws> ``. I copied
`mt_array(x, "deaths") / mt_array(x, "exposure")` from the getting-started
vignette for e0; for observed ASR I would have to copy the weight-matching
code of the causes vignette. Observed e0 has no confidence interval.

**Step 7. Inequality gap.** I adapted the 10-line `gap_paths()` function of
the getting-started vignette (`life_expectancy(fc$rates, fc$age)`, matching
`paste0(sx, ".5")` against `fc$strata$label`, `apply(..., quantile)`). It
worked, but I would not have written it myself.

**Step 8. Expected deaths with the office's population projection.**

```r
mt_deaths(fc, exposure = proj)          # data frame: "`exposure` must be [age, stratum] or [age, period, stratum]."
mt_deaths(fc, exposure = mt_data(population = proj, ...))   # "`exposure` does not match the draws."
mt_deaths(fc, exposure = mt_array(px, "population"))         # works, returns [age, period, stratum, draw]
```

The result is a 4-dimensional array without a summary. An exposure matrix
with the strata in reverse order changed the deaths by 24% on average
without any message: exposures are matched by their position in the array,
and the stratum names are ignored.
The help title says "Expected deaths", while the default `noise = TRUE`
returns Poisson draws of observable deaths.

**Step 9. Figures.** `mt_plot_life_expectancy(fc, at = 0)` draws all 10
strata as overlapping fans in one panel, without the observed years, and has
no `strata` or facet argument. `mt_plot_forecast()` plots the training data
after preprocessing as "observed": the 2020 rate at 60-64 for men in Q1 was
0.02148 in the data and is drawn as 0.01978 (interpolated), so the COVID
excess disappears from the figure. A stratum label that does not exist
(`strata = "Q1"`) returns a plot that fails when printed with the ggplot2
message `Faceting variables must have at least one value.`

**Step 10. Export and backtest.** `mt_export(fc, dir)` writes rate
quantiles only; `mt_export(life_expectancy(fc, at = 0), dir)` fails. The
backtest warned clearly that origin 2014 had too few years for Li-Lee.

## 4. Prioritised recommendations

### P1: blocking for an applied user

**P1.1 `mt_fill_suppressed()` must never return unconstrained cells as
estimates.** Evidence: Step 3 (432 death counts equal to populations,
`max_gap 0`, `mt_validate()` silent). Proposed change: keep cells that no
known total covers as `NA` and warn with their number; abort when none of
the hierarchy's totals exist in the data ("no stratum `IMD_quintile ==
"all"` in the data; add totals or use mt_aggregate_strata()"). Let
`mt_validate()` flag estimated cells whose rate exceeds the largest observed
rate at that age.

**P1.2 Support interval-censored suppression and age totals.** Evidence:
offices publish `"<3"`, `"1-2"` or `"<5"` together with all-age totals; the
package treats them as plain missing values and uses only stratum totals.
Proposed change: an argument such as `mt_data(censored = c("<3" = 2))` that
stores an upper bound, and use of the bound and of all-age totals
(`"Insgesamt"` rows) as constraints in `mt_fill_suppressed()`. A documented
simple fallback (`fill = "midpoint"`) would be enough for many users.

**P1.3 Name the real cause when cells are missing on the grid.** Evidence:
Step 4. Proposed change: `mt_as_exposure()` drops periods without deaths
from the grid (with a message), and the error in `mt_fit()` lists the
missing cells, for example "deaths missing for all ages in 2007 (population
only); restrict with mt_subset(periods = c(2008, 2023))".

**P1.4 Match exposures by stratum and age names.** Evidence: Step 8.
Proposed change: `mt_deaths()`, `mt_forecast(exposure =)`, `yll()`,
`mt_reconcile()` and `mt_aggregate_draws()` accept an `mt_data` object or a
long data frame of projected populations (with `type` as in
`mt_as_exposure()`), align by names, and abort on unnamed or mismatched
arrays.

**P1.5 Observed indicators and contrasts as functions.** Evidence: Steps 6
and 7; the vignettes need `fc$rates`, `fc$strata`, `$obs` and
`mortensor:::hierarchy_structure()` (regions vignette) for routine tasks.
Proposed change: methods of `life_expectancy()`, `asr()`, `prob_dying()` for
`mt_data` (with Chiang or Silcocks confidence intervals for e0 and Poisson or
gamma intervals for ASR), and one function for differences and ratios
between strata, path by path, for example
`mt_contrast(fc, quantity = "e0", var = "Quintil", ref = "5", type =
c("difference", "ratio"))`. Add the slope and relative index of inequality
(`sii()`, `rii()`), which public health reports on deprivation expect.

**P1.6 A vignette "Your own data, from CSV to report".** Evidence: every
vignette starts from a cleaned data set. Proposed content: reading a CSV with
`read.csv2()`, total rows, sex codes, suppression symbols, 31 December
populations, checking with `mt_validate()`, COVID years, sparse young ages
(wider groups such as 1-14, pooling, Whittaker smoothing from age 1),
forecast, observed and forecast e0/ASR in one table, gaps, expected deaths
with a projection, figures, export.

### P2: important

- **Sparse ages and smoothing.** The new plausibility warning is good;
  `life_expectancy()` should also report the share of affected draws.
  `mt_smooth()` and `step_smooth()` should leave age 0 unsmoothed and use the
  interval widths (Step 5), and the help should say so.
- **`mt_as_exposure()` default.** With `type = "midyear"` as default, a
  31 December table is silently treated as mid-year. Make `type` required.
- **Plots for reports.** Add `observed = TRUE` (raw data, with preprocessed
  values marked) to `mt_plot_life_expectancy()`, a `strata` argument and
  `facet` (for example by sex), and plots for ASR and contrasts. Check
  `strata` and `ages` against the available labels and list them in the
  error.
- **Units and metadata in outputs.** `asr()` returns no unit, standard or age
  range; with data from age 50 the values are per 100,000 of the standard
  population aged 50+, which a reader will compare with all-age ASRs. Add
  columns `standard`, `ages` and an argument `per = 1e5`; let
  `mt_summarise()` take `per` as well.
- **Export of report tables.** Let `mt_export()` accept any table returned by
  `life_expectancy()`, `asr()`, `mt_contrast()` and `mt_deaths()` summaries,
  with rounding and an optional wide layout (years as columns).
- **Standard populations.** Only `"esp2013"` exists. Add ESP 1976, WHO
  2000-2025 and a documented way to pass a national standard; an error for
  `asr(standard = "who")` should list the available names (currently "The
  standard population does not cover the grid ages").
- **Level labels.** Offer `mt_data(labels = list(Sex = c("1" = "male", "2" =
  "female")))` or `mt_relabel()`.
- **Defaults tied to example data.** `mt_plot_shares(var = "cause",
  total_level = "all")` should take both from the hierarchy stored in a
  reconciled forecast.
- **Stale help pages.** `man/mt_subset.Rd` lacked the `drop` argument that the
  regional vignette uses (regenerated during review). Run roxygen before
  each push that changes arguments.

### P3: nice to have

- Remove "mortensor 0.1" from help pages (`mt_cohere()`, `mt_disaggregate()`,
  `mt_smooth()`); new users do not know that version.
- `mt_forecast(h = 5, to = 2030)` silently uses `to`; warn when both are
  given.
- `life_expectancy()` on an `mt_data` object (until P1.5 exists) should say
  "compute rates with mt_array() first".
- Merge the help topics `mt_plot` and `mt_plot_more` into one plot family.
- Show the scale of backtest scores (log rate) in the `mt_score()` output and
  explain MASE above 1.

## 5. Proposed structure from a user's point of view

Organise the reference and the main help page by the steps of an analysis,
with the low-level layer in a separate section.

| Step | Functions (proposed changes in brackets) |
|---|---|
| 1. Read | `mt_data()` (with `labels`, `censored`), `mt_read_genesis()`, `mt_read_gbe()`, `mt_read_hmd()`, `as_mt_data()` |
| 2. Prepare | `mt_as_exposure()` (type required), `mt_subset()`, `mt_group_levels()`, `mt_relabel()` (new), `mt_fill_suppressed()`, `mt_validate()`, `mt_plot_rates()` |
| 3. Describe the past | `life_expectancy()`, `asr()`, `prob_dying()` for data (new methods), `mt_contrast()`, `sii()`, `rii()` (new) |
| 4. Model | constructors `rwd()`, `lee_carter()`, `li_lee()`, `tucker3()`, `apc_bayes()`, ...; `mt_prep()` with `step_*()`; `mt_fit()`, `mt_forecast()` |
| 5. Combine and make coherent | `mt_hierarchy()`, `mt_aggregate_strata()`, `mt_reconciled()`, `mt_ensemble()` |
| 6. Check | `mt_backtest()`, `mt_score()`, `mt_plot_scores()`, `mt_plot_coverage()` |
| 7. Report | `life_expectancy()`, `asr()`, `prob_dying()`, `yll()`, `mt_deaths()` (with summary), `mt_contrast()`, `mt_summarise()`, `mt_plot_*()`, `mt_export()` (any table) |
| Advanced | `mt_set_grid()`, `mt_regroup()`, `mt_ungroup()`, `mt_smooth()`, `mt_disaggregate()`, `mt_aggregate_draws()`, `mt_reconcile()`, `mt_cohere()`, `mt_tune()`, `mt_array()`, `mt_obs()`, `mt_grid()`, `mt_granularity()` |
| Internal (unexport or keyword internal) | `mt_operator()`, `mt_window()` (covered by `mt_subset(periods =)`), `mt_prep_apply()`, `new_mt_spec()` (move to an "Extending" article), axis helpers `axis_widths()`, `axis_midpoints()`, `is_open()` |

Naming:

- Six functions deal with grouping (`mt_regroup`, `mt_group_levels`,
  `mt_aggregate_strata`, `mt_aggregate_draws`, `mt_ungroup`,
  `mt_disaggregate`). Name them by what they act on: `mt_regroup_ages()` for
  age and period axes, `mt_combine_levels()` for levels of a stratum variable,
  `mt_add_totals()` for aggregate strata. Keep `mt_aggregate_draws()` and
  `mt_disaggregate()` in "Advanced".
- Pairs that differ by one letter (`mt_reconcile`/`mt_reconciled`,
  `mt_tune`/`mt_tuned`) are easy to confuse. Model wrappers are constructors
  and the package convention gives constructors noun names; name them
  `reconciled()` and `tuned()` (and `ensemble()`), after checking clashes with
  the packages listed in the export test.
- Derived quantities have no `mt_` prefix and mix abbreviations (`asr`,
  `yll`) with full names (`life_expectancy`). Either is acceptable if they
  are listed together in one help page "Indicators for reports" that states
  the unit of each.
- Use one argument name for selecting strata in all plots (`strata`), and
  replace `groups` in `mt_plot_shares()` with it.
- Readers named after sources (`mt_read_gbe()`, `mt_read_genesis()`) are fine
  as adapters; the general path (`mt_data()` on any CSV) should be the first
  thing a user sees in the README.

Functions tailored to particular data: `mt_read_gbe()` and
`mt_read_genesis()` (German sources, acceptable as adapters);
`mt_plot_shares()` defaults (`"cause"`, `"all"`) follow the example data;
`ons_cause_groups` and `de_cause_groups` are useful examples and should be
presented as templates for one's own cause groups. The modelling and
hierarchy functions are general.
