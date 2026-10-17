# Applied usability review 3: a regional public health user

Reviewer persona: Dr. A, epidemiologist at a regional public health institute
(rates, person-years, life tables, age standardisation, Poisson regression;
RStudio user who adapts vignette code; no data.table or S3 experience).
Goal: forecasts to 2035 of all-cause and cause-specific mortality by sex and
deprivation quintile, with report-ready life expectancy, age-standardised
rates (ASR), inequality gaps with intervals, and an explanation of which ages
and causes drive the gaps and their change (decomposition).

Version tested: branch `claude/vigilant-albattani-2ssw1x`, commit `c531ece`
(mortensor 0.1.9000, installed 2026-10-17), R 4.3.3, `library(mortensor)`.
Data: `ons_data` only. At most 300 forecast draws, two or three backtest
origins. Scripts and outputs (`common.R`, `s01_*.R` to `s14_*.R`, `*.out`,
`p*.png`) are in the `review3/` folder next to this report. I read
`dev/review-applied-1.md` (and skimmed `review-applied-2.md`) first; the
issues raised there that I re-met (observed indicators, contrasts, exposures
by name, stratum checks in plots) are fixed.

## 1. Summary verdict

The workflow for my goal now exists end to end on `ons_data`: one data
object, observed e0, ASR, gaps and SII with Poisson intervals, Li-Lee or
Lee-Carter forecasts to 2035, coherent cause forecasts, calibration, and a
decomposition that runs on observed data, on forecast draws and on both. All
quoted numbers that I re-ran in the deprivation and decomposition vignettes
reproduce (most exactly). `mt_decompose()` is fast, exact to 1e-14 in normal
cases, and its vignette explains the methods and the radix versus pooled
distinction better than most papers.

I would still publish wrong numbers in four places without any error or
warning: (1) a cause decomposition without `hierarchy` counts "all causes" as
a cause; (2) a subgroup decomposition on data with an aggregate level counts
the aggregate as a subgroup; (3) a backtest on data smoothed with
`mt_smooth()` scores against the smoothed values and makes smoothing look 33%
better; (4) decomposing "observed 2019 versus forecast 2035" attributes the
model's jump-off to ages, and the vignette's explanation of the narrowing gap
is mostly a jump-off effect. A fifth problem is substantive: the cause
approach the causes vignette recommends (bottom-up Lee-Carter per cause)
forecasts falling life expectancy at 50 for women in Q1 by 2035, driven by
the dementia trend, with some all-cause paths above a death rate of 1 at 90+
that no check reports. Smoothing, `anchored()` and `calibrated()` are well
explained where they are documented, but the narrowing of intervals by
smoothing is not mentioned on the help page of `step_smooth()`, in
`summary()` or in the model reference, and the vignette interprets the new
`dispersion` diagnostic opposite to its help page.

## 2. Findings

### P1: wrong results or blocking

**P1.1 `mt_decompose(cause = )` silently counts the all-cause level as a
cause when `hierarchy` is omitted.**
Code (`s03_decomp_causes.R`), with `xc` built exactly as in
`vignettes/decomposition.Rmd` (`mt_group_levels(..., total = "all", additive
= "rates")`):

```r
mt_decompose(xc, var = "quintile", ref = "Q1", cause = "cause")      # no hierarchy
#> female Q1 e50 25.36, Q5 31.71, gap 6.35   (no warning)
mt_decompose(xc, var = "quintile", ref = "Q1", cause = "cause",
             hierarchy = mt_hierarchy(cause = "all", additive = "rates"))
#> female Q1 e50 31.69, Q5 37.26, gap 5.57   (correct; life_expectancy() gives 31.69)
```

The rates are doubled, e50 is 6 years too low and the gap 0.8 years too
large. `?mt_decompose` calls `hierarchy` "Optional". `mt_validate(xc)`
already detects the problem ("total among levels ... cause = all ... equals
the sum of the other levels in 100% of cells"), and the object was created
with `total = "all"`. Without `cause =` at all, the function decomposes e50
of each cause stratum separately (values such as 323 years) without a
warning.
Fix: when `cause` (or `groups`) is given, run the total detection of
`mt_validate()` and abort with class `mt_error_hierarchy` if a level equals
the sum of the others and no hierarchy declares it; better, store the total
declared in `mt_group_levels(total = )` in the object and use it by default.
Warn when life expectancy is computed for strata of an `additive = "rates"`
variable (also in `life_expectancy()`).

**P1.2 `mt_decompose(groups = )` counts an aggregate stratum as a subgroup.**
Code: `xa <- mt_aggregate_strata(x, mt_hierarchy(sex = "both", quintile =
list(total = "England")))`, the object of the coherence vignette.

```r
mt_decompose(xa, periods = c(2001, 2019), groups = "quintile")
#> female: 80.533 -> 83.385, difference 2.852   (England counted as a 6th quintile)
mt_decompose(xa, periods = c(2001, 2019), groups = "quintile", hierarchy = h)
#> female: 80.418 -> 83.139, difference 2.721   (correct, as on x)
```

No warning. Fix: as in P1.1 (abort when a level of `groups` is the sum of the
others, or when `mt_aggregate_strata()` produced it; let
`mt_aggregate_strata()` store its hierarchy in the object).

**P1.3 Decomposing an observed year against a forecast year attributes the
jump-off to ages; the vignette's interpretation of the narrowing gap is
mostly a jump-off effect.**
Doc: `vignettes/decomposition.Rmd`, section "Forecasts": "Li-Lee narrows the
female gap from 6.5 to 4.6 years mainly at ages 60 to 84, where its
deviations of Q1 from the common trend return towards their average."
Code (`s04_forecast.R`, Li-Lee, 300 draws, seed 1):

```r
g <- mt_decompose(fc, var = "quintile", ref = "Q1", periods = c(2019, 2020, 2035))
```

| female gap Q5 - Q1 | 2019 observed | 2019 fitted (`mt_jump_off()`) | 2020 forecast | 2035 forecast |
|---|---|---|---|---|
| years | 6.46 | 5.69 | 5.12 | 4.57 |

Of the 1.89-year narrowing, 1.34 years (71%) occur between the last observed
year and the first forecast year, 0.77 of it because the Li-Lee fit of 2019
does not reproduce the observed gap. At ages 60-84, 0.84 of the 1.23 years
attributed to these ages is this first-year step. `anchored()` does not
change the 2035 decomposition (the offsets fade after 10 years). The
decomposition docs never mention the jump-off. A report would say "the gap
narrows because of ages 60-84" where most of the change is a model misfit in
the jump-off year.
Fix: (a) report a separate "jump-off" component when one period is observed
and the other forecast (observed to fitted in the last year), or offer
`base = c("observed", "fitted")`; (b) show `mt_jump_off()` next to the
decomposition in the vignette and rewrite the sentence; (c) mention in
`?mt_decompose` Details that observed-to-forecast comparisons include the
step from observed to fitted rates.

**P1.4 Backtests and fits on a pre-smoothed object are silently wrong.**
Doc: `vignettes/deprivation.Rmd` shows `xs <- mt_smooth(x, along = "period",
method = "whittaker")` and then backtests with `prep = smooth_t`; nothing
stops a user from passing `xs` itself.
Code (`s06_smooth.R`, Lee-Carter, origins 2010 and 2014, h = 5, 100 draws):

| data | CRPS | 95% coverage |
|---|---|---|
| `mt_backtest(lee_carter(), xs, ...)` (pre-smoothed) | 0.0461 | 0.81 |
| `mt_backtest(lee_carter(), x, ..., prep = step_smooth)` | 0.0699 | 0.73 |
| raw | 0.0686 | 0.81 |

The pre-smoothed backtest scores against smoothed "observations" (and the
smoother saw the test years), so smoothing looks 33% better than no
smoothing, while the honest comparison shows it slightly worse. No warning,
although `xs` carries the attribute `"smooth"`. `mt_fit(lee_carter(), xs)`
also drops the raw deaths (`fit$raw_D` is `NULL`), so the "smooth input, raw
residuals" behaviour announced in NEWS does not apply.
Fix: `mt_backtest()` aborts (and `mt_fit()` warns) when the data carry a
`"smooth"` attribute, naming `prep = mt_prep(step_smooth(...))`; or keep the
raw deaths in the `mt_data` returned by `mt_smooth()` and score against them.
Say in the vignette that `xs` is for checking only.

**P1.5 The recommended cause approach gives falling life expectancy to 2035,
and implausible paths are not reported.**
Doc: `vignettes/causes.Rmd` recommends the bottom-up sum of Lee-Carter models
per cause (better 5-year backtest for all causes; confirmed: CRPS 0.0377
against 0.0406) and forecasts with it. Code (`s08_causes_fc.R`, `s08e.R`,
ages 50+, to 2035, 300 draws):

| e50 in 2035, median (95% interval) | women Q1 | women Q5 |
|---|---|---|
| observed 2019 | 31.69 | 37.26 |
| Lee-Carter per cause, bottom-up | 30.84 (28.85-32.78) | 37.29 (34.79-40.73) |
| Lee-Carter on all causes | 33.10 (30.01-35.50) | 39.79 (37.94-41.92) |
| Li-Lee on all causes | 34.70 (30.35-37.76) | 38.84 (35.63-41.28) |

The cause model forecasts a fall of e50 for women in Q1 and a widening gap
(5.57 to 6.56 years), with dementia's contribution to the gap rising from 0.52
to 1.57 years. In 9 of 300 draws all-cause rates at 90+ exceed 1 (up to 1.41;
observed 0.18-0.24); the 97.5% quantile for women in Q1 is 0.82.
`mt_forecast()` gives no warning because `check_draws_plausible()` only
looks at closed age intervals; at 85-89 the same draws reach 0.43, below 1
but high enough that the life table caps q at 1 (l(90) = 0, e90 `NaN`; see
P2.3). The vignette discusses dementia
shares but never shows what the cause model implies for all-cause life
expectancy. This is the known effect that sums of cause-specific log-linear
extrapolations are dominated by the slowest-falling cause (Wilmoth 1995,
"Are mortality projections always more pessimistic when disaggregated by
cause of death?", Mathematical Population Studies 5: 293-319), which is not
cited.
Fix: in the causes vignette, compare all-cause e50 to 2035 for the three
approaches, cite Wilmoth (1995), and say that 5-year backtests do not support
16-year extrapolation of cause trends. Extend the plausibility warning to the
open interval (for example rates above 1, or life-table q capped at 1 in a
closed interval) and make sure it covers aggregates of reconciled forecasts.

### P2: important

**P2.1 Intervals of model-versus-model decompositions depend on the seeds.**
`mt_decompose(fc, fl, periods = 2035)` (Li-Lee against Lee-Carter, both
`seed = 1`) gives for women in Q1 a difference of -1.92 (-2.87 to -0.82);
with `seed = 2` for Lee-Carter -1.99 (-6.65 to 2.64). Draws are paired by
index, and identical seeds share innovations (correlation of the two e0
paths 0.995). The interval is an artefact, and the plot (`p5_mm.png`) shows
narrow bars for Q1 and wide ones for Q5. Fix: document that `y` is paired by
draw; warn when `x` and `y` come from different models; show only means (or
independent pairing) for such comparisons.

**P2.2 The default radix composition prints a value that is not the life
expectancy of the population.** `mt_decompose(x, periods = c(2001, 2019),
groups = "quintile")` prints female e0 80.418 to 83.139 (+2.72); the e0 of
the pooled rates is 80.649 to 83.631 (+2.98). The vignette explains this
well, but the printed object only says "radix composition". A report writer
will quote 80.4 as "life expectancy of England". Fix: label the totals
("radix-weighted mean of subgroup e0") in `print()` and `summary()`, or make
`composition` a required choice; consider `"pooled"` as default for
public-health use.

**P2.3 Decomposition residuals of up to half a year are not flagged.**
Cause decomposition of the bottom-up forecast (P1.5):
`max_residual` 0.29 and 0.40 years for the 2035 gaps, 0.49 for the change
2019-2035 of women in Q1. All draws are exact except a few with extreme
rates (draw 150: all-cause rate of women in Q1 at 85-89 of 0.43 and at 90+
of 1.13, so `life_table()` caps q at 1 at 85-89 and the Arriaga identity no
longer holds). `?mt_decompose`
calls Arriaga "Exact"; nothing explains a large residual. Fix: warn when
`max_residual` exceeds a tolerance and name the likely cause (implausible
rates), and say so in the help.

**P2.4 A fixed all-cause total zeroes cause-specific rates in most draws.**
`reconciled(lee_carter(), h, method = "wls", fixed = "total")` on the
cause data: `$reconciliation$share_set_to_zero` 0.5%, but 297 of 300 draws
have at least one cause rate set to 1e-12 (liver and alcohol 4,908 cells,
diabetes 2,927, dementia 2,376, other circulatory 1,935, mostly ages 60-79,
increasing towards 2035). The ASR of liver and alcohol for men in Q1 in 2035
falls from 127.9 (118.7-187.6) to 116.0 (85.7-177.2). The vignette says "the
other causes barely change", which holds for the female medians shown, not
for tails or for liver disease. No warning. Fix: warn when values are set to
zero (share of draws and strata), document the effect on cause-specific
intervals, and offer a proportional (top-down) allocation of the gap, which
is the usual choice when all-cause totals are fixed.

**P2.5 Smoothing: the narrowing of intervals is explained in some places
only, and the new diagnostic is explained in two contradictory ways.**
Measured (`s06_smooth.R`, Lee-Carter, e0 in 2035): 95% interval width for
women in Q1 5.6 years raw against 3.8 with `step_smooth(along = "period")`,
for men in Q5 2.9 against 1.6. `?mt_smooth` (Details), the deprivation
vignette and the FAQ explain this clearly; `?mt_prep` (`step_smooth`), the
"Data steps" table of `vignettes/models.Rmd`, the function map, `summary(fit)`
("preprocessing: step_smooth") and the forecast print do not.
`vignettes/deprivation.Rmd` reads dispersion 0.38-0.54 as "the smoother
removes noise and no visible signal", while `?mt_smooth_check` says "well
below 1 when the smoothing follows the noise". The 2D GAM gives dispersion
1.5-2.3 on the same data, which neither text interprets. Fix: one
interpretation in both places; a sentence on narrowing and calibration in
`?mt_prep`, in the models table and in `summary()` of fits with
`step_smooth(along = "period")` or `dims = 2`.

**P2.6 Model labels ignore preprocessing, so scores are silently pooled.**
`mt_backtest(lee_carter(), x, ...)` and the same with `prep = smooth` both
have `model = "Lee-Carter"`; `mt_score(rbind(b1, b2), by = "model")` returns
one row with n = 2000 (CRPS 0.0750 instead of 0.0743 and 0.0757). The
deprivation vignette works around this with `[, model := n]`. Plot titles of
smoothed forecasts also say only "Lee-Carter". Fix: add the preprocessing to
the label (for example "Lee-Carter + smooth(period)"), or warn in
`mt_score()` when one label covers different specifications.

**P2.7 No direct way to decompose the change of a gap.** My question "which
ages and causes drive the change of the gap from 2019 to 2035" needs the
difference of two decompositions. `summary(g, by = "age")` gives each year
separately; the change and its interval require matrix work on `g$paths`.
Fix: `summary(d, change = TRUE)` (path-by-path differences between periods
of the same comparison) or an example in the vignette; combine with the
jump-off component of P1.3.

**P2.8 Printing a filtered backtest fails.** `bt <- mt_backtest(...,
quantities = c("life_expectancy", "asr"))`; `bt[target ==
"life_expectancy"]` at the console gives `Error: The backtest holds no scores
for target "rate".` (the print method always scores rates). Fix: print the
targets present.

### P3: minor

- `vignettes/articles/function-map.Rmd`, "What is not covered": "life
  expectancy or age-standardised rates as targets are not yet available" is
  stale; the table on the same page lists them under `mt_backtest()`.
- `README.md` lists the vignettes but not "Decomposing differences in life
  expectancy".
- `mt_plot_decomposition()`: the title "Method: arriaga; bars: means; lines:
  95% intervals" overlaps the y-axis title, which is cut ("(years"); facet
  labels of model comparisons are truncated at width 7 ("ale, Q1: Lee-Carter
  vs Li-Lee (common by sex), 2"); with 12 causes, colorectal cancer and liver
  disease (brown) and IHD and other cancers (purple) are hard to tell apart.
  Put the method in the caption, wrap facet labels, use a palette with 12
  distinct hues or group small causes.
- `comparisons = c(4, 8)` (numbers of rows) is brittle; accept labels or a
  filter such as `to = "Q5"`. `comparisons = 99` returns a plot that fails
  when printed ("Faceting variables must have at least one value").
- `summary(d, by = "cause")` on a decomposition without causes silently
  returns the totals, while the plot stops with "The decomposition has no
  `cause`". Make both stop.
- Decomposition totals and bars are means, `life_expectancy()` reports
  medians: women Q1 2035 83.38 (mean) against 83.53 (median); change 3.55
  against 3.70. State this next to the tables, or add `centre =`.
- `mt_export(decomposition)` fails ("`x` must be an object of class
  <mt_draws>"); `mt_export(summary(d))` works. Accept the object.
- `vignettes/decomposition.Rmd`: "1.5 years of it at ages 60 and over";
  the vignette code gives 1.66 (500 draws, seed 1).
- `?calibrated`: "calibration then narrowed it (95% coverage 0.93 to 0.85 on
  ONS ...) (see the vignette ...)"; the vignette shows a different ensemble
  (Li-Lee with random walk, 0.864 to 0.850).
- FAQ "Why does the forecast start with a step?": "years = 1 ... was worse on
  the ONS data" refers to `li_lee(jump_off = "observed")`, not
  `anchored(years = 1)`; on England & Wales the observed jump-off was best
  (CRPS 0.0583 against 0.0648 anchored, `?anchored`), which the FAQ omits.
- `anchored()` with the default `decay = 10` has no effect on 2035, the year
  I report; say so in `?anchored` for long horizons. `summary()` of an
  anchored fit replaces the model description with "see the help page of its
  constructor" and does not show the size of the offsets.
- `mt_plot_coherence()` shows medians over ages and draws (within 3%),
  while single paths of the same forecast differ by -18% to +18%; add the
  path-level range so the figure does not suggest coherence.
- `mt_score(target = )` output has no unit: CRPS of rates and ASR are on the
  log scale, of life expectancy in years; backtest rows for `"asr"` show
  observed values such as -4.48 (log per person) while `asr()` reports 1,130
  per 100,000. Add a `unit` column.
- Observed indicators report the mean and quantiles of 500 Poisson draws,
  not the estimate from the observed rates (77.6684 against 77.6659); add an
  `estimate` column for reports.
- `vignettes/mortensor.Rmd` still computes observed e0 with `mt_array()` and
  `life_expectancy(m, mt_grid(x)$age)` although `life_expectancy(x)` exists.
- Terminology: `mt_contrast()` returns `level`/`ref`, `mt_decompose()`
  `from`/`to`; the FAQ writes `mt_hierarchy(cause = list(total = "all",
  additive = "rates"))`, the vignettes `mt_hierarchy(cause = "all", additive
  = "rates")` (both work; use one).

## 3. What works well

- **Observed indicators and gaps.** `life_expectancy(x)`, `asr(x)`,
  `mt_contrast(x, ...)` and `mt_inequality(x, ...)` give report tables with
  Poisson intervals in one call; every number of the deprivation vignette's
  "The past" section reproduced exactly.
- **Decomposition core.** `mt_decompose()` on observed data takes under a
  second, contributions add up to 1e-14, Arriaga, Horiuchi and stepwise agree
  to 0.005 years, and `quantity` accepts my own function (e65 by Horiuchi
  equals the built-in result to 1e-4). Kitagawa, radix and pooled
  composition are explained with a shares table that made the difference
  clear to me. Draw-by-draw intervals of forecast contributions, and the
  explanation of why single ages include zero while the total does not, are
  exactly what a report needs.
- **Error messages.** `ref` that is not a level, an unknown variable (with
  the list of variables), a year not in the data, `fixed` with
  `"bottom_up"`, an unknown label in `fixed`, `mt_smooth(along = "period")`
  with the default method, and an unknown stratum in `mt_plot_smooth()` all
  stop with a message that names the remedy.
- **`calibrated()`.** Deterministic, the deprivation calibration table
  reproduced to four decimals, the help page states the method, its
  evidence, its limits and the caveat on `by` for contrasts.
- **`anchored()` and `mt_jump_off()`.** The jump-off table is the clearest
  diagnostic I have seen for this problem; the forecast plot with observed
  points, fitted lines and dashed forecasts makes the step visible
  (`p8_le.png`).
- **Fixed totals.** `fixed = "total"` keeps England's e0 exactly (84.56) and
  the all-cause totals to 1e-8.
- **Backtests of derived quantities.** `mt_backtest(quantities = )` and
  `mt_score(target = )` worked first time; CRPS of e0 in years is easy to
  explain.
- **Honest vignettes.** The deprivation vignette shows where models fail
  (Q1), that smoothing narrowed coverage, and that the fitting period
  matters more than the intervals; the coherence vignette shows the cost of
  coherence. This is the right tone for a report.
- **Speed.** Li-Lee to 2035 with 300 draws and three decompositions: 5 s;
  130 cause strata bottom-up: 10 s; 2D GAM smoothing: 9 s for 10 strata.

## 4. Spot-checked numbers

"Exact" means identical at the printed precision; draws: vignette 500 or
1000, mine 300, seed 1.

| Location | Vignette | Mine | Status |
|---|---|---|---|
| mortensor.Rmd, observed e0 2001/2010/2019 (30 values) | e.g. female Q1 77.7/79.5/79.8 | same | exact |
| mortensor.Rmd, fitted 2019 above observed, women Q1 | about 0.4 | 0.43 | ok |
| deprivation.Rmd, e0 table 2001/2011/2019 (30 values) | | same | exact |
| deprivation.Rmd, ASR decline table (20 values) | | same | exact |
| deprivation.Rmd, gap and SII (12 rows) | female gap 5.17, 5.46, 6.46 | same | exact |
| deprivation.Rmd, calibration by quintile (6 rows) | Q1 1.0313, 0.3176 | same | exact |
| deprivation.Rmd, `mt_smooth_check()` (10 rows) | | same | exact |
| deprivation.Rmd, jump-off women Q1, Li-Lee / anchored | 80.891 / 80.374 | 80.899 / 80.382 | ok |
| deprivation.Rmd, "1.1 years above ... 0.4 from fit; anchored 0.5" | | 1.07; 0.43; 0.55 | ok |
| deprivation.Rmd, e0 2035, Lee-Carter / Li-Lee, women Q1 | 81.5 / 83.4 | 81.6 / 83.5 | ok |
| deprivation.Rmd, gap 2035 Li-Lee women | 4.6 (3.7-5.6) | 4.5 (3.7-5.6) | ok |
| deprivation.Rmd, Lee-Carter 2011-2019, men Q1 2035 | 75.8 (74.5-76.6) | same | exact |
| deprivation.Rmd, anchored and calibrated equals Li-Lee in 2035 | | 83.53 both | ok |
| decomposition.Rmd, gap 2019 women / men | 6.5 / 7.7 | 6.46 / 7.66 | ok |
| decomposition.Rmd, ages 40-84 | 5.2 / 6.3 | 5.21 / 6.31 | ok |
| decomposition.Rmd, below 40; infant (women) | 0.6 / 0.9; 0.2 | 0.55 / 0.93; 0.22 | ok |
| decomposition.Rmd, methods table (20 ages) | | same | exact |
| decomposition.Rmd, top named causes at 50 | respiratory, lung, IHD | women 0.89, 0.71, 0.65; men IHD 1.09, resp. 0.77, lung 0.71 | ok |
| decomposition.Rmd, England change radix | +2.7 / +3.7, composition negative | +2.72 / +3.74; -0.10 / -0.13 | ok |
| decomposition.Rmd, pooled higher by | 0.2-0.5 | 0.23-0.52 | ok |
| decomposition.Rmd, shares girls Q5 at birth; Q5 and Q1 at 75-79 | 18% to 16%; 19% to 24%; 19% to 14% | 18.4/15.6; 19.3/23.9; 18.8/14.3 | ok |
| decomposition.Rmd, Kitagawa rate effects | 140 to 360 | 138 to 357 | ok |
| decomposition.Rmd, men Q5 crude rate rose | yes | +58.5 | ok |
| decomposition.Rmd, gain women Q1 to 2035; ages 50-84 | 3.5; 2.5 | 3.55 (median 3.70); 2.54 | ok |
| decomposition.Rmd, intervals include zero at most ages | | 19 of 20 ages | ok |
| decomposition.Rmd, female gap 2019 to 2035 | 6.5 to 4.6 | 6.46 to 4.57 | ok (see P1.3) |
| decomposition.Rmd, Lee-Carter minus Li-Lee, women Q1; at 60+ | -1.9; 1.5 | -1.92; 1.66 | 60+ differs |
| coherence.Rmd, England against sum of parts, path quantiles | -15.8, 0.3, 20.1 | -17.5, 0.6, 18.1 | ok |
| coherence.Rmd, fixed total keeps England e0 | 84.55 | 84.56 base and fixed | ok |
| causes.Rmd, steepest and flattest gradients | | respiratory 3.56, lung 2.95, liver 2.76, diabetes 2.43; colorectal 1.18, breast/prostate 1.13 (women) | ok |
| causes.Rmd, bottom-up better for all causes than all-cause LC | | CRPS 0.0377 against 0.0406 (2 origins, 100 draws) | ok |
| causes.Rmd, fixed WLS: gap falls on dementia, others barely change | | dementia yes; liver and alcohol, men, median -9% to -15% | partly (P2.4) |
| `?mt_smooth`, 2D GAM about one second per stratum | | 8.6 s for 10 strata | ok |

Not re-run because of the time limit: the backtest tables of
`mortensor.Rmd`, `deprivation.Rmd` (seven origins) and `coherence.Rmd`, and
the cause backtest with five specifications.
