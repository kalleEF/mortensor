# Frequently asked questions

## Data

### How do I bring my own data in?

[`mt_data()`](../reference/mt_data.md) takes long data frames, one row
per cell, with a column for deaths, exposures or population and columns
for age, period and strata. Ages can be single years, labels such as
`"<1"`, `"01-04"`, `"85-89"`, `"90+"`, `"85 and over"`, or two columns
with lower and upper bounds:

``` r

parse_age_labels(c("<1", "01-04", "85-89", "90+", "85 and over"))
#>          label age_lo age_hi
#>         <char>  <num>  <num>
#> 1:          <1      0      1
#> 2:       01-04      1      5
#> 3:       85-89     85     90
#> 4:         90+     90    Inf
#> 5: 85 and over     85    Inf
d <- data.frame(region = "A", year = 2019, age = c("85-89", "90+"), deaths = c(120, 210),
                exposure = c(2500, 1400))
mt_data(deaths = d, exposure = d, by = "region")
#> <mt_data> 1 stratum by region, omega = 110
#>   deaths     5-year groups 85-89, 90+; periods 2019-2019 (width 1); 2 cells
#>   exposure   5-year groups 85-89, 90+; periods 2019-2019 (width 1); 2 cells
#>   grid       <mt_axis: age, 2 intervals> 85-89 90+
#>              <mt_axis: period, 1 interval> 2019
```

Deaths and exposures may come at different resolutions (for example
deaths in 5-year groups, exposures by single year of age). Tables in
wide format, with one value column per sex or cause, are read by listing
the columns
(`value = list(deaths = c(female = "Female", male = "Male"))`). Values
given as text with symbols for suppressed or missing cells (`"."`,
`"x"`, `"-"`) are recognised; other symbols are declared with `missing`
and `zero`. Readers exist for Human Mortality Database tables
([`mt_from_hmd()`](../reference/mt_from_hmd.md),
[`mt_read_hmd()`](../reference/mt_from_hmd.md)), GENESIS and GBE-Bund
exports from Destatis
([`mt_read_genesis()`](../reference/mt_read_genesis.md),
[`mt_read_gbe()`](../reference/mt_read_gbe.md)), `StMoMoData` and
`demogdata` objects ([`as_mt_data()`](../reference/as_mt_data.md)).

### I have populations on 31 December, not exposures. What now?

Read them as `population` and convert them with
[`mt_as_exposure()`](../reference/mt_as_exposure.md), whose `type` says
what the counts are: `"dec31"` for populations on 31 December (the
exposure of year t is the mean of the populations on 31 December of t -
1 and t), `"jan1"` for populations on 1 January, `"midyear"` for
mid-year estimates. Exposures are the denominator of the rates, so this
choice matters for the level of every forecast.

### What does `mt_validate()` check?

Cells with zero exposure, deaths without exposure, death rates above one
below age 100, missing values, strata without deaths or exposures,
levels that look like totals of other levels (for example `"Insgesamt"`
among the regions), grid cells without observations, and, with a
hierarchy, whether totals equal the sums of their parts (deaths and
exposures add up over regions or sexes; rates add up over causes,
`mt_hierarchy(cause = list(total = "all", additive = "rates"))`).

## Granularity and ungrouping

### Should I model single years of age or the age groups of my data?

Model at the resolution the question needs. Life expectancy,
probabilities of dying and age-standardised rates can be computed from
5-year groups (the life table uses the interval widths). Single years
are needed when results are reported by single year, when a model needs
them (cohort models such as APC need equal age and period steps), or
when data at different resolutions are combined. Backtests always score
against the deaths at their observed resolution, so ungrouping errors
count against the model, and
[`mt_backtest()`](../reference/mt_backtest.md) with and without
[`step_ungroup()`](../reference/mt_prep.md) answers the question for
your data.

### Which ungrouping method should I use?

The penalized composite link model
([`mt_ungroup()`](../reference/mt_ungroup.md), default) ungroups deaths
as rates with the exposures as offset and keeps group totals exactly. On
England & Wales data grouped into 5-year groups, it and the spline
methods (`"hyman_kannisto"`, `"natural_gompertz"`, `"gam"`) were almost
as close to the known single years as the Poisson noise allows; they
differed in the last, wide age group, where the PCLM was closest. Check
your own result with
[`mt_plot_ungroup()`](../reference/ungroup-diagnostics.md) and
[`mt_ungroup_check()`](../reference/ungroup-diagnostics.md); see
“Ungrouping age groups: methods and diagnostics”.

### Can I split the open interval 90+ into single years?

Only from assumptions: the data say nothing about the age pattern inside
the interval. [`mt_ungroup()`](../reference/mt_ungroup.md) warns when it
extrapolates. Prefer a grid whose open interval starts where the data’s
does, or data with a higher open age (for example exposures to 100+ with
deaths to 90+, where the exposures carry the age pattern).

## Models and uncertainty

### Which model should I use?

There is no single best model; a rolling-origin backtest on your own
data ([`mt_backtest()`](../reference/mt_backtest.md), horizons like the
ones you will report) decides. As a starting point:

| Situation | Try first | Why |
|----|----|----|
| One population, many deaths per cell, single years or 5-year groups | [`lee_carter()`](../reference/lee_carter.md), [`rwd()`](../reference/rwd.md) and `ensemble(lc = lee_carter(), rw = rwd())` | Lee-Carter is the standard; single models give too narrow intervals on real data, the ensemble is better calibrated |
| Several strata that should not diverge (sexes, regions, deprivation groups) | `li_lee(common_by = ...)` | a common trend with bounded deviations; pools information across strata |
| Small populations, many cells with few deaths | `li_lee(common_by = ...)`, wider age groups ([`mt_regroup()`](../reference/mt_regroup.md)), [`step_smooth()`](../reference/mt_prep.md); `rwd(observation_noise = TRUE)` in ensembles | pooling and smoothing reduce noise; the plain random walk accumulates it |
| Totals and parts must add up | `reconciled(model, hierarchy)` | bottom-up for joint models, WLS for separate models |
| Causes of death | separate or Li-Lee models per cause, `reconciled(..., method = "bottom_up")` with `cause = list(total = "all", additive = "rates")` | cause rates add up to the all-cause rate |
| Cohort effects (for example generations born around 1925-1945) | [`apc_bayes()`](../reference/apc_bayes.md), `stmomo("apc")` | age-period models attribute cohort patterns to period trends |
| Trend change in the recent past | `step_window(years)` with a backtest of several windows | a shorter fitting period follows the recent trend |

Check every fit with `summary(fit)`, `mt_parameters(fit)` and
`mt_plot_residuals(fit)`, and every forecast for warnings about
implausible rates.

### Why are my 95% intervals too narrow (or too wide)?

Intervals of most models describe the uncertainty of their own
assumptions, for example a constant drift. On real data, changes of
trend add uncertainty that these assumptions leave out, so intervals
tend to be too narrow over 5 to 10 years: in backtests on HMD and ONS
data the 95% coverage of Lee-Carter was 0.70 to 0.84. The plain random
walk is wider, because year-to-year changes of observed rates contain
observation noise twice; that extra width happened to give good 95%
coverage on real data. `rwd(observation_noise = TRUE)` removes the
double counting and is sharper (lower CRPS), with intervals that were
too narrow on real data. Mixtures of models
([`ensemble()`](../reference/ensemble.md)) widen the intervals where the
models disagree. [`recalibrated()`](../reference/recalibrated.md)
rescales the spread of any model’s draws to its forecast errors in the
last years of the training data (by default ten); on England & Wales men
it raised the 95% coverage of Lee-Carter from 0.75 to 0.88. Check
calibration with [`mt_backtest()`](../reference/mt_backtest.md) and
[`mt_plot_coverage()`](../reference/mt_plot.md).

### How do I choose smoothing and model parameters by forecast accuracy?

This is tuning: fit each candidate setting to the years up to `T - h`,
forecast the last `h` years, and keep the setting that forecast them
best, then refit it to all years. [`mt_tune()`](../reference/mt_tune.md)
does this for any grid of settings, including preprocessing, and
[`tuned()`](../reference/mt_tune.md) makes it part of a model, so that
an outer backtest repeats the choice at each origin without seeing the
test years:

``` r

grid <- expand.grid(k = c(10, 30), order = 1:6)
spec <- tuned(function(order, k) hyndman_ullah(order = order, smooth_args = list(k = k)), grid,
              inner_origins = 1, inner_h = 10, metric = "rmse")
```

`inner_origins = 1` is a single hold-out of the last ten years; several
origins (`inner_origins = 3, inner_h = 5`) average over more periods.
`metric = "rmse"` chooses by point accuracy, `"crps"` by the whole
forecast distribution, and `target = "life_expectancy"` by the quantity
of interest. [`recalibrated()`](../reference/recalibrated.md) answers a
different question: it keeps the model and rescales the width of its
intervals to the recent forecast errors. The two combine as
`recalibrated(tuned(...))`.

### Should I smooth the data before forecasting?

Smoothing over age can help where cells hold few deaths (rare causes,
small regions); check it in a backtest, as the applied vignettes do.
Smoothing over time removes the year-to-year noise from the fitted
trends and from the start of the forecast, but also the variation from
which models such as Lee-Carter estimate the uncertainty of their period
index: alone, it made intervals too narrow in backtests (95% coverage of
Lee-Carter 0.59 instead of 0.75 on England & Wales men). Most of the
gain there comes from [`recalibrated()`](../reference/recalibrated.md),
which re-estimates the spread against the raw observations: calibration
alone gave a CRPS of 0.0657 and a coverage of 0.88 (against 0.0687 and
0.75); Whittaker smoothing over time added nothing to it, the age-period
GAM a little (0.0646, coverage 0.91). On the ONS quintiles smoothing
over time made the forecasts worse (“Smoothing death rates”). Treat it
as optional and check it in a backtest. Use smoothing as a preprocessing
step (`mt_prep(step_smooth(...))`), never on the data before fitting
([`mt_backtest()`](../reference/mt_backtest.md) stops,
[`mt_fit()`](../reference/mt_fit.md) warns), so that backtests refit it
on each training window and models take their observation noise from the
raw rates. [`mt_smooth_check()`](../reference/mt_smooth_check.md) and
[`mt_plot_smooth()`](../reference/mt_smooth_check.md) show what the
smoothing removed.

### Why does the forecast start with a step?

Most models start from their fitted rates of the last observed year.
When the fit misses the level of the last years, the first forecast year
jumps away from the observations;
[`mt_jump_off()`](../reference/mt_jump_off.md) reports the step, and the
forecast plots show it as a gap between the solid (fitted) and the
dashed (forecast) line. `anchored(model, years = 3)` moves the start to
the observed level of the last three years, fading out over ten years
(so forecasts more than ten years ahead are those of the model); it
improved short-term accuracy on England & Wales men and on the ONS data.
Starting from the observed rates of the last year alone
(`lee_carter(jump_off = "observed")`, or
`anchored(years = 1, decay = Inf, shrink = FALSE)`) was best on England
& Wales men (CRPS 0.0583 against 0.0639 anchored to three years), whose
counts are large, and worse on the ONS quintiles, whose last year is
noisier. Decompositions of a change from an observed to a forecast year
include the step; [`mt_decompose()`](../reference/mt_decompose.md)
reports it (`jump_off`) and can start from the fitted rates
(`observed = "fitted"`).

### What is the difference between predictive draws and forecast rates?

For models of smooth rates (Lee-Carter, Li-Lee, Tucker3, FDM, StMoMo,
Bayesian APC, `rwd(observation_noise = TRUE)`),
[`mt_forecast()`](../reference/mt_forecast.md) returns draws of the
smooth death rates by default, and `predictive = TRUE` adds the
variation of observed rates around them (extra-Poisson and Poisson noise
given the exposures), which is what an observed future rate will show.
Models of observed rates ([`rwd()`](../reference/rwd.md),
[`arima_series()`](../reference/arima_series.md),
[`ets_series()`](../reference/arima_series.md)) forecast observed rates
directly, and `predictive` adds nothing. Backtests use predictive draws
because they are compared with observed rates. For life expectancy and
other summaries of the underlying mortality, use the default draws.

### How do I handle small populations and zero deaths?

Models that take logs of observed rates (the random walk, the SVD
version of Lee-Carter, Li-Lee) replace zero deaths by 0.5, which is
harmless for a few cells and distorting for many. The Poisson Lee-Carter
gives age groups without any death a constant rate and shrinks the age
pattern of change where the counts are too sparse, with a warning. With
many small cells:

- model at a coarser resolution (age groups instead of single years, or
  wider groups at young ages,
  [`mt_regroup()`](../reference/mt_regroup.md));
- pool information across strata, for example
  `li_lee(common_by = "sex")` for regions or deprivation groups, or
  borrow strength across ages and years within a stratum with the
  smoothing priors of [`apc_bayes()`](../reference/apc_bayes.md);
- smooth before fitting with [`step_smooth()`](../reference/mt_prep.md)
  (age 0 is left out by default);
- use `rwd(observation_noise = TRUE)` instead of
  [`rwd()`](../reference/rwd.md), alone or in an ensemble: the plain
  random walk treats the Poisson noise of observed rates as real change,
  and its paths over 10 or more years can reach impossible rates;
- combine rare causes of death into larger groups (see
  `ons_cause_groups` and
  [`mt_group_levels()`](../reference/mt_group_levels.md)).

[`mt_forecast()`](../reference/mt_forecast.md) warns when forecast death
rates exceed 1, the clearest sign that a model is unstable for the data.

### Which function aggregates what?

| Function | Combines | Example |
|----|----|----|
| [`mt_regroup()`](../reference/mt_regroup.md) | ages or periods | single years to 5-year groups, or to the groups of another table |
| [`mt_group_levels()`](../reference/mt_group_levels.md) | levels of one stratum variable, replacing them | deciles to quintiles, causes to cause groups |
| [`mt_aggregate_strata()`](../reference/mt_aggregate_strata.md) | strata along a hierarchy, adding totals and groups as new strata | regions to East, West and Germany, sexes to both sexes |
| [`mt_aggregate_draws()`](../reference/mt_aggregate_draws.md) | forecast draws along a hierarchy or to wider age groups | forecasts for Germany from forecasts by Land |

### How should I treat the COVID-19 years?

Either end the training data in 2019 (as in the ONS vignettes), or mark
the years as shocks so that the trend is estimated without them:

``` r

prep <- mt_prep(step_shocks(2020:2022, action = "interpolate"))
fit <- mt_fit(lee_carter(), data, prep = prep)
```

Whether mortality returns to the pre-pandemic trend is an assumption the
model does not test.

## Coherence

### How do I get forecasts that add up over regions, sexes or causes?

Describe the structure with
[`mt_hierarchy()`](../reference/mt_hierarchy.md) and either model the
bottom strata and sum them (`reconciled(..., method = "bottom_up")`), or
model every node and reconcile the draws (`method = "wls"`, `"mint"`,
…). Totals that are not in the data are added by
[`mt_aggregate_strata()`](../reference/mt_aggregate_strata.md). To keep
a trusted total, such as an all-cause forecast that the causes must add
up to, give it as `fixed` and use `method = "proportional"` (its parts
are scaled in proportion to their forecasts) or `"wls"` (the most
uncertain parts absorb the difference). Over long horizons, sums of
cause forecasts tend to be more pessimistic than forecasts of all causes
(Wilmoth 1995), so compare the two before choosing. Reconciliation of
deaths uses exposures of the forecast years; supply population
projections as `exposure` when they are available. See “Coherent
forecasts: totals and their parts” and the applied vignettes “Inequality
in England” and “Differences between the German Laender”.

### Will the gap between regions or deprivation groups keep growing?

That depends on the model, not only on the data. Separate models per
stratum extrapolate each stratum’s trend, so gaps can widen without
bound. [`li_lee()`](../reference/li_lee.md) and
[`product_ratio()`](../reference/li_lee.md) let strata deviate from a
common trend only through stationary processes, so gaps return towards
their past average. Backtests over a few years rarely separate the two;
the choice is an assumption to state.

### Which ages and causes explain a gap, and how much is composition?

[`mt_decompose()`](../reference/mt_decompose.md) splits a difference in
life expectancy (or in the probability of dying, an age-standardised or
crude rate, or any function of the rates) into contributions of ages and
causes: between strata in one year (`var`, `ref`), between years
(`periods`, which may combine an observed and a forecast year), or
between two forecasts (`y`). With `groups`, it also separates the
composition of a population of subgroups from the change of their
mortality: for groups defined at birth or for life, such as sex or
education, as in Riffe, Tursun-zade and Trias Llimós (2024)
(`composition = "radix"`); for areas such as deprivation quintiles,
whose populations move between groups, with `composition = "pooled"`.
Forecasts are decomposed draw by draw, so contributions have intervals;
see “Decomposing differences in life expectancy”.

## Practicalities

### Are forecasts reproducible?

Yes, with a `seed` (`mt_forecast(..., seed = 1)`,
`mt_backtest(..., seed = 1)`). Backtests and
[`recalibrated()`](../reference/recalibrated.md) seed every origin
separately, so their results do not depend on `cores`. Fits of several
strata in parallel (`mt_fit(..., cores = )`) run in forked processes;
for models with random steps in the fit, set `RNGkind("L'Ecuyer-CMRG")`
and a seed before the call to make them reproducible.

### How long does it take?

Lee-Carter, Li-Lee and random walks fit in seconds even for hundreds of
strata. One-dimensional PCLM ungrouping takes about two seconds per
stratum for 30 years of single ages; the two-dimensional version is much
slower. [`apc_bayes()`](../reference/apc_bayes.md) runs MCMC chains:
with `mcmc_control("fast")` from about ten seconds (8 age groups, 10
years) to one or two minutes per stratum (19 age groups, 14 years).

### Where do the validation numbers come from?

The `validation/` directory of the source repository runs benchmarks
with thresholds on every change (ungrouping accuracy, interval
calibration, coherence). The vignette “Comparison with demography, bamp
and StMoMo” reports accuracy and interval coverage of the models on four
data sets.

## References

Wilmoth JR (1995). Are mortality projections always more pessimistic
when disaggregated by cause of death? *Mathematical Population Studies*
5(4): 293-319. <doi:10.1080/08898489509525409>.
