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
as rates with the exposures as offset. Its estimates are the fitted
values of the model: they reproduce the total over all age groups of a
fit, and the diagnostics report the largest gap to a single group total
(`max_rel_diff`). `preserve_totals = TRUE` rescales the estimates to
every group total, a post-hoc rescaling of this package that creates
jumps at the group boundaries. `method = "monotone"` interpolates the
cumulative counts and reproduces every group total, without using the
exposures. Check your own result with
[`mt_plot_ungroup()`](../reference/ungroup-diagnostics.md) and
[`mt_ungroup_check()`](../reference/ungroup-diagnostics.md); see
“Ungrouping age groups: methods and diagnostics”.

### Can I split the open interval 90+ into single years?

Only from assumptions: the data say nothing about the age pattern inside
the interval. [`mt_ungroup()`](../reference/mt_ungroup.md) warns when it
extrapolates. Prefer a grid whose open interval starts where the data’s
does, or data with a higher open age (for example exposures to 100+ with
deaths to 90+, where the exposures carry the age pattern). When the grid
itself has an open interval (for example 100+),
[`mt_ungroup()`](../reference/mt_ungroup.md) splits its exposure over
single years up to the maximum age `omega` by the PCLM of the exposures
(`open_split = "pclm"`, default), as Rizzi, Gampe and Eilers (2015)
ungroup grouped exposures; `open_split = "stationary"` uses the weights
of a stationary population with the fitted rates, a method of this
package. On French data the rate of 100+ was too high by 0.12 to 0.15 on
the log scale with either split (see
[`?mt_ungroup`](../reference/mt_ungroup.md)).

## Models and uncertainty

### Which model should I use?

There is no single best model; a rolling-origin backtest on your own
data ([`mt_backtest()`](../reference/mt_backtest.md), horizons like the
ones you will report) decides. As a starting point:

| Situation | Try first | Why |
|----|----|----|
| One population, many deaths per cell, single years or 5-year groups | [`lee_carter()`](../reference/lee_carter.md), [`rwd()`](../reference/rwd.md) and `ensemble(lc = lee_carter(), rw = rwd())` | Lee-Carter is the standard; single models give too narrow intervals on real data, the ensemble is better calibrated |
| Several strata that should not diverge (sexes, regions, deprivation groups) | `li_lee(common_by = ...)` | a common trend with bounded deviations; pools information across strata |
| Small populations, many cells with few deaths | `li_lee(common_by = ...)`; adults: `saint(reference = ...)` against a national total; wider age groups ([`mt_regroup()`](../reference/mt_regroup.md)), [`step_smooth()`](../reference/mt_prep.md); `arima_series(order = c(0, 1, 1), drift = TRUE)` in ensembles | pooling and smoothing reduce noise; SAINT takes the trend from a large reference population; the plain random walk accumulates it, the random walk observed with noise separates it from the trend |
| Totals and parts must add up | `reconciled(model, hierarchy)` | bottom-up for joint models, WLS for separate models |
| Totals or gaps of strata forecast by separate models | `coupled(model, template = li_lee())`, `reconciled(model, hierarchy, dependence = "empirical")` or `"schaake"` (the dependence of past years) or a joint model such as [`cae()`](../reference/cae.md) ([`vignette("dependence")`](../articles/dependence.md)) | separate models give independent errors, so intervals of totals are too narrow and those of gaps too wide |
| More variation of the deaths than Poisson | `lee_carter(family = "nb")`, `stanmomo(family = "nb")` | negative binomial deaths (Delwarde, Denuit and Partrat 2007) widen the predictive intervals |
| Causes of death | separate or Li-Lee models per cause, `reconciled(..., method = "bottom_up")` with `cause = list(total = "all", additive = "rates")`; [`coda_causes()`](../reference/coda_causes.md) (experimental) | cause rates add up to the all-cause rate |
| Cohort effects (for example generations born around 1925-1945) | [`apc_bayes()`](../reference/apc_bayes.md), `stmomo("apc")`; at older ages `stmomo("m7")`, `stmomo("plat")` or [`apci()`](../reference/apci.md) | age-period models attribute cohort patterns to period trends |
| Cause-specific mortality such as cancer, in 5-year age groups | [`nordpred_apc()`](../reference/nordpred_apc.md) | the Nordpred age-period-cohort model (Møller et al. 2003) damps the trend in later periods; its intervals contain the uncertainty of the coefficients only, so use it for point projections or recalibrate it |
| Trend change in the recent past | `step_window(years)` with a backtest of several windows; [`p_splines()`](../reference/p_splines.md) | a shorter fitting period follows the recent trend; P-splines continue the trend of the last years of each age |
| Smooth age patterns, forecasts that keep the past age pattern and pace of change | [`cp_splines()`](../reference/cp_splines.md) | constrained P-splines (Camarda 2019) keep the forecast rates of ageing and of change within their past ranges; the intervals contain the uncertainty of the fit only |
| Long horizons, as mortality falls most at ever older ages | `lee_carter(rotate = TRUE)`, [`coda()`](../reference/coda.md), [`bohk_rau()`](../reference/bohk_rau.md) | the rotation moves the age pattern of decline to older ages as life expectancy rises (Li, Lee and Gerland 2013); CoDa forecasts the distribution of deaths by age (Oeppen 2008); the improvements of each age follow their own trend (Bohk-Ewald and Rau 2017) |

Check every fit with `summary(fit)`, `mt_parameters(fit)` and
`mt_plot_residuals(fit)`, and every forecast for warnings about
implausible rates.

### Why are my 95% intervals too narrow (or too wide)?

Intervals of most models describe the uncertainty of their own
assumptions, for example a constant drift, and their predictive draws
add the Poisson variation of the deaths only. On real data, changes of
trend and year effects add uncertainty that these assumptions leave out,
so intervals tend to be too narrow over 5 to 10 years: in the backtests
of the vignettes the 95% coverage of Lee-Carter was 0.79 on the ONS
quintiles (horizons 1 to 5, “Getting started with mortensor”) and 0.63
on England & Wales men (horizons 1 to 10, “Tuning and recalibration”).
The plain random walk is wider, because year-to-year changes of observed
rates contain observation noise twice; the local level model with drift,
`arima_series(order = c(0, 1, 1), drift = TRUE)`, separates the noise
from the trend (Harvey 1989). Some models leave out sources of
uncertainty by design: the draws of
[`cp_splines()`](../reference/cp_splines.md) contain the uncertainty of
the fit only (Camarda 2019 notes narrow intervals), those of
[`nordpred_apc()`](../reference/nordpred_apc.md) the uncertainty of the
coefficients only (95% coverage 0.11 in an example backtest on England &
Wales men), and the intervals of life expectancy of
[`p_splines()`](../reference/p_splines.md) are narrow (coverage 0.40 at
age 30 in the same backtest). The draws of
[`bohk_rau()`](../reference/bohk_rau.md) apply one quantile of the
improvements to all ages and years, as in Bohk-Ewald and Rau (2017);
`draws = "samples"` gives much narrower intervals of life expectancy.
Mixtures of models ([`ensemble()`](../reference/ensemble.md)) widen the
intervals where the models disagree, and the negative binomial family of
[`stanmomo()`](../reference/stanmomo.md) allows for overdispersion.
Death counts often vary more than Poisson counts:
`lee_carter(family = "nb", dispersion = "age")` estimates one dispersion
per age (Li, Hardy and Tan 2009) and widens the intervals where the
variation is large, and `p_splines(overdispersion = "age")` lets the
dispersion vary smoothly over age (Djeundje and Currie 2011); the
vignette “Overdispersion: extra variation in death counts” compares
them. [`recalibrated()`](../reference/recalibrated.md) rescales the
spread of any model’s draws to its forecast errors in the last years of
the training data (by default ten); on England & Wales men it raised the
95% coverage of Lee-Carter from 0.63 to 0.85. Check calibration with
[`mt_backtest()`](../reference/mt_backtest.md) and
[`mt_plot_coverage()`](../reference/mt_plot.md), and by horizon with
`mt_score_table(by = c("model", "horizon"))`. Where coverage falls with
the horizon, `recalibrated(scale = "per_horizon")` estimates one factor
per horizon, and `recalibrated(method = "errors")` adds resampled errors
of the model’s past forecasts; both are opt-in.

### Is one model really more accurate than another?

A backtest gives mean scores from a few origins, and two models often
differ by less than these scores vary between origins.
[`mt_dm_test()`](../reference/mt_dm_test.md) tests whether the
difference is larger than chance, horizon by horizon, with the
Diebold-Mariano test (Diebold and Mariano 1995). By default the long-run
variance of the loss differences is the equal-weighted cosine estimate
with fixed-smoothing critical values (Lazarus et al. 2018), which allows
for the overlap of forecasts from neighbouring origins and keeps the
size of the test nearly correct with few origins:

``` r

bt <- mt_compare(list(lc = lee_carter(), rw = rwd()), data, origins = 1995:2009, h = 5)
mt_dm_test(bt, ref = "lc")
```

With few origins the test has few degrees of freedom and finds only
large and consistent differences. More origins, or origins spaced `h`
years apart, whose forecasts do not overlap, give more informative
tests. `varestimator = "acf"` gives the test of
[`forecast::dm.test()`](https://pkg.robjhyndman.com/forecast/reference/dm.test.html)
with the small-sample correction of Harvey, Leybourne and Newbold
(1997); with 7 overlapping origins it rejected a true null hypothesis in
10% to 19% of simulated samples at the nominal 5% level (see
[`?mt_dm_test`](../reference/mt_dm_test.md)). `varestimator = "daniell"`
gives the fixed-smoothing test of Coroneo and Iacone (2020) with a
Daniell-kernel estimate: it has more degrees of freedom, and its size
was close to 5% for horizons up to 2 and too large beyond (15% at
`h = 5` with 7 origins). A non-significant difference means that the
backtest cannot tell the models apart; they may still differ.

### How should I weight the models of an ensemble?

[`ensemble()`](../reference/ensemble.md) gives equal weights by default.
Estimated weights come from an inner backtest on the training years.
`weights = "crps_stack"` is stacking of predictive distributions (Yao et
al. 2018) with the CRPS: the weights minimise the mean CRPS of the
mixture over the inner forecasts. `weights = "bias"` gives each model
the weight `exp(-|bias|)`, with the bias of its projected life
expectancy over the last 13 training years (Kontis et al. 2017).
`weights = "crps"` takes weights inverse to the mean CRPS, an extension
to distributions of the inverse-error weights of Shang (2012). Estimated
weights rest on a few inner origins from one period and can be unstable
(Yao et al. 2018). Compare them with equal weights in an outer backtest:

``` r

ens <- function(w) ensemble(lc = lee_carter(), rw = rwd(), fdm = hyndman_ullah(), weights = w)
bt <- mt_compare(list(equal = ens("equal"), stacking = ens("crps_stack"), bias = ens("bias")), data,
                 origins = 1996:2001, h = 10, quantities = "life_expectancy")
mt_score(bt)
```

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
index: on the ONS quintiles, Whittaker smoothing over time lowered the
95% coverage of Lee-Carter from 0.80 to 0.63 and raised its CRPS from
0.067 to 0.072 (“Smoothing death rates”).
[`recalibrated()`](../reference/recalibrated.md) re-estimates the spread
of a smoothed model against the raw observations. Treat smoothing as
optional and check it in a backtest. Use it as a preprocessing step
(`mt_prep(step_smooth(...))`), never on the data before fitting
([`mt_backtest()`](../reference/mt_backtest.md) stops,
[`mt_fit()`](../reference/mt_fit.md) warns), so that backtests refit it
on each training window and residuals and scores refer to the raw rates.
[`mt_smooth_check()`](../reference/mt_smooth_check.md) and
[`mt_plot_smooth()`](../reference/mt_smooth_check.md) show what the
smoothing removed.

### Why does the forecast start with a step?

Most models start from their fitted rates of the last observed year.
When the fit misses the level of the last years, the first forecast year
jumps away from the observations;
[`mt_jump_off()`](../reference/mt_jump_off.md) reports the step, and the
forecast plots show it as a gap between the solid (fitted) and the
dashed (forecast) line. `anchored(model)` starts any model from the
observed rates of the last training year and keeps this offset in every
forecast year, the jump-off adjustment of Lee and Miller (2001);
`lee_carter(jump_off = "observed")` and `li_lee(jump_off = "observed")`
do the same within the model. The observed rates of one year carry their
random variation into every forecast year, which matters where deaths
are few. Three published options reduce it within
[`lee_carter()`](../reference/lee_carter.md): the mean observed log
rates of the last years (`jump_off_years`; Stoeldraijer et al. 2018),
the Kalman filter estimate of the index in the last year
(`jump_off = "filtered"`; De Jong and Tickle 2006) and a fit constrained
to the observed rates of the last year (`jump_off = "constrained"`; Li,
Hardy and Tan 2009). Stoeldraijer et al. (2018) found the last observed
year most accurate and an average over several years most robust to data
updates; `mt_score(measures = "robustness")` gives their measure of
robustness for a backtest. The other settings of
[`anchored()`](../reference/anchored.md) are package-specific: the
offset of the last `years` years, shrunk where deaths are few
(`shrink = TRUE`) and fading out over `decay` years. “Where a forecast
starts: the jump-off” compares the starting points in a backtest.
Decompositions of a change from an observed to a forecast year include
the step; [`mt_decompose()`](../reference/mt_decompose.md) reports it
(`jump_off`) and can start from the fitted rates
(`observed = "fitted"`).

### What is the difference between predictive draws and forecast rates?

For models of expected rates (Lee-Carter, Li-Lee, CoDa, Tucker3, FDM,
P-splines and CP-splines, Bohk-Ewald and Rau, StMoMo, StanMoMo, Bayesian
APC, Nordpred), [`mt_forecast()`](../reference/mt_forecast.md) returns
draws of the expected death rates by default, and `predictive = TRUE`
draws the deaths given these rates and the exposures, which is what an
observed future rate will show. The deaths are Poisson, the observation
model of the Poisson Lee-Carter model (Brouhns, Denuit and Vermunt
2002), or negative binomial for `stanmomo(family = "nb")`; no other
observation noise is added. Models of observed rates
([`rwd()`](../reference/rwd.md),
[`arima_series()`](../reference/arima_series.md),
[`ets_series()`](../reference/arima_series.md)) forecast observed rates
directly, and `predictive` adds nothing unless they were fitted to
smoothed rates. Backtests use predictive draws because they are compared
with observed rates. For life expectancy and other summaries of the
underlying mortality, use the default draws.

### How do I handle small populations and zero deaths?

Models that take logs of observed rates (the random walk, the SVD
version of Lee-Carter, Li-Lee) replace zero deaths by 0.5, which is
harmless for a few cells and distorting for many. The Poisson Lee-Carter
gives age groups without any death a constant rate and shrinks the age
pattern of change where the counts are too sparse, with a warning (both
are safeguards of this package). With many small cells:

- model at a coarser resolution (age groups instead of single years, or
  wider groups at young ages,
  [`mt_regroup()`](../reference/mt_regroup.md));
- pool information across strata, for example
  `li_lee(common_by = "sex")` for regions or deprivation groups, or
  borrow strength across ages and years within a stratum with the
  smoothing priors of [`apc_bayes()`](../reference/apc_bayes.md);
- smooth before fitting with [`step_smooth()`](../reference/mt_prep.md)
  (age 0 is left out by default);
- use the local level model with drift, a random walk observed with
  noise (`arima_series(order = c(0, 1, 1), drift = TRUE)`; Harvey 1989),
  instead of [`rwd()`](../reference/rwd.md), alone or in an ensemble:
  the plain random walk treats the Poisson noise of observed rates as
  real change, and its paths over 10 or more years can reach impossible
  rates;
- combine rare causes of death into larger groups (see
  `ons_cause_groups` and
  [`mt_group_levels()`](../reference/mt_group_levels.md)).

[`mt_forecast()`](../reference/mt_forecast.md) warns when forecast death
rates exceed 1, the clearest sign that a model is unstable for the data.

In a backtest, a cell without deaths counts as half a death. The scores
of such cells can dominate the mean scores of the log rates;
[`mt_score()`](../reference/mt_score.md) reports their number
(`n_floored`). Score life expectancy or the age-standardised rate
(`quantities` in [`mt_backtest()`](../reference/mt_backtest.md)), which
sum the deaths over age first, or use coarser age groups
([`step_regroup()`](../reference/mt_prep.md)).

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
bottom strata and sum them (`reconciled(..., method = "bottom_up")`,
which fits the model to the bottom strata only), or model every node and
reconcile the draws (`method = "wls"`, `"mint"`, …). Totals that are not
in the data are added by
[`mt_aggregate_strata()`](../reference/mt_aggregate_strata.md). To keep
a trusted total, such as an all-cause forecast that the causes must add
up to, give it as `fixed` and use `method = "proportional"` (its parts
are scaled in proportion to their forecasts) or `"wls"` (the most
uncertain parts absorb the difference). `method = "energy_score"` learns
the reconciliation from the base model’s past forecasts (Panagiotelis et
al. 2023); in the package’s backtests it was less accurate than `"wls"`,
so compare it in a backtest before using it. Over long horizons, sums of
cause forecasts tend to be more pessimistic than forecasts of all causes
(Wilmoth 1995), so compare the two before choosing:
[`mt_pessimism()`](../reference/mt_pessimism.md) measures the difference
in life expectancy and the contribution of each cause. Reconciliation of
deaths uses exposures of the forecast years; supply population
projections as `exposure` when they are available. See “Coherent
forecasts: totals and their parts” and the applied vignettes “Inequality
in England” and “Differences between the German Laender”. “Causes of
death by deprivation: coherent and combined forecasts” reconciles causes
by sex and quintile in one crossed hierarchy,
`mt_hierarchy(sex = "both", quintile = "England", cause = list(total = "all", additive = "rates"))`,
and compares the methods and a per-cause ensemble with Lee-Carter in a
backtest.

### Will the gap between regions or deprivation groups keep growing?

That depends on the model, not only on the data. Separate models per
stratum extrapolate each stratum’s trend, so gaps can widen without
bound. [`li_lee()`](../reference/li_lee.md) and
[`coherent_fdm()`](../reference/coherent_fdm.md) let strata deviate from
a common trend only through stationary processes, so gaps return towards
their past level; [`coherent_fdm()`](../reference/coherent_fdm.md)
weights recent years more in its principal components (`beta = 0.05`, as
in Hyndman, Booth and Yasmeen 2013), so that level is weighted towards
recent years. In `coda(coherent = TRUE)` the point forecasts of the
deviations from the common trend level off (ARIMA(1,1,0) models;
Bergeron-Boucher et al. 2017): the point forecasts of the distributions
of deaths by age do not diverge, and their intervals widen with the
horizon. `cae(index = "coherent")` gives all strata one common drift and
lets their other period effects revert to zero (Enchev, Kleinow and
Cairns 2017), so the expected rates do not diverge.
[`relative_model()`](../reference/relative_model.md) gives
subpopulations, such as deprivation quintiles, the trend of their
national population plus their own drifting improvement differences, so
a gap that widened keeps widening; with `deviation = "constant"` the
differences stay at their last values (Villegas and Haberman 2014).
Backtests over a few years rarely tell these models apart; the choice is
an assumption to state.

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
whose populations move between groups, with `composition = "pooled"`, a
parametrisation of this package without a published source. Forecasts
are decomposed draw by draw, so contributions have intervals; see
“Decomposing differences in life expectancy”.

## Describing mortality

### How unequal are lifespans?

Life expectancy is the average age at death;
[`lifespan_variation()`](../reference/lifespan_variation.md) gives its
spread: life disparity (e-dagger), the life-table entropy, the Gini
coefficient, the standard deviation and the interquartile range of the
age at death, from age 0 or conditional on survival to an age `at`. It
works on observed data (with Poisson intervals) and on forecast draws.
All measures depend on the open age interval, e-dagger most: for French
women in 2006 it was 2.1 years higher with an open interval from 85 than
with single years of age to 100+. Use data with an open interval from 95
or higher, and compare measures only between tables with the same open
age. See “Lifespan variation”.

### How do I estimate excess deaths?

Fit a baseline model to the years before the shock and compare the
observed deaths with its predictive draws for the later years:

``` r

ex <- mt_excess(lee_carter(), data, baseline = c(2010, 2019), ages = c(0, 65, 80))
mt_plot_excess(ex, type = "p_score")
```

The intervals contain the uncertainty of the baseline trend and the
random variation of the death counts. The result depends on the baseline
model and years, which are assumptions to state; compare a few in the
same call. See “Excess mortality”.

## Practicalities

### Are forecasts reproducible?

Yes, with a `seed` (`mt_forecast(..., seed = 1)`,
`mt_backtest(..., seed = 1)`). Backtests and
[`recalibrated()`](../reference/recalibrated.md) seed every origin
separately, so their results do not depend on `cores`. Fits of several
strata (`mt_fit(..., cores = )`) give each stratum its own random number
stream, started from the seed of the model (`lee_carter(seed = )` for
the bootstrap of its parameters) or else from
[`set.seed()`](https://rdrr.io/r/base/Random.html), so models with
random steps in the fit give the same results with any number of cores.
The chains of [`apc_bayes()`](../reference/apc_bayes.md) are the
exception: with `mcmc_control(parallel = TRUE)`, the default, bamp runs
them in a socket cluster without seeded streams, so results differ
slightly between runs ([`?mcmc_control`](../reference/mcmc_control.md)).

### How long does it take?

Lee-Carter, Li-Lee and random walks fit in seconds even for hundreds of
strata. One-dimensional PCLM ungrouping takes about two seconds per
stratum for 30 years of single ages; the two-dimensional version, and
the two-dimensional P-spline smoothing of
[`mt_smooth()`](../reference/mt_smooth.md), about 10 seconds per stratum
([`mt_example()`](../reference/mt_example.md): 71 grid ages, 30 years).
[`apc_bayes()`](../reference/apc_bayes.md) runs MCMC chains: with
`mcmc_control("fast")` from about ten seconds (8 age groups, 10 years)
to one or two minutes per stratum (19 age groups, 14 years).
[`bohk_rau()`](../reference/bohk_rau.md) takes about 6 seconds per
stratum for 71 ages and 36 years.
[`cp_splines()`](../reference/cp_splines.md) refits the constrained
forecast for each of its `n_boot` bootstrap samples (100 by default), so
forecasts with many draws take time.

### Where do the accuracy numbers come from?

The accuracy and coverage numbers in the help pages and vignettes are
example results of backtests on Human Mortality Database data (England &
Wales men from StMoMo, France from demography) or on simulated data.
They hold for these data only; a backtest on your own data
([`mt_backtest()`](../reference/mt_backtest.md)) shows how a model
performs there. The vignette “Benchmarks: how well the pipeline works”
reports interval calibration on simulated data with known truth,
backtests of the standard models on four real data sets and the
behaviour on defective input tables. The vignette “Comparison with
demography, bamp and StMoMo” compares the models with those packages.
The unit tests of the package check the methods against published
results and reference implementations.

## References

Bergeron-Boucher M-P, Canudas-Romo V, Oeppen J, Vaupel JW (2017).
Coherent forecasts of mortality with compositional data analysis.
*Demographic Research* 37(17): 527-566. <doi:10.4054/DemRes.2017.37.17>.

Bohk-Ewald C, Rau R (2017). Probabilistic mortality forecasting with
varying age-specific survival improvements. *Genus* 73(1): 1.
<doi:10.1186/s41118-016-0017-8>.

Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear
regression approach to the construction of projected lifetables.
*Insurance: Mathematics and Economics* 31(3): 373-393.

Camarda CG (2019). Smooth constrained mortality forecasting.
*Demographic Research* 41(38): 1091-1130.
<doi:10.4054/DemRes.2019.41.38>.

Coroneo L, Iacone F (2020). Comparing predictive accuracy in small
samples using fixed-smoothing asymptotics. *Journal of Applied
Econometrics* 35(4): 391-409. <doi:10.1002/jae.2756>.

De Jong P, Tickle L (2006). Extending Lee-Carter mortality forecasting.
*Mathematical Population Studies* 13(1): 1-18.
<doi:10.1080/08898480500452109>.

Diebold FX, Mariano RS (1995). Comparing predictive accuracy. *Journal
of Business & Economic Statistics* 13(3): 253-263.

Djeundje VAB, Currie ID (2011). Smoothing dispersed counts with
applications to mortality data. *Annals of Actuarial Science* 5(1):
33-52. <doi:10.1017/S1748499510000047>.

Enchev V, Kleinow T, Cairns AJG (2017). Multi-population mortality
models: fitting, forecasting and comparisons. *Scandinavian Actuarial
Journal* 2017(4): 319-342. <doi:10.1080/03461238.2015.1133450>.

Harvey AC (1989). *Forecasting, Structural Time Series Models and the
Kalman Filter*. Cambridge University Press.

Harvey D, Leybourne S, Newbold P (1997). Testing the equality of
prediction mean squared errors. *International Journal of Forecasting*
13(2): 281-291.

Hyndman RJ, Booth H, Yasmeen F (2013). Coherent mortality forecasting:
the product-ratio method with functional time series models.
*Demography* 50(1): 261-283. <doi:10.1007/s13524-012-0145-5>.

Kontis V, Bennett JE, Mathers CD, Li G, Foreman K, Ezzati M (2017).
Future life expectancy in 35 industrialised countries: projections with
a Bayesian model ensemble. *The Lancet* 389(10076): 1323-1335.
<doi:10.1016/S0140-6736(16)32381-9>.

Lazarus E, Lewis DJ, Stock JH, Watson MW (2018). HAR inference:
recommendations for practice. *Journal of Business & Economic
Statistics* 36(4): 541-559.

Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
method for forecasting mortality. *Demography* 38(4): 537-549.

Li JS-H, Hardy MR, Tan KS (2009). Uncertainty in mortality forecasting:
an extension to the classical Lee-Carter approach. *ASTIN Bulletin*
39(1): 137-164. <doi:10.2143/AST.39.1.2038060>.

Li N, Lee R, Gerland P (2013). Extending the Lee-Carter method to model
the rotation of age patterns of mortality decline for long-term
projections. *Demography* 50(6): 2037-2051.
<doi:10.1007/s13524-013-0232-2>.

Møller B, Fekjær H, Hakulinen T, Sigvaldason H, Storm HH, Talbäck M,
Haldorsen T (2003). Prediction of cancer incidence in the Nordic
countries: empirical comparison of different approaches. *Statistics in
Medicine* 22(17): 2751-2766. <doi:10.1002/sim.1481>.

Oeppen J (2008). Coherent forecasting of multiple-decrement life tables:
a test using Japanese cause of death data. CODAWORK’08, Universitat de
Girona.

Panagiotelis A, Gamakumara P, Athanasopoulos G, Hyndman RJ (2023).
Probabilistic forecast reconciliation: properties, evaluation and score
optimisation. *European Journal of Operational Research* 306(2):
693-706. <doi:10.1016/j.ejor.2022.07.040>.

Riffe T, Tursun-zade R, Trias Llimós S (2024). Arriaga meets Kitagawa.
Life expectancy decomposition with population subgroups. MPIDR Working
Paper WP-2024-029. <doi:10.4054/MPIDR-WP-2024-029>.

Rizzi S, Gampe J, Eilers PHC (2015). Efficient estimation of smooth
distributions from coarsely grouped data. *American Journal of
Epidemiology* 182(2): 138-147. <doi:10.1093/aje/kwv020>.

Shang HL (2012). Point and interval forecasts of age-specific life
expectancies: a model averaging approach. *Demographic Research* 27(21):
593-644. <doi:10.4054/DemRes.2012.27.21>.

Stoeldraijer L, van Duin C, van Wissen L, Janssen F (2018). Comparing
strategies for matching mortality forecasts to the most recently
observed data: exploring the trade-off between accuracy and robustness.
*Genus* 74: 16. <doi:10.1186/s41118-018-0041-y>.

Villegas AM, Haberman S (2014). On the modeling and forecasting of
socioeconomic mortality differentials: an application to deprivation and
mortality in England. *North American Actuarial Journal* 18(1): 168-193.
<doi:10.1080/10920277.2013.866034>.

Wilmoth JR (1995). Are mortality projections always more pessimistic
when disaggregated by cause of death? *Mathematical Population Studies*
5(4): 293-319. <doi:10.1080/08898489509525409>.

Yao Y, Vehtari A, Simpson D, Gelman A (2018). Using stacking to average
Bayesian predictive distributions (with discussion). *Bayesian Analysis*
13(3): 917-1007. <doi:10.1214/17-BA1091>.
