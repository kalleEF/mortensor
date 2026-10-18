# Model reference

This vignette states, for every model of mortensor, the model equation,
the estimation method, how the forecast is generated, the defaults, and
which sources of uncertainty the forecast draws contain. It is meant as
the basis of a methods section.
[`summary()`](https://rdrr.io/r/base/summary.html) of a fitted model
prints the same information for the actual fit,
[`mt_parameters()`](../reference/mt_parameters.md) returns the
estimates.

``` r

library(mortensor)
x <- mt_example(deaths = "single", years = 1995:2019)
summary(mt_fit(lee_carter(), x))
#> Lee-Carter
#>   data: 2 strata; ages 30 to 100+; periods 1995-2019
#>   model and uncertainty:
#>     log m(x,t) = a(x) + b(x) k(t), estimated by Poisson maximum likelihood
#>     (Brouhns et al. 2002; Lee and Carter 1992). k(t) follows a random walk with
#>     drift; draws include its innovations and the uncertainty of the drift. a(x)
#>     and b(x) are fixed (parameter uncertainty with bootstrapped()). Jump-off:
#>     fitted rates.
#>   fit: Poisson deviance 3371.3 on 3550 cells (0.95 per cell); 3.9% of deviance residuals beyond +/-2
```

Notation: `m(x,t,i)` is the central death rate at age `x` in year `t` in
stratum `i` (deaths divided by exposure), `D` deaths, `E` exposures.

## Data steps that affect every model

| Step | What mortensor does | Function |
|----|----|----|
| Exposures | mid-year population, or the mean of two 1 January or 31 December populations | `mt_as_exposure(type = )` |
| Suppressed cells | raking (iterative proportional fitting; Deming and Stephan 1940) of the missing counts to published totals, within optional bounds | [`mt_fill_suppressed()`](../reference/mt_fill_suppressed.md) |
| Unknown ages | distributed over the age groups of the same stratum and year in proportion to the known deaths | [`mt_data()`](../reference/mt_data.md) |
| Ungrouping | penalised composite link model with exposures (Rizzi, Gampe and Eilers 2015) or spline methods | [`mt_ungroup()`](../reference/mt_ungroup.md), [`step_ungroup()`](../reference/mt_prep.md) |
| Smoothing | over age (Poisson P-splines; Whittaker; GAM), over time (Whittaker) or the age-period surface (P-splines, Currie, Durban and Eilers 2004; quasi-Poisson `te(age, year)` GAM, Wood 2017); age 0 is not smoothed; models are fitted to the smoothed rates, their observation noise, residuals and backtest scales come from the raw rates; smoothing over time narrows the intervals of models with a period index (check coverage, [`recalibrated()`](../reference/recalibrated.md)) | [`mt_smooth()`](../reference/mt_smooth.md), [`step_smooth()`](../reference/mt_prep.md), [`mt_smooth_check()`](../reference/mt_smooth_check.md) |
| Zero deaths | models on log rates use 0.5 deaths in a cell without deaths |  |

All preprocessing steps given as `prep` are refitted on every training
window of a backtest.

## Single-population models

### Lee-Carter: `lee_carter()`

`log m(x,t) = a(x) + b(x) k(t)`, with `sum b(x) = 1` and `sum k(t) = 0`
(Lee and Carter 1992).

- **Estimation**: Poisson maximum likelihood (default; Brouhns, Denuit
  and Vermunt 2002) or SVD of the log rates (`method = "svd"`, as
  [`demography::lca()`](https://pkg.robjhyndman.com/demography/reference/lca.html)),
  with `k(t)` re-estimated to match the observed total deaths
  (`adjust = "dt"`), the deaths by age (`"dxt"`, Booth, Maindonald and
  Smith 2002) or life expectancy (`"e0"`, Lee and Miller 2001).
  `choose_period = TRUE` fits the model to the longest recent period
  over which `k(t)` was close to linear (Booth, Maindonald and Smith
  2002). Age groups without any death get a constant rate of 0.5 deaths
  over the training period and `b(x) = 0`; when the Poisson likelihood
  has no finite maximum (very sparse counts), `b(x)` is shrunk towards
  equal shares by a ridge penalty. Both are reported as warnings.
- **Forecast**: `k(t)` follows a random walk with drift; the drift and
  the innovation variance are estimated from the year-to-year changes of
  `k(t)`.
- **Jump-off**: fitted rates of the last year (default) or observed
  rates (`jump_off = "observed"`).
- **Uncertainty in the draws**: innovations of `k(t)` and the
  uncertainty of its drift. `a(x)` and `b(x)` are fixed;
  `bootstrap = TRUE` or [`bootstrapped()`](../reference/bootstrapped.md)
  adds their uncertainty.
- **Rotation** (`rotate = TRUE`; Li, Lee and Gerland 2013): as life
  expectancy at birth rises from `e0_start` (80) to `e0_end` (102), the
  age pattern of decline moves from `b(x)` to an ultimate pattern with
  equal decline below 65 (the mean of `b` over ages 15-64) and the shape
  of `b(x)` above. The weight follows a sine curve with `power` 0.5. In
  each path and year the index is solved again so that life expectancy
  at birth stays that of the unrotated forecast; rotation therefore
  changes the age pattern of the forecast only. It needs ages from 0 and
  computes life tables for every path and year, which is about 20 times
  slower.

### Compositional data model: `coda()`

`clr(d(x,t) / a(x)) = b(x) k(t)` (Oeppen 2008; Bergeron-Boucher et al.
2017), where `d(x,t)` are the life-table deaths (radix 1), `a(x)` their
closed geometric mean over the years and `clr` the centred log-ratio
transformation. The model forecasts the distribution of deaths by age,
so the forecast deaths always form a life table.

- **Estimation**: first component of the SVD of the centred log-ratios,
  as `model.Oeppen()` of the R package MortalityForecast 0.8.1
  (Pascariu). The log rate of the open age interval is regressed on
  `k(t)`.
- **Forecast**: `k(t)` random walk with drift; deaths are transformed
  back to rates through the life table.
- **Jump-off**: observed deaths of the last year, changed by the model
  (default), or fitted deaths (`jump_off = "fitted"`).
- **Coherent version** (`coherent = TRUE` or `common_by`):
  `clr(d(x,t,i) / a(x,i)) = B(x) K(t) + b(x,i) k(t,i)`, with the common
  component from the average distribution of deaths of the strata;
  `K(t)` random walk with drift, `k(t,i)` AR(1) processes with mean, so
  strata stay close to the common trend (`model.OeppenC()`).
- **Uncertainty**: innovations and drift of `k(t)` (and `K(t)`),
  innovations of the AR(1) processes. `a` and `b` are fixed. All ages
  have equal weight on the clr scale, so the noise of ages with few
  deaths enters `k(t)`: intervals can be wider than those of Lee-Carter.

### Random walk with drift: `rwd()`

Each age: `log m(x,t) = log m(x,t-1) + mu(x) + e(x,t)`.

- **Estimation**: drift as the mean year-to-year change (over the last
  `drift_window` years); innovations are resampled vectors of observed
  changes across ages (`innovations = "bootstrap"`, keeps the
  correlation between ages) or normal (`"gaussian"`).
- **Uncertainty**: innovations and drift. The observed changes contain
  the Poisson noise of two years, so the intervals are wide where deaths
  are few, and over long horizons very wide.
- **`observation_noise = TRUE`**: log rate = level + noise, the level a
  random walk with drift (local level model). The noise variance is the
  Poisson part (one over the deaths) plus an extra part common to all
  ages; the level variance is estimated from the variance and first
  autocovariance of the changes and partly pooled across ages (ages with
  few deaths borrow from the median of the others). The level is Kalman
  filtered; draws include level, drift and, for predictive draws, noise.

### Functional demographic model: `hyndman_ullah()`

Calls
[`demography::fdm()`](https://pkg.robjhyndman.com/demography/reference/fdm.html)
(Hyndman and Ullah 2007): the rates of each year are smoothed over age
with
[`demography::smooth.demogdata()`](https://pkg.robjhyndman.com/demography/reference/smooth.demogdata.html)
(`smooth = TRUE`; settings in `smooth_args`), and the log rates are
decomposed into a mean function and functional principal components,
whose scores are forecast with time-series models. The number of
components and the smoothing settings can be chosen by forecast accuracy
with [`tuned()`](../reference/mt_tune.md). Draws are simulated from the
fitted model. Its intervals are those of the simulated paths, which are
somewhat narrower at short horizons than demography’s analytic intervals
(see the comparison vignette).

### Time series per age: `arima_series()`, `ets_series()`

Benchmarks that model the observed log rate of every age as a separate
time series with the forecast package (Hyndman and Khandakar 2008): an
ARIMA model chosen by
[`forecast::auto.arima()`](https://pkg.robjhyndman.com/forecast/reference/auto.arima.html)
(drift allowed) or an exponential smoothing model chosen by
[`forecast::ets()`](https://pkg.robjhyndman.com/forecast/reference/ets.html)
(damped trends allowed). **Uncertainty**: simulated future paths of the
fitted models, so the intervals grow with the horizon; ages are
simulated independently and parameter uncertainty is not included. The
series are observed rates, so their draws already contain observation
noise.

### StMoMo models: `stmomo()`

Calls [`StMoMo::fit()`](https://rdrr.io/pkg/StMoMo/man/fit.html) and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) (Villegas, Kaishev
and Millossovich 2018) for the Lee-Carter, CBD, APC, RH, M6, M7 and M8
families. M7 adds a quadratic age term and a cohort effect to CBD; M8 a
cohort effect that falls with age to zero at age `xc` (Cairns et al.
2009). The cohort effects of M6 to M8 follow an ARIMA model
(`gc_order`). Draws contain the period and cohort index uncertainty that
StMoMo simulates; parameter uncertainty is not included.

### Bayesian StMoMo models: `stanmomo()`

The same families (Lee-Carter, RH, APC, CBD, M6) fitted with Stan by
StanMoMo (Barigou et al. 2023), with a Poisson or negative binomial
likelihood. Period indexes are random walks with drift (CBD and M6:
correlated), cohort effects AR(2) processes.

- **Forecast**: StanMoMo forecasts inside the fit, up to `h_max` years.
- **Uncertainty**: posterior draws of all parameters and of the future
  indexes; predictive draws add Poisson and extra-Poisson noise.
- **Convergence**: the fit warns when the largest R-hat exceeds 1.05 or
  transitions diverged; models with cohort effects need longer chains.

### Bayesian age-period-cohort model: `apc_bayes()`

Calls
[`bamp::bamp()`](https://volkerschmid.github.io/bamp/reference/bamp.html)
(Schmid and Held 2007): logit death probability = age + period + cohort
effects with random-walk priors of order 1 or 2, binomial likelihood,
MCMC. Convergence is checked with the Gelman-Rubin statistic of the
identified quantities (reported by
[`mt_parameters()`](../reference/mt_parameters.md) as `rhat`); chains
that have not converged are refitted with longer runs. **Uncertainty**:
posterior samples of the forecast rates, including the uncertainty of
the parameters and of the random walks of period and cohort effects;
predictive draws add Poisson and extra-Poisson noise as for the other
models of smooth rates. Each stratum is fitted separately.
[`apc_select()`](../reference/apc_select.md) compares random-walk orders
by DIC.

## Models for several populations

### Li-Lee: `li_lee()`

`log m(x,t,i) = a(x,i) + B(x) K(t) + b(x,i) k(t,i)` (Li and Lee 2005).

- **Estimation**: the common factor `B(x) K(t)` is a Lee-Carter model
  (SVD, `k` refitted to deaths) of the pooled deaths and exposures of
  the strata in a group (`common_by`, for example one common factor per
  sex); the stratum factors come from an SVD of the residual log rates.
- **Forecast**: `K(t)` random walk with drift; each `k(t,i)` an AR(1)
  process without intercept, `|phi| <= 0.99`, so the deviations of a
  stratum from the common trend stay bounded (“coherent” in the sense
  that the strata do not diverge).
- **Uncertainty**: innovations and drift of `K(t)`, innovations of the
  AR(1) processes (independent across strata). `a`, `b`, `B` and `phi`
  are fixed.

### Product-ratio: `product_ratio()`, `coherent_fdm()`

The product-ratio method (Hyndman, Booth and Yasmeen 2013) models the
geometric mean of the strata (the product, which carries the common
trend) and the ratio of each stratum to it.
[`coherent_fdm()`](../reference/coherent_fdm.md) is the original method
with the functional data models of demography:

- **Estimation**: the rates of each stratum are smoothed over age; the
  product and the ratios are decomposed into principal components
  ([`demography::fdm()`](https://pkg.robjhyndman.com/demography/reference/fdm.html)).
  With `common_by`, each group has its own product; in a hierarchy only
  the bottom strata enter the product.
- **Forecast**: the scores of the product follow ARIMA models, those of
  the ratios stationary ARFIMA models, so the ratios return towards
  their recent level and the strata do not diverge.
- **Uncertainty**: simulated scores and resampled model errors of the
  product and of each ratio; the ratios of each path are rescaled to
  multiply to one, as in
  [`demography::forecast.fdmpr()`](https://pkg.robjhyndman.com/demography/reference/forecast.fdmpr.html).

[`product_ratio()`](../reference/li_lee.md) uses Lee-Carter components
in place of functional data models.

### Tucker3 and cause shares (experimental)

[`tucker3()`](../reference/tucker3.md) decomposes the array of log rates
`[age, period, stratum]` (or of the log-ratio coordinates of life-table
deaths, `transform = "coda"`) into factor matrices and a core; the
period factors are forecast as a multivariate random walk with drift or
with ARIMA/ETS models, and the in-sample residual variation is added to
the draws. [`cause_shares()`](../reference/cause_shares.md) forecasts
all-cause mortality with any model and the shares of causes with a
Tucker model of their isometric log-ratio coordinates, so the causes add
up to all causes in every draw. Neither has a published evaluation for
mortality; compare them with the models above in a backtest.

[`coda_causes()`](../reference/coda_causes.md) is the multiple-decrement
form of [`coda()`](../reference/coda.md) that Oeppen (2008) proposed for
causes of death: the deaths of a multiple-decrement life table by age
and cause form one composition with one index `k(t)`, so the causes add
up to all causes in every draw. Only the abstract of Oeppen (2008) is
published, so the implementation follows its natural reading; in a
backtest on simulated causes it was less accurate than
[`cause_shares()`](../reference/cause_shares.md).

## Jump-off: `anchored()`

Forecasts start from the fitted rates of the last year for most models.
[`anchored()`](../reference/anchored.md) moves the start towards the
observed level of the last `years` years for any model: the log ratio of
observed to fitted deaths per age and stratum, shrunk towards the
stratum mean where deaths are few, multiplies the draws with a weight
that fades out over `decay` years.
[`mt_jump_off()`](../reference/mt_jump_off.md) reports the step between
the last observed (and fitted) year and the first forecast year.

## Combining models

### Ensembles: `ensemble()`

A mixture of the predictive distributions of its components: each draw
comes from one component, for all strata at once, chosen with
probabilities equal to the weights (equal, inverse CRPS of an inner
backtest, or stacking). The mixture contains the uncertainty of each
component and the disagreement between them.

### Different models for groups of strata: `stratified()`

One model per level of a stratum variable, for example Lee-Carter for
all-cause mortality and Li-Lee for the causes of death. The groups are
fitted and forecast separately and their draws combined; inside
[`reconciled()`](../reference/reconciled.md) they can be made coherent,
for example with the causes scaled to the all-cause forecast.

### Reconciled forecasts: `reconciled()`

Forecasts of all nodes of a hierarchy (bottom strata and their totals)
are made to add up draw by draw: bottom-up sums, or projections of the
base forecasts with weights from their variances (WLS), the shrunk
covariance (MinT; Wickramasuriya, Athanasopoulos and Hyndman 2019), the
structure (structural scaling; Athanasopoulos et al. 2017) or equal
weights (OLS; Hyndman et al. 2011). Parts that a projection would make
negative are held at zero and the rest projected again, so that `fixed`
nodes (a trusted total) are kept; with fixed nodes, `"proportional"`
instead scales their parts in proportion to their forecasts (top-down
allocation). Deaths add up over every hierarchy variable; exposures add
up over variables that divide the population and are shared by the
levels of variables such as causes of death (`additive = "rates"`). The
exposures of the forecast years are those of the last observed year
unless projections are given.

Score-optimal reconciliation (`method = "energy_score"`; Panagiotelis et
al. 2023) learns the reconciliation from the base model’s own past
forecasts: the reconciled draws are `S (d + G y)`, where `y` are the
base draws of all nodes, and `d` and `G` (per age) minimise the energy
score of forecasts from `n_origins` inner origins in the last `years` of
the training data, by stochastic gradient descent from OLS. It has many
parameters and few training forecasts; in the package’s backtests it was
less accurate than WLS.

### Calibration on recent years: `recalibrated()`

Refits a model at several origins in the last `years` (default 10) of
the training data, compares its forecasts with the observed rates of
that window, and rescales the spread of the draws around their median
per cell by a factor `s(h) = s1 h^gamma` that minimises the CRPS (or the
interval score of one interval), with the Poisson and extra-Poisson
variation of the observations as in a backtest (minimum-CRPS estimation;
Gneiting et al. 2005). Optionally a drift correction is estimated as
well. Ranks within each cell are kept, so the draws remain paths.
Calibration widens models that were persistently too narrow (Lee-Carter
on England & Wales: 95% coverage 0.75 to 0.88 in backtests) and can
narrow mixtures that were already wide. `ages` gives each age band its
own factors; `quantity = "life_expectancy"` or `"asr"` calibrates the
intervals of these quantities (which also depend on how the ages move
together); `n_boot` reports bootstrap intervals of the factors. Beyond
the longest horizon of the window the factors stay constant. Calibration
keeps the model and its settings; choosing the settings themselves is
tuning (below).

### Bootstrap and tuning: `bootstrapped()`, `tuned()`

[`bootstrapped()`](../reference/bootstrapped.md) refits a model to
Poisson resamples of the deaths (parameter uncertainty, including
preprocessing with `type = "data"`).
[`tuned()`](../reference/mt_tune.md) chooses settings by how well they
forecast the last years of the training data only
([`mt_tune()`](../reference/mt_tune.md) shows the scores): model
settings such as the number of components or the rank, and with
`spec_fun` returning `list(spec = , prep = )` preprocessing settings
such as the smoothing parameter or the fitting window.
`inner_origins = 1, inner_h = 10` is the hold-out of the last ten years
often used with
[`demography::fdm()`](https://pkg.robjhyndman.com/demography/reference/fdm.html);
several origins make the choice depend less on one period. The criterion
can be the CRPS, the interval score or point accuracy
(`metric = "rmse"`) of rates, life expectancy or ASR (`target`);
`select = "one_se"` takes the simplest setting within one standard error
of the best, and `by` tunes groups of strata separately. The fitting
period matters for Lee-Carter forecasts when the pace of decline changed
(Booth, Maindonald and Smith 2002); `step_window(years)` makes it a
preprocessing step that tuning compares across windows.

## Sources of uncertainty at a glance

| Model | Index innovations | Drift | Parameters | Poisson noise | Extra noise |
|----|----|----|----|----|----|
| [`lee_carter()`](../reference/lee_carter.md) | yes | yes | with `bootstrap = TRUE` | predictive draws | predictive draws |
| [`rwd()`](../reference/rwd.md) | yes (observed changes) | yes | no | contained in the changes | contained in the changes |
| `rwd(observation_noise = TRUE)` | yes (level) | yes | no | predictive draws | predictive draws |
| [`arima_series()`](../reference/arima_series.md), [`ets_series()`](../reference/arima_series.md) | yes | as the fitted model | no | contained in the series | contained in the series |
| [`li_lee()`](../reference/li_lee.md), [`product_ratio()`](../reference/li_lee.md) | yes (common and stratum) | common factor | no | predictive draws | predictive draws |
| [`coherent_fdm()`](../reference/coherent_fdm.md) | yes (product and ratios) | as fdm (product) | no | predictive draws | predictive draws |
| [`tucker3()`](../reference/tucker3.md) | yes (period factors) | yes | no | predictive draws | residual variation and predictive draws |
| [`cause_shares()`](../reference/cause_shares.md) | of the total model and the share factors | yes | no | predictive draws | residual variation and predictive draws |
| [`coda()`](../reference/coda.md), [`coda_causes()`](../reference/coda_causes.md) | yes | yes | no | predictive draws | predictive draws |
| [`hyndman_ullah()`](../reference/hyndman_ullah.md) | yes | as fdm | no | predictive draws | predictive draws |
| [`stmomo()`](../reference/stmomo.md) | yes | yes | no | predictive draws | predictive draws |
| [`stanmomo()`](../reference/stanmomo.md) | yes | posterior | yes (posterior) | predictive draws | predictive draws |
| [`apc_bayes()`](../reference/apc_bayes.md) | yes | posterior | yes (posterior) | predictive draws | predictive draws |
| [`ensemble()`](../reference/ensemble.md) | of the components | of the components | of the components, plus model choice | of the components | of the components |
| [`stratified()`](../reference/stratified.md) | of the model of each group | of the model of each group | of the model of each group | as the models | as the models |
| [`bootstrapped()`](../reference/bootstrapped.md) | of the base model | of the base model | yes | as the base model | as the base model |
| [`reconciled()`](../reference/reconciled.md) | of the base model, projected |  |  | bottom strata, summed | bottom strata, summed |
| [`recalibrated()`](../reference/recalibrated.md) | of the base model, rescaled | rescaled (optional correction) |  | as the base model | as the base model |
| [`anchored()`](../reference/anchored.md) | of the base model | of the base model |  | as the base model | as the base model |

Predictive draws (`predictive = TRUE`, always used in backtests)
describe observable rates. For models of smooth rates they add Poisson
deaths given the exposures and an extra-Poisson part estimated from the
in-sample residual variance. For models of observed rates
([`rwd()`](../reference/rwd.md) without `observation_noise`,
[`arima_series()`](../reference/arima_series.md),
[`ets_series()`](../reference/arima_series.md)) the draws already
describe observable rates and nothing is added. Forecast draws of models
of smooth rates without `predictive` describe the underlying rates.

## Evaluation

[`mt_backtest()`](../reference/mt_backtest.md) refits preprocessing and
model at each origin and scores predictive draws against the observed
log rates of the test years with the continuous ranked probability score
(CRPS; Gneiting and Raftery 2007), the weighted interval score (Bracher
et al. 2021), interval coverage, the probability integral transform and
the absolute error of the median. Cells with estimated (suppressed)
deaths are not scored unless `score_estimated = TRUE`.
[`mt_dm_test()`](../reference/mt_dm_test.md) tests differences in
accuracy between two models of the same backtest with the
Diebold-Mariano test (Diebold and Mariano 1995) and the correction of
Harvey, Leybourne and Newbold (1997).

## References

Athanasopoulos G, Hyndman RJ, Kourentzes N, Petropoulos F (2017).
Forecasting with temporal hierarchies. *European Journal of Operational
Research* 262(1): 60-74.

Barigou K, Goffard P-O, Loisel S, Salhi Y (2023). Bayesian model
averaging for mortality forecasting using leave-future-out validation.
*International Journal of Forecasting* 39(2): 674-690.

Bergeron-Boucher M-P, Canudas-Romo V, Oeppen J, Vaupel JW (2017).
Coherent forecasts of mortality with compositional data analysis.
*Demographic Research* 37(17): 527-566. <doi:10.4054/DemRes.2017.37.17>.

Booth H, Maindonald J, Smith L (2002). Applying Lee-Carter under
conditions of variable mortality decline. *Population Studies* 56(3):
325-336.

Bracher J, Ray EL, Gneiting T, Reich NG (2021). Evaluating epidemic
forecasts in an interval format. *PLoS Computational Biology* 17(2):
e1008618.

Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear
regression approach to the construction of projected lifetables.
*Insurance: Mathematics and Economics* 31(3): 373-393.

Cairns AJG, Blake D, Dowd K, Coughlan GD, Epstein D, Ong A, Balevich I
(2009). A quantitative comparison of stochastic mortality models using
data from England and Wales and the United States. *North American
Actuarial Journal* 13(1): 1-35.

Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
mortality rates. *Statistical Modelling* 4(4): 279-298.

Deming WE, Stephan FF (1940). On a least squares adjustment of a sampled
frequency table when the expected marginal totals are known. *Annals of
Mathematical Statistics* 11(4): 427-444.

Diebold FX, Mariano RS (1995). Comparing predictive accuracy. *Journal
of Business & Economic Statistics* 13(3): 253-263.

Gneiting T, Raftery AE, Westveld AH, Goldman T (2005). Calibrated
probabilistic forecasting using ensemble model output statistics and
minimum CRPS estimation. *Monthly Weather Review* 133(5): 1098-1118.

Gneiting T, Raftery AE (2007). Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical
Association* 102(477): 359-378.

Harvey D, Leybourne S, Newbold P (1997). Testing the equality of
prediction mean squared errors. *International Journal of Forecasting*
13(2): 281-291.

Hyndman RJ, Ahmed RA, Athanasopoulos G, Shang HL (2011). Optimal
combination forecasts for hierarchical time series. *Computational
Statistics & Data Analysis* 55(9): 2579-2589.

Hyndman RJ, Khandakar Y (2008). Automatic time series forecasting: the
forecast package for R. *Journal of Statistical Software* 27(3): 1-22.

Hyndman RJ, Booth H, Yasmeen F (2013). Coherent mortality forecasting:
the product-ratio method with functional time series models.
*Demography* 50(1): 261-283.

Hyndman RJ, Ullah MS (2007). Robust forecasting of mortality and
fertility rates: a functional data approach. *Computational Statistics &
Data Analysis* 51(10): 4942-4956.

Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
*Journal of the American Statistical Association* 87(419): 659-671.

Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
method for forecasting mortality. *Demography* 38(4): 537-549.

Li N, Lee R (2005). Coherent mortality forecasts for a group of
populations: an extension of the Lee-Carter method. *Demography* 42(3):
575-594.

Li N, Lee R, Gerland P (2013). Extending the Lee-Carter method to model
the rotation of age patterns of mortality decline for long-term
projections. *Demography* 50(6): 2037-2051.
<doi:10.1007/s13524-013-0232-2>.

Oeppen J (2008). Coherent forecasting of multiple-decrement life tables:
a test using Japanese cause of death data. CODAWORK’08, Universitat de
Girona.

Panagiotelis A, Gamakumara P, Athanasopoulos G, Hyndman RJ (2023).
Probabilistic forecast reconciliation: properties, evaluation and score
optimisation. *European Journal of Operational Research* 306(2):
693-706. <doi:10.1016/j.ejor.2022.07.040>.

Pascariu MD. MortalityForecast: Standard tools to compare and evaluate
various mortality forecasting methods. R package version 0.8.1.
<https://github.com/mpascariu/MortalityForecast>

Rizzi S, Gampe J, Eilers PHC (2015). Efficient estimation of smooth
distributions from coarsely grouped data. *American Journal of
Epidemiology* 182(2): 138-147.

Schmid VJ, Held L (2007). Bayesian age-period-cohort modeling and
prediction: BAMP. *Journal of Statistical Software* 21(8): 1-15.

Villegas AM, Kaishev VK, Millossovich P (2018). StMoMo: an R package for
stochastic mortality modeling. *Journal of Statistical Software* 84(3):
1-38.

Wood SN (2017). *Generalized Additive Models: An Introduction with R*
(2nd ed.). Chapman and Hall/CRC.

Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast
reconciliation for hierarchical and grouped time series through trace
minimization. *Journal of the American Statistical Association*
114(526): 804-819.
