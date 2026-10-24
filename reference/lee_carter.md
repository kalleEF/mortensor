# Lee-Carter model

`log m[x, t] = a[x] + b[x] k[t]` with the constraints `sum(b) = 1` and
`sum(k) = 0`, fitted per stratum:

## Usage

``` r
lee_carter(
  method = c("poisson", "svd"),
  jump_off = c("fitted", "observed", "filtered", "constrained"),
  bootstrap = FALSE,
  n_boot = 100,
  adjust = c("dt", "dxt", "e0", "none"),
  choose_period = FALSE,
  min_years = 20,
  rotate = FALSE,
  rotation = c(e0_start = 80, e0_end = 102, power = 0.5),
  family = c("poisson", "nb"),
  dispersion = c("global", "age"),
  jump_off_years = 1
)
```

## Arguments

- method:

  `"poisson"` (maximum likelihood with the distribution of `family`) or
  `"svd"`.

- jump_off:

  `"fitted"` (forecasts start from the fitted rates of the last year),
  `"observed"` (from the observed rates of the last year; Lee and Miller
  2001; with `jump_off_years > 1` from their mean over the last years,
  Stoeldraijer et al. 2018), `"filtered"` (from the rates at the Kalman
  filter estimate of the index in the last year; De Jong and
  Tickle 2006) or `"constrained"` (the maximum likelihood fit is
  constrained to the observed rates of the last year, the forecast
  origin correction of Li, Hardy and Tan 2009; `method = "poisson"`
  only). See "Jump-off" above.

- bootstrap:

  Semi-parametric bootstrap of the parameters.

- n_boot:

  Number of bootstrap fits.

- adjust:

  Re-estimation of `k[t]` for `method = "svd"`: `"dt"`, `"dxt"`, `"e0"`
  or `"none"`.

- choose_period:

  Choose the fitting period (Booth, Maindonald and Smith 2002).

- min_years:

  `choose_period` considers fitting periods of more than `min_years`
  years (`minperiod` of
  [`demography::lca()`](https://pkg.robjhyndman.com/demography/reference/lca.html)).

- rotate:

  Rotate the age pattern of decline `b` towards an ultimate pattern as
  life expectancy at birth rises (Li, Lee and Gerland 2013).

- rotation:

  Settings of the rotation, a named vector: life expectancy at birth at
  which the rotation starts (`e0_start`, default 80) and is complete
  (`e0_end`, default 102), and the exponent of the weights (`power`,
  default 0.5). The defaults are those of Li, Lee and Gerland (2013,
  eq. 2) and MortCast; elements left out keep their defaults.

- family:

  Distribution of the deaths for `method = "poisson"`: `"poisson"`
  (Brouhns, Denuit and Vermunt 2002) or `"nb"` (negative binomial;
  Delwarde, Denuit and Partrat 2007).

- dispersion:

  Dispersion of the negative binomial deaths (`family = "nb"`):
  `"global"` (one for all ages; Delwarde, Denuit and Partrat 2007) or
  `"age"` (one per age group of the data; Li, Hardy and Tan 2009).

- jump_off_years:

  For `jump_off = "observed"`: the number of last observed years whose
  mean log rates give the jump-off rates (default 1, the last year; 2 to
  5 in Stoeldraijer et al. 2018).

## Value

A model specification.

## Details

- `method = "poisson"`: Poisson log-bilinear model (Brouhns, Denuit and
  Vermunt 2002), maximum likelihood by Newton updates;

- `method = "svd"`: original method (Lee and Carter 1992): singular
  value decomposition of centred log rates, with `k[t]` re-estimated in
  a second stage (`adjust`, as in
  [`demography::lca()`](https://pkg.robjhyndman.com/demography/reference/lca.html)):
  so that the fitted deaths match the observed total deaths of each year
  (`"dt"`, Lee and Carter 1992), by Poisson maximum likelihood of the
  deaths by age given `a` and `b` (`"dxt"`, Booth, Maindonald and Smith
  2002), so that the fitted life expectancy at the first age matches the
  observed one (`"e0"`, Lee and Miller 2001), or not at all (`"none"`).

`choose_period = TRUE` fits the model to the last years only, starting
in the year that minimises the ratio of the mean deviance of a fit with
a linear `k[t]` to the mean deviance of the Lee-Carter fit (Booth,
Maindonald and Smith 2002): the longest recent period over which the
decline was close to linear
(`demography::lca(chooseperiod = TRUE, breakmethod = "bms")` gives the
same period). Fitted rates of earlier years are missing.

Compared with
[`demography::lca()`](https://pkg.robjhyndman.com/demography/reference/lca.html),
`k[t]` is centred again after the second stage (a constant shift of `k`
and `a`, with identical fitted rates), and life tables use `a_x = n/2`
at age 0 (as [`life_table()`](life_table.md)), so `adjust = "e0"`
matches a slightly different life expectancy than demography's life
table where age 0 is in the data.

The period index follows a random walk with drift. Draws include the
innovations and the uncertainty of the drift. With `bootstrap = TRUE`
the parameters are re-estimated on `n_boot` Poisson resamples of the
fitted deaths of the fitting period (semi-parametric bootstrap with
fitted deaths of Renshaw and Haberman 2008), which adds the uncertainty
of `a` and `b`. The bootstrap of Brouhns, Denuit and Van Keilegom (2005)
resamples the observed deaths instead:
`bootstrapped(lee_carter(), type = "data")`
([`bootstrapped()`](bootstrapped.md)).

Age groups without deaths in the fitting period get a constant rate (0.5
deaths over the period) and no trend; with `adjust = "e0"` they enter
the life tables at that rate.

**Safeguards of this package for sparse counts.** Where the Poisson
likelihood has no finite maximum (`method = "poisson"` with very few
deaths), the fit falls back to two package-specific estimators
(reference unknown), each with a warning (class `mt_warning_model_fit`;
within [`mt_fit()`](mt_fit.md) it is reported with the stratum): first
to the Poisson fit with a ridge penalty that shrinks `b` towards equal
shares `1 / n_ages`, and if that also fails to the SVD estimates
(`method = "svd"` with `adjust = "dt"`). Their forecasts are not those
of the Poisson Lee-Carter model; wider age groups
([`mt_regroup()`](mt_regroup.md)) or pooling ([`li_lee()`](li_lee.md))
avoid them.

Draws describe the expected rates. Predictive draws
(`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
given the rates and the exposures, the observation model of the Poisson
Lee-Carter model (Brouhns, Denuit and Vermunt 2002); no other
observation noise is added. If intervals are too narrow on real data,
that is a property of the model: see
[`recalibrated()`](recalibrated.md), [`ensemble()`](ensemble.md), or
overdispersed models such as `family = "nb"` below or
`stanmomo("lc", family = "nb")` ([`stanmomo()`](stanmomo.md)).

**Jump-off.** The jump-off rates are the rates of the last year `T` of
the fitting period, from which the forecast starts. Stoeldraijer et al.
(2018) compare six options for the Lee-Carter model; their names are in
quotes:

- `jump_off = "fitted"`: "model values", `a[x] + b[x] k[T]` (Lee and
  Carter 1992);

- `jump_off = "observed"`: "last observed values", the observed log
  rates of year `T` (Lee and Miller 2001): `a[x]` is replaced by them
  and `k` is 0 in year `T`;

- `jump_off = "observed", jump_off_years = n` (`n` = 2 to 5 in the
  paper): "average n years observed". `a[x]` is the mean of the observed
  log rates of the last `n` years and `k` is 0 at the midpoint of these
  years; the random walk with drift carries the mean from the midpoint
  to year `T`, so the jump-off rates are that mean plus
  `b[x] d (n - 1) / 2`, with `d` the drift of `k`. The paper does not
  say how the index at the midpoint relates to the fitted index of year
  `T`; this reading (the forecast does not depend on `k[T]`) reproduces
  its results for French men (below), where the mean of the fitted index
  over the `n` years does not.

Stoeldraijer et al. call a forecast accurate if its out-of-sample errors
in historical data are small (mean absolute error of life expectancy at
65 one and five years ahead), and robust if it changes little when one
more year of data is added (standard deviation of the changes of the
forecast life expectancy at 65 in 2020 between successive fitting
periods). In eight Western European countries the last observed values
were the most accurate for most countries and an average of observed
years the most robust; more years averaged gave more robust and less
accurate forecasts. For French men (Human Mortality Database in
[`StanMoMo::FRMaleData`](https://rdrr.io/pkg/StanMoMo/man/FRMaleData.html),
a later release than the paper's; ages 0-94 and 95+; fitting periods
from 1960 to 2005, ..., 2014; `method = "svd"`), the mean absolute
errors and the standard deviations of the last observed values and of
the averages are within 0.03 years of Table 8a of the paper. Those of
the model values differ by up to 0.08 years: the mean error of the
fitted log rates at ages 65 and over in 2005-2014 is -0.002 with these
data and 0.005 in the paper (Table 9a), and a shift of the log rates at
these ages by 0.007 moves life expectancy at 65 by 0.05 years.
[`anchored()`](anchored.md) moves the jump-off of any model.

**Filtered jump-off** (`jump_off = "filtered"`; De Jong and Tickle
2006). The paper writes the Lee-Carter model in state-space form (eq.
6): `y[t] = a + b alpha[t] + eps[t]` for the vector of log rates of year
`t`, and `alpha[t + 1] = delta + alpha[t] + lambda eta[t]`, with
independent normal `eps[t]` (variance `sigma^2` at every age) and
`eta[t]` (variance `sigma^2`), `b'b = 1` and the signal-to-noise ratio
`lambda`. The forecast starts from
`a + b E(alpha[T] | y[1], ..., y[T])`, the Kalman filter estimate of the
index in the jump-off year (eq. 15). The paper estimates `b` and
`lambda` by maximum likelihood with the Kalman filter, then `a`, `delta`
and `sigma` given them (Section 5), for a model with B-spline age
effects (eq. 2; eq. 6 is its case without splines, the model fitted
here). Here `a` and `b` are those of the Lee-Carter fit (`method`), in
the scale `sum(b) = 1` of `k`; given them, `sigma`, `lambda`, `delta`
and the level of the index in the first year are estimated by maximum
likelihood of the paper's model for the observed log rates. The filtered
index is computed in closed form (generalised least squares for the
level and the drift and the best linear prediction of the index), which
equals the Kalman filter with a diffuse initial state. Without state
noise (`lambda = 0`, which is the estimate in some short or noisy
series) it is the value in year `T` of the linear trend fitted to the
per-year indices `z[t] = b'(y[t] - a) / b'b`; with little observation
noise (large `lambda`) it approaches `z[T]`, the projection of the
observed log rates of year `T` on the age pattern, which is `k[T]` for
`method = "svd", adjust = "none"`. Only the start changes: the forecast
keeps the drift and the innovations of `k` (the paper's forecast uses
the maximum likelihood drift and adds the variance of the filtered
index). Package-specific (reference unknown): cells without deaths enter
with 0.5 deaths, as for `jump_off = "observed"`; age groups without
deaths in the fitting period are left out; with `bootstrap = TRUE` each
bootstrap fit filters the observed rates with its own `a` and `b`.
[`mt_parameters()`](mt_parameters.md) reports the filtered index
(`k_filtered`), its standard error (`se_k_filtered`, which includes the
uncertainty of the level and the drift) and the estimates
`filter_sigma`, `filter_lambda` and `filter_drift` (per year, in the
units of `k`). The paper's data (Australian women, 1921-2000) are not in
an R package; tests check the limiting cases, a Kalman filter and
recovery on simulated data.

**Negative binomial version** (`family = "nb"`; Delwarde, Denuit and
Partrat 2007). Following the paper:

- the deaths are negative binomial with mean
  `mu = E exp(a[x] + b[x] k[t])` and variance `mu + phi mu^2`, a Poisson
  distribution mixed over a gamma distributed random effect of each cell
  (Section 2.2). The dispersion `phi` (the paper's `k`) is one parameter
  for all ages and years;

- `a`, `b`, `k` and `phi` are estimated by maximum likelihood (eq. 5):
  uni-dimensional Newton updates of `a`, `b` and `k` (eq. 6-8), each
  given the other parameters, alternated with the estimate of `phi`
  given the fitted deaths;

- `k[t]` follows a random walk with drift (Section 3.5);

- predictive draws (`mt_forecast(predictive = TRUE)`) are negative
  binomial deaths with the estimated dispersion (size `1 / phi`), the
  observation model of Section 2.2.

Package-specific (reference unknown): the Newton updates start from the
Poisson estimates (the paper: `a = 0`, `b = 1`, `k = 0`, or random
values); `phi` is maximised by Newton steps on `log(phi)` with step
halving (the paper: Newton steps on `1 / phi`, eq. 9); the estimate is 0
(the Poisson model) where the derivative of the log-likelihood at
`phi = 0` is not positive; the log-likelihood uses gamma functions, so
deaths need not be whole numbers. The paper states no identification
constraints; `sum(b) = 1` and `sum(k) = 0` of the Poisson model
(Brouhns, Denuit and Vermunt 2002) agree with its Figures 2 and 3. As
for the Poisson model, draws include the innovations of `k` and the
uncertainty of its drift; the paper reports point forecasts only. With
`bootstrap = TRUE` the resampled deaths are negative binomial with the
estimated dispersion, and each draw has the dispersion of its bootstrap
fit (package-specific). The ridge and SVD fallbacks are not used: a fit
that does not converge fails. Like other maximum likelihood estimates of
a variance, `phi` is biased downwards when many parameters are estimated
relative to the number of cells. For France, ages 0-99, the fit to the
Human Mortality Database series in
[`demography::fr.mort`](https://pkg.robjhyndman.com/demography/reference/fr.mort.html)
reproduces the paper's dispersions (Section 3.3), AIC (Section 3.4) and
drifts and standard deviations of `k` (Sections 3.5 and 3.6) within 1%;
the paper's column labelled `sigma^2` holds the standard deviation
(divisor `T - 1`).

**Dispersion by age** (`family = "nb"` with `dispersion = "age"`; Li,
Hardy and Tan 2009). Following the paper:

- each cell consists of clusters whose deaths are Poisson given a gamma
  distributed frailty with mean 1 and a variance that depends on age.
  The deaths of the cell are then negative binomial with mean
  `mu = E exp(a[x] + b[x] k[t])` and variance `mu + alpha[x] mu^2`
  (Section 3, eq. 5-12), with one dispersion `alpha[x]` per age group of
  the data; `alpha[x]` near 0 gives the Poisson model;

- `a`, `b`, `k` and `alpha` are estimated by maximum likelihood (eq.
  13): the Newton updates of `a`, `b` and `k` of Appendix 1, alternated
  with the estimate of each `alpha[x]` given the fitted deaths of its
  age group;

- `jump_off = "constrained"` applies the forecast origin correction
  (Section 3 and Appendix 2): `a[x]` is the observed log rate of the
  last year `T` (the paper writes `a[x] = m[x, T]`; `a` is on the log
  scale of eq. 1), `k[T] = 0`, and `b` and `k` are estimated by maximum
  likelihood under these constraints. The fitted rates of the last year
  are the observed ones, and the log-likelihood is lower than without
  the constraints (footnote 3). With `family = "poisson"` it gives the
  original model with the correction, which the paper fits for
  comparison (Section 4.2);

- interval forecasts come from the parametric bootstrap of Section 3:
  `bootstrap = TRUE` simulates negative binomial deaths from the fitted
  model, re-estimates `a`, `b` and `k` (with `jump_off = "constrained"`,
  `a` comes from the simulated deaths of the last year, as the
  correction is part of the estimation in Section 3), and simulates `k`
  from the time series model of each refit. The dispersions are those of
  the fit: step 2(a) re-estimates `a`, `b` and `k`. The paper's method
  is
  `lee_carter(family = "nb", dispersion = "age", jump_off = "constrained", bootstrap = TRUE)`;

- predictive draws (`mt_forecast(predictive = TRUE)`) are negative
  binomial deaths with the dispersion of their age group (size
  `1 / alpha[x]`), the observation model of eq. (10). The paper's
  intervals describe the central death rates, which are the draws
  without `predictive`.

Package-specific (reference unknown): as for one dispersion, the updates
start from the Poisson estimates (Appendix 1 recommends them), and each
`alpha[x]` is maximised by Newton steps on `log(alpha[x])` with step
halving (the paper: Newton steps on `alpha[x]`). Its estimate is 0 where
the derivative of the log-likelihood at 0 is not positive. The paper
fits an ARIMA model to `k` without naming the model it chose; here `k`
follows a random walk with drift, and the draws of each bootstrap fit
include the uncertainty of its drift, as all draws of `lee_carter()`.
With `jump_off = "constrained"`, an age group without deaths in the last
year starts from the rate of 0.5 deaths. Age groups without deaths in
all years keep their constant rate and have Poisson deaths (dispersion
0). As with one dispersion, a fit that does not converge fails;
`jump_off = "constrained"` has no fallbacks either. Each `alpha[x]` is
estimated from cells that also estimate `a[x]`, `b[x]` and `k`. Like
other maximum likelihood estimates of a variance it is biased downwards:
by about `2 / T` with `T` years, and more at the ages that carry most of
the information on `k` (many deaths and a small dispersion). The paper
uses data of the Human Mortality Database for the United States and
Canada, ages 0-99, 1950-2004 (Section 4.1). The packages demography and
StMoMo do not hold these series, so its estimates are not reproduced.

`rotate = TRUE` rotates the age pattern of decline as life expectancy
rises (Li, Lee and Gerland 2013; status: published, with the deviations
listed below). In long-lived populations the decline of mortality slows
down at young ages and speeds up at old ages, which a fixed `b[x]` does
not capture over long horizons. For each path and year `t`, the life
expectancy at birth `e0` of the unrotated forecast gives the weight
`w = (0.5 (1 + sin(pi / 2 (2 u - 1))))^power` with
`u = (e0 - e0_start) / (e0_end - e0_start)`; `w` is 0 below `e0_start`
and 1 from `e0_end` (eq. 2 of the paper). The age pattern of that year
is `B[x, t] = (1 - w) b[x] + w bu[x]` (eq. 3). The defaults are those of
the paper: the rotation starts at 80 (about the median life expectancy
of the 20 populations with the highest life expectancy in 2005-2010), is
complete at 102 (the value with which the ratio of the death rates at
age 0 and at ages 15-19 approaches its estimated long-run level of 7.7),
and the exponent is 0.5. The ultimate pattern `bu` has the mean of `b`
over ages 15-64 at all ages below 65, as in the paper. From age 65 it
keeps the shape of `b`, scaled so that the first age group from 65
equals that mean, and it sums to 1 like `b`: the paper only notes that
lowering `b` at young ages raises it at older ages because `b` sums to
1, and this scaling is that of MortCast (Sevcikova et al. 2016). The
index is then solved again, path by path and year by year, so that the
rotated rates have the same life expectancy at birth as the unrotated
ones (eq. 4): rotation changes the age pattern of the forecast and keeps
its life expectancy at birth. Decline thus moves from young ages to the
ages where `bu` exceeds `b`. Where `b` at 65 is already above its mean
over ages 15-64 (as for French women), it moves mostly to ages 15-64.
Fitted rates, and paths with `e0` below `e0_start`, are unchanged.
Rotation needs ages from 0 (life expectancy at birth) with age groups
starting at 15 and 65. With a closed last age group, a warning notes
that life expectancy is then truncated at the highest age. Rotated
forecasts compute life tables for every path and year and are much
slower (about 20 times for single years of age).

Deviations of the rotation from Li, Lee and Gerland (2013):

- the paper forecasts `b` and `k` for both sexes combined (Li and
  Lee 2005) and takes the weights from the life expectancy of both sexes
  combined; `lee_carter()` fits each stratum separately, so the weights
  come from the life expectancy of the stratum (MortCast takes the path
  of life expectancy from a separate model);

- the paper illustrates the method with the median path of life
  expectancy; here every draw is rotated with its own path of life
  expectancy;

- rates change from the jump-off year `T`,
  `log m[x, t] = log m[x, T] + B[x, t] (k*[t] - k[T])`. Eq. (4) writes
  `a[x] + B[x, t] K[t]` with `a` the last observed log rates (as
  `jump_off = "observed"`) or their mean over the years; with
  `jump_off = "fitted"` (or `"filtered"`) the fitted (or filtered) rates
  of year `T` are the starting point (package-specific), which keeps the
  jump-off rates when `B` departs from `b`;

- the paper extends the open age group to age 130 with a logistic model
  before fitting (footnote 3); here the age groups of the data are used,
  and life expectancy closes the open interval with `L = l / m`
  ([`life_table()`](life_table.md)).

## References

Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
*Journal of the American Statistical Association* 87(419): 659-671.

Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear
regression approach to the construction of projected lifetables.
*Insurance: Mathematics and Economics* 31(3): 373-393.
[doi:10.1016/S0167-6687(02)00185-3](https://doi.org/10.1016/S0167-6687%2802%2900185-3)

Booth H, Maindonald J, Smith L (2002). Applying Lee-Carter under
conditions of variable mortality decline. *Population Studies* 56(3):
325-336.
[doi:10.1080/00324720215935](https://doi.org/10.1080/00324720215935)

Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
method for forecasting mortality. *Demography* 38(4): 537-549.
[doi:10.1353/dem.2001.0036](https://doi.org/10.1353/dem.2001.0036)

Stoeldraijer L, van Duin C, van Wissen L, Janssen F (2018). Comparing
strategies for matching mortality forecasts to the most recently
observed data: exploring the trade-off between accuracy and robustness.
*Genus* 74: 16.
[doi:10.1186/s41118-018-0041-y](https://doi.org/10.1186/s41118-018-0041-y)

De Jong P, Tickle L (2006). Extending Lee-Carter mortality forecasting.
*Mathematical Population Studies* 13(1): 1-18.
[doi:10.1080/08898480500452109](https://doi.org/10.1080/08898480500452109)

Delwarde A, Denuit M, Partrat C (2007). Negative binomial version of the
Lee-Carter model for mortality forecasting. *Applied Stochastic Models
in Business and Industry* 23(5): 385-401.
[doi:10.1002/asmb.679](https://doi.org/10.1002/asmb.679)

Li JS-H, Hardy MR, Tan KS (2009). Uncertainty in mortality forecasting:
an extension to the classical Lee-Carter approach. *ASTIN Bulletin*
39(1): 137-164.
[doi:10.2143/AST.39.1.2038060](https://doi.org/10.2143/AST.39.1.2038060)

Brouhns N, Denuit M, Van Keilegom I (2005). Bootstrapping the Poisson
log-bilinear model for mortality forecasting. *Scandinavian Actuarial
Journal* 2005(3): 212-224.

Renshaw AE, Haberman S (2008). On simulation-based approaches to risk
measurement in mortality with specific reference to Poisson Lee-Carter
modelling. *Insurance: Mathematics and Economics* 42(2): 797-816.

Li N, Lee R, Gerland P (2013). Extending the Lee-Carter method to model
the rotation of age patterns of mortality decline for long-term
projections. *Demography* 50(6): 2037-2051.
[doi:10.1007/s13524-013-0232-2](https://doi.org/10.1007/s13524-013-0232-2)

Li N, Lee R (2005). Coherent mortality forecasts for a group of
populations: an extension of the Lee-Carter method. *Demography* 42(3):
575-594.
[doi:10.1353/dem.2005.0021](https://doi.org/10.1353/dem.2005.0021)

Sevcikova H, Li N, Kantorova V, Gerland P, Raftery AE (2016).
Age-specific mortality and fertility rates for probabilistic population
projections. In Schoen R (ed.) *Dynamic Demographic Analysis*. Springer,
285-310.
[doi:10.1007/978-3-319-26603-9_15](https://doi.org/10.1007/978-3-319-26603-9_15)

## Examples

``` r
x <- mt_example(deaths = "single")
fit <- mt_fit(lee_carter(), x)
fc <- mt_forecast(fit, h = 10, n_draws = 200, seed = 1)
# \donttest{
# rotation needs ages from 0: French women (HMD, via demography)
if (requireNamespace("demography", quietly = TRUE)) {
  fr <- as_mt_data(demography::fr.mort, series = "female", years = 1970:2006)
  fr <- mt_regroup(fr, age = c(0, 1, seq(5, 100, 5)))
  fc_rot <- mt_forecast(mt_fit(lee_carter(rotate = TRUE), fr), h = 50, n_draws = 100, seed = 1)
  life_expectancy(fc_rot)
}
#>     stratum    sex  label period     mean   q0.025     q0.1     q0.5     q0.9
#>       <int> <char> <char>  <num>    <num>    <num>    <num>    <num>    <num>
#>  1:       1 female female   2007 84.45832 83.97909 84.14176 84.42570 84.84040
#>  2:       1 female female   2008 84.64582 83.98058 84.14838 84.59950 85.21669
#>  3:       1 female female   2009 84.81986 84.07303 84.20940 84.77582 85.50512
#>  4:       1 female female   2010 84.95914 84.04569 84.24819 84.94668 85.64804
#>  5:       1 female female   2011 85.17228 84.19406 84.39411 85.14099 86.06945
#>  6:       1 female female   2012 85.36369 84.29550 84.57521 85.32251 86.25054
#>  7:       1 female female   2013 85.52862 84.25485 84.64479 85.50286 86.57259
#>  8:       1 female female   2014 85.71318 84.34690 84.67043 85.68552 86.81340
#>  9:       1 female female   2015 85.88541 84.46973 84.64564 85.91419 87.05366
#> 10:       1 female female   2016 86.04160 84.41473 84.82095 86.08765 87.20145
#> 11:       1 female female   2017 86.21922 84.48068 84.94156 86.28038 87.39625
#> 12:       1 female female   2018 86.35767 84.40915 85.00913 86.40913 87.58079
#> 13:       1 female female   2019 86.55185 84.48270 85.21895 86.64661 87.86014
#> 14:       1 female female   2020 86.73739 84.59386 85.34588 86.88488 88.21593
#> 15:       1 female female   2021 86.86154 84.42512 85.50065 86.93892 88.28391
#> 16:       1 female female   2022 87.02331 84.72102 85.58790 87.17942 88.46582
#> 17:       1 female female   2023 87.17882 84.52893 85.72703 87.32831 88.81129
#> 18:       1 female female   2024 87.34424 84.46988 85.89195 87.45908 88.97988
#> 19:       1 female female   2025 87.50461 84.85015 86.03837 87.58242 89.08385
#> 20:       1 female female   2026 87.64083 84.90858 86.10474 87.72627 89.32925
#> 21:       1 female female   2027 87.80483 85.23695 86.25216 87.89847 89.50748
#> 22:       1 female female   2028 87.99057 85.60346 86.35903 87.97488 89.43914
#> 23:       1 female female   2029 88.14179 85.67133 86.42451 88.09798 89.69419
#> 24:       1 female female   2030 88.27834 85.71681 86.63274 88.27427 89.95481
#> 25:       1 female female   2031 88.44130 85.77055 86.79347 88.51256 90.18044
#> 26:       1 female female   2032 88.58123 85.83978 86.80131 88.55663 90.17031
#> 27:       1 female female   2033 88.75500 85.89025 86.91354 88.85676 90.58709
#> 28:       1 female female   2034 88.92816 85.86039 87.07892 88.94462 90.73743
#> 29:       1 female female   2035 89.06049 85.90906 87.27126 89.17776 90.87124
#> 30:       1 female female   2036 89.21073 86.11898 87.51358 89.36751 91.02247
#> 31:       1 female female   2037 89.35301 86.21238 87.67122 89.44646 91.30539
#> 32:       1 female female   2038 89.52548 86.20872 87.73986 89.49441 91.54192
#> 33:       1 female female   2039 89.70539 86.27487 87.96464 89.78436 91.90117
#> 34:       1 female female   2040 89.84683 86.38395 88.11563 89.88735 92.10813
#> 35:       1 female female   2041 89.99395 86.85480 88.16030 89.96020 92.17889
#> 36:       1 female female   2042 90.15294 87.11016 88.25911 90.18968 92.35955
#> 37:       1 female female   2043 90.30115 87.20544 88.46409 90.35192 92.40212
#> 38:       1 female female   2044 90.42379 87.30900 88.34859 90.37858 92.47518
#> 39:       1 female female   2045 90.56684 87.39265 88.53738 90.59868 92.69714
#> 40:       1 female female   2046 90.74251 87.68452 88.77240 90.80616 92.93426
#> 41:       1 female female   2047 90.92607 87.92870 88.88210 90.97651 93.05826
#> 42:       1 female female   2048 91.03234 87.83601 88.92326 91.06034 93.36926
#> 43:       1 female female   2049 91.16372 87.79822 89.18376 91.10862 93.48960
#> 44:       1 female female   2050 91.32277 88.06035 89.31103 91.36373 93.77092
#> 45:       1 female female   2051 91.46829 88.03537 89.38917 91.48872 93.95286
#> 46:       1 female female   2052 91.58919 87.91411 89.43010 91.60156 94.11008
#> 47:       1 female female   2053 91.71784 88.32247 89.39964 91.74799 94.15596
#> 48:       1 female female   2054 91.84898 88.25525 89.64994 91.90272 94.28004
#> 49:       1 female female   2055 91.95132 88.26744 89.70716 91.92643 94.40703
#> 50:       1 female female   2056 92.04926 88.41287 89.67235 92.07860 94.61926
#>     stratum    sex  label period     mean   q0.025     q0.1     q0.5     q0.9
#>       <int> <char> <char>  <num>    <num>    <num>    <num>    <num>    <num>
#>       q0.975    at
#>        <num> <num>
#>  1: 84.97797     0
#>  2: 85.37398     0
#>  3: 85.71163     0
#>  4: 85.86549     0
#>  5: 86.23777     0
#>  6: 86.46993     0
#>  7: 86.85593     0
#>  8: 87.11083     0
#>  9: 87.46653     0
#> 10: 87.94733     0
#> 11: 88.06047     0
#> 12: 88.30347     0
#> 13: 88.54289     0
#> 14: 88.87291     0
#> 15: 89.02873     0
#> 16: 89.12558     0
#> 17: 89.42112     0
#> 18: 89.53384     0
#> 19: 89.76067     0
#> 20: 89.93365     0
#> 21: 90.07433     0
#> 22: 90.20059     0
#> 23: 90.47060     0
#> 24: 90.38945     0
#> 25: 90.48946     0
#> 26: 90.92321     0
#> 27: 91.15902     0
#> 28: 91.50720     0
#> 29: 91.64974     0
#> 30: 91.85898     0
#> 31: 92.02428     0
#> 32: 92.12758     0
#> 33: 92.53931     0
#> 34: 92.53557     0
#> 35: 92.65133     0
#> 36: 92.83783     0
#> 37: 92.97059     0
#> 38: 93.15846     0
#> 39: 93.36540     0
#> 40: 93.44503     0
#> 41: 93.76265     0
#> 42: 94.03238     0
#> 43: 94.21914     0
#> 44: 94.28079     0
#> 45: 94.42662     0
#> 46: 94.70687     0
#> 47: 94.90896     0
#> 48: 95.02335     0
#> 49: 95.14042     0
#> 50: 95.19589     0
#>       q0.975    at
#>        <num> <num>
# }
# negative binomial deaths (Delwarde, Denuit and Partrat 2007)
fit_nb <- mt_fit(lee_carter(family = "nb"), x)
mt_parameters(fit_nb)[parameter == "nb_dispersion"]
#>    stratum     parameter  index        value
#>     <char>        <char> <char>        <num>
#> 1:  female nb_dispersion        2.913877e-05
#> 2:    male nb_dispersion        0.000000e+00
# dispersion by age with the forecast origin correction and the parametric
# bootstrap (Li, Hardy and Tan 2009)
fit_lht <- mt_fit(lee_carter(family = "nb", dispersion = "age", jump_off = "constrained",
                             bootstrap = TRUE, n_boot = 20), mt_subset(x, sex == "female"))
mt_parameters(fit_lht)[parameter == "nb_dispersion"]
#>     stratum     parameter  index        value
#>      <char>        <char> <char>        <num>
#>  1:  female nb_dispersion     30 1.332502e-02
#>  2:  female nb_dispersion     31 4.772923e-03
#>  3:  female nb_dispersion     32 4.053318e-03
#>  4:  female nb_dispersion     33 1.236417e-03
#>  5:  female nb_dispersion     34 0.000000e+00
#>  6:  female nb_dispersion     35 0.000000e+00
#>  7:  female nb_dispersion     36 0.000000e+00
#>  8:  female nb_dispersion     37 8.228632e-04
#>  9:  female nb_dispersion     38 4.278345e-03
#> 10:  female nb_dispersion     39 0.000000e+00
#> 11:  female nb_dispersion     40 0.000000e+00
#> 12:  female nb_dispersion     41 0.000000e+00
#> 13:  female nb_dispersion     42 2.691890e-03
#> 14:  female nb_dispersion     43 0.000000e+00
#> 15:  female nb_dispersion     44 3.367764e-05
#> 16:  female nb_dispersion     45 6.100974e-03
#> 17:  female nb_dispersion     46 1.888541e-03
#> 18:  female nb_dispersion     47 1.450254e-03
#> 19:  female nb_dispersion     48 0.000000e+00
#> 20:  female nb_dispersion     49 0.000000e+00
#> 21:  female nb_dispersion     50 0.000000e+00
#> 22:  female nb_dispersion     51 2.742291e-03
#> 23:  female nb_dispersion     52 3.108202e-04
#> 24:  female nb_dispersion     53 0.000000e+00
#> 25:  female nb_dispersion     54 1.143663e-03
#> 26:  female nb_dispersion     55 5.628666e-03
#> 27:  female nb_dispersion     56 0.000000e+00
#> 28:  female nb_dispersion     57 8.045295e-04
#> 29:  female nb_dispersion     58 0.000000e+00
#> 30:  female nb_dispersion     59 0.000000e+00
#> 31:  female nb_dispersion     60 0.000000e+00
#> 32:  female nb_dispersion     61 0.000000e+00
#> 33:  female nb_dispersion     62 8.316641e-05
#> 34:  female nb_dispersion     63 8.191533e-04
#> 35:  female nb_dispersion     64 1.580029e-04
#> 36:  female nb_dispersion     65 0.000000e+00
#> 37:  female nb_dispersion     66 3.316936e-04
#> 38:  female nb_dispersion     67 5.760218e-04
#> 39:  female nb_dispersion     68 0.000000e+00
#> 40:  female nb_dispersion     69 0.000000e+00
#> 41:  female nb_dispersion     70 1.959614e-04
#> 42:  female nb_dispersion     71 0.000000e+00
#> 43:  female nb_dispersion     72 0.000000e+00
#> 44:  female nb_dispersion     73 0.000000e+00
#> 45:  female nb_dispersion     74 3.327211e-04
#> 46:  female nb_dispersion     75 0.000000e+00
#> 47:  female nb_dispersion     76 2.120416e-04
#> 48:  female nb_dispersion     77 0.000000e+00
#> 49:  female nb_dispersion     78 0.000000e+00
#> 50:  female nb_dispersion     79 0.000000e+00
#> 51:  female nb_dispersion     80 0.000000e+00
#> 52:  female nb_dispersion     81 4.490368e-04
#> 53:  female nb_dispersion     82 1.416356e-04
#> 54:  female nb_dispersion     83 0.000000e+00
#> 55:  female nb_dispersion     84 3.430378e-05
#> 56:  female nb_dispersion     85 9.622345e-05
#> 57:  female nb_dispersion     86 0.000000e+00
#> 58:  female nb_dispersion     87 0.000000e+00
#> 59:  female nb_dispersion     88 0.000000e+00
#> 60:  female nb_dispersion     89 0.000000e+00
#> 61:  female nb_dispersion     90 0.000000e+00
#> 62:  female nb_dispersion     91 3.060295e-05
#> 63:  female nb_dispersion     92 1.258580e-05
#> 64:  female nb_dispersion     93 0.000000e+00
#> 65:  female nb_dispersion     94 2.127695e-05
#> 66:  female nb_dispersion     95 0.000000e+00
#> 67:  female nb_dispersion     96 7.320695e-05
#> 68:  female nb_dispersion     97 0.000000e+00
#> 69:  female nb_dispersion     98 2.175004e-04
#> 70:  female nb_dispersion     99 1.609601e-04
#> 71:  female nb_dispersion   100+ 6.562962e-04
#>     stratum     parameter  index        value
#>      <char>        <char> <char>        <num>
# jump-off: mean of the last three observed years (Stoeldraijer et al.
# 2018) and the Kalman filter estimate (De Jong and Tickle 2006)
fit_av <- mt_fit(lee_carter(jump_off = "observed", jump_off_years = 3), x)
fit_kf <- mt_fit(lee_carter(jump_off = "filtered"), x)
mt_parameters(fit_kf)[grepl("filter", parameter)]
#>     stratum     parameter  index        value
#>      <char>        <char> <char>        <num>
#>  1:  female    k_filtered   2019 -13.31314483
#>  2:  female se_k_filtered   2019   0.39234808
#>  3:  female  filter_sigma          0.06480888
#>  4:  female filter_lambda          0.81877306
#>  5:  female  filter_drift         -0.87267723
#>  6:    male    k_filtered   2019 -12.32094921
#>  7:    male se_k_filtered   2019   0.38782858
#>  8:    male  filter_sigma          0.05381630
#>  9:    male filter_lambda          1.79159067
#> 10:    male  filter_drift         -0.94477645
```
