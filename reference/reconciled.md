# Models with reconciled forecasts

Wraps a model specification so that its forecasts are reconciled over a
hierarchy every time it is forecast, also inside
[`mt_backtest()`](mt_backtest.md). If the training data contain only the
bottom strata, the aggregate strata are added with
[`mt_aggregate_strata()`](mt_aggregate_strata.md) before fitting, so
that the forecasts cover every node; exposures for the forecast
(`exposure` of [`mt_forecast()`](mt_forecast.md)) can then be given for
the bottom strata only, and those of the aggregates are their sums.
Predictive draws (used by backtests) are coherent as well: given the
reconciled rates, the Poisson variation of the deaths of the bottom
strata is drawn independently per cell and draw (none for models of
observed rates such as [`rwd()`](rwd.md); for
[`ensemble()`](ensemble.md) and [`stratified()`](stratified.md) models
as for the component or group each draw comes from; negative binomial
for `stanmomo(family = "nb")`), and the aggregates are the sums of the
bottom strata.

## Usage

``` r
reconciled(
  spec,
  hierarchy,
  method = c("wls", "mint", "structural", "ols", "bottom_up", "proportional",
    "energy_score"),
  shrink = "auto",
  dependence = c("none", "ecc", "schaake", "empirical", "auto"),
  fixed = NULL,
  years = 10,
  n_origins = 5,
  horizon = 5,
  control = list(),
  seed = 1,
  template = NULL,
  covariance = c("draws", "errors")
)
```

## Arguments

- spec:

  Base model specification (fitted to all strata, or to the bottom
  strata with `method = "bottom_up"`, see Details).

- hierarchy:

  An [`mt_hierarchy()`](mt_hierarchy.md). When the data hold the
  aggregate strata, their exposures are compared with those the
  hierarchy implies, with a warning (class
  `mt_warning_hierarchy_exposure`) when they differ, for example causes
  of death declared without `additive = "rates"`.

- method:

  A method of [`mt_reconcile()`](mt_reconcile.md), or `"energy_score"`
  (score-optimal reconciliation, see Details).

- shrink:

  See [`mt_reconcile()`](mt_reconcile.md); with `covariance = "errors"`
  the shrinkage intensity for the errors.

- dependence:

  `"none"` (default: the draws are reconciled as they are), `"ecc"`,
  `"schaake"`, `"empirical"` or `"auto"` (the same as `"none"`, kept for
  compatibility). With `"ecc"` the draws of the base model take the
  dependence of the draws of `template` before reconciliation (ensemble
  copula coupling, Schefzik, Thorarinsdottir and Gneiting 2013): the
  same as `reconciled(coupled(spec, template), ...)`, see
  [`coupled()`](coupled.md). With `"schaake"` the draws of every age,
  forecast year and stratum take the ranks of the observed changes of
  distinct training years (the Schaake shuffle of Clark et al. 2004),
  which also replaces the dependence between the ages and years of each
  stratum. With `"empirical"` the whole paths of each stratum take the
  ranks of the observed annual changes of its log age-standardised rate
  in distinct training years (the Schaake shuffle with the stratum as
  margin): every stratum keeps its paths, and only their pairing across
  strata changes. [`mt_reconcile()`](mt_reconcile.md) describes both.
  Their templates come from the training data of each fit, also inside
  [`mt_backtest()`](mt_backtest.md). Draws of joint models such as
  [`li_lee()`](li_lee.md) and [`tucker3()`](tucker3.md), which carry the
  dependence between strata, are never reordered by `"schaake"` or
  `"empirical"` (they then give a warning and are ignored). The
  observation noise of predictive draws is independent across cells
  either way.

- fixed:

  Strata whose forecasts are kept (see
  [`mt_reconcile()`](mt_reconcile.md)); not with `"energy_score"`.

- years, n_origins, horizon:

  Training of `method = "energy_score"` and of `covariance = "errors"`:
  the length of the training window (the last `years` observed years),
  the number of inner origins in it (as in
  [`recalibrated()`](recalibrated.md)) and the number of years each
  inner forecast covers (horizons 1 to `horizon`, cut at the last
  training year).

- control:

  Settings of the training for `method = "energy_score"` (`n_draws` also
  for `covariance = "errors"`): a list with any of `n_draws` (draws per
  inner forecast, default 100), `draws_per_step` (draws per training
  forecast in each step, 50), `max_iter` (number of steps, 500), `eta`
  (learning rate, 0.001), `beta1` (0.9), `beta2` (0.999), `tol` (1e-4)
  and `epsilon` (1e-8), as in `ProbReco::scoreopt.control()`. The
  optimisation stops early when no parameter moves by more than `tol` in
  a step, so `tol` must be below `eta` (0 runs all steps).
  `max_iter = 0` gives OLS reconciliation.

- seed:

  Random seed of the training for `method = "energy_score"` and
  `covariance = "errors"` (inner forecasts, observation noise and the
  draws of each step), so that fits are reproducible. The Schaake
  shuffles of `dependence` use the `seed` of
  [`mt_forecast()`](mt_forecast.md).

- template:

  With `dependence = "ecc"`: the model specification whose draws give
  the dependence, usually a joint model such as [`cae()`](cae.md) or
  [`li_lee()`](li_lee.md).

- covariance:

  `"draws"` (default: the weights of `"wls"` and `"mint"` from the base
  draws of each age and forecast year, see
  [`mt_reconcile()`](mt_reconcile.md)) or `"errors"` (from the errors of
  inner forecasts in the training data, see Details).

## Value

A model specification. With `method = "energy_score"` its fit holds
`$reconcile$scoreopt`: `d` `[bottom stratum, age]`, `G`
`[bottom stratum, stratum, age]` on the scale of the reconciliation, the
inner `origins`, the `control` settings and a `table` per age with the
number of training forecasts, the scale of the bottom strata, the number
of steps, and the mean in-sample energy score with OLS (`energy_ols`)
and with the estimated \\d\\ and \\G\\ (`energy_opt`). With
`covariance = "errors"` its fit holds `$reconcile$errors`: the errors
`[age, horizon, stratum, origin]`, the number of error vectors with
observed values for every stratum per age and horizon (`n`) and the
inner `origins`. With `method = "bottom_up"`, `$reconcile$bottom_only`
tells whether the base model was fitted to the bottom strata only.
Forecasts with `dependence = "schaake"` report the shuffle in
`$schaake`, those with `dependence = "empirical"` in `$empirical` (see
[`mt_reconcile()`](mt_reconcile.md)).

## Details

**Bottom-up.** With `method = "bottom_up"`, the base model is fitted to
the bottom strata only, and the aggregates are the sums of their draws
(of deaths, or of rates for variables with `additive = "rates"`).
Wrappers inside the base model ([`ensemble()`](ensemble.md) weights,
[`tuned()`](mt_tune.md), [`recalibrated()`](recalibrated.md)) then run
their inner backtests on the bottom strata. The fit describes the
aggregates by sums as well: their deaths and exposures in the training
years, and their fitted rates (fitted deaths of the bottom strata over
the exposures). Levels of [`stratified()`](stratified.md) models that
are aggregates of the hierarchy, such as
`stratified("cause", all = ...)`, have no strata in this fit and are
skipped. Joint models that pool the bottom strata and give each
aggregate parameters of its own ([`li_lee()`](li_lee.md),
[`cae()`](cae.md), [`coherent_fdm()`](coherent_fdm.md),
`coda(coherent = TRUE)`) estimate the same parameters of the bottom
strata as when fitted to all strata; their draws of the bottom strata
have the same distribution and differ by the random numbers drawn. The
other methods fit the base model to all strata, and so does bottom-up
for models whose fit of the bottom strata uses the aggregates:
[`tucker3()`](tucker3.md) and [`coda_multiway()`](coda_multiway.md)
decompose all strata together, [`saint()`](saint.md) takes its trend
from reference strata, and [`coda_causes()`](coda_causes.md) needs the
all-cause stratum. Observation noise of predictive draws is drawn before
the reordering of `dependence`, for the base draws, and moves with their
values; reordering and noise are independent, so this gives the
distribution of noise drawn after the reordering, and with the same seed
every bottom stratum keeps exactly the predictive draws of
`dependence = "none"`.

Whether coherence costs accuracy, and which method and base model work
best, depends on the data: compare them in a backtest at every level of
the hierarchy (`mt_score_table(..., hierarchy = )`), as in the vignette
"Coherent forecasts: totals and their parts".

**Score-optimal reconciliation.** `method = "energy_score"` follows
Panagiotelis et al. (2023). The reconciled draws are \\\tilde y = S (d +
G \hat y)\\, where \\\hat y\\ is a base draw of all nodes, \\S\\ sums
the bottom strata to the nodes, and \\d\\ and \\G\\ minimise the total
energy score of past forecasts. Each time the model is fitted (also at
every origin of an outer [`mt_backtest()`](mt_backtest.md)), the base
model is refitted at `n_origins` origins in the last `years` of the
training data, as in [`recalibrated()`](recalibrated.md), and forecasts
the next `horizon` years. Every origin, horizon and age with observed
values for all strata is one training forecast. \\d\\ and \\G\\ are
estimated separately for each age, from the training forecasts of all
origins and horizons, on the scale of the other methods (deaths, or
rates when the levels share exposures). The estimation is stochastic
gradient descent (Adam) started from OLS reconciliation, as in the
package ProbReco. The score is that of the predictive distribution the
forecast gives: reconciled draws plus the observation noise of the
bottom strata. In the forecast, bottom values below zero are set to zero
and the aggregates are their sums. The estimates and a table of the
training are in the fit (see Value).

As in the paper (Section 6, Algorithm 1), \\d\\ and \\G\\ are free (no
constraint such as \\G S = I\\), the objective is the total energy score
with exponent 1 estimated from pairs of base draws (their Eq. 5), and
Adam starts from OLS reconciliation with their settings (`eta = 0.001`,
`beta1 = 0.9`, `beta2 = 0.999`, `epsilon = 1e-8`, their Section 7). The
gradient is written out analytically; ProbReco differentiates the same
estimate automatically; with the same draws both follow the same
optimisation path. The departures from the paper:

- training forecasts: Algorithm 1 refits the base models on rolling
  windows of the length of the final fit at consecutive origins (250
  one-step forecasts in their simulations), and the objective (their
  Eq. 4) is that of one horizon. Here the inner models use all training
  years up to each of `n_origins` origins (expanding windows, shorter
  than the final fit), and the forecasts of horizons 1 to `horizon`
  enter one objective, so that one \\d\\ and \\G\\ serve all horizons;

- one \\d\\ and \\G\\ per age: every age is a hierarchy of its own;

- the score includes the observation noise of the bottom strata, so it
  is the score of the predictive distribution of the deaths or rates;

- `draws_per_step` is 50 by default (the paper uses 250 draws per score
  in its simulations), and the draws come from a fixed set of inner
  draws (see below);

- stopping: the paper stops when the change in all gradients is less
  than 10% of `eta`; ProbReco and mortensor stop when no parameter moves
  by more than `tol` (1e-4, 10% of the default `eta`) in a step, or
  after `max_iter` steps (500, the default of ProbReco);

- the rescaling before the optimisation and the bottom values set to
  zero in the forecast are package-specific (reference unknown). The
  paper does not treat negative values; the non-negative quadratic
  program of the projection methods (see
  [`mt_reconcile()`](mt_reconcile.md)) does not apply to the map \\S
  (d + G \hat y)\\.

Two details differ from ProbReco. (1) Each step estimates the energy
score from two samples of `draws_per_step` base draws of every training
forecast, which ProbReco simulates independently. Here they are taken
from the `n_draws` inner draws without replacement when there are at
least `2 * draws_per_step` of them, and otherwise from two disjoint
random halves, so that a draw is never paired with itself (such pairs
would contribute a distance of zero to the score). (2) Before the
optimisation, the bottom values are divided by their mean observed value
and every node, as an input, by its own mean observed value; the
estimates are transformed back afterwards. This leaves the optimum
unchanged, but Adam moves every parameter by steps of about `eta`, so
the rescaling changes the path of the optimisation, and since it stops
after `max_iter` steps or when the steps fall below `tol`, also the
estimates, relative to ProbReco, which optimises on the original scale.
The rescaling is package-specific; it makes the steps comparable across
ages and between small and large nodes.

The method has many parameters and few training forecasts. \\G\\ has one
parameter per bottom stratum and node at each age (20 x 54 for 2 sexes x
2 regions x 5 quintiles with all aggregates); the defaults give about 20
training forecasts per age, which overlap. The only regularisation is
the limited number of steps from the OLS start (`control$max_iter`,
`control$eta`). \\d\\ and \\G\\ then fit the training forecasts much
better than OLS, mostly through level corrections of single strata.
These help in new years only if the errors of the training years recur:
they do not after a change of trend, and the inner models are fitted to
shorter series than the final model. In the package's backtests on
simulated regional data and on the ONS data by sex and deprivation
quintile, the method was less accurate than `"wls"` at every level of
the hierarchy, also with fewer and smaller steps. The applications in
the paper had longer training windows. Compare the method with `"wls"`
or `"mint"` in a backtest before using it.

The training refits the base model at `n_origins` origins; the
optimisation takes about 3 seconds per age for 54 strata with the
defaults. `control$max_iter = 0` gives OLS reconciliation at the cost of
the inner forecasts.

**Weights from past errors.** `covariance = "errors"` with
`method = "wls"` or `"mint"`. MinT weights the base forecasts by the
covariance matrix of their `h`-step forecast errors (Wickramasuriya,
Turlach and Hyndman 2020, Section 2.1, after Wickramasuriya et al.
2019); Hyndman et al. (2011, Section 3) note that reconciled intervals
need this matrix. With `covariance = "errors"` it is estimated from
errors of past forecasts. Each time the model is fitted (also at every
origin of an outer [`mt_backtest()`](mt_backtest.md)), the base model is
refitted at `n_origins` inner origins in the last `years` of the
training data, as for `"energy_score"`, and forecasts the next `horizon`
years (cut at the last training year) with `control$n_draws` draws and
the observed exposures. The error of a stratum is its observed value
minus the mean of its draws, on the scale of the reconciliation (deaths,
or rates when the levels share exposures). For each age and inner
horizon `h`, `"mint"` takes the shrinkage estimator of MinT: the sample
covariance of the `n` error vectors `e`, \\\hat W = n^{-1} \sum e e'\\
(errors of mean zero), shrunk towards its diagonal, \\\lambda \\
diag(\hat W) + (1 - \lambda) \hat W\\, with \\\lambda = \sum\_{i \ne j}
\widehat{Var}(r\_{ij}) / \sum\_{i \ne j} r\_{ij}^2\\ from the
correlations \\r\_{ij}\\ of the errors (Schäfer and Strimmer 2005;
Wickramasuriya et al. 2019; Wickramasuriya et al. 2020, Section 5: the
estimator "shrinks the sample correlation matrix towards an identity
matrix"). Panagiotelis et al. (2023, Table 1) give the same form and
write \\\lambda\\ with the covariances. `"wls"` takes the diagonal, the
mean squared errors. `shrink` fixes \\\lambda\\; `shrink = 1` gives the
result of `"wls"`. The draws are then projected as with the default
weights, as Panagiotelis et al. (2023, Section 7.3, MinTShr) do for
probabilistic forecasts. The errors are reported in the fit (see Value).

Deviations: the papers estimate the weights from the in-sample
one-step-ahead residuals of the base models and take them as
proportional for all horizons (Wickramasuriya et al. 2019). Here the
`h`-step errors of inner forecasts estimate the weights of each horizon,
one set per age (every age is a hierarchy of its own). Package-specific
(reference unknown): forecast years beyond `horizon` take the weights of
the longest inner horizon (a projection does not change when its weights
are multiplied by a constant, so this is the papers' proportionality
beyond `horizon`); an inner horizon with fewer than 2 error vectors with
observed values for every stratum takes the weights of the next shorter
one; ages without any take the weights from the draws, with a warning.
Each age and horizon has at most `n_origins` errors, from overlapping
forecasts; with many strata the shrinkage intensity is then often high,
and `"mint"` near `"wls"` with the error variances.

## References

Panagiotelis A, Gamakumara P, Athanasopoulos G, Hyndman RJ (2023).
Probabilistic forecast reconciliation: properties, evaluation and score
optimisation. *European Journal of Operational Research* 306(2):
693-706.
[doi:10.1016/j.ejor.2022.07.040](https://doi.org/10.1016/j.ejor.2022.07.040)

Schefzik R, Thorarinsdottir TL, Gneiting T (2013). Uncertainty
quantification in complex simulation models using ensemble copula
coupling. *Statistical Science* 28(4): 616-640.
[doi:10.1214/13-STS443](https://doi.org/10.1214/13-STS443)

Clark M, Gangopadhyay S, Hay L, Rajagopalan B, Wilby R (2004). The
Schaake shuffle: a method for reconstructing space-time variability in
forecasted precipitation and temperature fields. *Journal of
Hydrometeorology* 5(1): 243-262.

Hyndman RJ, Ahmed RA, Athanasopoulos G, Shang HL (2011). Optimal
combination forecasts for hierarchical time series. *Computational
Statistics & Data Analysis* 55(9): 2579-2589.

Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast
reconciliation for hierarchical and grouped time series through trace
minimization. *Journal of the American Statistical Association*
114(526): 804-819.
[doi:10.1080/01621459.2018.1448825](https://doi.org/10.1080/01621459.2018.1448825)

Wickramasuriya SL, Turlach BA, Hyndman RJ (2020). Optimal non-negative
forecast reconciliation. *Statistics and Computing* 30: 1167-1182.
[doi:10.1007/s11222-020-09930-0](https://doi.org/10.1007/s11222-020-09930-0)

Schäfer J, Strimmer K (2005). A shrinkage approach to large-scale
covariance matrix estimation and implications for functional genomics.
*Statistical Applications in Genetics and Molecular Biology* 4(1):
Article 32.
[doi:10.2202/1544-6115.1175](https://doi.org/10.2202/1544-6115.1175)

## Examples

``` r
x <- mt_example_regions(years = 2000:2019, regions = 3)
h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
fit <- mt_fit(reconciled(lee_carter(), h, method = "wls"), mt_aggregate_strata(x, h))
fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
# \donttest{
# the Schaake shuffle of the separate Lee-Carter forecasts, then MinT with
# the covariance of the errors of forecasts from 4 inner origins
sm <- reconciled(lee_carter(), h, method = "mint", dependence = "schaake", covariance = "errors",
                 years = 8, n_origins = 4, horizon = 3)
fit_sm <- mt_fit(sm, mt_aggregate_strata(x, h))
fit_sm$reconcile$errors$n[1:3, ]
#>        horizon
#> age     1 2 3
#>   40-44 4 4 3
#>   45-49 4 4 3
#>   50-54 4 4 3
fc_sm <- mt_forecast(fit_sm, h = 5, n_draws = 50, seed = 1)
fc_sm$schaake$dates
#> [1] 15
# score-optimal reconciliation, trained on forecasts of 2012-2019 from 4 origins
es <- reconciled(lee_carter(), h, method = "energy_score", years = 8, n_origins = 4, horizon = 3,
                 control = list(max_iter = 200))
fit_es <- mt_fit(es, mt_aggregate_strata(x, h))
fit_es$reconcile$scoreopt$table
#>        age forecasts      scale iterations energy_ols energy_opt
#>     <char>     <num>      <num>      <num>      <num>      <num>
#>  1:  40-44        11   80.74848        200   110.6156   92.85656
#>  2:  45-49        11  120.94848        200   183.2950  132.70090
#>  3:  50-54        11  187.08182        200   222.4106  166.41359
#>  4:  55-59        11  293.97273        200   249.2886  202.12100
#>  5:  60-64        11  461.76970         22   399.3858  404.22064
#>  6:  65-69        11  711.45455         17   562.5127  650.40601
#>  7:  70-74        11 1068.51515        200   898.4196  700.22405
#>  8:  75-79        11 1529.67879        200   882.9746  721.68237
#>  9:  80-84        11 2004.54545         36  1168.2083 1128.61859
#> 10:  85-89        11 2283.83333         13  1395.2534 2110.01867
#> 11:    90+        11 4201.22727         13  2289.7388 3204.27096
fc_es <- mt_forecast(fit_es, h = 5, n_draws = 50, seed = 1)
# }
```
