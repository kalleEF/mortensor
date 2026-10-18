# Models with reconciled forecasts

Wraps a model specification so that its forecasts are reconciled over a
hierarchy every time it is forecast, also inside
[`mt_backtest()`](mt_backtest.md). If the training data contain only the
bottom strata, the aggregate strata are added with
[`mt_aggregate_strata()`](mt_aggregate_strata.md) before fitting, so
that the base model forecasts every node; exposures for the forecast
(`exposure` of [`mt_forecast()`](mt_forecast.md)) can then be given for
the bottom strata only, and those of the aggregates are their sums.
Predictive draws (used by backtests) are coherent as well: the
observation noise of the base model (Poisson and extra-Poisson variation
for models of smooth rates, none for models of observed rates such as
[`rwd()`](rwd.md); for [`ensemble()`](ensemble.md) and
[`stratified()`](stratified.md) models that of the component or group
each draw comes from) is drawn for the bottom strata only, independently
per cell, and the aggregates are their sums.

## Usage

``` r
reconciled(
  spec,
  hierarchy,
  method = c("wls", "mint", "structural", "ols", "bottom_up", "proportional",
    "energy_score"),
  shrink = "auto",
  dependence = c("auto", "empirical", "none"),
  fixed = NULL,
  years = 10,
  n_origins = 5,
  horizon = 5,
  control = list(),
  seed = 1
)
```

## Arguments

- spec:

  Base model specification (fitted to all strata).

- hierarchy:

  An [`mt_hierarchy()`](mt_hierarchy.md).

- method:

  A method of [`mt_reconcile()`](mt_reconcile.md), or `"energy_score"`
  (score-optimal reconciliation, see Details).

- shrink:

  See [`mt_reconcile()`](mt_reconcile.md).

- dependence:

  `"auto"` (default: `"none"` for models fitted jointly across strata,
  whose draws already carry their dependence, such as
  [`li_lee()`](li_lee.md) and [`tucker3()`](tucker3.md); `"empirical"`
  for models fitted per stratum), `"empirical"` or `"none"`. With
  `"empirical"` the paths of the strata are reordered before
  reconciliation (see [`mt_reconcile()`](mt_reconcile.md)), so that
  aggregates keep the observed co-movement of the strata; the
  observation noise of predictive draws is independent across cells, as
  for models without a hierarchy.

- fixed:

  Strata whose forecasts are kept (see
  [`mt_reconcile()`](mt_reconcile.md)); not with `"energy_score"`.

- years, n_origins, horizon:

  Training of `method = "energy_score"`: the length of the training
  window (the last `years` observed years), the number of inner origins
  in it (as in [`recalibrated()`](recalibrated.md)) and the number of
  years each inner forecast covers (horizons 1 to `horizon`, cut at the
  last training year).

- control:

  Settings of the training for `method = "energy_score"`: a list with
  any of `n_draws` (draws per inner forecast, default 100),
  `draws_per_step` (draws per training forecast in each step, 50),
  `max_iter` (number of steps, 500), `eta` (learning rate, 0.001),
  `beta1` (0.9), `beta2` (0.999), `tol` (1e-4) and `epsilon` (1e-8), as
  in `ProbReco::scoreopt.control()`. The optimisation stops early when
  no parameter moves by more than `tol` in a step, so `tol` must be
  below `eta` (0 runs all steps). `max_iter = 0` gives OLS
  reconciliation.

- seed:

  Random seed of the training for `method = "energy_score"` (inner
  forecasts, observation noise and the draws of each step), so that fits
  are reproducible.

## Value

A model specification. With `method = "energy_score"` its fit holds
`$reconcile$scoreopt`: `d` `[bottom stratum, age]`, `G`
`[bottom stratum, stratum, age]` on the scale of the reconciliation, the
inner `origins`, the `control` settings and a `table` per age with the
number of training forecasts, the scale of the bottom strata, the number
of steps, and the mean in-sample energy score with OLS (`energy_ols`)
and with the estimated \\d\\ and \\G\\ (`energy_opt`).

## Details

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

## References

Panagiotelis A, Gamakumara P, Athanasopoulos G, Hyndman RJ (2023).
Probabilistic forecast reconciliation: properties, evaluation and score
optimisation. *European Journal of Operational Research* 306(2):
693-706.
[doi:10.1016/j.ejor.2022.07.040](https://doi.org/10.1016/j.ejor.2022.07.040)

## Examples

``` r
x <- mt_example_regions(years = 2000:2019, regions = 3)
h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
fit <- mt_fit(reconciled(lee_carter(), h, method = "wls"), mt_aggregate_strata(x, h))
fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
# \donttest{
# score-optimal reconciliation, trained on forecasts of 2012-2019 from 4 origins
es <- reconciled(lee_carter(), h, method = "energy_score", years = 8, n_origins = 4, horizon = 3,
                 control = list(max_iter = 200))
fit_es <- mt_fit(es, mt_aggregate_strata(x, h))
fit_es$reconcile$scoreopt$table
#>        age forecasts      scale iterations energy_ols energy_opt
#>     <char>     <num>      <num>      <num>      <num>      <num>
#>  1:  40-44        11   80.74848        200   112.7909   97.18129
#>  2:  45-49        11  120.94848        200   181.4128  136.10143
#>  3:  50-54        11  187.08182        200   222.5451  174.70273
#>  4:  55-59        11  293.97273        200   260.7576  218.00871
#>  5:  60-64        11  461.76970         30   396.7081  362.29530
#>  6:  65-69        11  711.45455         18   581.3779  647.50836
#>  7:  70-74        11 1068.51515        200   881.1605  720.61784
#>  8:  75-79        11 1529.67879         38   906.8957  893.64292
#>  9:  80-84        11 2004.54545         35  1202.3287 1175.47593
#> 10:  85-89        11 2283.83333         18  1418.6084 1838.65283
#> 11:    90+        11 4201.22727         13  2269.7176 3327.87922
fc_es <- mt_forecast(fit_es, h = 5, n_draws = 50, seed = 1)
# }
```
