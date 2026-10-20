# Reconcile forecast draws over a hierarchy

Makes forecasts of nested or crossed strata add up, draw by draw. For
each age, year and draw, the base forecasts `y` of all strata in the
hierarchy (bottom series and aggregates) are projected onto the coherent
subspace, `y - W U (U' W U)^-1 U' y`, where `U' y = 0` are the
aggregation constraints (each aggregate equals the sum of its bottom
series). With `additive = "deaths"` in the hierarchy, the projection
works on expected deaths (rates times exposures; exposures of aggregates
are the sums of their bottom series), and reconciled rates are
reconciled deaths divided by exposures. Methods choose `W`:

## Usage

``` r
mt_reconcile(
  draws,
  hierarchy = NULL,
  method = c("wls", "mint", "structural", "ols", "bottom_up", "proportional"),
  shrink = "auto",
  exposure = NULL,
  dependence = c("none", "schaake", "empirical"),
  seed = NULL,
  fixed = NULL,
  var = NULL,
  total = NULL,
  parts = NULL
)
```

## Arguments

- draws:

  An [mt_draws](mt_draws.md) object with the strata of the hierarchy.

- hierarchy:

  An [`mt_hierarchy()`](mt_hierarchy.md) (or, for one variable, `var`
  and `total`).

- method:

  `"wls"` (default), `"mint"`, `"structural"`, `"ols"`, `"bottom_up"` or
  `"proportional"` (with `fixed`). The vignette "Coherent forecasts:
  totals and their parts" compares them. Score-optimal reconciliation
  (`"energy_score"`) is learnt from past forecasts and is available in
  [`reconciled()`](reconciled.md) only.

- shrink:

  Shrinkage intensity for `"mint"` (number between 0 and 1), or
  `"auto"`.

- exposure:

  Exposures `[age, period, stratum]` or `[age, stratum]` for the
  forecast years (default: last training year). Only bottom strata are
  used, and they are enough; aggregates are summed.

- dependence:

  `"none"` (default: keep the draws as they are), `"schaake"` (the
  Schaake shuffle of the draws of separately fitted strata before the
  reconciliation, see Details) or `"empirical"` (deprecated). Draws of
  joint models are not reordered: `"schaake"` then gives a warning and
  is ignored. The shuffle needs the training data of the draws (as
  [`mt_forecast()`](mt_forecast.md) returns them), with at least two
  more training years than forecast years.

- seed:

  Random seed for the Schaake shuffle (order of the draws, choice of
  dates, ties).

- fixed:

  Strata whose forecasts are kept unchanged (labels, or `"total"` for
  the grand total of the hierarchy), for example a trusted all-cause
  forecast that the causes must add up to. The projection gives these
  nodes near-zero variance, so the other strata absorb the whole gap;
  `$reconciliation$max_relative_change_fixed` reports how closely they
  were kept; they move, with a warning, when their own forecasts do not
  add up (a node fixed together with all its parts) or when every bottom
  stratum below them would have to be negative. Needs a projection
  method or `"proportional"`.

- var, total, parts:

  One-variable form: stratum variable and its total level (`parts` is
  ignored: all other levels are parts).

## Value

The reconciled `mt_draws`, with diagnostics in `$reconciliation`
(method, largest relative gap before reconciliation, share of bottom
values held at zero, mean shrinkage intensity). After the Schaake
shuffle, `$schaake` holds the number of dates (`dates`), their start
years (`years`), the start year of the date of each draw (`draw_year`),
the block of each draw (`block`), the number of blocks (`blocks`), the
number of shuffled margins (`margins`) and the share of margins with
tied template values (`tied`); `$component` is then dropped, since the
values of one draw come from different base draws.

## Details

- `"ols"`: identity (Hyndman, Ahmed, Athanasopoulos and Shang 2011);

- `"structural"`: diagonal with the number of bottom series in each node
  (structural scaling, Athanasopoulos, Hyndman, Kourentzes and
  Petropoulos 2017);

- `"wls"`: diagonal with the variance of each stratum's base draws
  (weighted least squares, Hyndman, Lee and Wang 2016; see below);

- `"mint"`: covariance of the base draws across strata, shrunk towards
  its diagonal with the intensity of Schäfer and Strimmer (2005) by
  default (minimum-trace reconciliation, Wickramasuriya, Athanasopoulos
  and Hyndman 2019; see below);

- `"bottom_up"`: aggregates are the sums of the bottom forecasts.

**`W` of `"wls"` and `"mint"`.** In the papers, `W` holds the variances
(Hyndman, Lee and Wang 2016) or the shrunk covariance matrix
(Wickramasuriya et al. 2019) of the in-sample one-step-ahead residuals
of the base models. Here `W` is the variance (`"wls"`) or the
Schäfer-Strimmer shrunk covariance (`"mint"`) of the base forecast draws
across draws, computed separately in each age and forecast year cell:
the variance of the model's own `h`-step forecast of that cell (of
expected deaths, or of rates when the levels share exposures; without
the Poisson variation of the deaths unless predictive draws are
reconciled). It changes with the horizon and follows the uncertainty the
model reports, not its past errors: a stratum whose model is too narrow
is changed too little. For strata fitted by separate models the base
draws are independent across strata, so the shrunk covariance is close
to its diagonal and `"mint"` gives nearly the result of `"wls"`; it
differs where a joint model (for example [`li_lee()`](li_lee.md))
simulates correlated strata, or after the Schaake shuffle.
[`reconciled()`](reconciled.md) with `covariance = "errors"` estimates
`W` from the errors of past forecasts instead, the input of the papers.

**Non-negative forecasts.** Projections can make bottom values negative
where the base forecasts disagree strongly. For the draws where this
happens, the bottom values `b` are those of the non-negative
reconciliation of Wickramasuriya, Turlach and Hyndman (2020): they
minimise the same weighted distance `(y - S b)' W^-1 (y - S b)` subject
to `b >= 0` (their Eq. 3), and the aggregates are their sums. The
quadratic program is solved exactly by their block principal pivoting
algorithm (their Algorithm 3.1: all infeasible values exchanged while
their number falls, at most 3 times without a fall, then one value at a
time), started from the unconstrained solution, with the projection
above for each set of bottom values held at zero. Infeasibilities below
1e-10 times the largest absolute base forecast of the draw count as zero
(the paper uses 1e-12 in absolute terms). Draws without negative values
keep the unconstrained solution, which is then the optimum. A warning
reports when more than 5% of the draws have bottom values at zero,
because the distributions of the affected strata (often small causes)
are then truncated at zero. Earlier versions of mortensor held negative
values at zero and projected the rest again, a heuristic without a
guaranteed optimum; it was removed.

**Fixed nodes.** `fixed` keeps the forecasts of some nodes, as in the
reconciliation with immutable forecasts of Zhang, Kang, Panagiotelis and
Li (2023). The package gives the fixed nodes a near-zero variance in `W`
(1e-10 times the largest variance); with `"mint"` the other nodes get
their covariance conditional on the fixed nodes. The projection is then
the weighted least-squares solution with the fixed nodes as equality
constraints, up to a relative change of the fixed nodes of about 1e-9,
and the other strata absorb the whole gap. With negative bottom values,
the same weights enter the non-negative quadratic program, so the fixed
nodes are kept as equality constraints there as well, as long as
non-negative parts can add up to them. When they cannot (a fixed node
below zero, or a node fixed together with parts whose forecasts do not
add up to it), the fixed nodes move, with a warning.

**Proportional.** `"proportional"` keeps the `fixed` nodes and scales
the bottom strata below each of them by one factor per draw, so that
they add up to the fixed node; nodes between them are recomputed as
sums. This is a single-level top-down method (Gross and Sohl 1990) with
forecast proportions (Athanasopoulos, Ahmed and Hyndman 2009): the
proportions are those of the bottom forecasts below the fixed node in
the same draw, while Athanasopoulos et al. (2009) derive them level by
level down the hierarchy (the same for two levels). It never sets a part
to zero. Fixed nodes must not share parts.

**Dependence between strata.** Strata forecast by separate models have
independent draws. With `dependence = "none"` (default), the draws are
reconciled as they are, so sums of separately fitted parts have narrower
intervals than they would if the parts moved together. Joint models (for
example [`li_lee()`](li_lee.md), [`coherent_fdm()`](coherent_fdm.md) and
[`coda()`](coda.md) with `coherent = TRUE`) carry the dependence between
strata in the model itself; their draws are never reordered (`"schaake"`
then gives a warning and is ignored). Two methods give separately fitted
strata a dependence before the reconciliation: the Schaake shuffle takes
it from the observed changes of the training years
(`dependence = "schaake"`, below), ensemble copula coupling (Schefzik,
Thorarinsdottir and Gneiting 2013) from the draws of a joint model
([`mt_couple()`](coupled.md) for draws before `mt_reconcile()`, and
[`coupled()`](coupled.md) or `reconciled(dependence = "ecc")` for
models).

**The Schaake shuffle** (`dependence = "schaake"`; Clark, Gangopadhyay,
Hay, Rajagopalan and Wilby 2004, Section 3b). Status: published, adapted
to mortality forecasts as listed below. In the paper, the ensemble
members of each station and variable are reordered so that their ranks
follow those of observations from as many dates of the historical
record, the same dates for all stations and variables (their eqs. 5-12);
the later lead times take the following days of the same dates, which
keeps the persistence over lead times. Here a margin is an age, a
forecast year and a stratum. The template value of a date, a training
year `tau`, for the margin of forecast year `h` (counted from the last
training year) is the observed change of the log death rate from `tau`
to `tau + h`: one observed trajectory per date, over consecutive years.
The dates are the `n` training years whose next `H` years are observed
(`H` the number of forecast years, so `n` is the number of training
years minus `H`), from the training data of the draws only. In every
margin the draws receive the sorted values of that margin in the rank
order of the template values of their dates. Every margin keeps exactly
its values, and the draws take the rank correlations of the observed
changes between strata, between ages and between forecast years. Ties
among template values are ranked at random, as in the paper (Section
4c); margins with a missing template value keep the order of their
draws, as the paper's code skips vectors with missing values (appendix),
with a warning. Schefzik et al. (2013, Section 3.4) describe the shuffle
as the empirical copula of the historical record applied to samples of
the margins; ensemble copula coupling takes the copula from the draws of
a model instead.

Package-specific (reference unknown): (1) the template holds changes of
log rates from the start year of each date, where the paper takes the
observed values themselves (from dates within seven days of the forecast
date), so that the ranks follow the variation of past changes around the
trend of mortality; (2) deaths of zero count as 0.5 in the template; (3)
the paper has as many dates as ensemble members. With more draws than
dates, the draws are taken in a random order in blocks of `n`, and each
block is shuffled with all `n` dates; a last block with fewer draws uses
a random subset of distinct dates, as do fewer draws than dates (the
paper selects its dates at random). The draws of a block take the ranks
of the template among themselves, and the same dates recur in every
block.

The shuffle assumes that the dependence of the training years holds in
the forecast years (Clark et al. 2004, Section 5). It sets the
dependence between the ages and years of a stratum as well, so measures
over ages, such as the life expectancy of a single stratum, change:
Lee-Carter draws move all ages of a stratum together, and after the
shuffle they move together as much as the observed changes did, which
narrows the intervals of life expectancy where these are less correlated
between ages. The observed changes carry the Poisson variation of the
deaths, which is independent between strata and ages; where deaths are
few, it weakens the dependence the template gives. Templates are short:
with 19 training years and 5 forecast years there are 14 dates, with 16
forecast years 3.

`dependence = "empirical"` is deprecated (it warns): a package-specific
reordering of whole paths of each stratum by one-year changes of log
rates, averaged over ages, of training years drawn with replacement. It
keeps the dependence between the ages and years of each stratum. Use
`"schaake"`, the published method.

The one-variable interface with `var` and `total` is still accepted. It
reconciles deaths when the exposures of the parts add up to that of the
total (sexes, regions) and rates when each part has the exposure of the
total (causes of death) or there are no exposures; give a hierarchy with
`additive` to choose.

## References

Hyndman RJ, Ahmed RA, Athanasopoulos G, Shang HL (2011). Optimal
combination forecasts for hierarchical time series. *Computational
Statistics & Data Analysis* 55(9): 2579-2589.

Hyndman RJ, Lee AJ, Wang E (2016). Fast computation of reconciled
forecasts for hierarchical and grouped time series. *Computational
Statistics & Data Analysis* 97.

Athanasopoulos G, Hyndman RJ, Kourentzes N, Petropoulos F (2017).
Forecasting with temporal hierarchies. *European Journal of Operational
Research* 262(1): 60-74.

Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast
reconciliation for hierarchical and grouped time series through trace
minimization. *Journal of the American Statistical Association*
114(526): 804-819.
[doi:10.1080/01621459.2018.1448825](https://doi.org/10.1080/01621459.2018.1448825)

Wickramasuriya SL, Turlach BA, Hyndman RJ (2020). Optimal non-negative
forecast reconciliation. *Statistics and Computing* 30: 1167-1182.
[doi:10.1007/s11222-020-09930-0](https://doi.org/10.1007/s11222-020-09930-0)

Zhang B, Kang Y, Panagiotelis A, Li F (2023). Optimal reconciliation
with immutable forecasts. *European Journal of Operational Research*
308(2): 650-660.
[doi:10.1016/j.ejor.2022.11.035](https://doi.org/10.1016/j.ejor.2022.11.035)

Gross CW, Sohl JE (1990). Disaggregation methods to expedite product
line forecasting. *Journal of Forecasting* 9(3): 233-254.

Athanasopoulos G, Ahmed RA, Hyndman RJ (2009). Hierarchical forecasts
for Australian domestic tourism. *International Journal of Forecasting*
25: 146-166.

Clark M, Gangopadhyay S, Hay L, Rajagopalan B, Wilby R (2004). The
Schaake shuffle: a method for reconstructing space-time variability in
forecasted precipitation and temperature fields. *Journal of
Hydrometeorology* 5(1): 243-262.

Schefzik R, Thorarinsdottir TL, Gneiting T (2013). Uncertainty
quantification in complex simulation models using ensemble copula
coupling. *Statistical Science* 28(4): 616-640.
[doi:10.1214/13-STS443](https://doi.org/10.1214/13-STS443)

Schäfer J, Strimmer K (2005). A shrinkage approach to large-scale
covariance matrix estimation and implications for functional genomics.
*Statistical Applications in Genetics and Molecular Biology* 4(1):
Article 32.
[doi:10.2202/1544-6115.1175](https://doi.org/10.2202/1544-6115.1175)

## Examples

``` r
x <- mt_example_regions(years = 2005:2019, regions = 3)
h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
xa <- mt_aggregate_strata(x, h)
fc <- mt_forecast(mt_fit(lee_carter(), xa), h = 3, n_draws = 50, seed = 1)
rc <- mt_reconcile(fc, h, method = "wls")
rc$reconciliation
#> $method
#> [1] "wls"
#> 
#> $max_relative_gap_before
#> [1] 0.9210924
#> 
#> $share_set_to_zero
#> [1] 0
#> 
#> $mean_shrinkage
#> [1] NA
#> 
#> $share_draws_with_zeros
#> [1] 0
#> 
#> $max_relative_change_fixed
#> [1] NA
#> 
# the strata take the dependence of the observed changes of 12 training years
rs <- mt_reconcile(fc, h, method = "wls", dependence = "schaake", seed = 1)
rs$schaake$years
#>  [1] 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016
```
