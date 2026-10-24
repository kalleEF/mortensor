# Relative model of subpopulations and a reference population

Joint model of subpopulations (for example deprivation quintiles) and
the population they belong to (for example England): the relative model
of Villegas and Haberman (2014, Section 2.4). The reference population
has a Lee-Carter model, with a cohort effect if `cohort = TRUE`; each
subpopulation's rates are the reference rates times an age pattern of
level differences and an index of improvement differences. Section and
equation numbers refer to that paper.

## Usage

``` r
relative_model(
  reference,
  by = NULL,
  cohort = FALSE,
  deviation = c("mrwd", "rwd", "constant"),
  trend_from = NULL
)
```

## Arguments

- reference:

  The reference population: stratum labels, a named list of values of
  stratum variables (for example `list(quintile = "England")`), or an
  [mt_data](mt_data.md) object (see Details).

- by:

  Stratum variable(s) with a model each, for example `"sex"`. `NULL`:
  one model for all strata.

- cohort:

  Whether the reference model has cohort effects (eq. 16; the paper's
  "APC reference"). `FALSE`: Lee-Carter reference.

- deviation:

  Forecast of the subpopulation indexes: `"mrwd"`, `"rwd"` or
  `"constant"` (see Details).

- trend_from:

  `NULL` (all years), or the first year of the reference index used for
  its random walk with drift.

## Value

A model specification (joint over strata). The fit (`fit$fits$models`,
one per group) holds the reference parameters (`ref`: `alpha`, `beta`,
`kappa`, `gamma` and the time series parameters), the subpopulation
parameters (`alpha`, `beta`, `kappa`, `drift`, `sigma`) and the
log-likelihoods, effective numbers of parameters and BIC of both parts.

## Model (Section 2.4, eqs. 16-18)

The reference rates are
`log mu'[x, t] = alpha'[x] + beta'[x] kappa'[t] + gamma'[t - x]` (eq.
16), the age-period-cohort model of Renshaw and Haberman (2006) with
cohort effects `gamma'`; with `cohort = FALSE` the cohort effects are
left out and the reference has the Lee-Carter model (the paper's "LC
reference"). The rates of subpopulation `g` are
`log mu[x, t, g] = log mubar[x, t] + alpha[x, g] + beta[x] kappa[t, g]`
(eq. 17). `exp(alpha[x, g])` is the average ratio of the mortality of
`g` to the reference (level differences); `kappa[t, g]` is the index of
the improvements of `g` relative to the reference (a falling index:
faster improvement), with the age pattern `beta[x]` common to all
subpopulations. With `beta[x] > 0`, subpopulations keep their order as
long as their indexes keep it (Section 2.4).

`mubar[x, t]` is the geometric mean of the reference rates over the
single ages of the subpopulations' age group (eq. 18). The reference can
cover more ages and years than the subpopulations (Section 2.4): the
paper used England and Wales at ages 10-99 in 1961-2009 as the reference
for English deprivation quintiles at ages 50-84 in 1981-2007.

## Reference

`reference` is either

- strata of the data: stratum labels or a named list of values of
  stratum variables, as in [`saint()`](saint.md), for example
  `list(quintile = "England")` after
  [`mt_aggregate_strata()`](mt_aggregate_strata.md). The reference
  strata of a group of `by` are pooled (sum of deaths and exposures)
  into one reference population. A single reference stratum of a group
  is forecast by the reference model; all other strata of the group are
  subpopulations. Reference and subpopulations then share ages and
  years, and `mubar` is the reference rate of the same age group;

- an [mt_data](mt_data.md) object with the reference population, which
  can cover more ages and years (Section 2.4). With `by`, its strata are
  matched to the groups of `by` by the values of these variables; the
  strata that match a group are pooled. Without `by`, all its strata are
  pooled into one reference. Every age group of the subpopulations must
  be a union of age groups of the reference, and every year of the
  subpopulations a year of the reference. All strata of the data are
  subpopulations. Preprocessing (`mt_fit(prep = )`) applies to the data,
  not to the reference.

Package-specific (reference unknown): where the reference has age groups
wider than one year inside an age group of the subpopulations, `mubar`
is their geometric mean weighted by the widths (eq. 18 for single ages);
years of the reference after the last year of the data are not used, so
that backtests use no later data (the paper's reference ends two years
after its subpopulation data, Section 3).

## Estimation (Sections 2.5 and 3.1)

Deaths of the reference and of the subpopulations are independent
Poisson counts, and all parameters maximise the total log-likelihood
(eq. 25), so the subpopulations also inform the reference parameters.
The paper uses "straightforward extensions" of the Newton-Raphson
procedure of Brouhns, Denuit and Vermunt (2002). Here each iteration
updates `alpha'`, `kappa'`, `beta'`, `gamma'`, `alpha`, `kappa` and
`beta` in turn by one-dimensional Newton steps of the total
log-likelihood, started from the Poisson Lee-Carter fit of the reference
and a singular value decomposition of the subpopulations' log ratios to
it. `alpha'` takes its steps from the reference deaths only:
`alpha[x, g]` absorbs any change of `alpha'` in the subpopulation part,
so its score there vanishes at the maximum and the steps reach the same
maximum (package-specific). Package-specific (reference unknown): with
cohort effects these iterations converge slowly, as known for this model
(Section 2.5), and their fixed point is found with Anderson acceleration
(Anderson 1965; Walker and Ni 2011), which extrapolates from the last
ten iterations. The iterations of the Lee-Carter reference are not
accelerated: they converge within a few hundred iterations, and
extrapolations can lead them to a fixed point of lower likelihood. The
iterations stop when one iteration changes no fitted log rate by more
than `1e-9` (at most 5,000 iterations; a warning of class
`mt_warning_relative_model` otherwise). The fit reproduces
`StMoMo::rh(approxConst = TRUE)` for a reference without subpopulations.

Constraints: `sum(beta') = 1`, `sum(kappa') = 0`, `sum(gamma') = 0` over
the cohorts (eqs. 19-21), `sum(beta) = 1` and `sum(kappa[, g]) = 0` for
every `g` (eqs. 22-23). With `cohort = TRUE`, the cohort effects also
have no linear trend, `sum((c - cbar) gamma'[c]) = 0` (eq. 26), the
approximate constraint of Hunt and Villegas (2015). It is imposed after
every iteration by their approximately invariant transformation, as in
`StMoMo::rh(approxConst = TRUE)` (Villegas, Kaishev and Millossovich
2018), the package of the first author: `alpha'`, `beta'`, `kappa'` and
`gamma'` change so that the cohort effects lose their linear trend and
the fitted rates stay the same where `kappa'` is linear. The
transformation is skipped where `kappa'` has no trend (package-specific
safeguard).

[`mt_parameters()`](mt_parameters.md) reports, per part, the
log-likelihood (with the `log d!` terms of eqs. 24-25), the effective
number of parameters and the BIC of Section 3.1, `L - 0.5 nu log N` with
`N` the number of cells (higher is better). The subpopulation part has
`nu = k m + k + n m - 1 - m` parameters for `k` age groups, `n` years
and `m` subpopulations, which gives the 171 of Table 1; the Lee-Carter
reference has `2 k' + n' - 2` (227 in Table 1), the age-period-cohort
reference `2 k' + n' + (k' + n' - 1) - 4`, counting every cohort (363
for the dimensions of Table 1, which reports 353).

## Forecasts (Sections 2.4 and 3.3)

Forecasts start from the fitted rates of the last year (eq. 32). The
reference index `kappa'` follows a random walk with drift (Section 2.4);
`trend_from` estimates it from the years from `trend_from` only, as the
paper did with the years after 1980, where the index was close to linear
(Section 3.3, following Haberman and Renshaw 2009). The subpopulation
indexes follow, by `deviation`:

- `"mrwd"`: a multivariate random walk with drift with correlated normal
  innovations (eq. 4), the paper's model (Sections 2.4 and 3.3);

- `"rwd"`: independent random walks with drift (eq. 3). Section 2.1
  notes that each subpopulation has the same forecast distribution under
  both processes, and that differences between subpopulations get wider
  intervals from independent walks; Section 3.1 compares both in a
  backtest (for the independent Lee-Carter model). With the same seed,
  the draws of the first subpopulation of a group are identical for
  both;

- `"constant"`: the indexes stay at their values of the last year, so
  all subpopulations improve as the reference (Scenario 2 of Section
  3.4, eq. 33: level differences only).

The drifts are the mean year-to-year changes and the innovation
(co)variances those of the changes around them, with divisor the number
of changes minus one (as
[`StMoMo::mrwd()`](https://rdrr.io/pkg/StMoMo/man/mrwd.html)). Draws
include the innovations only: parameters, drifts and covariances are
fixed, as in the paper, whose intervals "ignore any provision for
parameter uncertainty" (Section 3.1). The reference index and the
subpopulation indexes are independent, separate time series models as in
the paper. Each group of `by` (for example each sex) is a separate
model, as the paper fits women and men separately. Predictive draws
(`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
given the rates.

The drifting indexes let gaps between subpopulations widen or narrow
without limit; [`li_lee()`](li_lee.md),
[`coherent_fdm()`](coherent_fdm.md) and [`saint()`](saint.md) keep the
deviations from a common trend stationary.

Cohort effects: forecasts need no cohort beyond the fitted ones while
the horizon is at most the distance between the youngest age of the
subpopulations and the youngest age of the reference (the paper's case,
Section 2.4). Later cohorts follow an ARIMA(1,1,0) model with drift of
the fitted cohort effects, drawn with its innovations (package-specific:
the paper does not extrapolate cohorts; this is the default of
[`StMoMo::forecast.fitStMoMo()`](https://rdrr.io/pkg/StMoMo/man/forecast.fitStMoMo.html)).

## Data

The reference needs single-year consecutive ages without an open
interval for `cohort = TRUE` (the paper uses single ages, Section 3),
positive exposures, and deaths at every age (and in every cohort).
Subpopulations need deaths in every age group and every year.

## References

Villegas AM, Haberman S (2014). On the modeling and forecasting of
socioeconomic mortality differentials: an application to deprivation and
mortality in England. *North American Actuarial Journal* 18(1): 168-193.
[doi:10.1080/10920277.2013.866034](https://doi.org/10.1080/10920277.2013.866034)

Renshaw AE, Haberman S (2006). A cohort-based extension to the
Lee-Carter model for mortality reduction factors. *Insurance:
Mathematics and Economics* 38(3): 556-570.
[doi:10.1016/j.insmatheco.2005.12.001](https://doi.org/10.1016/j.insmatheco.2005.12.001)

Hunt A, Villegas AM (2015). Robustness and convergence in the Lee-Carter
model with cohort effects. *Insurance: Mathematics and Economics* 64:
186-202.
[doi:10.1016/j.insmatheco.2015.05.004](https://doi.org/10.1016/j.insmatheco.2015.05.004)

Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear
regression approach to the construction of projected lifetables.
*Insurance: Mathematics and Economics* 31(3): 373-393.
[doi:10.1016/S0167-6687(02)00185-3](https://doi.org/10.1016/S0167-6687%2802%2900185-3)

Haberman S, Renshaw A (2009). On age-period-cohort parametric mortality
rate projections. *Insurance: Mathematics and Economics* 45(2): 255-270.
[doi:10.1016/j.insmatheco.2009.07.006](https://doi.org/10.1016/j.insmatheco.2009.07.006)

Villegas AM, Kaishev VK, Millossovich P (2018). StMoMo: an R package for
stochastic mortality modeling. *Journal of Statistical Software* 84(3):
1-38. [doi:10.18637/jss.v084.i03](https://doi.org/10.18637/jss.v084.i03)

Anderson DG (1965). Iterative procedures for nonlinear integral
equations. *Journal of the ACM* 12(4): 547-560.
[doi:10.1145/321296.321305](https://doi.org/10.1145/321296.321305)

Walker HF, Ni P (2011). Anderson acceleration for fixed-point
iterations. *SIAM Journal on Numerical Analysis* 49(4): 1715-1735.
[doi:10.1137/10078356X](https://doi.org/10.1137/10078356X)

## See also

[`saint()`](saint.md), another model of populations relative to a
reference; [`cae()`](cae.md), whose populations also keep drifting
indexes.

## Examples

``` r
x <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile"),
             age = "age_group", value = c(deaths = "all_causes", population = "population"))
x <- mt_as_exposure(x, type = "midyear")
x <- mt_group_levels(x, "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)),
                     name = "quintile")
x <- mt_subset(x, ages = c(50, 85))
xa <- mt_aggregate_strata(x, mt_hierarchy(quintile = list(total = "England")))
fit <- mt_fit(relative_model(reference = list(quintile = "England"), by = "sex"), xa)
p <- mt_parameters(fit)
p[parameter == "bic"]
#>                  stratum parameter  index     value
#>                   <char>    <char> <char>     <num>
#> 1:      reference female       bic        -1062.329
#> 2: subpopulations female       bic        -4076.166
#> 3:        reference male       bic        -1197.945
#> 4:   subpopulations male       bic        -4268.884
fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
```
