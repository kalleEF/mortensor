# SAINT model: trend of a reference population and stationary spreads

The SAINT model (Spread Adjusted InterNational Trend) of Jarner and
Kryger (2011) in the revised form of Jarner and Jallbjørn (2022). It
forecasts the adult mortality of small populations as deviations
(spreads) from the trend of a large reference population. The trend
carries the long-term development. The spreads are stationary, so the
populations do not diverge from the reference.

## Usage

``` r
saint(
  reference,
  by = NULL,
  ecm = NULL,
  spread_ar = "diagonal",
  background = TRUE
)
```

## Arguments

- reference:

  Reference strata: stratum labels, or a named list of values of stratum
  variables, for example `list(region = "Deutschland")`.

- by:

  Stratum variable(s) with a trend each, for example `"sex"`. `NULL`:
  one trend for all strata.

- ecm:

  `NULL`, or the group of `by` whose mean changes give the common drift
  of the error-correction model that couples the level and slope of two
  trends (Section 5.2.1), for example `"female"`. Needs `by` with two
  groups.

- spread_ar:

  `A` of the VAR(1) of the spread: `"diagonal"` (diagonal, estimated,
  `0 <= A_ii <= 0.99`), `"full"` (unrestricted, eq. 5.12), or a number
  `a` in `[0, 1)` for `A = a I` (the paper's choice for Denmark: 0.99).

- background:

  Whether the trend has the background mortality `exp(zeta_t)` of eq.
  5.3 (default). `FALSE` leaves it out, a model of the class of Section
  4.2, which allows empty components. Use it where the data do not
  identify a background, for example rates smoothed over age; a warning
  says when the estimate of the background reaches its lower bound,
  which makes its random walk and the forecasts at young ages unstable.

## Value

A model specification (joint over strata).

## Details

The reference is one or more strata of the data (`reference`), for
example the national total or several countries. Their deaths and
exposures are pooled into one reference population (Jarner and Kryger
2011, Section 3.1). With `by` (for example `"sex"`), each group gets its
own trend; both papers fit women and men separately. All other strata of
a group are targets with a spread from the group's trend. A group with a
single reference stratum forecasts that stratum by the trend. In a group
with several reference strata, each of them is also a target, as Jarner
and Kryger (2011, concluding remarks) propose for a group of related
populations.

## Version

This is the version of Jarner and Jallbjørn (2022), the model in use at
the Danish pension fund ATP. Its authors replaced three components of
the original model: the deterministic parametric trend by a stochastic
frailty model, so that forecasts describe the uncertainty of the trend
(Section 2.2); the trend's age pattern by one that fits high ages better
(Section 2.4); and the three spread regressors (level, slope, curvature)
by five (Section 5.4). Section and equation numbers below refer to
Jarner and Jallbjørn (2022) unless the 2011 paper is named.

## Trend (Sections 4 and 5.1)

The reference rates are
`mu(t, x) = exp(-sigma2 M(t, x)) mu0(t, x) + exp(zeta_t)` (eqs.
5.1-5.3), with baseline
`log mu0(t, x) = alpha_t + beta_t (x - 75) + kappa_t (x - 75) 1{x < 75}`
(eq. 5.2), background mortality `exp(zeta_t)`, and Gamma-distributed
frailty with mean 1 and variance `sigma2` (Example 4.1).
`exp(-sigma2 M)` is the mean frailty of the survivors of a cohort, where
`M(t, x)` cumulates the observed rates of the cohort, background
subtracted, from the first age of the data; years before the first year
take the rates of the first year (eqs. 4.7-4.8). Deaths are Poisson (eq.
4.9). For a given `sigma2`, the indexes of every year are the estimates
of the EM algorithm of Appendix B; `sigma2` maximises the profile
pseudo-likelihood (eq. 4.14). `kappa` and `zeta` follow random walks
without drift (eq. 5.4). `alpha` and `beta` follow a bivariate random
walk with drift (Section 4.5) or, with `ecm`, the error-correction model
of Section 5.2.1 for two groups (eqs. 5.5-5.6): the equilibrium
differences are the differences at the jump-off year and the common
drift is the mean change of the group named by `ecm` (women in the
paper). The median forecast is then a random walk with this drift for
both groups, and in every path the differences between the groups stay
stationary. Rates are forecast by eqs. 4.15-4.17: the mean frailty
follows from the baseline cumulated along each cohort. Selection by
frailty makes the rates of improvement at high ages rise over time
(Section 3.3).

## Spread (Section 5.4)

A target's rates are the trend times `exp(y_t' r(x))` (eq. 5.10), with
regressors `r_i(x) = min(1, max(0, x_i - x) / 20)` and knots `x_i` = 40,
60, 80, 100, 120 (eq. 5.11). The spread parameters `y_t` of each year
are Poisson maximum likelihood estimates given the fitted trend (eq.
5.9); the trend is estimated first and kept fixed, the two-stage
estimation of Jarner and Kryger (2011, concluding remarks). They follow
a VAR(1) model without mean, `y_t = A y_{t-1} + u_t` with normal
innovations (eq. 5.12), so the spread reverts towards zero. `spread_ar`
sets `A`; the paper recommends a diagonal `A` with `0 <= A_ii < 1` and
imposes `A = 0.99 I` for Denmark, a deviation half-life of about 69
years (Section 5.4.1). Each target has a VAR(1) of its own, "a separate
VAR model" per target population (Section 5.4.1).

## Forecasts and uncertainty

Draws include the innovations of the trend indexes and of the spread.
Trend and spread are independent (Jarner and Kryger 2011, eq. 5).
Parameters are fixed, as in Jarner and Kryger (2011, Section 2.3.3), who
leave parameter uncertainty aside. Forecasts start from the fitted rates
of the last year; the paper uses no jump-off correction (2022, Section
5.3). Predictive draws (`mt_forecast(predictive = TRUE)`) add Poisson
variation of the deaths given the rates, the observation model of the
papers.

## Ages

Both papers use single ages 20-100 and model adult mortality only (2011,
Section 1.2). `saint()` needs consecutive single-year age groups from
age 20 without an open interval: the frailty term follows cohorts year
by year. Ages below 20 are refused. Ages above 100 are accepted; the
parametric trend extends to them (2022, Section 5.1.1). Drop an open
interval with `mt_subset(x, ages = c(20, 100))`; estimate single years
from age groups with [`mt_ungroup()`](mt_ungroup.md) or
[`step_ungroup()`](mt_prep.md).

## Corrections and package-specific choices

Three formulas of Jarner and Jallbjørn (2022) are read as evident
misprints. In (B14) and (B16) the denominator `c F - G` is the rate
`c F + G` of (B7)-(B8). In the second line of (4.16) the baseline of
year `t_max` is that of year `t - 1`. The common drift of Section 5.2.1
divides the sum of changes by `t_max - t_min + 1`; the package divides
by the number of changes, the empirical average that the text describes.

Package-specific (reference unknown):

- the fixed point of the EM algorithm is computed by Newton-Raphson
  sweeps (Fisher scoring) over the parameters of each year with the mean
  frailty held fixed, alternated with updates of the mean frailty;
  Section 4.4 names Newton-Raphson as the alternative to EM, and the
  tests check that both give the same estimates;

- `sigma2` is searched from 0 by doubling (0.025, 0.05, ..., 3.2) and
  then by [`stats::optimize()`](https://rdrr.io/r/stats/optimize.html)
  around the best value; a warning (class `mt_warning_saint`) is given
  at the upper end;

- starting values, and a lower bound of `1e-10` on the background rate,
  with a warning when it binds (see `background`);

- on age ranges narrower than 20-100: `kappa` is set to 0 without ages
  on both sides of 75, and spread regressors that are 0 at all ages, or
  1 at all ages like the level, are left out;

- the innovations of `kappa` and `zeta` are independent of each other
  and of `alpha` and `beta` (the paper gives their marginal
  distributions), and the spread innovations of different targets are
  independent; variances are maximum likelihood estimates;

- the adjustment coefficients and the innovation covariance of the
  error-correction model are Gaussian maximum likelihood estimates given
  the equilibrium differences and the drift; a warning is given when the
  model does not correct deviations (eigenvalues of eq. 5.8 outside (-1,
  1));

- `spread_ar = "diagonal"`: each `A_ii` is the least-squares coefficient
  of its component, bounded to `[0, 0.99]` (a warning names the strata
  where 0.99 binds); `"full"`: conditional maximum likelihood (the 2011
  paper uses Yule-Walker estimates of the same model), with a warning
  when `A` is not stationary.

`background = FALSE` is an opt-in variant of this package: the SAINT
specification of Section 5.1 has a background, the class of models of
Section 4.2 allows it to be left out.

Not implemented: the deterministic trend and the three regressors of
2011, the error-correction model with free equilibrium and drift, and
frailty distributions other than Gamma (Appendix A).

## References

Jarner SF, Kryger EM (2011). Modelling adult mortality in small
populations: the SAINT model. *ASTIN Bulletin* 41(2): 377-418.
[doi:10.2143/AST.41.2.2136982](https://doi.org/10.2143/AST.41.2.2136982)

Jarner SF, Jallbjørn S (2022). The SAINT model: a decade later. *ASTIN
Bulletin* 52(2): 483-517.
[doi:10.1017/asb.2021.37](https://doi.org/10.1017/asb.2021.37)

## See also

[`li_lee()`](li_lee.md) and [`coherent_fdm()`](coherent_fdm.md), other
models whose strata keep bounded distances to a common trend.

## Examples

``` r
x <- mt_example_regions(years = 1995:2019, regions = 3, ages = "single")
h <- mt_hierarchy(region = "all", quintile = "all")
x <- mt_aggregate_strata(mt_subset(x, ages = c(40, 100)), h)
x <- mt_subset(x, quintile == "all", drop = TRUE)
fit <- mt_fit(saint(reference = list(region = "all"), by = "sex"), x)
fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
p <- mt_parameters(fit)
p[parameter == "sigma2"]
#>         stratum parameter  index       value
#>          <char>    <char> <char>       <num>
#> 1: trend female    sigma2        0.000000000
#> 2:   trend male    sigma2        0.002233368
```
