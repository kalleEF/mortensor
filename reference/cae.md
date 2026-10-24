# Common age effect model for several populations

Joint model of the strata as populations that share their age effects
and keep their own period effects: the common age effect (CAE) model of
Kleinow (2015). Status: published. The model and its estimation follow
Kleinow (2015), and section and equation numbers refer to that paper
unless another is named; `index = "coherent"` forecasts with the model
of Enchev, Kleinow and Cairns (2017); the default forecast and the other
choices listed below are package-specific.

## Usage

``` r
cae(order = 1, common_by = NULL, index = c("rwd", "coherent"))
```

## Arguments

- order:

  Number of common age effects (`p` in the paper).

- common_by:

  Stratum variable(s) defining groups of strata with their own common
  age effects (package-specific, see Details).

- index:

  Time-series model of the period effects: `"rwd"` (one multivariate
  random walk with a drift for every population and order,
  package-specific) or `"coherent"` (a common drift for the first period
  effects and mean-reverting autoregressions for the others; Enchev,
  Kleinow and Cairns 2017). See Details.

## Value

A model specification. The fit (`fit$fits`) holds the common age effects
`beta` (one matrix per group), the centring `a` `[age, stratum]`, the
period effects `kappa` `[order, year, stratum]`, the drift `drift`, the
autoregressive coefficient `phi` (1 for a random walk) and the
innovation covariance `sigma` of the stacked period effects (order
within stratum), the common drift of each group `common_drift`
(`index = "coherent"`), `mse` and `bic`.

## Details

**Model (Sections 2 and 3.1, eq. 1-3).** The log death rates
`y[x, t, i] = log(D[x, t, i] / E[x, t, i])` of each population `i` are
centred over the years: `m[x, t, i] = y[x, t, i] - a[x, i]`, with
`a[x, i]` the mean over the years. The centred rates follow
`m[, , i] = beta kappa[, , i] + e[, , i]`. The age effects `beta` (ages
x `order`) are common to all populations; the period effects
`kappa[, , i]` (`order` x years) belong to each population. With one
population the model is the Lee-Carter model with `order` age-period
terms, estimated by a singular value decomposition; for `order = 1` it
gives the fitted rates and forecasts of
`lee_carter(method = "svd", adjust = "none")`.

**Identification (Section 2).** The columns of `beta` have length 1 and
are orthogonal to each other. The period effects of each population sum
to zero over the years, because the rates are centred. The sign of each
column of `beta` is chosen so that its sum is positive
(package-specific).

**Estimation (Section 3.2).** With `Q[i] = m[, , i] t(m[, , i])`, `beta`
is the orthogonal matrix that minimises the sum over populations of the
squared off-diagonal elements of `t(beta) Q[i] beta`: the paper's
least-squares modification of the common principal component analysis of
Flury (1984). Every population has the same weight whatever its size, as
in the paper, which prefers this to maximum likelihood because large
populations would dominate the age effects. The first `order` columns
are kept, and `kappa[, , i] = t(beta) m[, , i]`. The paper requires more
years than age groups (Section 2), and so does `cae()`: with fewer
years, use wider age groups ([`mt_regroup()`](mt_regroup.md)) or fewer
ages ([`mt_subset()`](mt_subset.md)); Kleinow fitted ages 18-52 and
53-87 separately.

**Forecasts.** Kleinow (2015) does not forecast. `index` chooses the
time-series model of the period effects. Both models start from the
fitted rates of the last year, as in Enchev, Kleinow and Cairns (2017,
Section 6.3) ([`anchored()`](anchored.md) starts them from the observed
rates), and keep `a` and `beta` fixed. Predictive draws
(`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
given the rates and the exposures (Kleinow makes no assumption about the
distribution of the errors).

- `index = "rwd"` (default, package-specific): the period effects of all
  populations and orders follow one multivariate random walk with drift,
  `kappa[t] = kappa[t - 1] + mu + e[t]`, with normal innovations whose
  covariance matrix is estimated from the year-to-year changes (divided
  by the number of changes minus one), so the innovations are correlated
  between populations and between orders. This is the random walk with
  drift of Lee and Carter (1992) for every index, and the multivariate
  random walk with drift that StMoMo uses for several period indexes
  (Villegas, Kaishev and Millossovich 2018), here over all populations
  together. Draws include the innovations and, as in
  [`lee_carter()`](lee_carter.md), the uncertainty of the drift (normal,
  with the innovation covariance divided by the number of changes;
  StMoMo keeps the drift fixed). Each population has its own drift, so
  the forecasts of different populations can diverge. Enchev, Kleinow
  and Cairns (2017, Section 6.2, Table 10) use this process for the
  population-specific period effects of two variants of the Li-Lee model
  and note that it lets the populations diverge.

- `index = "coherent"`: the forecasting model for the CAE model of
  Enchev, Kleinow and Cairns (2017, Section 6.2, Table 10, Model 3). The
  period effects of order 1 follow a multivariate random walk with one
  drift `d` common to all populations,
  `kappa1[t] = kappa1[t - 1] + d + e1[t]`, so that the populations do
  not diverge quickly. The period effects of order 2 follow a vector
  autoregression of order 1 with a diagonal matrix `Phi`, which reverts
  to 0, `kappa2[t] = Phi kappa2[t - 1] + e2[t]`. The expected log rates
  of two populations then keep the difference of their first period
  effects in the last year, and the difference of their second period
  effects decays to 0, the mean over the fitting years. The variance of
  the difference of the first period effects still grows with the
  horizon unless their innovations are perfectly correlated. With one
  population and `order = 1`, the draws equal those of `index = "rwd"`.

**Choices for `index = "coherent"` (package-specific; the paper is
silent on them).**

- `d` is the mean year-to-year change of the first period effects over
  the years and the populations, each population with the same weight
  (as in the estimation of `beta`). With `common_by`, each group has its
  own drift, from the strata of the group (the bottom strata of a
  hierarchy); the aggregates of a hierarchy take the drift of their
  group.

- Each diagonal element of `Phi` is the least-squares estimate through
  the origin for its population, bounded to `|phi| <= 0.99` as in
  [`li_lee()`](li_lee.md), because the paper's process is stationary; a
  warning of class `mt_warning_cae` names the strata where the bound
  binds. Orders 3 and higher follow the same process as order 2 (the
  paper has `order = 2`).

- One covariance matrix of the innovations of all populations and orders
  is estimated from the residuals (divided by the number of changes
  minus one), as for `"rwd"`. The paper gives each process a volatility
  matrix, does not say whether the innovations of the two processes are
  correlated, and reports correlations below 0.3 between them for its
  data (Section 6.1).

- Draws include the innovations and the uncertainty of `d` (normal, with
  the variance of the mean change), and keep `Phi` fixed.

**Estimation in Enchev, Kleinow and Cairns (2017).** The paper fits the
CAE model with `order = 2` by Poisson maximum likelihood (Section 4.3),
with age effects that sum to 1, period effects that sum to 0 over the
years for each population, and second period effects that sum to 0 over
the populations in each year (Table 4), a constraint that lowers the
likelihood (Section 4.2). `cae()` keeps the least-squares estimates of
Kleinow (2015) for both values of `index`: the period effects are
identified by the common principal components, and `index = "coherent"`
applies the paper's forecasting model to them. Rescaling a column of
`beta` leaves the forecasts unchanged. Rotating the age effects changes
them, so forecasts from the paper's estimates can differ. The paper's
data (men aged 60-89 in six countries, 1961-2010) are not in the Human
Mortality Database extracts of demography and StMoMo, and the paper
reports its forecasts as figures.

**Other package-specific choices.**

- The criterion is minimised by Jacobi rotations (Cardoso and Souloumiac
  1996), each of which minimises it in one plane, started from the
  eigenvectors of the sum of the `Q[i]`, until no rotation of more than
  1e-12 radians lowers the criterion by more than 1e-20 times the sum of
  squares of the matrices (at most 500 sweeps over all pairs of columns;
  a warning of class `mt_warning_cae` reports when this limit is
  reached). The paper uses the F-G algorithm (Flury and
  Constantine 1985) with the modification of Clarkson (1988). The
  criterion can have several local minima; the estimate is the one
  reached from this start.

- The columns of `beta` are ordered by the sum over populations of the
  diagonal elements of `t(beta) Q[i] beta`. The paper orders them by the
  diagonal of one population and notes that the order can differ between
  populations, with no general solution (Section 3.2).

- Zero deaths are replaced by 0.5 before taking logarithms.

- `common_by` (for example `"sex"`) gives each group of strata its own
  common age effects, estimated from the strata of the group. The paper
  names this extension and leaves it for further research (Section 3.1).
  The innovations of the period effects of all strata have one
  covariance matrix, with correlations across groups.

- If the data carry a hierarchy (see
  [`mt_aggregate_strata()`](mt_aggregate_strata.md)), `beta` is
  estimated from the bottom strata only, and the aggregates get their
  period effects by projection on it, as the other strata (as in
  [`li_lee()`](li_lee.md)).

**Choice of `order`.** The fit holds the mean squared error of the
centred log rates (eq. 4) and the approximate BIC of eq. (5),
`N log(MSE) + k log(N)` with `N` the number of rates and
`k = order (K T + n)` parameters for `K` populations, `T` years and `n`
ages (one `beta` per group of `common_by`);
[`mt_parameters()`](mt_parameters.md) reports them (rows `mse` and
`bic`). In the paper's application (men in ten countries), `order = 2`
had the lowest BIC in both age ranges (Table 2). Compare orders also in
a backtest.

## References

Kleinow T (2015). A common age effect model for the mortality of
multiple populations. *Insurance: Mathematics and Economics* 63:
147-152.
[doi:10.1016/j.insmatheco.2015.03.023](https://doi.org/10.1016/j.insmatheco.2015.03.023)

Flury BN (1984). Common principal components in k groups. *Journal of
the American Statistical Association* 79(388): 892-898.
[doi:10.1080/01621459.1984.10477108](https://doi.org/10.1080/01621459.1984.10477108)

Flury BN, Constantine G (1985). Algorithm AS 211: the F-G
diagonalization algorithm. *Applied Statistics* 34(2): 177-183.
[doi:10.2307/2347375](https://doi.org/10.2307/2347375)

Clarkson DB (1988). Remark AS R71: a remark on algorithm AS 211. The F-G
diagonalization algorithm. *Applied Statistics* 37(1): 147-151.
[doi:10.2307/2347513](https://doi.org/10.2307/2347513)

Cardoso J-F, Souloumiac A (1996). Jacobi angles for simultaneous
diagonalization. *SIAM Journal on Matrix Analysis and Applications*
17(1): 161-164.
[doi:10.1137/S0895479893259546](https://doi.org/10.1137/S0895479893259546)

Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
*Journal of the American Statistical Association* 87(419): 659-671.

Enchev V, Kleinow T, Cairns AJG (2017). Multi-population mortality
models: fitting, forecasting and comparisons. *Scandinavian Actuarial
Journal* 2017(4): 319-342.
[doi:10.1080/03461238.2015.1133450](https://doi.org/10.1080/03461238.2015.1133450)

Villegas AM, Kaishev VK, Millossovich P (2018). StMoMo: an R package for
stochastic mortality modeling. *Journal of Statistical Software* 84(3):
1-38. [doi:10.18637/jss.v084.i03](https://doi.org/10.18637/jss.v084.i03)

## See also

[`li_lee()`](li_lee.md) and [`coherent_fdm()`](coherent_fdm.md) for
coherent joint models; [`coupled()`](coupled.md) to give separate
forecasts the dependence of a joint model.

## Examples

``` r
# more years (30) than age groups (9): ages 50 and over in 5-year groups
x <- mt_subset(mt_regroup(mt_example(), age = "deaths"), ages = c(50, Inf))
fit <- mt_fit(cae(order = 2), x)
p <- mt_parameters(fit)
p[parameter %in% c("mse", "bic")]
#>    stratum parameter  index         value
#>     <char>    <char> <char>         <num>
#> 1:   model       mse         3.047876e-04
#> 2:   model       bic        -3.503547e+03
fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
# common drift and mean-reverting second period effects
x2 <- mt_subset(mt_example_regions(years = 1990:2019, regions = 2), ages = c(60, Inf))
fit2 <- mt_fit(cae(order = 2, index = "coherent"), x2)
fit2$fits$common_drift
#> [1] -0.03004355
fc2 <- mt_forecast(fit2, h = 10, n_draws = 100, seed = 1)
```
