# Common age effect model for several populations

Joint model of the strata as populations that share their age effects
and keep their own period effects: the common age effect (CAE) model of
Kleinow (2015). Status: published, with the package-specific forecast
and choices listed below; section and equation numbers refer to that
paper.

## Usage

``` r
cae(order = 1, common_by = NULL)
```

## Arguments

- order:

  Number of common age effects (`p` in the paper).

- common_by:

  Stratum variable(s) defining groups of strata with their own common
  age effects (package-specific, see Details).

## Value

A model specification. The fit (`fit$fits`) holds the common age effects
`beta` (one matrix per group), the centring `a` `[age, stratum]`, the
period effects `kappa` `[order, year, stratum]`, the drift `drift` and
innovation covariance `sigma` of the stacked period effects (order
within stratum), `mse` and `bic`.

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

**Forecasts (package-specific).** Kleinow (2015) does not forecast. Here
the period effects of all populations and orders follow one multivariate
random walk with drift, `kappa[t] = kappa[t - 1] + mu + e[t]`, with
normal innovations whose covariance matrix is estimated from the
year-to-year changes (divided by the number of changes minus one), so
the innovations are correlated between populations and between orders.
This is the random walk with drift of Lee and Carter (1992) for every
index, and the multivariate random walk with drift that StMoMo uses for
several period indexes (Villegas, Kaishev and Millossovich 2018), here
over all populations together. Draws include the innovations and, as in
[`lee_carter()`](lee_carter.md), the uncertainty of the drift (normal,
with the innovation covariance divided by the number of changes; StMoMo
keeps the drift fixed). `a` and `beta` are fixed. Forecasts start from
the fitted rates of the last year ([`anchored()`](anchored.md) starts
them from the observed rates). Each population has its own drift, so the
forecasts of different populations can diverge; [`li_lee()`](li_lee.md)
and [`coherent_fdm()`](coherent_fdm.md) keep them together. Predictive
draws (`mt_forecast(predictive = TRUE)`) add Poisson variation of the
deaths given the rates and the exposures (the paper makes no assumption
about the distribution of the errors).

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
  The period effects of all strata still form one random walk, with
  innovations correlated across groups.

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
```
