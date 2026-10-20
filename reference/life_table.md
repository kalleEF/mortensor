# Life tables from death rates on any age axis

Period life table for single or grouped ages with an open last interval.
Interval widths come from the age axis, so abridged and single-year
schedules are handled by the same code. The relations are the standard
ones of Preston, Heuveline and Guillot (2001):

## Usage

``` r
life_table(
  mx,
  age = NULL,
  ax = NULL,
  open = c("reciprocal", "gompertz", "kannisto", "makeham", "weibull", "beard",
    "logquad", "logistic", "perks", "lynch_brown", "aic"),
  k = 10,
  radix = 1e+05,
  omega = 110,
  exposure = NULL
)
```

## Arguments

- mx:

  Death rates, one per interval of `age`.

- age:

  An [mt_axis](mt_axis.md) (default: single ages from 0 with an open
  last interval).

- ax:

  Optional average years lived in the interval by those who die.

- open:

  `"reciprocal"`, a hazard law (`"gompertz"`, `"kannisto"`, `"makeham"`,
  `"weibull"`, `"beard"`, `"logquad"`, `"logistic"`, `"perks"`,
  `"lynch_brown"`) or `"aic"`.

- k:

  Number of ages used for the fit of a hazard law.

- radix:

  Starting number of survivors.

- omega:

  Age up to which a fitted hazard is applied (default 110, the default
  `omega` of [`mt_data()`](mt_data.md)); ignored with
  `open = "reciprocal"`.

- exposure:

  Optional person-years of the age intervals, which weight the
  maximum-likelihood fit of a hazard law with the observed numbers at
  risk; by default the numbers of the life table. Needed for
  `open = "aic"`.

## Value

A `data.table` with `age_lo`, `age_hi`, `n`, `mx`, `ax`, `qx`, `lx`,
`dx`, `Lx`, `Tx`, `ex`. With a Gompertz or Kannisto closure, the
attribute `"tail"` holds the fitted `a` and `b`, the ages of the fit and
whether it converged. With another law it holds `alpha`, `beta`, `gamma`
and `delta` (`NA` where the law has fewer parameters), the
log-likelihood `loglik` and `aic` (with `exposure` only), the ages and
`converged`; with `"aic"` one such row per law, with `delta_aic` and the
`selected` law.

## Details

- `q_x = n m / (1 + (n - a) m)` for a closed interval of width `n`;
  where this exceeds 1 (rates above about `2 / n`, far beyond observed
  values), the constant-hazard value `1 - exp(-n m)` is used, so that
  survivors stay positive;

- `a_x` defaults to `n / 2` for every closed interval, the simple
  convention of deaths spread evenly over the interval (supply `ax` to
  override). At age 0, `n / 2` overstates the average time lived by
  infants who die; the Human Mortality Database uses the age-0 values of
  Andreev and Kingkade (2015), which depend on sex and on the infant
  death rate and are not implemented here. The effect on e0 at low
  infant mortality is below 0.01 years;

- open interval: `"reciprocal"` (default) uses `L = l / m` (Preston,
  Heuveline and Guillot 2001). `"gompertz"` and `"kannisto"` close the
  table with the Gompertz hazard `mu(x) = a exp(b x)` or the Kannisto
  hazard `mu(x) = a exp(b x) / (1 + a exp(b x))`, whose logit is linear
  in age and which levels off at old ages (Thatcher, Kannisto and Vaupel
  1998, chapter 2). They are fitted as in Thatcher, Kannisto and Vaupel
  (1998, chapter 6), by maximum likelihood with binomial deaths: of the
  `N` persons who reach age `x`, `d` die before `x + 1` with probability
  `q = 1 - exp(-H)`, where `H` is the integral of the hazard over the
  year. The fit uses the last `k` closed single-year ages. `N` and `d`
  are the survivors and deaths of the life table itself, the numbers of
  a stationary population with these rates, unless `exposure` gives the
  observed person-years `E` (then `d = m E` and `N = E + (1 - a) d`).
  Thatcher et al. fit observed numbers, at ages 80 to 98; the life-table
  numbers are a package-specific choice for schedules without counts
  (reference unknown). For the cohorts born in 1871-1880 the life-table
  numbers are the observed ones, and the fit reproduces the published
  parameters (their Appendix D) from their probabilities of dying (Table
  6.1). Earlier versions of mortensor fitted least squares to the log or
  logit rates. The fitted hazard is applied year by year from the start
  of the open interval to `omega` (survival `exp(-H)` and person-years
  `l (1 - exp(-H)) / H` in each year), and the survivors at `omega` are
  closed with `L = l / H` at the integrated hazard of the following
  year;

- `"makeham"`, `"weibull"`, `"beard"`, `"logquad"`, `"logistic"`,
  `"perks"` and `"lynch_brown"` close the table with the other hazard
  laws that Feehan (2018, Table 1 and Appendix E) compares, for age `x`
  in years:

  - Makeham: `gamma + alpha exp(beta x)`;

  - Weibull: `alpha x^(beta - 1)` with `x` the age, as Thatcher et
    al. (1998) fit `a x^b` (Feehan's fits count age 80 as `z = 1`,
    Appendix E.9, which changes this law only);

  - Beard (gamma-Gompertz):
    `alpha exp(beta x) / (1 + delta exp(beta x))`;

  - log-quadratic: `exp(alpha + beta x + gamma x^2)` (the quadratic
    model of Coale and Kisker in Thatcher et al. 1998);

  - logistic: `gamma + alpha exp(beta x) / (1 + delta exp(beta x))` (the
    logistic model of Thatcher et al. 1998 in another parameterisation);

  - Perks (gamma-Makeham):
    `(gamma + alpha exp(beta x)) / (1 + delta exp(beta x))`;

  - Lynch-Brown: `alpha + beta atan(gamma (x - delta))`.

  The parameters are positive, except `beta` of Makeham, all parameters
  of the log-quadratic law and `alpha` and `delta` of Lynch-Brown, which
  are real (Feehan's Appendix E). They are fitted by the same binomial
  maximum likelihood as the Kannisto closure, with the integrated
  hazards of Feehan's Appendix E in closed form (for the log-quadratic
  law by Gauss-Legendre quadrature with 20 nodes per year, as its closed
  form needs the error function of an imaginary argument when
  `gamma > 0`). Feehan maximises with BFGS from one heuristic and ten
  random starts; the package takes Levenberg-Marquardt steps on the
  Fisher information from several deterministic starts and keeps the
  highest likelihood (package-specific, reference unknown);

- `"aic"` fits all nine laws and closes each schedule with the law of
  the lowest Akaike information criterion, `AIC = -2 loglik + 2 k` with
  `k` parameters, the criterion by which Feehan (2018, section 4.2)
  compares the models. It needs `exposure`: with the numbers of the life
  table (radix 1) the likelihood has no sample size;

- a rate of 0 in the open interval (no deaths, as in small populations)
  would give an infinite life expectancy: the last positive rate of the
  closed intervals is used for the open interval instead (column `mx`
  keeps the given rate). Only a schedule without any positive rate gives
  an infinite value.

## References

Andreev EM, Kingkade WW (2015). Average age at death in infancy and
infant mortality level: reconsidering the Coale-Demeny formulas at
current levels of low mortality. *Demographic Research* 33: 363-390.

Feehan DM (2018). Separating the signal from the noise: evidence for
deceleration in old-age death rates. *Demography* 55(6): 2025-2044.
[doi:10.1007/s13524-018-0728-x](https://doi.org/10.1007/s13524-018-0728-x)

Preston SH, Heuveline P, Guillot M (2001). *Demography: Measuring and
Modeling Population Processes*. Blackwell.

Thatcher AR, Kannisto V, Vaupel JW (1998). *The Force of Mortality at
Ages 80 to 120*. Odense University Press.

## See also

[`mt_open_check()`](mt_open_check.md) to compare the laws out of sample.

## Examples

``` r
m <- 3e-4 + 2e-5 * exp(0.1 * (30:100))
lt <- life_table(m, single_ages(30, 100))
lt[age_lo %in% c(30, 65, 100)]
#>    age_lo age_hi     n           mx       ax           qx         lx         dx
#>     <num>  <num> <num>        <num>    <num>        <num>      <num>      <num>
#> 1:     30     31     1 0.0007017107 0.500000 0.0007014646 100000.000   70.14646
#> 2:     65     66     1 0.0136028327 0.500000 0.0135109391  87531.735 1182.63594
#> 3:    100    Inf   Inf 0.4408293159 2.268452 1.0000000000   1460.141 1460.14118
#>          Lx         Tx        ex
#>       <num>      <num>     <num>
#> 1: 99964.93 4973452.39 49.734524
#> 2: 86940.42 1596839.74 18.242981
#> 3:  3312.26    3312.26  2.268452
# Kannisto closure fitted by maximum likelihood to ages 90-99
lk <- life_table(m, single_ages(30, 100), open = "kannisto")
attr(lk, "tail")
#>         law            a         b age_from age_to converged
#>      <char>        <num>     <num>    <num>  <num>    <lgcl>
#> 1: kannisto 1.045127e-06 0.1338894       90    100      TRUE
# the log-quadratic law of Feehan (2018)
attr(life_table(m, single_ages(30, 100), open = "logquad"), "tail")
#>        law     alpha       beta        gamma delta loglik   aic age_from age_to
#>     <char>     <num>      <num>        <num> <num>  <num> <num>    <num>  <num>
#> 1: logquad -9.903961 0.07871943 0.0001176787    NA     NA    NA       90    100
#>    converged
#>       <lgcl>
#> 1:      TRUE
```
