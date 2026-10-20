# Check old-age mortality models out of sample

Fits the hazard laws of [`life_table()`](life_table.md) to the deaths at
ages `from` to `cut - 1` and scores their extrapolation against the
observed deaths at ages `cut` and above, the design of Thatcher,
Kannisto and Vaupel (1998, chapter 6): they fit six models at ages 80-98
and compare the predicted probabilities of dying with the observed ones
at ages 99 and over. The predictions are based only on the data below
the cut, so they are independent of the observations they are compared
with.

## Usage

``` r
mt_open_check(
  data,
  cut = 99,
  from = 80,
  to = NULL,
  laws = c("gompertz", "kannisto", "makeham", "weibull", "beard", "logquad", "logistic",
    "perks", "lynch_brown"),
  by_age = FALSE
)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object with deaths and exposures in single
  years of age from `from` to `to`.

- cut:

  First age of the test ages; the laws are fitted to the ages `from` to
  `cut - 1` (Thatcher et al. 1998: 99, fitting ages 80-98).

- from:

  First age of the fit (default 80, the lowest age of the data of
  Thatcher et al.).

- to:

  Last test age (default: the last single year of age of the data).

- laws:

  Hazard laws to check, see [`life_table()`](life_table.md) (default:
  all nine models of Feehan 2018).

- by_age:

  Return one row per test age instead of the sums.

## Value

A `data.table` with the stratum columns, `period` and `law`, `n_par`
(number of parameters), `loglik`, `aic` and `delta_aic` of the fit, the
out-of-sample sums `nll` and `chisq` over `n_ages` test ages, `rank` (by
`nll` within each stratum and period, 1 = best) and `converged`. With
`by_age = TRUE`, one row per test age with `age`, `N`, `d`, the observed
`q_obs = d / N`, the predicted `q_fit`, `nll` and `chisq`. No
thresholds: the check shows which law extrapolated best in these data.

## Details

Every law is fitted by maximum likelihood with binomial deaths: of the
`N_x` persons who reach age `x`, `d_x` die before `x + 1` with
probability `q_x = 1 - exp(-H_x)`, where `H_x` is the integral of the
hazard over the year (Thatcher et al. 1998, chapter 6; Feehan 2018,
equation 4). The deaths are those of the data and `N = E + d / 2`, the
person-years `E` plus half the deaths, as for
[`life_table()`](life_table.md) with `exposure` (Thatcher et al. use the
numbers reaching each age directly). At the test ages, with the
predicted `q`, Thatcher et al. compare the models by

- `nll`: minus the binomial log-likelihood of the observed deaths,
  `-(d log q + (N - d) log(1 - q))`, without the binomial coefficients,
  which are the same for all models (their Table 6.3); the model with
  the lowest sum fits best on the likelihood criterion (`rank`);

- `chisq`: `(d - N q)^2 / (N q (1 - q))`, distributed like chi-squared
  with one degree of freedom per age for a correct model and accurate
  data (their Table 6.4), so that the sum over `n_ages` ages has
  expected value `n_ages`.

The columns `aic` and `delta_aic` give Akaike's information criterion of
the fit below the cut, `-2 loglik + 2 k` with `k` parameters, and its
difference to the lowest value of the schedule, as Feehan (2018, section
4.2) compares the models; without the binomial coefficients, which
cancel in the differences. Thatcher et al. pool the data of decades and
13 countries: pool periods with [`mt_regroup()`](mt_regroup.md)
(`period = `) before the check. Cells whose deaths exceed the numbers at
risk (`d > N`, possible for period data with tiny exposures) are left
out (package-specific safeguard).

## References

Feehan DM (2018). Separating the signal from the noise: evidence for
deceleration in old-age death rates. *Demography* 55(6): 2025-2044.
[doi:10.1007/s13524-018-0728-x](https://doi.org/10.1007/s13524-018-0728-x)

Thatcher AR, Kannisto V, Vaupel JW (1998). *The Force of Mortality at
Ages 80 to 120*. Odense University Press.

## See also

[`life_table()`](life_table.md) for the laws and their use as a closure
of the open age interval.

## Examples

``` r
# expected deaths of a cohort under a Kannisto hazard: the laws fitted at
# 80-94 are scored at ages 95-104
x <- 80:104
H <- log1p(5e-5 * exp(0.11 * (x + 1))) / 0.11 - log1p(5e-5 * exp(0.11 * x)) / 0.11
N <- 1e5 * exp(-cumsum(c(0, H[-length(H)])))
d <- round(N * (1 - exp(-H)))
sim <- mt_data(deaths = data.frame(age = x, year = 2000, deaths = d),
               exposure = data.frame(age = x, year = 2000, exposure = N - d / 2))
chk <- mt_open_check(sim, cut = 95, laws = c("gompertz", "kannisto", "logquad", "beard"))
chk[, list(law, nll, chisq, n_ages, rank)]
#>         law      nll    chisq n_ages  rank
#>      <char>    <num>    <num>  <num> <int>
#> 1: gompertz 212.7159 8.625589      9     4
#> 2: kannisto 209.1797 1.061330      9     1
#> 3:  logquad 209.1924 1.087074      9     3
#> 4:    beard 209.1810 1.063967      9     2
```
