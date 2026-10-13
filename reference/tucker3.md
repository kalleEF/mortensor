# Tucker3 tensor model across strata

Joint model of the array of log rates `[age, period, stratum]` (or of
centred log-ratio transformed life-table deaths, `transform = "coda"`)
through a Tucker3 decomposition `X = G x1 A x2 B x3 C` with ranks
`ranks = c(age, period, stratum)`, estimated by higher-order orthogonal
iteration (HOOI; De Lathauwer, De Moor and Vandewalle 2000). The period
factors (columns of `B`) are forecast as a multivariate random walk with
drift (bootstrap of joint year-to-year changes, which keeps their
correlation) or with ARIMA/ETS models per column (`ts = "arima"`,
`"ets"`, package forecast).

## Usage

``` r
tucker3(
  ranks = c(3, 3, 2),
  transform = c("log", "coda"),
  ts = c("rwd", "arima", "ets"),
  residual_noise = TRUE,
  maxit = 200,
  tol = 1e-08
)
```

## Arguments

- ranks:

  Tucker ranks for age, period and stratum, or `"bic"` / `"aic"` to
  choose them by an information criterion on the training data
  (`n log(RSS / n)` plus the penalty times the number of parameters of
  the factor matrices and the core; ranks up to 6 for age and period and
  up to 4 for strata). For selection by forecast accuracy use
  [`tuned()`](mt_tune.md).

- transform:

  `"log"` or `"coda"`.

- ts:

  Forecast method for the period factors: `"rwd"`, `"arima"`, `"ets"`.

- residual_noise:

  Add the in-sample residual variation (per age and stratum, normal on
  the transformed scale) to the draws. It represents structure that the
  low-rank model does not capture. On the calibration benchmark with
  independent period trends in 10 strata, 95% predictive intervals
  covered 0.95 on average with it and 0.87 without it.

- maxit, tol:

  Settings of the HOOI iterations.

## Value

A model specification.

## Details

The model is experimental: tensor decompositions have been proposed for
mortality of several populations, but this specification has no
published evaluation; compare it with established models such as
[`li_lee()`](li_lee.md) in a backtest before using it.

With `transform = "coda"`, life-table deaths (from age `x0`, the first
grid age) form a composition over the age cells including the open
interval; the open interval's death rate is forecast as an auxiliary
random walk with drift on its log, so the open age group is kept.

Ranks are reduced to the array dimensions where needed, so the model
also runs with one stratum.

## References

De Lathauwer L, De Moor B, Vandewalle J (2000). On the best rank-1 and
rank-(R1, R2, ..., RN) approximation of higher-order tensors. *SIAM
Journal on Matrix Analysis and Applications* 21(4): 1324-1342.
[doi:10.1137/S0895479898346995](https://doi.org/10.1137/S0895479898346995)

Kolda TG, Bader BW (2009). Tensor decompositions and applications. *SIAM
Review* 51(3): 455-500.
[doi:10.1137/07070111X](https://doi.org/10.1137/07070111X)

## Examples

``` r
x <- mt_example(deaths = "single")
fit <- mt_fit(tucker3(ranks = c(3, 3, 2)), x)
fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
```
