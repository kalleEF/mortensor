# Generalised age-period-cohort models (StMoMo)

Wrapper around StMoMo (Villegas, Kaishev and Millossovich 2018):
age-period-cohort (`"apc"`), Cairns-Blake-Dowd (`"cbd"`, logit link),
Renshaw-Haberman (`"rh"`), M6, M7 and M8 (`"m6"`, `"m7"`, `"m8"`, logit
link) and Lee-Carter (`"lc"`). Models are fitted per stratum with real
ages (the lower bounds of single-year intervals; models with a cohort
effect need single years of age), deaths and central exposures; for the
logit-link models initial exposures are approximated as central
exposures plus half the deaths, and forecast death probabilities are
converted to central rates (`m = -log(1 - q)`). Cohorts observed in
fewer than `clip` cells get zero weight; their fitted rates use the
effect of the nearest estimated cohort. Draws come from
[`simulate()`](https://rdrr.io/r/stats/simulate.html) of the fitted
model (period indexes simulated as a multivariate random walk with
drift, the cohort index as an ARIMA process).

## Usage

``` r
stmomo(
  model = c("apc", "cbd", "rh", "m6", "m7", "m8", "lc"),
  clip = 3,
  gc_order = NULL,
  xc = NULL
)
```

## Arguments

- model:

  `"apc"`, `"cbd"`, `"rh"`, `"m6"`, `"m7"`, `"m8"` or `"lc"`.

- clip:

  Cohorts with fewer observations get zero weight.

- gc_order:

  ARIMA order for the cohort index. Default: `c(2, 0, 0)` for M7,
  `c(1, 1, 0)` for the other cohort models.

- xc:

  Age at which the cohort effect of M8 vanishes (`model = "m8"` only).
  Cairns et al. (2009) estimate it; StMoMo's `m8()` needs a fixed value,
  and its example uses the highest age of the data (89 for ages 55 to
  89). Default (`NULL`): the highest age of the data of each fit.
  Forecasts can depend strongly on it; compare values in a backtest
  ([`mt_compare()`](mt_compare.md)) or choose one with
  [`tuned()`](mt_tune.md).

## Value

A model specification.

## Details

M6, M7 and M8 are the cohort extensions of the CBD model of Cairns et
al. (2009), with `x` the age, `xbar` the mean age of the data and `s2`
the mean of `(x - xbar)^2`:

- M6: `logit q(x,t) = k1(t) + (x - xbar) k2(t) + g(t - x)`.

- M7: M6 plus a quadratic age term `((x - xbar)^2 - s2) k3(t)`; the
  cohort effect has no linear or quadratic trend (constraints of
  StMoMo's `m7()`).

- M8: `logit q(x,t) = k1(t) + (x - xbar) k2(t) + (xc - x) g(t - x)`: the
  cohort effect fades with age and vanishes at age `xc`.

The CBD family is meant for older ages (StMoMo's examples use ages 55 to
89). For M7 the cohort index is forecast with ARIMA(2,0,0) with a
non-zero mean, as in the StMoMo vignette, since its constraints remove
any trend; the other cohort models use ARIMA(1,1,0) with drift.

## References

Cairns AJG, Blake D, Dowd K, Coughlan GD, Epstein D, Ong A, Balevich I
(2009). A quantitative comparison of stochastic mortality models using
data from England and Wales and the United States. *North American
Actuarial Journal* 13(1): 1-35.

Villegas AM, Kaishev VK, Millossovich P (2018). StMoMo: an R package for
stochastic mortality modeling. *Journal of Statistical Software* 84(3):
1-38. [doi:10.18637/jss.v084.i03](https://doi.org/10.18637/jss.v084.i03)

## Examples

``` r
if (requireNamespace("StMoMo", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 1995:2019)
  fc <- mt_forecast(mt_fit(stmomo("cbd"), x), h = 5, n_draws = 50, seed = 1)
  old <- mt_regroup(mt_set_grid(x, age = single_ages(55, 89, open = FALSE)),
                    age = single_ages(55, 89, open = FALSE))
  fc7 <- mt_forecast(mt_fit(stmomo("m7"), old), h = 5, n_draws = 50, seed = 1)
}
```
