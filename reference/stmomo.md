# Generalised age-period-cohort models (StMoMo)

Wrapper around StMoMo (Villegas, Kaishev and Millossovich 2018):
age-period-cohort (`"apc"`), Cairns-Blake-Dowd (`"cbd"`, logit link),
Renshaw-Haberman (`"rh"`), M6 (`"m6"`, logit link) and Lee-Carter
(`"lc"`). Models are fitted per stratum with real ages (the lower bounds
of single-year intervals; models with a cohort effect need single years
of age), deaths and central exposures; for the logit-link models initial
exposures are approximated as central exposures plus half the deaths,
and forecast death probabilities are converted to central rates
(`m = -log(1 - q)`). Cohorts observed in fewer than `clip` cells get
zero weight; their fitted rates use the effect of the nearest estimated
cohort. Draws come from
[`simulate()`](https://rdrr.io/r/stats/simulate.html) of the fitted
model (period and cohort indexes simulated as multivariate random walks
with drift and ARIMA, respectively).

## Usage

``` r
stmomo(
  model = c("apc", "cbd", "rh", "m6", "lc"),
  clip = 3,
  gc_order = c(1, 1, 0)
)
```

## Arguments

- model:

  `"apc"`, `"cbd"`, `"rh"`, `"m6"` or `"lc"`.

- clip:

  Cohorts with fewer observations get zero weight.

- gc_order:

  ARIMA order for the cohort index.

## Value

A model specification.

## References

Villegas AM, Kaishev VK, Millossovich P (2018). StMoMo: an R package for
stochastic mortality modeling. *Journal of Statistical Software* 84(3):
1-38. [doi:10.18637/jss.v084.i03](https://doi.org/10.18637/jss.v084.i03)

## Examples

``` r
if (requireNamespace("StMoMo", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 1995:2019)
  fc <- mt_forecast(mt_fit(stmomo("cbd"), x), h = 5, n_draws = 50, seed = 1)
}
```
