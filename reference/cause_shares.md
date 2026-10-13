# Coherent multi-cause forecasts with cause shares

Forecasts all-cause mortality with any model (`total`) and the shares of
the causes in all-cause mortality with a compositional tensor model, so
that cause-specific forecasts add up to the all-cause forecast in every
draw (a Tucker4 decomposition of isometric log-ratio coordinates; an
experimental method without a published application to mortality, so
compare it with [`reconciled()`](reconciled.md) separate models in a
backtest):

## Usage

``` r
cause_shares(
  total = lee_carter(),
  var = "cause",
  total_level = "all",
  ranks = c(3, 3, 2, 2),
  residual = "auto",
  residual_noise = TRUE
)
```

## Arguments

- total:

  Model specification for all-cause mortality.

- var:

  Stratum variable with the causes.

- total_level:

  Level of `var` for all causes.

- ranks:

  Tucker ranks for age, year, group and ILR coordinate.

- residual:

  Add a residual cause so that listed causes may cover less than all
  deaths: `TRUE`, `FALSE`, or `"auto"` (added when the listed causes
  leave more than 0.1% of all-cause mortality unexplained in the median
  cell; a residual that is zero almost everywhere would only add a
  degenerate coordinate).

- residual_noise:

  Add in-sample residual variation of the ILR coordinates to the draws.

## Value

A model specification.

## Details

1.  shares of each cause (and of a residual cause, the part of all-cause
    mortality not covered by the listed causes) per age, year and group;

2.  isometric log-ratio (ILR) coordinates of the shares (orthonormal
    Helmert basis);

3.  Tucker decomposition of the array `[age, year, group, coordinate]`
    (HOOI) with ranks `ranks`; the year factors are forecast as a
    multivariate random walk with drift (bootstrap of joint changes) and
    the in-sample residual variation is added to the draws;

4.  inverse ILR gives share draws, which multiply the all-cause draws.

Shares of rare causes are noisy; smoothing the rates first
(`prep = mt_prep(step_smooth())`) helps. In backtests on simulated cause
data ([`mt_example_causes()`](mt_example_causes.md)), the 95% coverage
was 0.99 without smoothing and 0.94 to 0.97 with it.

The data need a stratum variable `var` with the all-cause level
`total_level` and the cause levels; other stratum variables (for example
sex) define the groups.

## References

Egozcue JJ, Pawlowsky-Glahn V, Mateu-Figueras G, Barcelo-Vidal C (2003).
Isometric logratio transformations for compositional data analysis.
*Mathematical Geology* 35(3): 279-300.

## Examples

``` r
x <- mt_example_causes(years = 2010:2019)
fit <- mt_fit(cause_shares(total = lee_carter(), ranks = c(2, 2, 2, 2)), x)
fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
```
