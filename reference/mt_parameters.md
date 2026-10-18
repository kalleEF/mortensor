# Parameters of a fitted model

Returns the estimated parameters of a fit in one long table, for every
model type that has interpretable parameters:

## Usage

``` r
mt_parameters(fit)
```

## Arguments

- fit:

  An `mt_fit` object.

## Value

A `data.table` with `stratum`, `parameter`, `index` (age label, year or
position), `value` and, for Bayesian fits, `lower` and `upper`; for
ensembles also `component`, for fits by groups also `group`.

## Details

- [`lee_carter()`](lee_carter.md): `a` (age pattern), `b` (age pattern
  of change), `k` (period index), and the drift and innovation standard
  deviation of `k`; with `rotate = TRUE` also the ultimate pattern
  `b_ultimate`;

- [`li_lee()`](li_lee.md), [`product_ratio()`](product_ratio.md): the
  common factor of each group (`B`, `K`, its drift and standard
  deviation), and per stratum the deviations `a`, `b`, `k` and their
  AR(1) coefficient `phi` and standard deviation;

- [`coda()`](coda.md): the centre `a` (geometric mean distribution of
  life-table deaths), the loadings `b`, the index `k` with its drift and
  standard deviation, the share of variation explained by the first
  component and the parameters of the open age interval; the coherent
  model gives per group of `common_by` the common `a`, `B`, `K` and per
  stratum `a`, `b`, `k` and their AR(1) parameters.
  [`coda_causes()`](coda_causes.md): per group, as [`coda()`](coda.md),
  with parts by age and cause;

- [`rwd()`](rwd.md): drift and innovation standard deviation per age
  (and, with `observation_noise = TRUE`, the variances of level and
  noise);

- [`apc_bayes()`](apc_bayes.md): posterior medians and 90% intervals of
  the age, period and cohort effects (logit scale, bamp's
  identification), and the Gelman-Rubin statistics (`rhat`, values below
  1.1 indicate convergence) of the identified quantities.

Fits of wrapped packages (StMoMo, demography) keep the original fitted
object in `fit$fits[[i]]`, where the methods of those packages apply.

Wrappers return the parameters of the model they wrap
([`recalibrated()`](recalibrated.md) adds its factors). For an
[`ensemble()`](ensemble.md), the rows of each component carry its name
in column `component`, and rows with `parameter = "weight"` give the
weights; components without interpretable parameters (wrapped packages)
contribute their weight only. For [`stratified()`](stratified.md) and
`tuned(by = )`, the rows of each group carry its name in column `group`.

## See also

[`summary.mt_fit()`](fit-inspection.md),
[`residuals.mt_fit()`](fit-inspection.md).

## Examples

``` r
x <- mt_example(deaths = "single")
p <- mt_parameters(mt_fit(lee_carter(), x))
head(p)
#>    stratum parameter  index     value
#>     <char>    <char> <char>     <num>
#> 1:  female         a     30 -7.530513
#> 2:  female         a     31 -7.438225
#> 3:  female         a     32 -7.453512
#> 4:  female         a     33 -7.353122
#> 5:  female         a     34 -7.325712
#> 6:  female         a     35 -7.275708
```
