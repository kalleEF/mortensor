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
  life-table deaths), the loadings `b`, the index `k` with the order
  (`order_k`), coefficients (`arima_k`) and innovation standard
  deviation (`sd_k`) of its ARIMA model, the share of variation
  explained by the component, the parameters of the open age interval
  (`a_open`, `b_open`) or of the Kannisto model of each year
  (`kannisto_a`, `kannisto_b`); the coherent model gives per group of
  `common_by` the common `a`, `B`, `K` (with `order_K`, `arima_K`,
  `sd_K`) and per stratum `a`, `b`, `k` and the ARIMA model of `k`.
  [`coda_causes()`](coda_causes.md): per group, as [`coda()`](coda.md),
  with parts by age and cause (`b1`, `k1`, `b2`, `k2`, ... with several
  components);

- [`rwd()`](rwd.md): drift and innovation standard deviation per age
  (and, with `observation_noise = TRUE`, the variances of level and
  noise);

- [`apc_bayes()`](apc_bayes.md): posterior medians and 90% intervals of
  the age, period and cohort effects (logit scale, bamp's
  identification), and the Gelman-Rubin statistics (`rhat`, values below
  1.1 indicate convergence) of the identified quantities;

- [`p_splines()`](p_splines.md), [`cp_splines()`](cp_splines.md): the
  smoothing parameters over age (`lambda_age`) and years
  (`lambda_year`), the effective dimension (`ed`), the criterion (`bic`
  or `aic`) and the dispersion (deviance divided by the number of cells
  minus the effective dimension); for [`cp_splines()`](cp_splines.md)
  also the bounds of the rates of ageing and of change per age
  (`ageing_lower`, `ageing_upper`, `change_lower`, `change_upper`).

- [`bohk_rau()`](bohk_rau.md): posterior medians and 90% intervals of
  the curve parameters `b1` and `b2` by age and of the variance and
  hyperparameters, the least-squares values `theta1`, `theta2`
  (exponential core model) and the Gelman-Rubin statistics (`rhat`);

- [`nordpred_apc()`](nordpred_apc.md): estimates and 90% Wald intervals
  of the age, period and cohort effects (cohorts by central year of
  birth) and of the drift, the projected drift, the number of periods
  fitted, the choice of the recent trend and the p-values of the
  goodness of fit and of the test for the recent trend.

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
