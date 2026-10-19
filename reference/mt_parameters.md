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
  `b_ultimate`; with `family = "nb"` also the dispersion `nb_dispersion`
  (variance `mu + phi mu^2`);

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

- [`coda_multiway()`](coda_multiway.md): the common age patterns `b`,
  time indexes `k` (with `order_k`, `arima_k`, `sd_k`), population
  markers `gamma` and core `g` (index `p,q,r`), numbered with several
  components, and the share of variation explained; per stratum the
  centre `a`, the parameters of the open age interval or of the Kannisto
  model of each year;

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

- [`saint()`](saint.md): per trend (`"trend <group>"`) the indexes
  `alpha`, `beta`, `kappa`, `zeta` by year, the frailty variance
  `sigma2`, the pseudo-log-likelihood `loglik`, the `drift` and
  covariance `cov` of `alpha` and `beta` (or, with `ecm`, the parameters
  of the error-correction model: `a11`, `a21`, `a32`, `a42`, `C_D`,
  `C_delta`, `lambda`, `Omega`) and the standard deviations `sd_kappa`,
  `sd_zeta`; per stratum with a spread the spread parameters `y1` ...
  `y5` by year and their VAR(1) (`ar` or `A`, and `Omega`);

- [`bohk_rau()`](bohk_rau.md): posterior medians and 90% intervals of
  the curve parameters `b1` and `b2` by age and of the variance and
  hyperparameters, the least-squares values `theta1`, `theta2`
  (exponential core model) and the Gelman-Rubin statistics (`rhat`);

- [`apci()`](apci.md): `alpha`, `beta` (by age), `kappa` (by year),
  `gamma` (by year of birth; 0 for the cohorts with four or fewer
  observations), the mean year `ybar`, the smoothing parameters
  (`lambda_alpha`, `lambda_beta`), the effective dimension `ed` and the
  `bic`, and the ARIMA models of `kappa` and `gamma` (`order_kappa`,
  `arima_kappa`, `sd_kappa`, `order_gamma`, `arima_gamma`, `sd_gamma`);

- [`nordpred_apc()`](nordpred_apc.md): estimates and 90% Wald intervals
  of the age, period and cohort effects (cohorts by central year of
  birth) and of the drift, the projected drift, the number of periods
  fitted, the choice of the recent trend and the p-values of the
  goodness of fit and of the test for the recent trend;

- [`cae()`](cae.md): the common age effects `beta1`, `beta2`, ... of
  each group (rows labelled `common <g> (<strata>)`), and per stratum
  the centring `a`, the period effects `kappa1`, `kappa2`, ..., their
  drift (`drift_kappa1`, ...), the standard deviation of their
  innovations (`sd_kappa1`, ...) and the correlation of their
  innovations with those of the same order of every stratum
  (`cor_kappa1`, ..., indexed by the other stratum); rows of stratum
  `"model"` give the mean squared error of the centred log rates (`mse`)
  and the approximate BIC of Kleinow (2015, eq. 5; `bic`).

Fits of wrapped packages (StMoMo, demography) keep the original fitted
object in `fit$fits[[i]]`, where the methods of those packages apply.

Wrappers return the parameters of the model they wrap
([`recalibrated()`](recalibrated.md) adds its factors;
[`coupled()`](coupled.md) adds those of the template, with column
`component` `"base"` or `"template"`). For an
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
