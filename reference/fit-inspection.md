# Fitted rates, residuals and summaries of fitted models

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
in-sample fitted death rates as an array `[age, period, stratum]` (`NA`
for models without fitted values, such as the plain random walk).
[`residuals()`](https://rdrr.io/r/stats/residuals.html) returns deviance
residuals (default) or Pearson residuals of the death counts, assuming
Poisson deaths given the fitted rates, as a long table; values far
beyond 2 in absolute value point to structure that the model misses (for
example cohort effects) or to overdispersion.
[`summary()`](https://rdrr.io/r/base/summary.html) describes the model,
its sources of uncertainty and the fit, and the fitting period where
`lee_carter(choose_period = TRUE)` chose a shorter one (as
[`print()`](https://rdrr.io/r/base/print.html) does);
[`mt_plot_residuals()`](mt_plot_residuals.md) shows the residuals as an
age-period heat map.

## Usage

``` r
# S3 method for class 'mt_fit'
fitted(object, ...)

# S3 method for class 'mt_fit'
residuals(object, type = c("deviance", "pearson"), ...)

# S3 method for class 'mt_fit'
summary(object, ...)
```

## Arguments

- object:

  An `mt_fit` object.

- ...:

  Not used.

- type:

  `"deviance"` or `"pearson"`.

## Value

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html): an array;
[`residuals()`](https://rdrr.io/r/stats/residuals.html): a `data.table`
with `stratum`, `age`, `period`, `deaths`, `fitted_deaths`, `residual`;
[`summary()`](https://rdrr.io/r/base/summary.html): an object of class
`summary.mt_fit`, printed.

## Examples

``` r
x <- mt_example(deaths = "single")
fit <- mt_fit(lee_carter(), x)
summary(fit)
#> Lee-Carter
#>   data: 2 strata; ages 30 to 100+; periods 1990-2019
#>   model and uncertainty:
#>     log m(x,t) = a(x) + b(x) k(t), estimated by Poisson maximum likelihood
#>     (Brouhns et al. 2002). k(t) follows a random walk with drift; draws include
#>     its innovations and the uncertainty of the drift. a(x) and b(x) are fixed
#>     (parameter uncertainty with bootstrapped()). Predictive draws add Poisson
#>     variation of the deaths given the rates. Jump-off: fitted rates.
#>   fit: Poisson deviance 4105.7 on 4260 cells (0.96 per cell); 4.1% of deviance residuals beyond +/-2
head(residuals(fit))
#>    stratum    age period deaths fitted_deaths   residual
#>     <char> <char>  <num>  <num>         <num>      <num>
#> 1:  female     30   1990     84      79.27054  0.5260417
#> 2:  female     31   1990    100      88.53964  1.1929981
#> 3:  female     32   1990     93      91.57006  0.1490452
#> 4:  female     33   1990    104      99.82076  0.4154301
#> 5:  female     34   1990    105     109.53592 -0.4364418
#> 6:  female     35   1990     96     107.46697 -1.1267445
```
