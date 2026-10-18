# Choose the random-walk orders of a Bayesian APC model by DIC

Fits [`apc_bayes()`](apc_bayes.md) with every combination of first- and
second-order random walks for age, period and cohort (and, optionally,
with and without overdispersion) and ranks the variants by the deviance
information criterion (Spiegelhalter et al. 2002) that bamp reports,
summed over strata. DIC measures in-sample fit; to choose by forecast
accuracy, compare the candidates with [`mt_backtest()`](mt_backtest.md)
instead. Selection on all data followed by a backtest of the selected
model reuses the test years once.

## Usage

``` r
apc_select(
  data,
  prep = NULL,
  orders = c("rw1", "rw2"),
  overdispersion = FALSE,
  mcmc = mcmc_control("fast"),
  ...
)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- prep:

  Optional [`mt_prep()`](mt_prep.md) pipeline.

- orders:

  Random-walk orders to try.

- overdispersion:

  Values of `overdispersion` to try.

- mcmc:

  MCMC settings, see [`mcmc_control()`](mcmc_control.md).

- ...:

  Further arguments of [`apc_bayes()`](apc_bayes.md).

## Value

A list of class `mt_apc_selection` with `best` (the chosen
specification) and `table` (DIC, effective number of parameters `pd` and
deviance at the posterior mean per candidate, sorted by DIC; DIC =
deviance at the mean + 2 pd).

## References

Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A (2002). Bayesian
measures of model complexity and fit. *Journal of the Royal Statistical
Society B* 64(4): 583-639.

## Examples

``` r
# \donttest{
if (requireNamespace("bamp", quietly = TRUE)) {
  x <- mt_example(years = 2005:2019)
  sel <- apc_select(x, prep = mt_prep(step_ungroup()), orders = "rw1", mcmc = mcmc_control("test"))
  sel$table
}
#>       age period cohort overdispersion     dic       pd deviance_at_mean
#>    <char> <char> <char>         <lgcl>   <num>    <num>            <num>
#> 1:    rw1    rw1    rw1          FALSE 1029.48 131.4061         766.6683
#>    converged
#>       <lgcl>
#> 1:      TRUE
# }
```
