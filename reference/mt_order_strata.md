# Enforce an ordering between strata in forecast draws

For each draw, age and year where the rate of the `lower` stratum
exceeds that of the `upper` stratum (for example female above male
mortality), both rates are replaced so that their mean is kept and their
ratio equals the mean observed ratio of the last `years` training years
(bounded to 0.02-0.98), applied path by path and reported as the share
of adjusted cells.

## Usage

``` r
mt_order_strata(draws, var, lower, upper, years = 5)
```

## Arguments

- draws:

  An [mt_draws](mt_draws.md) object.

- var:

  Stratum variable that distinguishes the pair (for example `"sex"`).

- lower, upper:

  Values of `var` that should satisfy lower \<= upper.

- years:

  Years used for the observed ratio.

## Value

The adjusted `mt_draws`, with the share of adjusted cells in
`$coherence`.

## Details

This is a post-hoc adjustment of this package without a published source
(reference unknown): it changes the draws and leaves the fitted model as
it is, so the adjusted draws no longer follow the model's distribution.
Models that keep strata together by construction are
[`li_lee()`](li_lee.md), [`coherent_fdm()`](coherent_fdm.md) and
`coda(coherent = TRUE)` ([`coda()`](coda.md));
[`reconciled()`](reconciled.md) makes totals and their parts coherent.

## Examples

``` r
x <- mt_example(deaths = "single")
fc <- mt_forecast(mt_fit(rwd(), x), h = 10, n_draws = 50, seed = 1)
fc2 <- mt_order_strata(fc, "sex", lower = "female", upper = "male")
fc2$coherence
#>     lower  upper share_adjusted
#>    <char> <char>          <num>
#> 1: female   male     0.07346479
```
