# Different models for different groups of strata

Fits one model specification to each group of strata, defined by the
levels of a stratum variable, and combines their forecasts into one set
of draws. A typical use is a trusted all-cause forecast from one model
and cause-specific forecasts from a model that pools sparse causes,
reconciled so that the causes add up to all causes:

## Usage

``` r
stratified(var, ..., .default = NULL)
```

## Arguments

- var:

  Stratum variable that defines the groups.

- ...:

  Named model specifications: each name is a level of `var` whose strata
  the model is fitted to.

- .default:

  Model for the strata of all other levels (`NULL`: every level must be
  named).

## Value

A model specification.

## Details

`reconciled(stratified("cause", all = lee_carter(), .default = li_lee(common_by = c("sex", "cause"))), h, method = "proportional", fixed = all_cause_labels)`.

Joint models (for example [`li_lee()`](li_lee.md)) are joint within
their group only; the draws of different groups are independent, so sums
over strata of different groups have narrower intervals than they would
if the groups moved together. [`reconciled()`](reconciled.md) keeps the
draws as they are by default (`dependence = "none"`). Its
`dependence = "empirical"` reorders the paths of every stratum so that
their ranks follow the changes of a resampled training year: a
package-specific heuristic (see [`mt_reconcile()`](mt_reconcile.md))
that imposes a rank dependence taken from the training data without a
model of it, and that also replaces the dependence within joint groups.
Each group is fitted and forecast as if it were the whole data, so
backtests, [`anchored()`](anchored.md),
[`recalibrated()`](recalibrated.md) and [`reconciled()`](reconciled.md)
work as for any model.

## Examples

``` r
x <- mt_example(deaths = "single", years = 1990:2019)
spec <- stratified("sex", female = lee_carter(), male = rwd())
fc <- mt_forecast(mt_fit(spec, x), h = 5, n_draws = 50, seed = 1)
fc$model
#> [1] "By sex (female: Lee-Carter; male: Random walk with drift)"
```
