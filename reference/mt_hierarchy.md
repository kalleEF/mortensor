# Hierarchies of strata

Describes how strata add up, for example women and men to both sexes,
regions to the country and deprivation quintiles to all quintiles.
Several variables give a crossed (grouped) structure: with `sex`,
`region` and `quintile` every combination of a level or a total is a
node, from the bottom series (one sex, one region, one quintile) to the
national total for both sexes. Intermediate groups of levels, such as
federal states grouped into East and West, are given with `groups`.

## Usage

``` r
mt_hierarchy(..., additive = c("deaths", "rates"))
```

## Arguments

- ...:

  One argument per aggregated stratum variable: the label of its total
  level (for example `region = "Germany"`), or a list with `total`,
  optionally `groups` (a named list of intermediate groups, each a
  vector of bottom levels) and `additive`.

- additive:

  Default additivity of the variables: `"deaths"` or `"rates"`.

## Value

An object of class `mt_hierarchy`.

## Details

Deaths add up over every hierarchy variable. Exposures add up over
variables that divide the population (sex, regions, deprivation
quintiles: `additive = "deaths"`, the default) and are shared by the
levels of variables that do not (causes of death: `additive = "rates"`,
whose rates add up to the all-cause rate). Variables of both kinds can
be crossed, for example regions and causes: deaths are summed over both,
exposures over regions only, and forecasts are reconciled on the deaths
scale with these exposures.

## See also

[`mt_aggregate_strata()`](mt_aggregate_strata.md),
[`mt_validate()`](mt_validate.md), [`mt_reconcile()`](mt_reconcile.md),
[`reconciled()`](reconciled.md).

## Examples

``` r
h <- mt_hierarchy(sex = "both", region = list(total = "DE", groups = list(East = c("BE", "SN"))),
                  quintile = "all")
h
#> <mt_hierarchy>
#>   sex: total "both"; deaths and exposures add up
#>   region: total "DE"; groups East; deaths and exposures add up
#>   quintile: total "all"; deaths and exposures add up
# regions crossed with causes of death, which share the exposure
mt_hierarchy(region = "DE", cause = list(total = "all", additive = "rates"))
#> <mt_hierarchy>
#>   region: total "DE"; deaths and exposures add up
#>   cause: total "all"; deaths add up, exposures are shared
```
