# Add aggregate strata to data

Sums deaths, exposures and populations of the bottom strata into every
node of the hierarchy (all combinations of levels, groups and totals of
the hierarchy variables). Cells are aggregated where all member strata
observe the same age and period interval; other cells are dropped with a
warning. Aggregate strata that are already in the data are kept as
observed unless `replace = TRUE`; use [`mt_validate()`](mt_validate.md)
with the hierarchy to check that they agree with the sums. Members of
`groups` that are not levels of the data give a warning; a group without
any is an error.

## Usage

``` r
mt_aggregate_strata(data, hierarchy, replace = FALSE)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object with (at least) the bottom strata.

- hierarchy:

  An [`mt_hierarchy()`](mt_hierarchy.md).

- replace:

  Replace aggregate strata that are already in the data.

## Value

An `mt_data` object with all nodes as strata; the hierarchy is stored in
`$hierarchy`.

## Examples

``` r
x <- mt_example_regions(years = 2010:2014)
xa <- mt_aggregate_strata(x, mt_hierarchy(sex = "both", region = "all", quintile = "all"))
nrow(mt_strata(xa))
#> [1] 90
```
