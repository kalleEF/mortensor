# Combine levels of a stratum variable

Adds an aggregation layer by combining levels of one stratum variable
into groups: causes of death into cause groups, deprivation deciles into
quintiles, districts into regions, single nationalities into larger
groups. Death counts are summed. Exposures and populations are summed
when the variable divides the population (regions, deciles,
`additive = "deaths"`) and kept when every level shares the same
population (causes of death, `additive = "rates"`). A group with one
member renames the level.

## Usage

``` r
mt_group_levels(
  data,
  var,
  groups,
  other = NULL,
  total = NULL,
  keep = character(),
  additive = c("deaths", "rates"),
  name = NULL
)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- var:

  Stratum variable whose levels are combined.

- groups:

  Named list (the levels of `var` in each group), or a function that
  maps a vector of levels to group names (for example
  `function(d) paste0("Q", ceiling(as.numeric(d) / 2))` for deciles to
  quintiles).

- other:

  Name of a residual group, or `NULL` to keep unlisted levels as they
  are.

- total:

  Level of `var` that is the total over all levels (for example all
  causes), or `NULL`. It is kept, and with `other` the residual group is
  the total minus the other groups.

- keep:

  Levels kept as they are when `other` is given.

- additive:

  `"deaths"` (exposures are summed) or `"rates"` (exposures are shared
  by the levels and kept).

- name:

  New name of the stratum variable (for example `"imd_quintile"`), or
  `NULL` to keep `var`.

## Value

The `mt_data` object with the combined levels. Cells are summed where
every member level is observed; others are dropped with a warning.

## Details

A residual group ("all other causes") is defined in one of two ways.
With `other` alone it is the sum of the levels not listed in `groups`.
With `other` and `total` it is the total minus all other groups, so it
also holds deaths that no level covers (for example when the source
lists selected causes and all causes). Negative remainders (the listed
levels exceed the total) are set to zero with a warning.

To add totals and groups as extra strata while keeping the original
levels, use [`mt_aggregate_strata()`](mt_aggregate_strata.md) instead.

## Examples

``` r
x <- mt_example_causes(years = 2015:2019)
y <- mt_group_levels(x, "cause", list(noncancer = c("cvd", "other")), additive = "rates")
mt_strata(y)
#>    stratum    sex     cause            label
#>      <int> <char>    <char>           <char>
#> 1:       1 female noncancer female.noncancer
#> 2:       2 female       all       female.all
#> 3:       3 female    cancer    female.cancer
#> 4:       4   male noncancer   male.noncancer
#> 5:       5   male       all         male.all
#> 6:       6   male    cancer      male.cancer
# deciles to quintiles, with a new name for the aggregation layer
d <- expand.grid(year = 2020, age = "90+", decile = 1:10)
d$deaths <- 10 + d$decile
d$exposure <- 1000
z <- mt_data(deaths = d, exposure = d, by = "decile")
z <- mt_group_levels(z, "decile", function(l) paste0("Q", ceiling(as.numeric(l) / 2)),
                     name = "quintile")
mt_strata(z)
#>    stratum quintile  label
#>      <int>   <char> <char>
#> 1:       1       Q1     Q1
#> 2:       2       Q2     Q2
#> 3:       3       Q3     Q3
#> 4:       4       Q4     Q4
#> 5:       5       Q5     Q5
```
