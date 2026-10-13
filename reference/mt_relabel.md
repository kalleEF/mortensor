# Rename levels of a stratum variable

Gives levels of a stratum variable new labels, for example sex codes
`"1"` and `"2"` as `"male"` and `"female"`. Levels mapped to the same
new label are combined (see [`mt_group_levels()`](mt_group_levels.md)).

## Usage

``` r
mt_relabel(data, var, labels, additive = c("deaths", "rates"))
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- var:

  Stratum variable.

- labels:

  Named character vector, `c(old = "new")`.

- additive:

  `"deaths"` (exposures are summed) or `"rates"` (exposures are shared
  by the levels and kept).

## Value

The `mt_data` object with the new labels.

## Examples

``` r
d <- data.frame(Sex = rep(1:2, each = 2), year = 2020:2021, age = "90+", deaths = 5:8,
                exposure = 100)
x <- mt_data(deaths = d, exposure = d, by = "Sex")
mt_strata(mt_relabel(x, "Sex", c("1" = "male", "2" = "female")))
#>    stratum    Sex  label
#>      <int> <char> <char>
#> 1:       1   male   male
#> 2:       2 female female
```
