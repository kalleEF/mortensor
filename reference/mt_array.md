# Arrays on the modelling grid

Returns a measure as an array `[age, period, stratum]` on the grid.
Cells observed at grid resolution are filled; cells that are not
observed are `NA`. Observed cells that are coarser than the grid raise
an error unless `coarse = "na"`.

## Usage

``` r
mt_array(data, measure = "deaths", coarse = c("error", "na"))
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- measure:

  Measure name.

- coarse:

  `"error"` or `"na"`.

## Value

A numeric array with dimnames.

## Examples

``` r
x <- mt_example()
dim(mt_array(x, "exposure"))
#> [1] 71 30  2
```
