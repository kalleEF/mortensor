# Aggregation operator from grid cells to observed cells

A sparse matrix with one row per observed cell of one stratum and one
column per grid cell (age varies fastest), with entry 1 where the grid
cell lies inside the observed cell. It links latent quantities on the
grid to what is observed, for example in the penalized composite link
model.

## Usage

``` r
mt_operator(data, measure, stratum)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- measure:

  Measure name.

- stratum:

  Stratum id.

## Value

A `dgCMatrix`. Attribute `rows` gives the observed cells.

## Examples

``` r
x <- mt_example()
C <- mt_operator(x, "deaths", 1)
dim(C)
#> [1]  390 2130
```
