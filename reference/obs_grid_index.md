# Map observed cells to grid indices

Map observed cells to grid indices

## Usage

``` r
obs_grid_index(obs, grid, measure = "observations", strict = TRUE)
```

## Value

A `data.table` with `a1`, `a2`, `p1`, `p2` (first and last grid cell
covered on each axis), aligned with the rows of `obs`. Cells outside the
grid get `NA`.
