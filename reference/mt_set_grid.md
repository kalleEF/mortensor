# Set the modelling grid

The grid is the resolution on which models work. Every observed cell
must be a union of grid cells, so the grid is at least as fine as the
finest table. Without a grid set, the grid is the coarsest one on which
every table is exact: the age groups of the data when all tables have
the same groups, single years when one table (for example populations)
has single years. Coarser grids are obtained with
[`mt_regroup()`](mt_regroup.md), single years from grouped deaths with
[`mt_ungroup()`](mt_ungroup.md).

## Usage

``` r
mt_set_grid(data, age = NULL, period = NULL)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- age, period:

  [mt_axis](mt_axis.md) objects (`NULL` keeps the current axis).

## Value

The updated `mt_data` object.

## Examples

``` r
x <- mt_example()
x <- mt_set_grid(x, age = single_ages(30, 100))
mt_grid(x)$age
#> <mt_axis: age, 71 intervals> 30 31 32 33 ... 98 99 100+
```
