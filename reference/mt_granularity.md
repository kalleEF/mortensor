# Describe the granularity of each table

Describe the granularity of each table

## Usage

``` r
mt_granularity(data)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

## Value

A `data.table` with one row per measure and stratum: number of cells,
age widths, start of the open interval, period widths and range.

## Examples

``` r
mt_granularity(mt_example())
#>     measure stratum  label cells                     ages open_from   periods
#>      <char>   <int> <char> <int>                   <char>    <char>    <char>
#> 1:   deaths       1 female   390 5-year groups 30-89, 90+        90 1990-2019
#> 2:   deaths       2   male   390 5-year groups 30-89, 90+        90 1990-2019
#> 3: exposure       1 female  2130 single years 30-99, 100+       100 1990-2019
#> 4: exposure       2   male  2130 single years 30-99, 100+       100 1990-2019
#>    period_width
#>          <char>
#> 1:            1
#> 2:            1
#> 3:            1
#> 4:            1
```
