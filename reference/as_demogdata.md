# Convert data to demography and StMoMo objects

Exports observed deaths and exposures for cross-checks with the original
packages: `as_demogdata()` returns a `demogdata` object of the
demography package with one series per stratum (rates and exposures as
population), and `as_stmomodata()` a `StMoMoData` object of the StMoMo
package for one stratum (central exposures). The data must be on the
grid (deaths and exposures in the same age groups and single years).
[`as_mt_data()`](as_mt_data.md) converts in the other direction.

## Usage

``` r
as_demogdata(data, name = "mortensor")

as_stmomodata(data, stratum = NULL, name = "mortensor")
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- name:

  Name of the population.

- stratum:

  Label of the stratum for `as_stmomodata()` (default: the first).

## Value

A `demogdata` or `StMoMoData` object.

## Examples

``` r
if (requireNamespace("demography", quietly = TRUE)) {
  x <- mt_example(deaths = "single", years = 2001:2010)
  dd <- as_demogdata(x)
  names(dd$rate)
}
#> Registered S3 methods overwritten by 'demography':
#>   method      from 
#>   print.lca   e1071
#>   summary.lca e1071
#> [1] "female" "male"  
```
