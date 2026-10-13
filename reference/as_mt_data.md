# Convert other mortality data formats

Methods for long data frames, `StMoMoData` objects from StMoMo (for
example
[`StMoMo::EWMaleData`](https://rdrr.io/pkg/StMoMo/man/EWMaleData.html),
Human Mortality Database data for England and Wales) and `demogdata`
objects from demography (for example
[`demography::fr.mort`](https://pkg.robjhyndman.com/demography/reference/fr.mort.html),
HMD data for France). For `demogdata`, deaths are reconstructed as rate
times exposure and rounded.

## Usage

``` r
as_mt_data(x, ...)

# S3 method for class 'data.frame'
as_mt_data(x, deaths = "deaths", exposure = "exposure", ...)

# S3 method for class 'StMoMoData'
as_mt_data(x, ...)

# S3 method for class 'demogdata'
as_mt_data(x, series = NULL, years = NULL, ...)
```

## Arguments

- x:

  Object to convert.

- ...:

  Passed to methods; for data frames, arguments of
  [`mt_data()`](mt_data.md).

- deaths, exposure:

  Column names of deaths and exposures in a data frame.

- series:

  Series of a `demogdata` object to keep (default all except `"total"`).

- years:

  Optional subset of years.

## Value

An [mt_data](mt_data.md) object.

## Examples

``` r
if (requireNamespace("StMoMo", quietly = TRUE)) {
  ew <- as_mt_data(StMoMo::EWMaleData)
  ew
}
#> <mt_data> 1 stratum by sex, omega = 110
#>   deaths     single years 0-100; periods 1961-2011 (width 1); 5151 cells
#>   exposure   single years 0-100; periods 1961-2011 (width 1); 5151 cells
#>   grid       <mt_axis: age, 101 intervals> 0 1 2 3 ... 98 99 100
#>              <mt_axis: period, 51 intervals> 1961 1962 1963 1964 ... 2009 2010 2011
```
