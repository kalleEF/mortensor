# Example stratified mortality data

Simulated mortality data with a realistic Gompertz-Makeham age
structure, sex differentials, and deprivation gradients. Designed for
demonstrating the `mortensor` smoothing, decomposition, and forecasting
pipeline.

## Usage

``` r
example_mortality
```

## Format

A `mort_data` object with components:

- mx:

  3D array `[13 x 20 x 10]` of age-specific mortality rates

- pop:

  3D array `[13 x 20 x 10]` of population exposure

- ages:

  Integer vector of 13 age-group start values (30, 35, ..., 90)

- years:

  Integer vector of 20 calendar years (2001–2020)

- groups:

  data.frame with 10 rows and columns `sex` (Male/Female) and `quintile`
  (Q1–Q5)

## Source

Synthetically generated. See `data-raw/create_example_data.R`.

## Details

The dataset contains 13 five-year age groups (30–34 to 90+), 20 years
(2001–2020), and 10 groups (2 sexes x 5 deprivation quintiles).
Mortality rates follow a Gompertz law with period improvement trends and
a COVID-like shock in 2020. Population counts are synthetic but
plausible.

## Examples

``` r
data(example_mortality)
print(example_mortality)
#> <mort_data>
#>   Cause:   All causes 
#>   Ages:    13 groups ( 30-34  ..  90+ )
#>   Years:   2001 - 2020  ( 20  years)
#>   Groups:  10  ( sex x quintile )
#>   Pop:     yes 
#>   Smoothed: FALSE 

# Compute life expectancy
le <- life_expectancy(example_mortality)
head(le)
#>     year   group       ex    sex quintile
#>    <int>  <char>    <num> <char>   <char>
#> 1:  2001 Male_Q1 56.73900   Male       Q1
#> 2:  2002 Male_Q1 56.88602   Male       Q1
#> 3:  2003 Male_Q1 57.47525   Male       Q1
#> 4:  2004 Male_Q1 57.79061   Male       Q1
#> 5:  2005 Male_Q1 57.46339   Male       Q1
#> 6:  2006 Male_Q1 57.89281   Male       Q1
```
