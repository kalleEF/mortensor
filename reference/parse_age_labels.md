# Parse age-group labels

Understands the label styles of common sources, for example `"<1"`,
`"Under 1"`, `"01-14"`, `"1 to 4"`, `"30-34"`, `"30-34 years"`,
`"1 and under 5"`, `"90+"`, `"85 and over"`, `"85 years and older"`,
`"90 & over"`, `"Age 85+"`, `"0"`, `"110+"`, and the German labels of
the Federal Statistical Office (`"unter 1 Jahr"`,
`"15 bis unter 20 Jahre"`, `"1 bis 4 Jahre"`, `"85 Jahre und mehr"`, and
single years written as `"68-Jaehrige"` with an umlaut).

## Usage

``` r
parse_age_labels(x)
```

## Arguments

- x:

  Character vector of labels.

## Value

A `data.table` with `label`, `age_lo`, `age_hi` (`Inf` for open
intervals), in the order of `x`.

## Examples

``` r
parse_age_labels(c("<1", "01-14", "15-19", "90+"))
#>     label age_lo age_hi
#>    <char>  <num>  <num>
#> 1:     <1      0      1
#> 2:  01-14      1     15
#> 3:  15-19     15     20
#> 4:    90+     90    Inf
parse_age_labels(c("Under 1", "1 to 4", "85 years and over"))
#>                label age_lo age_hi
#>               <char>  <num>  <num>
#> 1:           Under 1      0      1
#> 2:            1 to 4      1      5
#> 3: 85 years and over     85    Inf
```
