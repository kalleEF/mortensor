# Standard populations

- `"esp2013"`: the revised European Standard Population 2013 (Eurostat
  2013, Revision of the European Standard Population, report
  KS-RA-13-028), in 5-year groups with 0, 1-4 and 95+, total 100,000;

- `"esp1976"`: the European Standard Population of 1976 (0, 1-4, 5-year
  groups, 85+), total 100,000;

- `"who2000"`: the WHO World Standard Population, world average
  2000-2025 (Ahmad et al. 2001), 5-year groups to 100+; the published
  weights sum to 100,035 and are used as published (rates are divided by
  the sum of the weights).

## Usage

``` r
standard_population(name = c("esp2013", "esp1976", "who2000"))
```

## Arguments

- name:

  Name of the standard.

## Value

A `data.table` with `age_lo`, `age_hi`, `weight`; the name is in
attribute `"name"`.

## Details

A national standard is used by passing a data frame with `age_lo`,
`age_hi` (`Inf` for the open group) and `weight` to the `standard`
argument of [`asr()`](derived.md).

## References

Ahmad OB, Boschi-Pinto C, Lopez AD, Murray CJL, Lozano R, Inoue M
(2001). Age standardization of rates: a new WHO standard. GPE Discussion
Paper Series No. 31. World Health Organization, Geneva.

## Examples

``` r
standard_population()
#>     age_lo age_hi weight
#>      <num>  <num>  <num>
#>  1:      0      1   1000
#>  2:      1      5   4000
#>  3:      5     10   5500
#>  4:     10     15   5500
#>  5:     15     20   5500
#>  6:     20     25   6000
#>  7:     25     30   6000
#>  8:     30     35   6500
#>  9:     35     40   7000
#> 10:     40     45   7000
#> 11:     45     50   7000
#> 12:     50     55   7000
#> 13:     55     60   6500
#> 14:     60     65   6000
#> 15:     65     70   5500
#> 16:     70     75   5000
#> 17:     75     80   4000
#> 18:     80     85   2500
#> 19:     85     90   1500
#> 20:     90     95    800
#> 21:     95    Inf    200
#>     age_lo age_hi weight
#>      <num>  <num>  <num>
standard_population("who2000")
#>     age_lo age_hi weight
#>      <num>  <num>  <num>
#>  1:      0      5   8860
#>  2:      5     10   8690
#>  3:     10     15   8600
#>  4:     15     20   8470
#>  5:     20     25   8220
#>  6:     25     30   7930
#>  7:     30     35   7610
#>  8:     35     40   7150
#>  9:     40     45   6590
#> 10:     45     50   6040
#> 11:     50     55   5370
#> 12:     55     60   4550
#> 13:     60     65   3720
#> 14:     65     70   2960
#> 15:     70     75   2210
#> 16:     75     80   1520
#> 17:     80     85    910
#> 18:     85     90    440
#> 19:     90     95    150
#> 20:     95    100     40
#> 21:    100    Inf      5
#>     age_lo age_hi weight
#>      <num>  <num>  <num>
```
