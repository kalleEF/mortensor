# Read deaths and populations from GENESIS (Destatis) flat files

Reads tables of the German Federal Statistical Office in the "flat" CSV
format of GENESIS-Online and the regional database (columns `time`,
`1_variable_code`, `1_variable_attribute_label`, ..., `value`): deaths
(for example table 23211-0002) and populations (for example table
12411-0013, on 31 December).

## Usage

``` r
mt_read_genesis(deaths, population = NULL, years = NULL)
```

## Arguments

- deaths:

  Path of a flat CSV file with deaths, or a data frame read from one.

- population:

  Path or data frame with populations, or `NULL`.

- years:

  Optional subset of years of the deaths (populations of the year before
  are kept for exposures from 31 December populations).

## Value

An [mt_data](mt_data.md) object with deaths and populations.

## Details

Variables are recognised by their codes: `GES` (sex), codes starting
with `ALT` (age, German labels such as `"unter 1 Jahr"`,
`"15 bis unter 20 Jahre"`, `"85 Jahre und mehr"`, `"68-Jaehrige"` with
umlaut), `DLAND`, `KREISE` or `REGBEZ` (region), codes starting with
`TODUR` (cause of death; the total "Insgesamt" becomes level `"all"`,
the other levels are the GENESIS codes, their labels are in attribute
`"cause_labels"`), and `DINSG` (Germany as a whole, dropped). Other
variables are kept as stratum variables named by their code. Values
`"-"` (nothing) are zero; `"."`, `"..."`, `"x"` and `"/"` are missing.
Stratum variables that only the population table has (for example
regions, when the deaths are national) are summed out of the
populations.

The reader makes no analysis decisions: populations stay populations,
and the periods are those of the source (the year of the reference
date). The next steps are usually
[`mt_as_exposure()`](mt_as_exposure.md) with `type = "dec31"`,
[`mt_regroup()`](mt_regroup.md) with `age = "deaths"` (single-year
populations to the age groups of the deaths),
[`mt_fill_suppressed()`](mt_fill_suppressed.md) and
[`mt_group_levels()`](mt_group_levels.md) (cause groups, with
`total = "all"` for a residual group).

## Examples

``` r
d <- data.frame(time = 2020, `1_variable_code` = "GES", `1_variable_attribute_label` = "weiblich",
                `2_variable_code` = "ALT027", `2_variable_attribute_label` = "85 Jahre und mehr",
                value = "120", check.names = FALSE)
p <- data.frame(time = c("2019-12-31", "2020-12-31"), `1_variable_code` = "GES",
                `1_variable_attribute_label` = "weiblich", `2_variable_code` = "ALT103",
                `2_variable_attribute_label` = "85 Jahre und mehr", value = c(900, 950),
                check.names = FALSE)
x <- mt_read_genesis(d, p)
x <- mt_as_exposure(x, type = "dec31")
mt_obs(x, "exposure")
#>    stratum    sex  label age_lo age_hi period_lo period_hi value
#>      <int> <char> <char>  <num>  <num>     <num>     <num> <num>
#> 1:       1 female female     85    Inf      2020      2021   925
#>          provenance
#>              <char>
#> 1: population_dec31
```
