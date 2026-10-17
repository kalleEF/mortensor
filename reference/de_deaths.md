# Deaths and populations by Land, Germany

Deaths by year (1998-2024), Land, sex, age group and ICD-10 chapter, and
populations on 31 December (1997-2024) by Land, sex and single year of
age. The vignette "Differences between the German Laender" reads them
with [`mt_data()`](mt_data.md) and prepares them with
[`mt_as_exposure()`](mt_as_exposure.md),
[`mt_regroup()`](mt_regroup.md), [`mt_subset()`](mt_subset.md) and
[`mt_fill_suppressed()`](mt_fill_suppressed.md).

## Usage

``` r
de_deaths

de_population

de_cause_groups
```

## Format

`de_deaths`: a data frame with columns `year`, `region`, `sex`
(`"male"`, `"female"`), `age`, `cause`, `deaths`. `de_population`: a
data frame with columns `year` (population on 31 December), `region`,
`sex`, `age` (`"unter 1 Jahr"`, single years, `"90 Jahre und mehr"`),
`population`. `de_cause_groups`: a named list of chapter codes.

An object of class `data.frame` with 81536 rows and 5 columns.

An object of class `list` of length 7.

## Source

Deaths: Statistisches Bundesamt (Destatis), Todesursachenstatistik,
table "Sterbefaelle, Sterbeziffern (ab 1998)" of the
Gesundheitsberichterstattung des Bundes, <https://www.gbe-bund.de>,
retrieved 13 October 2026. Population: Statistisches Bundesamt
(Destatis), Genesis-Online, table 12411-0013, retrieved 13 October 2026.
Datenlizenz Deutschland - Namensnennung - Version 2.0
(<https://www.govdata.de/dl-de/by-2-0>); reshaped, suppressed cells kept
as missing (own processing, `data-raw/de_data.R`).

## Details

`de_deaths` holds the national total (`region = "Deutschland"`) and the
16 Laender, 18 age groups (column `age`: `"Unter 1 Jahr"`,
`"1 Jahr bis unter 15 Jahre"`, 5-year groups, `"90 Jahre und aelter"`
with umlaut), all causes (`cause = "all"`) and the 18 ICD-10 chapters of
the source table (codes such as `"C00-D48"`; the labels are in attribute
`"cause_labels"`). Deaths that the source suppresses for confidentiality
(`"."`) are `NA`: 13% of the cells, mostly rare chapters in small
Laender; 33 all-cause cells (Bremen, Saarland, Hamburg, Sachsen-Anhalt,
ages 1-34). Deaths from pregnancy and childbirth in men (`"X"`, not
applicable) are 0. The causes add up to all causes and the Laender to
Germany wherever all cells are published; 2025 is left out because its
chapters did not add up. Injuries and poisonings (S00-T98) are the
chapter of external causes.

`de_cause_groups` lists cause groups for cause-specific analyses at ages
50 and over: neoplasms, diseases of the circulatory, respiratory and
digestive systems, mental and nervous system disorders combined
(dementia is coded in both chapters), endocrine diseases and external
causes; the remaining chapters form `"rest"` with
[`mt_group_levels()`](mt_group_levels.md).

## Examples

``` r
head(de_deaths)
#>   year            region    sex          age   cause deaths
#> 1 1998 Baden-Württemberg female Unter 1 Jahr A00-B99     NA
#> 2 1998 Baden-Württemberg female Unter 1 Jahr C00-D48      0
#> 3 1998 Baden-Württemberg female Unter 1 Jahr D50-D90      0
#> 4 1998 Baden-Württemberg female Unter 1 Jahr E00-E90      4
#> 5 1998 Baden-Württemberg female Unter 1 Jahr F00-F99      0
#> 6 1998 Baden-Württemberg female Unter 1 Jahr G00-H95     NA
x <- mt_data(deaths = subset(de_deaths, cause == "all" & year >= 2018),
             population = subset(de_population, year >= 2017),
             by = c("region", "sex"))
#> Warning: 10 missing deaths value(s) left out; suppressed counts can be estimated with mt_fill_suppressed().
x <- mt_as_exposure(x, type = "dec31")
x
#> <mt_data> 34 strata by region, sex, omega = 110
#>   deaths     mixed widths (1,5,14) 0-89, 90+; periods 2018-2024 (width 1); 4274 cells
#>   exposure   single years 0-89, 90+; periods 2018-2024 (width 1); 20384 cells
#>   population single years 0-89, 90+; periods 2017-2024 (width 1); 23296 cells
#>   grid       <mt_axis: age, 91 intervals> 0 1 2 3 ... 88 89 90+
#>              <mt_axis: period, 7 intervals> 2018 2019 2020 2021 2022 2023 2024
```
