# Deaths and populations by deprivation decile, England 2001-2019

Death registrations and mid-year populations by sex, age group and
decile of the English Index of Multiple Deprivation (IMD), with deaths
by underlying cause (ICD-10), from an ONS ad hoc release. Years 2001 to
2019 (before the COVID-19 pandemic) of the first workbook; the column of
deaths involving COVID-19 is dropped (zero in these years). Built by
`data-raw/ons_data.R`.

## Usage

``` r
ons_data

ons_cause_groups
```

## Format

`ons_data`: a data frame with 7,600 rows (19 years x 2 sexes x 20 age
groups x 10 deciles) and columns `year`, `sex` (`"female"`, `"male"`),
`age_group` (`"<1"`, `"01-04"`, 5-year groups, `"90+"`), `imd_decile`
(1-10), `imd_version`, `population`, `all_causes` and 29 cause columns.
`ons_cause_groups`: a named list of cause columns.

An object of class `list` of length 11.

## Source

Office for National Statistics (2025). Death registrations and
populations by Index of Multiple Deprivation, England. Ad hoc release
3064, 13 October 2025.
<https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/deaths/adhocs/3064deathregistrationsandpopulationsbyindexofmultipledeprivationengland>.
Source: Office for National Statistics licensed under the Open
Government Licence v3.0.

## Details

IMD decile 1 is the most deprived. Each year uses the IMD version
current at the time (`imd_version`: 2004 for 2001-2005, 2007 for
2006-2008, 2010 for 2009-2012, 2015 for 2013-2017, 2019 for 2018-2019),
so the deciles are reclassified four times. Deaths are by year of
registration for usual residents of England. The causes are mutually
exclusive (see the ICD-10 lookup of the workbook); neonatal deaths have
no cause and count only in `all_causes`.

`ons_cause_groups` lists the cause groups used in the vignettes for ages
50 and over, for [`mt_group_levels()`](mt_group_levels.md): lung,
colorectal, breast or prostate (sex-specific cancers; for men including
the few male breast cancers), other cancers, ischaemic heart disease,
stroke, other circulatory (heart failure, hypertension, atrial
fibrillation), dementia, chronic lower respiratory (COPD, asthma),
diabetes, liver disease and alcohol-specific causes. All other deaths,
including suicide and the causes not listed in the table, form the
remainder, all causes minus these groups
(`other = "rest", total = "all_causes"`). Small causes are combined
within their ICD-10 chapter so that cells by sex, quintile, 5-year age
group and year have enough deaths.

## Examples

``` r
head(ons_data[, 1:8])
#>   year    sex age_group imd_decile imd_version population all_causes
#> 1 2001 female        <1          1        2004      35977        307
#> 2 2001 female        <1          2        2004      32042        234
#> 3 2001 female        <1          3        2004      28759        169
#> 4 2001 female        <1          4        2004      27016        152
#> 5 2001 female        <1          5        2004      25451         89
#> 6 2001 female        <1          6        2004      24814         94
#>   breast_cancer
#> 1             0
#> 2             0
#> 3             0
#> 4             0
#> 5             0
#> 6             0
# all causes by sex and deprivation quintile
x <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile"),
             age = "age_group", value = c(deaths = "all_causes", population = "population"))
x <- mt_as_exposure(x, type = "midyear")
x <- mt_group_levels(x, "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)),
                     name = "imd_quintile")
x
#> <mt_data> 10 strata by sex, imd_quintile, omega = 110
#>   deaths     mixed widths (1,4,5) 0-89, 90+; periods 2001-2019 (width 1); 3800 cells
#>   exposure   mixed widths (1,4,5) 0-89, 90+; periods 2001-2019 (width 1); 3800 cells
#>   population mixed widths (1,4,5) 0-89, 90+; periods 2001-2019 (width 1); 3800 cells
#>   grid       <mt_axis: age, 20 intervals> 0 1-4 5-9 10-14 ... 80-84 85-89 90+
#>              <mt_axis: period, 19 intervals> 2001 2002 2003 2004 ... 2017 2018 2019
# cause groups at ages 50 and over: one column per cause
causes <- names(ons_data)[-(1:6)]
y <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile", "cause"),
             age = "age_group", value = list(deaths = causes, population = "population"))
y <- mt_subset(mt_as_exposure(y, type = "midyear"), ages = c(50, Inf))
y <- mt_group_levels(y, "cause", ons_cause_groups, other = "rest", total = "all_causes",
                     additive = "rates")
unique(mt_strata(y)$cause)
#>  [1] "lung_cancer"       "colorectal_cancer" "breast_prostate"  
#>  [4] "other_cancers"     "ihd"               "stroke"           
#>  [7] "other_circulatory" "dementia"          "respiratory"      
#> [10] "diabetes"          "liver_alcohol"     "all_causes"       
#> [13] "rest"             
```
