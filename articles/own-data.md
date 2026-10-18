# Your own data, from CSV to report

This vignette goes from files as a statistical office publishes them to
a report with forecasts, observed indicators and inequality measures.
The files are simulated, with the peculiarities of real exports: German
column names and age labels, sexes coded 1 and 2, small counts
suppressed as `"."`, zeros written as `"-"`, rows with totals over all
ages, and populations on 31 December.

``` r

library(mortensor)
library(data.table)
```

## The files

The code below writes two CSV files (semicolon-separated, as in many
European exports). With your own data, you start from your files.

``` r

set.seed(42)
ages <- c("unter 1 Jahr", "1 bis unter 5 Jahre", paste(seq(5, 85, 5), "bis unter", seq(10, 90, 5), "Jahre"),
          "90 Jahre und mehr")
lo <- c(0, 1, seq(5, 90, 5))
g <- CJ(Jahr = 2005:2022, Geschlecht = 1:2, Quintil = 1:5, k = seq_along(ages))
g[, Altersgruppe := ages[k]]
g[, bev := round(4e4 * exp(-lo[k] / 45) * (1 + 0.1 * (Quintil == 1)))]
g[, rate := exp(-10.1 + 0.09 * (lo[k] + 2) + 0.3 * (Geschlecht == 1) + 0.12 * (5 - Quintil) - 0.015 * (Jahr - 2005))]
g[k == 1, rate := 0.004]
g[Jahr %in% 2020:2022, rate := rate * 1.12]                 # pandemic years
g[, Gestorbene := rpois(.N, bev * pmin(rate, 0.7))]
deaths <- g[, .(Jahr, Geschlecht, Quintil, Altersgruppe,
                Gestorbene = ifelse(Gestorbene == 0, "-", ifelse(Gestorbene < 3, ".", Gestorbene)))]
totals <- g[, .(Altersgruppe = "Insgesamt", Gestorbene = as.character(sum(Gestorbene))),
            by = .(Jahr, Geschlecht, Quintil)]
deaths <- rbind(deaths, totals)
pop <- g[, .(Stichtag = Jahr, Geschlecht, Quintil, Altersgruppe, Bevoelkerung = bev)]
pop <- rbind(pop[Stichtag == 2005][, Stichtag := 2004], pop)   # 31 December of the year before
f_deaths <- file.path(tempdir(), "gestorbene.csv")
f_pop <- file.path(tempdir(), "bevoelkerung.csv")
fwrite(deaths, f_deaths, sep = ";")
fwrite(pop, f_pop, sep = ";")
```

## Reading

[`read.csv2()`](https://rdrr.io/r/utils/read.table.html) reads
semicolon-separated files; reading every column as text keeps symbols
such as `"."` intact. [`mt_data()`](../reference/mt_data.md) then gets
the names of the columns: stratum variables (`by`), age, period and the
value of each table.

``` r

d <- read.csv2(f_deaths, colClasses = "character")
p <- read.csv2(f_pop, colClasses = "character")
head(d)
#>   Jahr Geschlecht Quintil          Altersgruppe Gestorbene
#> 1 2005          1       1          unter 1 Jahr        194
#> 2 2005          1       1   1 bis unter 5 Jahre          4
#> 3 2005          1       1  5 bis unter 10 Jahre          9
#> 4 2005          1       1 10 bis unter 15 Jahre         10
#> 5 2005          1       1 15 bis unter 20 Jahre         13
#> 6 2005          1       1 20 bis unter 25 Jahre         13
```

The population table names the year `Stichtag`, so it is added with
[`mt_add_population()`](../reference/mt_add.md) and its own column
names:

``` r

x <- mt_data(deaths = d, by = c("Geschlecht", "Quintil"), age = "Altersgruppe", period = "Jahr",
             value = "Gestorbene")
#> Warning in ingest_obs(x, by, age, period, value, measure, age_open, missing =
#> missing, : 180 row(s) of the deaths table with age or period labelled as a
#> total ("Insgesamt") are left out.
#> Warning in ingest_obs(x, by, age, period, value, measure, age_open, missing =
#> missing, : 142 missing deaths value(s) left out; suppressed counts can be
#> estimated with mt_fill_suppressed().
x <- mt_add_population(x, p, age = "Altersgruppe", period = "Stichtag", value = "Bevoelkerung")
x
#> <mt_data> 10 strata by Geschlecht, Quintil, omega = 110
#>   deaths     mixed widths (1,4,5) 0-89, 90+; periods 2005-2022 (width 1); 3458 cells
#>   population mixed widths (1,4,5) 0-89, 90+; periods 2004-2022 (width 1); 3800 cells
#>   grid       <mt_axis: age, 20 intervals> 0 1-4 5-9 10-14 ... 80-84 85-89 90+
#>              <mt_axis: period, 19 intervals> 2004 2005 2006 2007 ... 2020 2021 2022
```

The warnings report what [`mt_data()`](../reference/mt_data.md) did with
the peculiarities: rows labelled `"Insgesamt"` (totals over all ages)
are left out, `"-"` is read as zero and the suppressed cells (`"."`) as
missing.

## Preparing

Each step is a decision; the functions make it explicit.

``` r

x <- mt_relabel(x, "Geschlecht", c("1" = "male", "2" = "female"))
x <- mt_relabel(x, "Quintil", c("1" = "Q1", "2" = "Q2", "3" = "Q3", "4" = "Q4", "5" = "Q5"))
x <- mt_as_exposure(x, type = "dec31")      # populations on 31 December
mt_validate(x)
#> <mt_validation> 1 finding(s): 0 error, 0 warning, 1 info
#>   [info] grid coverage (deaths, all, n = 142): grid cells not covered by observations
```

Without published totals, suppressed counts are estimated within the
range that the suppression rule implies (here 1 or 2 deaths), as
expected deaths at the pooled rate of the age group and year:

``` r

x <- mt_fill_suppressed(x, lower = 1, upper = 2)
attr(x, "fill")
#>    filled missing inconsistent iterations max_gap
#>     <int>   <int>        <int>      <int>   <num>
#> 1:    142       0            0          0       0
```

The pandemic years 2020-2022 are left out of the fitting period:

``` r

x <- mt_subset(x, periods = c(2005, 2019))
```

## The past: observed indicators

Life expectancy, age-standardised rates and the gap between quintiles
are computed from the observed rates; the intervals reflect Poisson
variation of the death counts.

``` r

e0 <- life_expectancy(x)
e0[period == 2019, .(Geschlecht, Quintil, q0.025, q0.5, q0.975)]
#>     Geschlecht Quintil   q0.025     q0.5   q0.975
#>         <char>  <char>    <num>    <num>    <num>
#>  1:       male      Q1 72.54543 73.04346 73.51449
#>  2:       male      Q2 74.09485 74.56660 75.04868
#>  3:       male      Q3 74.85565 75.37452 75.93295
#>  4:       male      Q4 76.61884 77.15492 77.65594
#>  5:       male      Q5 77.76862 78.30108 78.77808
#>  6:     female      Q1 76.03602 76.49766 76.98199
#>  7:     female      Q2 77.32803 77.79842 78.31247
#>  8:     female      Q3 79.20649 79.71905 80.29379
#>  9:     female      Q4 79.45129 80.06467 80.62082
#> 10:     female      Q5 81.77045 82.31385 82.87126
gap <- mt_contrast(x, var = "Quintil", ref = "Q5")
gap[period == 2019 & level == "Q1"]
#>    Geschlecht  level    ref period      mean    q0.025      q0.5    q0.975
#>        <char> <char> <char>  <num>     <num>     <num>     <num>     <num>
#> 1:       male     Q1     Q5   2019 -5.259102 -5.987421 -5.245696 -4.541609
#> 2:     female     Q1     Q5   2019 -5.807940 -6.510433 -5.804283 -5.086910
#>    p_positive        quantity       type
#>         <num>          <char>     <char>
#> 1:          0 life_expectancy difference
#> 2:          0 life_expectancy difference
mt_inequality(x, "Quintil", order = paste0("Q", 5:1), quantity = "life_expectancy")[period == 2019]
#>    Geschlecht  index period       mean     q0.025       q0.5     q0.975
#>        <char> <char>  <num>      <num>      <num>      <num>      <num>
#> 1:       male    sii   2019 -6.6226995 -7.4169824 -6.6084044 -5.7942876
#> 2:       male    rii   2019  0.9161272  0.9064890  0.9163064  0.9262289
#> 3:     female    sii   2019 -7.0064219 -7.7924671 -7.0167733 -6.1357340
#> 4:     female    rii   2019  0.9153267  0.9062257  0.9151347  0.9253908
#>           quantity
#>             <char>
#> 1: life_expectancy
#> 2: life_expectancy
#> 3: life_expectancy
#> 4: life_expectancy
```

## Choosing a model

Three candidates are compared in a backtest with three forecast origins
and five years ahead: Lee-Carter per stratum; Li-Lee with a common trend
per sex, which keeps the quintiles from diverging; and an ensemble of
Lee-Carter and the random walk with drift per age.

``` r

cand <- list("Lee-Carter" = lee_carter(), "Li-Lee" = li_lee(common_by = "Geschlecht"),
             "Lee-Carter + random walk" = ensemble(lc = lee_carter(), rw = rwd()))
bt <- mt_compare(cand, x, origins = 2012:2014, h = 5, n_draws = 200)
#> Warning in check_draws_plausible(arr, fit$age, fit$strata, spec$label):
#> Lee-Carter: forecast death rates above 1 in 21.5% of the draws (strata
#> female.Q3; ages 10-14). The model is unstable here, often because of sparse
#> counts or a steep trend extrapolated far ahead; consider wider age groups
#> (mt_regroup()), pooling (li_lee(common_by = )), smoothing (step_smooth()) or an
#> ensemble.
#> Warning in check_draws_plausible(arr, fit$age, fit$strata, spec$label):
#> Lee-Carter: forecast death rates above 1 in 4.5% of the draws (strata
#> female.Q3; ages 10-14). The model is unstable here, often because of sparse
#> counts or a steep trend extrapolated far ahead; consider wider age groups
#> (mt_regroup()), pooling (li_lee(common_by = )), smoothing (step_smooth()) or an
#> ensemble.
#> Warning in report_fit_problems(diagnostics, spec$label, classes): Li-Lee
#> (common by Geschlecht): Li-Lee: the AR(1) coefficient of the stratum index is
#> bounded at |phi| = 0.99 in 1 stratum/strata (male.Q5): their deviations from
#> the common trend are close to a unit root and revert very slowly; coherence may
#> not suit them. (diagnostics of 10 stratum/strata of the fit).
#> Warning in check_draws_plausible(arr, fit$age, fit$strata, spec$label):
#> Lee-Carter: forecast death rates above 1 in 21.5% of the draws (strata
#> female.Q3; ages 10-14). The model is unstable here, often because of sparse
#> counts or a steep trend extrapolated far ahead; consider wider age groups
#> (mt_regroup()), pooling (li_lee(common_by = )), smoothing (step_smooth()) or an
#> ensemble.
#> Warning in check_draws_plausible(arr, fit$age, fit$strata, spec$label): Random
#> walk with drift: forecast death rates above 1 in 3.5% of the draws (strata
#> female.Q3, female.Q5; ages 10-14, 15-19). The model is unstable here, often
#> because of sparse counts or a steep trend extrapolated far ahead; consider
#> wider age groups (mt_regroup()), pooling (li_lee(common_by = )), smoothing
#> (step_smooth()) or an ensemble.
#> Warning in check_draws_plausible(arr, fit$age, fit$strata, spec$label):
#> Lee-Carter: forecast death rates above 1 in 4.5% of the draws (strata
#> female.Q3; ages 10-14). The model is unstable here, often because of sparse
#> counts or a steep trend extrapolated far ahead; consider wider age groups
#> (mt_regroup()), pooling (li_lee(common_by = )), smoothing (step_smooth()) or an
#> ensemble.
#> Warning in check_draws_plausible(arr, fit$age, fit$strata, spec$label): Random
#> walk with drift: forecast death rates above 1 in 2.5% of the draws (strata
#> female.Q3, female.Q5; ages 10-14, 15-19). The model is unstable here, often
#> because of sparse counts or a steep trend extrapolated far ahead; consider
#> wider age groups (mt_regroup()), pooling (li_lee(common_by = )), smoothing
#> (step_smooth()) or an ensemble.
scores <- mt_score_table(bt, metrics = c("crps", "cov95"))
scores
#>                       model  crps cov95
#>                      <fctr> <num> <num>
#> 1:               Lee-Carter 0.148 0.909
#> 2:                   Li-Lee 0.175 0.880
#> 3: Lee-Carter + random walk 0.156 0.968
```

Warnings about implausible rates in some backtest fits point to young
age groups with few deaths. The model with the lowest CRPS (the rows
follow the order of the list) is used for the forecast:

``` r

best <- cand[[which.min(scores$crps)]]
best$label
#> [1] "Lee-Carter"
```

## Forecasting

``` r

fit <- mt_fit(best, x)
summary(fit)
#> Lee-Carter
#>   data: 10 strata; ages 0 to 90+; periods 2005-2019
#>   model and uncertainty:
#>     log m(x,t) = a(x) + b(x) k(t), estimated by Poisson maximum likelihood
#>     (Brouhns et al. 2002). k(t) follows a random walk with drift; draws include
#>     its innovations and the uncertainty of the drift. a(x) and b(x) are fixed
#>     (parameter uncertainty with bootstrapped()). Predictive draws add Poisson
#>     variation of the deaths given the rates. Jump-off: fitted rates.
#>   fit: Poisson deviance 2481.2 on 3000 cells (0.83 per cell); 2.7% of deviance residuals beyond +/-2
fc <- mt_forecast(fit, to = 2035, n_draws = 500, seed = 1)
e0_fc <- life_expectancy(fc)
e0_fc[period == 2035, .(Geschlecht, Quintil, q0.025, q0.5, q0.975)]
#>     Geschlecht Quintil   q0.025     q0.5   q0.975
#>         <char>  <char>    <num>    <num>    <num>
#>  1:       male      Q1 73.57397 76.04022 78.39923
#>  2:       male      Q2 74.60027 77.18007 80.09615
#>  3:       male      Q3 75.40221 77.73909 80.69551
#>  4:       male      Q4 76.89530 79.84707 83.14944
#>  5:       male      Q5 77.46679 80.63718 83.65660
#>  6:     female      Q1 77.59210 79.20813 80.75888
#>  7:     female      Q2 78.12215 80.56010 83.13295
#>  8:     female      Q3 77.48239 83.08839 90.03644
#>  9:     female      Q4 78.86362 82.58275 87.21109
#> 10:     female      Q5 82.17392 86.09262 90.35157
mt_contrast(fc, var = "Quintil", ref = "Q5")[period == 2035 & level == "Q1"]
#>    Geschlecht  level    ref period      mean     q0.025      q0.5     q0.975
#>        <char> <char> <char>  <num>     <num>      <num>     <num>      <num>
#> 1:       male     Q1     Q5   2035 -4.646715  -8.513676 -4.663307 -0.7107595
#> 2:     female     Q1     Q5   2035 -6.912613 -11.376681 -6.909357 -2.6093638
#>    p_positive        quantity       type
#>         <num>          <char>     <char>
#> 1:      0.012 life_expectancy difference
#> 2:      0.000 life_expectancy difference
```

``` r

mt_plot_life_expectancy(fc, facet = "Geschlecht", strata = c("male.Q1", "male.Q5", "female.Q1", "female.Q5"))
```

![](own-data_files/figure-html/unnamed-chunk-12-1.png)

## Report tables

[`mt_export()`](../reference/mt_export.md) writes any of these tables,
or the forecast draws, as CSV:

``` r

f <- mt_export(e0_fc, tempdir(), file = "e0_forecast.csv", digits = 4)
head(read.csv(f), 3)
#>   stratum Geschlecht Quintil   label period  mean q0.025  q0.1  q0.5  q0.9
#> 1       1       male      Q1 male.Q1   2020 73.25  72.83 72.96 73.26 73.54
#> 2       1       male      Q1 male.Q1   2021 73.44  72.84 73.01 73.45 73.88
#> 3       1       male      Q1 male.Q1   2022 73.61  72.86 73.09 73.60 74.14
#>   q0.975 at
#> 1  73.68  0
#> 2  74.12  0
#> 3  74.41  0
```

## What to check in your own data

- the warnings of [`mt_data()`](../reference/mt_data.md): totals,
  unknown ages, symbols that are not declared, levels that differ
  between tables (for example `"F"` and `"female"`);
- [`mt_validate()`](../reference/mt_validate.md): impossible values,
  missing cells, levels that are totals of other levels;
- the type of the population counts in
  [`mt_as_exposure()`](../reference/mt_as_exposure.md);
- warnings of [`mt_fit()`](../reference/mt_fit.md) and
  [`mt_forecast()`](../reference/mt_forecast.md) about sparse counts and
  implausible rates, and `summary(fit)`.
