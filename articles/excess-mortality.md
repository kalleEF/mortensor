# Excess mortality

Excess deaths are the deaths observed in a period of crisis (a heatwave,
an epidemic, a war) minus the deaths expected without it. The expected
deaths are not observed; they come from a model of the years before the
crisis, the *baseline*. [`mt_excess()`](../reference/mt_excess.md) fits
any mortality model of the package to the baseline years, forecasts the
later years with the observed exposures, and compares the observed
deaths with each simulated path of expected deaths. The intervals
therefore contain the uncertainty of the baseline trend and the random
variation of the death counts. The WHO estimates of the excess mortality
of the COVID-19 pandemic follow this approach with their own baseline
model (Msemburi et al. 2023).

The results are

- the **excess**: observed minus expected deaths;
- the **P-score**: the excess in percent of the expected deaths;
- `p_positive`: the share of paths in which the observed deaths exceed
  the expected deaths;

by year, and optionally by age group, stratum, and summed over years.

## The heatwave of August 2003 in France

Deaths and exposures of France (Human Mortality Database, via the
demography package), single years of age to 100+. The baseline is a
Lee-Carter model of the years 1990-2002:

``` r

library(mortensor)
library(data.table)
fr <- as_mt_data(demography::fr.mort)
x <- mt_regroup(mt_subset(fr, periods = c(1990, 2006)), age = 0:100)
fit <- mt_fit(lee_carter(), mt_subset(x, periods = c(1990, 2002)))
ex <- mt_excess(fit, x, by = character(), seed = 1)
ex[, list(period, observed, expected_q0.5, excess_q0.025, excess_q0.5, excess_q0.975, p_score_q0.5)]
#>    period observed expected_q0.5 excess_q0.025 excess_q0.5 excess_q0.975
#>     <num>    <num>         <num>         <num>       <num>         <num>
#> 1:   2003   552339      537269.0      5606.375     15070.0     25445.525
#> 2:   2004   509429      539842.5    -44524.750    -30413.5    -16279.675
#> 3:   2005   527535      543595.0    -33115.700    -16060.0       798.525
#> 4:   2006   516416      548697.5    -52128.375    -32281.5    -13218.850
#>    p_score_q0.5
#>           <num>
#> 1:     2.804926
#> 2:    -5.633773
#> 3:    -2.954405
#> 4:    -5.883296
```

`by = character()` sums women and men. In 2003, about 15,000 more people
died than expected (95% interval 5,600 to 25,400), 2.8% more than
expected. Fouillet et al. (2006) counted about 15,000 excess deaths
between 1 and 20 August 2003 against the mortality of 2000-2002, so the
rest of the year was close to the baseline.

``` r

mt_plot_excess(ex)
```

![](excess-mortality_files/figure-html/plot-expected-1.png)

The years after 2003 lie below the baseline: 30,000 fewer deaths than
expected in 2004 and as many in 2006. Part of the deficit of 2004 may be
deaths brought forward by the heatwave, but the deficit continues: from
2003 mortality fell faster than the trend of 1990-2002. A baseline
extrapolated over several years becomes uncertain, and so do excess
deaths far from the baseline period.

By sex and age group:

``` r

mt_excess(fit, x, periods = 2003, seed = 1)[, list(sex, observed, excess_q0.025, excess_q0.5, excess_q0.975)]
#>       sex observed excess_q0.025 excess_q0.5 excess_q0.975
#>    <char>    <num>         <num>       <num>         <num>
#> 1: female   273315      3448.800     11671.5     19341.575
#> 2:   male   279024     -1337.075      3531.0      8953.625
ex_age <- mt_excess(fit, x, periods = 2003, by = character(), ages = c(0, 45, 65, 75, 85), seed = 1)
ex_age[, list(age, observed, excess_q0.5, p_score_q0.025, p_score_q0.5, p_score_q0.975)]
#>       age observed excess_q0.5 p_score_q0.025 p_score_q0.5 p_score_q0.975
#>    <char>    <num>       <num>          <num>        <num>          <num>
#> 1:   0-44    29115      -439.0     -4.4707459    -1.485417      1.3511059
#> 2:  45-64    81495      2498.5      1.5390342     3.162798      4.8867019
#> 3:  65-74    89304     -1400.5     -3.5247795    -1.544025      0.3337944
#> 4:  75-84   167633      5619.0      0.5302193     3.468219      6.3495506
#> 5:    85+   184792      8982.0      3.3076824     5.108924      6.9832477
```

``` r

mt_plot_excess(ex_age, type = "p_score")
```

![](excess-mortality_files/figure-html/plot-age-1.png)

Three quarters of the excess were women, and most were 75 and older; the
P-score was highest at 85+. Fouillet et al. (2006) also found the excess
concentrated among women and at old ages. `ages` gives the start of each
age group; the deaths are summed path by path before the quantiles are
taken.

The same comparison on the scale of life expectancy: the observed life
expectancy at birth minus that of each path of expected rates.

``` r

le <- mt_excess(fit, x, periods = 2003, quantity = "life_expectancy", seed = 1)
le[, list(sex, observed, expected_q0.5, excess_q0.025, excess_q0.5, excess_q0.975)]
#>       sex observed expected_q0.5 excess_q0.025 excess_q0.5 excess_q0.975
#>    <char>    <num>         <num>         <num>       <num>         <num>
#> 1: female 82.95594      83.20279    -0.5685339  -0.2468486     0.1022864
#> 2:   male 75.87930      76.00291    -0.3828217  -0.1236114     0.0994448
```

Life expectancy of women was 0.25 years lower than expected; the
interval includes zero, because an annual life expectancy hides a shock
of two weeks.

## The choice of baseline

The baseline is an assumption. A few alternatives, each fitted to the
years before 2003:

``` r

base <- list("Lee-Carter 1990-2002" = fit,
             "Lee-Carter 1993-2002" = mt_fit(lee_carter(), mt_subset(x, periods = c(1993, 2002))),
             "random walk 1990-2002" = mt_fit(rwd(), mt_subset(x, periods = c(1990, 2002))),
             "Li-Lee 1990-2002" = mt_fit(li_lee(), mt_subset(x, periods = c(1990, 2002))))
rbindlist(lapply(names(base), function(b) {
  mt_excess(base[[b]], x, periods = 2003, by = character(), seed = 1)[, baseline := b]
}))[, list(baseline, excess_q0.025, excess_q0.5, excess_q0.975, p_score_q0.5)]
#>                 baseline excess_q0.025 excess_q0.5 excess_q0.975 p_score_q0.5
#>                   <char>         <num>       <num>         <num>        <num>
#> 1:  Lee-Carter 1990-2002      4662.100    15269.50      24719.80     2.843114
#> 2:  Lee-Carter 1993-2002      5804.200    15274.50      23506.88     2.844072
#> 3: random walk 1990-2002      6604.257    15194.38      26120.47     2.828732
#> 4:      Li-Lee 1990-2002      -258.550    14395.00      28768.22     2.675929
```

The medians agree within 900 deaths. Li-Lee, whose common factor ties
the trends of women and men together, gives the widest interval of the
sum: when the trend errors of the parts are correlated, they do not
cancel in the total. Report the baseline model and years with the
results, and compare a few.

## A large shock: France 1918

The influenza pandemic of 1918 against a baseline of 1900-1913 (ages
0-94 and 95+, because some exposures at the oldest ages are zero in
these years):

``` r

y <- mt_regroup(mt_subset(fr, periods = c(1900, 1919)), age = 0:95)
e18 <- mt_excess(lee_carter(), y, baseline = c(1900, 1913), seed = 1)
e18[, list(sex, period, excess_q0.025, excess_q0.5, excess_q0.975, p_score_q0.5)]
#>        sex period excess_q0.025 excess_q0.5 excess_q0.975 p_score_q0.5
#>     <char>  <num>         <num>       <num>         <num>        <num>
#>  1: female   1914    -11119.275     17744.5      47573.07     6.274322
#>  2: female   1915    -21114.825     23712.0      62769.37     8.547329
#>  3: female   1916    -37333.100     18283.5      66942.35     6.812711
#>  4: female   1917    -40205.525     28591.5      77379.35    10.843196
#>  5: female   1918     37314.225    112526.5     172727.88    42.546285
#>  6: female   1919    -31413.825     45546.0     109408.95    17.512169
#>  7:   male   1914    278445.150    305674.5     329615.65    98.833110
#>  8:   male   1915    316801.950    351930.5     381905.85   117.832756
#>  9:   male   1916    202813.575    239698.5     274085.35    83.982005
#> 10:   male   1917    111154.650    159475.0     199189.42    57.155606
#> 11:   male   1918    220695.450    275881.0     322360.60    98.098397
#> 12:   male   1919     -8753.025     54239.5     103360.67    19.707532
```

``` r

mt_plot_excess(e18, type = "p_score")
```

![](excess-mortality_files/figure-html/plot-1918-1.png)

The excess of men starts in 1914, with the First World War, so for men
it mixes war deaths and the pandemic. That of women is concentrated in
1918: 113,000 deaths (95% interval 37,000 to 173,000), 43% more than
expected. The intervals widen from 1914 on, because the expected deaths
of 1918 are forecast five years ahead of the baseline. By age, for 1918:

``` r

e18a <- mt_excess(lee_carter(), y, baseline = c(1900, 1913), periods = 1918, ages = c(0, 15, 45, 65), seed = 1)
e18a[, list(sex, age, observed, excess_q0.5, p_score_q0.5)]
#>       sex    age observed excess_q0.5 p_score_q0.5
#>    <char> <char>    <num>       <num>        <num>
#> 1: female   0-14    56670     25833.5    83.775761
#> 2: female  15-44   102306     55792.0   119.946683
#> 3: female  45-64    71760     14789.5    25.959927
#> 4: female    65+   146271     15674.0    12.001807
#> 5:   male   0-14    60876     25979.0    74.444812
#> 6:   male  15-44   275811    226843.0   463.248014
#> 7:   male  45-64    94349     15915.0    20.290945
#> 8:   male    65+   126074      7537.5     6.358801
```

Women aged 15-44 died at more than twice the expected rate, a mark of
the 1918 pandemic, which was most lethal for young adults.

## Germany, 2020-2024, by Land

The German data of the package (Destatis) give deaths by Land, sex and
age group. A few small cells are suppressed;
[`mt_fill_suppressed()`](../reference/mt_fill_suppressed.md) estimates
them from the totals of Germany:

``` r

library(mortensor)
library(data.table)
de <- mt_data(deaths = de_deaths[de_deaths$cause == "all", ], population = de_population, by = c("region", "sex"))
de <- mt_as_exposure(de, type = "dec31")
de <- mt_regroup(mt_subset(de, periods = c(2010, 2024)), age = "deaths")
de <- mt_fill_suppressed(de, mt_hierarchy(region = "Deutschland"))
laender <- mt_subset(de, region != "Deutschland", drop = TRUE)
```

The populations are given by Land, so the national total is formed by
summing the Länder
([`mt_aggregate_strata()`](../reference/mt_aggregate_strata.md)). A
Lee-Carter model of the national total, both sexes, fitted to 2010-2019:

``` r

h <- mt_hierarchy(region = "Deutschland", sex = "both")
germany <- mt_subset(mt_aggregate_strata(laender, h), region == "Deutschland" & sex == "both", drop = TRUE)
ex_de <- mt_excess(lee_carter(), germany, baseline = c(2010, 2019), by = character(), seed = 1)
ex_de[, list(period, observed, excess_q0.025, excess_q0.5, excess_q0.975, p_score_q0.5, p_positive)]
#>    period observed excess_q0.025 excess_q0.5 excess_q0.975 p_score_q0.5
#>     <num>    <num>         <num>       <num>         <num>        <num>
#> 1:   2020   985572      -963.875     26354.5      54663.78     2.747500
#> 2:   2021  1023687     13959.925     57289.5      97750.07     5.928151
#> 3:   2022  1066341     47014.275    101341.5     145831.95    10.501715
#> 4:   2023  1028206      2275.975     66324.5     124618.85     6.895288
#> 5:   2024  1007758    -32997.600     30685.5      95076.57     3.140555
#>    p_positive
#>         <num>
#> 1:      0.970
#> 2:      0.996
#> 3:      1.000
#> 4:      0.978
#> 5:      0.828
mt_excess(lee_carter(), germany, baseline = c(2010, 2019), by = character(), cumulative = TRUE,
          seed = 1)[, list(from, period, excess_q0.025, excess_q0.5, excess_q0.975)]
#>     from period excess_q0.025 excess_q0.5 excess_q0.975
#>    <num>  <num>         <num>       <num>         <num>
#> 1:  2020   2020      -963.875     26354.5      54663.78
#> 2:  2020   2021     21745.725     85125.0     145962.10
#> 3:  2020   2022     72705.150    182873.5     284930.35
#> 4:  2020   2023     85289.000    250826.0     401462.32
#> 5:  2020   2024     61848.500    282996.5     483027.62
```

The excess rose from about 26,000 in 2020 to 101,000 in 2022 and fell
afterwards; from 2020 to 2022 it summed to 183,000, with a wide interval
(73,000 to 285,000). The interval is wide because the expected deaths of
2022 are a forecast three years beyond a baseline of ten years.

Summing separate models of the 32 strata (Land by sex) instead gives a
much narrower interval for the same total:

``` r

mt_excess(lee_carter(), laender, baseline = c(2010, 2019), by = character(),
          seed = 1)[, list(period, excess_q0.025, excess_q0.5, excess_q0.975)]
#>    period excess_q0.025 excess_q0.5 excess_q0.975
#>     <num>         <num>       <num>         <num>
#> 1:   2020     11068.050     19094.0      28479.55
#> 2:   2021     35794.700     48118.5      60482.52
#> 3:   2022     76635.350     91800.0     106571.65
#> 4:   2023     36244.150     54795.5      73615.62
#> 5:   2024      1853.475     22284.0      41176.32
```

Separate models draw their trend errors independently, so these cancel
in the sum, while in reality the Länder share much of their trend. For a
total, use a model of the total, or a joint model of its parts such as
[`li_lee()`](../reference/li_lee.md), whose common factor carries the
shared trend. The warning suppressed here concerns the oldest men of
Bremen, a small population.

P-scores of the Länder in 2021, from a Li-Lee model with one common
factor per sex:

``` r

ex_l <- mt_excess(li_lee(common_by = "sex"), laender, baseline = c(2010, 2019), by = "region", periods = 2021,
                  seed = 1)
ex_l[order(-p_score_q0.5), list(region, p_score_q0.025, p_score_q0.5, p_score_q0.975, cells, skipped)]
#>                     region p_score_q0.025 p_score_q0.5 p_score_q0.975 cells
#>                     <char>          <num>        <num>          <num> <int>
#>  1:              Thüringen      7.4532384   15.1557232      21.292852    36
#>  2:                Sachsen      6.3765040   13.6719614      19.846265    36
#>  3:         Sachsen-Anhalt      4.1907960   10.8027397      17.254436    36
#>  4:            Brandenburg      0.7964173    8.1966265      14.434801    36
#>  5: Mecklenburg-Vorpommern     -0.2892798    6.3543089      12.241962    36
#>  6:                 Bayern     -1.2681683    5.5483044      11.024289    36
#>  7:                 Hessen     -0.9902274    5.2028681      10.884268    36
#>  8:               Saarland     -1.7159657    4.5899119      10.228890    36
#>  9:                Hamburg     -2.7470839    4.0671515       9.742152    36
#> 10:                 Bremen     -2.4111913    3.8805974       9.832265    34
#> 11:                 Berlin     -2.9159538    3.8432410       9.440658    36
#> 12:      Baden-Württemberg     -2.9372918    3.6408126       9.007041    36
#> 13:    Nordrhein-Westfalen     -3.2942193    3.0946965       8.314086    36
#> 14:        Rheinland-Pfalz     -3.8786642    2.6646230       8.372199    36
#> 15:          Niedersachsen     -4.4369139    2.0238411       7.120832    36
#> 16:     Schleswig-Holstein     -5.9877643    0.2301982       5.290194    36
#>     skipped
#>       <int>
#>  1:       0
#>  2:       0
#>  3:       0
#>  4:       0
#>  5:       0
#>  6:       0
#>  7:       0
#>  8:       0
#>  9:       0
#> 10:       2
#> 11:       0
#> 12:       0
#> 13:       0
#> 14:       0
#> 15:       0
#> 16:       0
```

``` r

mt_plot_excess(ex_l, type = "p_score")
```

![](excess-mortality_files/figure-html/plot-laender-1.png)

In 2021 the excess was highest in Thüringen, Sachsen and Sachsen-Anhalt
(P-scores 11-15%); in most western Länder the interval contains zero.
`cells` and `skipped` count the observed cells used and left out: cells
whose deaths were estimated (suppressed cells filled above) are not
compared unless `include_estimated = TRUE`.

## Do the intervals hold?

When the baseline model describes the years without a shock, the 95%
interval of the excess should contain the true excess in 95% of cases.
The validation suite of the package checks this on simulated data
([`mt_example()`](../reference/mt_example.md), a Lee-Carter baseline,
deaths at ages 70 and older raised by 10% in two years, 20 replicates):
the intervals contained zero where there was no shock in 97.5% of the
cases, and the injected deaths in 98% of the cases where there was one.
On real data, the coverage depends on how well the baseline model
forecasts the years of the shock; a backtest of the baseline model on
years without a shock ([`mt_backtest()`](../reference/mt_backtest.md))
shows how wide its intervals need to be.

## References

Acosta RJ, Irizarry RA (2022). A flexible statistical framework for
estimating excess mortality. *Epidemiology* 33(3): 346-353.
<doi:10.1097/EDE.0000000000001445>.

Fouillet A, Rey G, Laurent F, Pavillon G, Bellec S, Guihenneuc-Jouyaux
C, Clavel J, Jougla E, Hémon D (2006). Excess mortality related to the
August 2003 heat wave in France. *International Archives of Occupational
and Environmental Health* 80(1): 16-24. <doi:10.1007/s00420-006-0089-4>.

Msemburi W, Karlinsky A, Knutson V, Aleshin-Guendel S, Chatterji S,
Wakefield J (2023). The WHO estimates of excess mortality associated
with the COVID-19 pandemic. *Nature* 613: 130-137.
<doi:10.1038/s41586-022-05522-2>.
