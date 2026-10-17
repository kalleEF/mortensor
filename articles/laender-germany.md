# Differences between the German Laender

*This vignette is precomputed (`dev/precompute-vignettes.R`) because its
backtests take about half an hour.*

How will life expectancy in the 16 German Laender develop until 2035,
and will the remaining gap between East and West Germany close? This
vignette answers with one complete pipeline: official tables with
suppressed cells, a hierarchy of Laender, East and West and Germany, a
model chosen in a backtest of both the Laender and the East-West gap,
and forecasts that add up at every level. Each step refers to the
vignette that explains it.

``` r

library(mortensor)
library(data.table)
```

## Data

`de_deaths` holds deaths by year, Land, sex, age group and ICD-10
chapter (Destatis, published by the Federal Health Monitoring,
GBE-Bund); `de_population` the populations on 31 December by single year
of age (GENESIS-Online, table 12411-0013). Both are licensed under the
Datenlizenz Deutschland - Namensnennung - Version 2.0.

``` r

x <- mt_data(deaths = de_deaths, population = de_population, by = c("region", "sex", "cause"))
#> Warning in ingest_obs(x, by, age, period, value, measure, age_open, missing =
#> missing, : 41396 missing deaths value(s) left out; suppressed counts can be
#> estimated with mt_fill_suppressed().
x <- mt_as_exposure(x, type = "dec31")
x <- mt_subset(x, periods = c(1998, 2019))
x <- mt_regroup(x, age = "deaths")
```

- **Exposures**: the person-years of year *t* are the mean of the
  populations on 31 December of *t* - 1 and *t* (`type = "dec31"`).
- **Years**: 1998-2019, before the COVID-19 pandemic.
- **Ages**: the age groups of the deaths (0, 1-14, 15-19, …, 90+).

**Suppressed cells.** The warning says that the source suppresses small
counts. The published totals constrain them: the Laender add up to
Germany and the chapters to all causes.
[`mt_fill_suppressed()`](../reference/mt_fill_suppressed.md) estimates
the suppressed counts by raking them to these totals (“Your own data”):

``` r

h_fill <- mt_hierarchy(region = "Deutschland", cause = list(total = "all", additive = "rates"))
x <- mt_fill_suppressed(x, h_fill)
attr(x, "fill")
#>    filled missing inconsistent iterations      max_gap
#> 1:  33225       0            0       1959 7.975615e-07
```

All totals are matched (`max_gap`). For all causes only a few cells at
young ages were suppressed, so the estimates matter little here.

## The hierarchy

The 16 Laender by sex are the bottom series. East (with Berlin) and West
Germany, Germany, and both sexes are their sums. All forecasts below are
made for the 32 bottom series and summed
(`reconciled(method = "bottom_up")`, “Coherent forecasts: totals and
their parts”):

``` r

xa <- mt_subset(x, cause == "all", region != "Deutschland", drop = TRUE)
laender <- unique(mt_strata(xa)$region)
east <- c("Berlin", "Brandenburg", "Mecklenburg-Vorpommern", "Sachsen", "Sachsen-Anhalt",
          grep("^Th", laender, value = TRUE))
h <- mt_hierarchy(region = list(total = "Deutschland", groups = list(West = setdiff(laender, east), East = east)),
                  sex = "both")
xh <- mt_aggregate_strata(xa, h)
mt_strata(xh, h)[, .N, by = level]
#>        level  N
#> 1:    bottom 32
#> 2: aggregate 24
#> 3:     total  1
bu <- function(model) reconciled(model, h, method = "bottom_up")
```

The observed difference in life expectancy at birth between East and
West:

``` r

gap_obs <- mt_contrast(xh, var = "region", ref = "West")
mt_plot_contrast(observed = gap_obs, levels = "East")
```

![](laender-germany-plot-gap-observed-1.png)

Women in the East caught up with the West by about 2014 and lived 0.1
years longer in 2019. Men in the East still live about 1.2 years less
(1.7 years in 1998).

## Model: a backtest of candidates

Small Laender such as Bremen and Saarland have few deaths per age group,
so their rates are noisy. The candidates differ in how the Laender may
move relative to each other (“Choosing and combining models”):

- **Lee-Carter** for each Land and sex: each Land continues its own
  trend;
- **Li-Lee** and the **coherent fdm** (product-ratio method with
  functional data models): the Laender of each sex share a common trend
  and deviate from it with stationary factors, so small Laender borrow
  strength from the others and gaps return towards their past level;
- the **ensemble** of Lee-Carter and the coherent fdm.

The East-West gap forecast by three of them:

``` r

cand <- list("Lee-Carter" = bu(lee_carter()), "Li-Lee" = bu(li_lee(common_by = "sex")),
             "coherent fdm" = bu(coherent_fdm(common_by = "sex")),
             "Lee-Carter + coherent fdm" = bu(ensemble(lc = lee_carter(), cf = coherent_fdm(common_by = "sex"))))
gaps <- lapply(cand[c("Lee-Carter", "Li-Lee", "coherent fdm")], function(m) {
  mt_contrast(mt_forecast(mt_fit(m, xh), to = 2035, n_draws = 300, seed = 1), var = "region", ref = "West", at = 0)
})
mt_plot_contrast(gaps, observed = gap_obs, levels = "East")
```

![](laender-germany-plot-gap-models-1.png)

The backtest fits each candidate to the years up to 2008, 2011 and 2014,
forecasts the following five years, and scores the rates and life
expectancy of every series and the gap between each region and the West
(`contrast`, “Backtesting and scores”):

``` r

bt <- mt_compare(cand, xh, origins = c(2008, 2011, 2014), h = 5, n_draws = 200, quantities = "life_expectancy",
                 contrast = c(var = "region", ref = "West"))
mt_score_table(bt, metrics = c("crps", "cov95"), targets = c("rate", "life_expectancy"), hierarchy = h)
#>                         model     level   crps cov95 crps_le cov95_le
#>  1:                Lee-Carter    bottom 0.0915 0.814   0.374    0.776
#>  2:                Lee-Carter aggregate 0.0665 0.692   0.361    0.533
#>  3:                Lee-Carter     total 0.0477 0.437   0.360    0.267
#>  4:                    Li-Lee    bottom 0.0894 0.851   0.392    0.786
#>  5:                    Li-Lee aggregate 0.0647 0.764   0.362    0.700
#>  6:                    Li-Lee     total 0.0502 0.600   0.343    0.533
#>  7:              coherent fdm    bottom 0.0817 0.839   0.280    0.682
#>  8:              coherent fdm aggregate 0.0591 0.727   0.250    0.553
#>  9:              coherent fdm     total 0.0458 0.470   0.236    0.333
#> 10: Lee-Carter + coherent fdm    bottom 0.0815 0.872   0.300    0.795
#> 11: Lee-Carter + coherent fdm aggregate 0.0583 0.754   0.281    0.583
#> 12: Lee-Carter + coherent fdm     total 0.0443 0.481   0.290    0.333
```

``` r

mt_score_table(bt[region == "East" & sex != "both"], metrics = c("crps", "cov95"), targets = "life_expectancy_gap")
#>                        model crps_le_gap cov95_le_gap
#> 1:                Lee-Carter      0.0761        1.000
#> 2:                    Li-Lee      0.1550        0.400
#> 3:              coherent fdm      0.1520        0.467
#> 4: Lee-Carter + coherent fdm      0.0754        1.000
```

``` r

mt_plot_scores(bt[region %in% laender & sex != "both"], by = "model", target = "life_expectancy")
```

![](laender-germany-plot-scores-1.png)

Lee-Carter continues the trends of 1998-2019: the gap of men keeps
narrowing, and women in the East pull further ahead. Li-Lee and the
coherent fdm hold the gaps near their recent level, with women in the
East slightly behind again. The backtest decides between these
assumptions:

- For the Laender, the coherent fdm and the ensemble are the most
  accurate (CRPS of the rates 0.082, against 0.089 for Li-Lee and 0.092
  for Lee-Carter); the ensemble’s intervals cover the rates best (95%
  coverage 0.87). The intervals of life expectancy are too narrow for
  every model (95% coverage 0.68 to 0.80).
- For the totals (East, West and Germany), the ensemble forecasts the
  rates most accurately, but all intervals are too narrow there.
- The East-West gap is forecast best by the ensemble and Lee-Carter
  (CRPS 0.075 and 0.076), whose intervals covered every observed gap;
  those of Li-Lee and the coherent fdm covered it in only 40% and 47% of
  the cases.

The ensemble of Lee-Carter and the coherent fdm is the most accurate for
the Laender and, together with Lee-Carter, for the gap, so the pipeline
continues with it: combining models pays here.

### Jump-off

Starting the forecasts at the observed level of the last three years
([`anchored()`](../reference/anchored.md), “Where a forecast starts: the
jump-off”), checked in the same backtest:

``` r

ens <- ensemble(lc = lee_carter(), cf = coherent_fdm(common_by = "sex"))
bt_anchor <- mt_compare(list("+ anchored" = bu(anchored(ens, years = 3))), xh, origins = c(2008, 2011, 2014),
                        h = 5, n_draws = 200, quantities = "life_expectancy", contrast = c(var = "region", ref = "West"))
bt_steps <- rbind(bt[model == "Lee-Carter + coherent fdm"], bt_anchor)
mt_score_table(bt_steps[sex != "both"], metrics = c("crps", "cov95"))
#>                        model   crps cov95 crps_le cov95_le crps_le_gap
#> 1: Lee-Carter + coherent fdm 0.0766 0.835   0.301    0.730       0.108
#> 2:                + anchored 0.0759 0.839   0.293    0.738       0.104
#>    cov95_le_gap
#> 1:        0.996
#> 2:        1.000
```

Anchoring made every quantity slightly more accurate (CRPS of the rates
0.0759 against 0.0766), so the pipeline uses it. Smoothing the rates
over age made the forecasts of this ensemble slightly worse in a
separate run and is not used.

## Forecasts to 2035

``` r

fc <- mt_forecast(mt_fit(bu(anchored(ens, years = 3)), xh), to = 2035, n_draws = 1000, seed = 1)
```

Life expectancy at birth in 2035, both sexes, with 95% intervals:

``` r

e35 <- life_expectancy(fc)[period == 2035 & sex == "both"]
mt_plot_compare(list(forecast = fc), periods = 2035, strata = rev(e35[order(q0.5), label]))
```

![](laender-germany-plot-laender-1.png)

The intervals of small Laender are the widest. The paths of the Laender
add up to those of East, West and Germany, so these intervals are
consistent with each other.

## The East-West gap

The difference in life expectancy at birth between East and West
Germany, path by path, with the observed differences:

``` r

gap <- mt_contrast(fc, var = "region", ref = "West", at = 0)
mt_plot_contrast(gap, observed = gap_obs, levels = "East")
```

![](laender-germany-plot-gap-1.png)

``` r

gap[level == "East" & period %in% c(2025, 2035),
    .(sex, period, gap = sprintf("%.2f (%.2f to %.2f)", q0.5, q0.025, q0.975))]
#>       sex period                    gap
#> 1: female   2025   0.10 (-0.29 to 0.66)
#> 2: female   2035   0.21 (-0.45 to 1.13)
#> 3:   male   2025 -1.09 (-1.54 to -0.58)
#> 4:   male   2035 -0.91 (-1.68 to -0.01)
#> 5:   both   2025 -0.54 (-0.85 to -0.08)
#> 6:   both   2035  -0.40 (-0.92 to 0.42)
```

In 2035, women in the East are forecast to live 0.2 years longer than
those in the West (95% interval -0.45 to 1.13), men 0.9 years less
(-1.68 to -0.01). The forecast continues the slow narrowing of the male
gap, but the gap most likely does not close by 2035.

Which ages make up the gap of men in 2019 and in 2035 (“Decomposing
differences in life expectancy”):

``` r

dg <- mt_decompose(fc, var = "region", ref = "West", periods = c(2019, 2035))
mt_plot_decomposition(dg, by = "age", comparisons = list(to = "East", sex = "male"))
```

![](laender-germany-plot-decomp-1.png)

The male gap arises at ages 45 to 74, in 2019 as in 2035; at older ages
the contributions are small and uncertain.

## Summary of the pipeline

| Step | Decision here | Evidence | Function | Explained in |
|----|----|----|----|----|
| Exposures | mean of 31 December populations |  | `mt_as_exposure(type = "dec31")` | “Your own data” |
| Suppressed cells | raking to Land and national totals | all totals matched | [`mt_fill_suppressed()`](../reference/mt_fill_suppressed.md) | “Your own data” |
| Hierarchy | Laender, East and West, Germany; sexes |  | [`mt_hierarchy()`](../reference/mt_hierarchy.md), [`mt_aggregate_strata()`](../reference/mt_aggregate_strata.md) | “Coherent forecasts” |
| Coherence | bottom-up sums of the Laender |  | `reconciled(method = "bottom_up")` | “Coherent forecasts” |
| Model | ensemble of Lee-Carter and coherent fdm | most accurate for the Laender and for the gap | [`ensemble()`](../reference/ensemble.md), [`coherent_fdm()`](../reference/coherent_fdm.md), `mt_compare(contrast = )` | “Choosing and combining models” |
| Jump-off | observed level of the last three years | slightly better scores | [`anchored()`](../reference/anchored.md) | “Where a forecast starts” |
| Gap | East - West, path by path |  | [`mt_contrast()`](../reference/mt_contrast.md), [`mt_decompose()`](../reference/mt_decompose.md) | “Decomposing differences” |

## References

Hyndman RJ, Booth H, Yasmeen F (2013). Coherent mortality forecasting:
the product-ratio method with functional time series models.
*Demography* 50(1): 261-283.
