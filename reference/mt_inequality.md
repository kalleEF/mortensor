# Slope and relative index of inequality, path by path

Summarises a gradient over the ordered levels of a stratum variable (for
example deprivation quintiles) by a linear regression of the quantity on
the relative rank of each level in the population, weighted by
population size (Mackenbach and Kunst 1997). Moreno-Betancur et al.
(2015) define the slope and the relative index of inequality as
parameters of a regression of the outcome on the relative rank: the
difference and the ratio of the expected outcomes at the two ends of the
rank scale. Here, with the linear regression (identity link), the slope
index of inequality (SII) is the fitted value at the top (rank 1) minus
that at the bottom (rank 0), the slope; the relative index of inequality
(RII) is the ratio of the fitted values of the same linear regression at
rank 1 and rank 0, `(b0 + b1) / b0` with intercept `b0` and slope `b1`
(not the exponentiated slope of a regression with a log link). The
relative rank of a level is the midpoint of its cumulative population
share, with the levels ordered as in `order`. Computed path by path,
with the population shares of the exposures of each year for observed
rates and of the last observed year for forecasts.

## Usage

``` r
mt_inequality(
  x,
  var,
  order = NULL,
  quantity = c("asr", "life_expectancy", "prob_dying", "rate", "lifespan_variation"),
  probs = c(0.025, 0.5, 0.975),
  ...
)
```

## Arguments

- x:

  An [mt_draws](mt_draws.md) or [mt_data](mt_data.md) object.

- var:

  Ordered stratum variable (for example `"quintile"`).

- order:

  Levels of `var` from rank 0 to rank 1 (default: their order in the
  data, without the total and the groups of a hierarchy stored in `x` by
  [`mt_aggregate_strata()`](mt_aggregate_strata.md),
  [`reconciled()`](reconciled.md) or `mt_group_levels(total = )`). For
  deprivation, list the least deprived level first so that a positive
  SII means higher values in more deprived groups. Levels not listed,
  such as the total of a hierarchy, are left out.

- quantity:

  `"life_expectancy"`, `"asr"`, `"prob_dying"`, `"rate"` or
  `"lifespan_variation"`.

- probs:

  Quantile levels.

- ...:

  Arguments of the quantity: `at` for life expectancy, `standard`,
  `per`, `exposure` and `age_groups` for [`asr()`](derived.md), `from`
  and `to` for [`prob_dying()`](derived.md), `age` (an age label) for
  `"rate"`, `measure` (one index, default `"edagger"`) and `at` for
  [`lifespan_variation()`](lifespan_variation.md).

## Value

A `data.table` with the other stratum variables, `index` (`"sii"` or
`"rii"`), `period`, `mean` and quantiles.

## References

Mackenbach JP, Kunst AE (1997). Measuring the magnitude of
socio-economic inequalities in health: an overview of available measures
illustrated with two examples from Europe. *Social Science & Medicine*
44(6): 757-771.

Moreno-Betancur M, Latouche A, Menvielle G, Kunst AE, Rey G (2015).
Relative index of inequality and slope index of inequality: a structured
regression framework for estimation. *Epidemiology* 26(4): 518-527.

## Examples

``` r
x <- mt_example_regions(years = 2010:2019, regions = 2)
mt_inequality(x, "quintile", order = paste0("Q", 5:1), quantity = "asr")
#>        sex region  index period          mean        q0.025          q0.5
#>     <char> <char> <char>  <num>         <num>         <num>         <num>
#>  1: female     R1    sii   2010  -813.1256544  -845.4876158  -813.9808547
#>  2: female     R1    sii   2011  -823.4020314  -855.0685687  -823.1995401
#>  3: female     R1    sii   2012  -846.2667159  -877.9198136  -845.8553194
#>  4: female     R1    sii   2013  -779.0186946  -811.5352950  -779.2561265
#>  5: female     R1    sii   2014  -796.2281897  -833.2032070  -796.4263111
#>  6: female     R1    sii   2015  -782.1403142  -811.6118573  -782.2393719
#>  7: female     R1    sii   2016  -809.7202770  -843.0385963  -808.7848328
#>  8: female     R1    sii   2017  -833.4852068  -866.6315060  -833.4392795
#>  9: female     R1    sii   2018  -846.4930033  -876.8483770  -846.1142497
#> 10: female     R1    sii   2019  -863.3234237  -891.6693884  -864.4822846
#> 11: female     R1    rii   2010     0.5461505     0.5322274     0.5456264
#> 12: female     R1    rii   2011     0.5358881     0.5231633     0.5360637
#> 13: female     R1    rii   2012     0.5220938     0.5088161     0.5222567
#> 14: female     R1    rii   2013     0.5342469     0.5199507     0.5342163
#> 15: female     R1    rii   2014     0.5326151     0.5167819     0.5324500
#> 16: female     R1    rii   2015     0.5365736     0.5234315     0.5366370
#> 17: female     R1    rii   2016     0.5254045     0.5108581     0.5255549
#> 18: female     R1    rii   2017     0.5177664     0.5031492     0.5173148
#> 19: female     R1    rii   2018     0.5114061     0.4985839     0.5112631
#> 20: female     R1    rii   2019     0.5076578     0.4957194     0.5074201
#> 21: female     R2    sii   2010  -860.9319614  -930.0490379  -863.5406460
#> 22: female     R2    sii   2011  -851.7015469  -919.0138171  -850.9937453
#> 23: female     R2    sii   2012  -888.7949409  -955.3337850  -889.2619785
#> 24: female     R2    sii   2013  -872.5227755  -938.1700840  -870.6554141
#> 25: female     R2    sii   2014  -905.6683995  -962.3913300  -906.8963030
#> 26: female     R2    sii   2015  -891.7330673  -958.8089958  -892.5342184
#> 27: female     R2    sii   2016  -840.9410327  -906.6258574  -840.8722113
#> 28: female     R2    sii   2017  -861.8172792  -923.4897951  -861.3518968
#> 29: female     R2    sii   2018  -895.8191786  -961.0201528  -896.4476401
#> 30: female     R2    sii   2019  -817.5199249  -881.6330471  -817.0246879
#> 31: female     R2    rii   2010     0.5491475     0.5211052     0.5480671
#> 32: female     R2    rii   2011     0.5491745     0.5237685     0.5496910
#> 33: female     R2    rii   2012     0.5328002     0.5064995     0.5322818
#> 34: female     R2    rii   2013     0.5430576     0.5176165     0.5436813
#> 35: female     R2    rii   2014     0.5105880     0.4875032     0.5101321
#> 36: female     R2    rii   2015     0.5149813     0.4877047     0.5143841
#> 37: female     R2    rii   2016     0.5441732     0.5181653     0.5441188
#> 38: female     R2    rii   2017     0.5256042     0.5010630     0.5255880
#> 39: female     R2    rii   2018     0.5093281     0.4842775     0.5092825
#> 40: female     R2    rii   2019     0.5260878     0.4991492     0.5259044
#> 41:   male     R1    sii   2010 -1340.7585020 -1381.9928478 -1341.8739374
#> 42:   male     R1    sii   2011 -1370.8057999 -1419.2184686 -1371.4053685
#> 43:   male     R1    sii   2012 -1372.4596921 -1416.2128780 -1372.8447232
#> 44:   male     R1    sii   2013 -1383.9255423 -1430.2325960 -1383.7989540
#> 45:   male     R1    sii   2014 -1325.6162722 -1368.5785984 -1325.7925874
#> 46:   male     R1    sii   2015 -1270.2046582 -1319.0288937 -1268.7124348
#> 47:   male     R1    sii   2016 -1268.4991509 -1316.8621501 -1267.9912588
#> 48:   male     R1    sii   2017 -1328.4187987 -1367.9534374 -1328.8939011
#> 49:   male     R1    sii   2018 -1293.1476623 -1332.2741537 -1294.3387607
#> 50:   male     R1    sii   2019 -1304.0622636 -1346.8911245 -1304.7918892
#> 51:   male     R1    rii   2010     0.5355597     0.5245694     0.5353503
#> 52:   male     R1    rii   2011     0.5387278     0.5269436     0.5386644
#> 53:   male     R1    rii   2012     0.5348445     0.5238705     0.5349149
#> 54:   male     R1    rii   2013     0.5281488     0.5164459     0.5280926
#> 55:   male     R1    rii   2014     0.5307790     0.5190121     0.5306158
#> 56:   male     R1    rii   2015     0.5292821     0.5160826     0.5294863
#> 57:   male     R1    rii   2016     0.5412943     0.5270927     0.5413442
#> 58:   male     R1    rii   2017     0.5100232     0.4993929     0.5100441
#> 59:   male     R1    rii   2018     0.5125362     0.5012702     0.5122770
#> 60:   male     R1    rii   2019     0.5075545     0.4958366     0.5071810
#> 61:   male     R2    sii   2010 -1451.6965718 -1547.2666672 -1453.7116663
#> 62:   male     R2    sii   2011 -1393.8954402 -1492.3513617 -1394.6621612
#> 63:   male     R2    sii   2012 -1453.2759451 -1551.9522023 -1454.1223787
#> 64:   male     R2    sii   2013 -1433.8389560 -1529.3652534 -1431.7481827
#> 65:   male     R2    sii   2014 -1390.8207941 -1491.2665848 -1387.8416111
#> 66:   male     R2    sii   2015 -1479.2440084 -1570.3759459 -1479.9637190
#> 67:   male     R2    sii   2016 -1305.4177315 -1396.8460874 -1306.1229150
#> 68:   male     R2    sii   2017 -1429.9325147 -1526.3932870 -1429.7293780
#> 69:   male     R2    sii   2018 -1367.2259156 -1456.8708419 -1365.8872436
#> 70:   male     R2    sii   2019 -1552.7698282 -1641.9578798 -1552.7519569
#> 71:   male     R2    rii   2010     0.5297382     0.5066733     0.5295156
#> 72:   male     R2    rii   2011     0.5389256     0.5147631     0.5395640
#> 73:   male     R2    rii   2012     0.5253062     0.5020753     0.5255032
#> 74:   male     R2    rii   2013     0.5272505     0.5032857     0.5271738
#> 75:   male     R2    rii   2014     0.5265515     0.5007919     0.5271275
#> 76:   male     R2    rii   2015     0.5081126     0.4856063     0.5076150
#> 77:   male     R2    rii   2016     0.5395917     0.5155057     0.5396213
#> 78:   male     R2    rii   2017     0.5173961     0.4924136     0.5175333
#> 79:   male     R2    rii   2018     0.5338159     0.5090549     0.5342850
#> 80:   male     R2    rii   2019     0.4788511     0.4579897     0.4789011
#>        sex region  index period          mean        q0.025          q0.5
#>     <char> <char> <char>  <num>         <num>         <num>         <num>
#>            q0.975 quantity
#>             <num>   <char>
#>  1:  -778.3241430      asr
#>  2:  -789.5933766      asr
#>  3:  -814.8970092      asr
#>  4:  -748.4308464      asr
#>  5:  -763.1550175      asr
#>  6:  -753.4625835      asr
#>  7:  -780.8311241      asr
#>  8:  -800.3484721      asr
#>  9:  -814.4136063      asr
#> 10:  -829.0282828      asr
#> 11:     0.5605077      asr
#> 12:     0.5501327      asr
#> 13:     0.5359602      asr
#> 14:     0.5474949      asr
#> 15:     0.5471287      asr
#> 16:     0.5498628      asr
#> 17:     0.5391299      asr
#> 18:     0.5330180      asr
#> 19:     0.5253130      asr
#> 20:     0.5221797      asr
#> 21:  -788.8112392      asr
#> 22:  -783.3913085      asr
#> 23:  -823.3340965      asr
#> 24:  -807.8896900      asr
#> 25:  -841.3252061      asr
#> 26:  -825.5026222      asr
#> 27:  -771.6379164      asr
#> 28:  -798.3948705      asr
#> 29:  -830.5884501      asr
#> 30:  -750.5694728      asr
#> 31:     0.5787146      asr
#> 32:     0.5754077      asr
#> 33:     0.5578405      asr
#> 34:     0.5701303      asr
#> 35:     0.5357334      asr
#> 36:     0.5433206      asr
#> 37:     0.5736705      asr
#> 38:     0.5518164      asr
#> 39:     0.5374529      asr
#> 40:     0.5552831      asr
#> 41: -1290.1699494      asr
#> 42: -1320.5334680      asr
#> 43: -1323.6218473      asr
#> 44: -1337.6360511      asr
#> 45: -1276.3168849      asr
#> 46: -1220.4927956      asr
#> 47: -1223.6709785      asr
#> 48: -1285.4083966      asr
#> 49: -1250.7219389      asr
#> 50: -1265.3713375      asr
#> 51:     0.5490372      asr
#> 52:     0.5514010      asr
#> 53:     0.5469309      asr
#> 54:     0.5407971      asr
#> 55:     0.5437287      asr
#> 56:     0.5427190      asr
#> 57:     0.5535703      asr
#> 58:     0.5213386      asr
#> 59:     0.5245542      asr
#> 60:     0.5192530      asr
#> 61: -1358.3626570      asr
#> 62: -1298.1725472      asr
#> 63: -1353.7038810      asr
#> 64: -1337.9113877      asr
#> 65: -1299.5342761      asr
#> 66: -1384.6798314      asr
#> 67: -1213.5882463      asr
#> 68: -1333.9936999      asr
#> 69: -1271.2379780      asr
#> 70: -1457.4973606      asr
#> 71:     0.5527591      asr
#> 72:     0.5618838      asr
#> 73:     0.5504849      asr
#> 74:     0.5527491      asr
#> 75:     0.5524176      asr
#> 76:     0.5318742      asr
#> 77:     0.5644056      asr
#> 78:     0.5419475      asr
#> 79:     0.5579053      asr
#> 80:     0.5017499      asr
#>            q0.975 quantity
#>             <num>   <char>
```
