# Lee-Carter model

`log m[x, t] = a[x] + b[x] k[t]` with the constraints `sum(b) = 1` and
`sum(k) = 0`, fitted per stratum:

## Usage

``` r
lee_carter(
  method = c("poisson", "svd"),
  jump_off = c("fitted", "observed"),
  bootstrap = FALSE,
  n_boot = 100,
  adjust = c("dt", "dxt", "e0", "none"),
  choose_period = FALSE,
  min_years = 20,
  rotate = FALSE,
  rotation = c(e0_start = 80, e0_end = 102, power = 0.5),
  family = c("poisson", "nb")
)
```

## Arguments

- method:

  `"poisson"` (maximum likelihood with the distribution of `family`) or
  `"svd"`.

- jump_off:

  `"fitted"` (forecasts start from the fitted rates of the last year) or
  `"observed"` (from the observed rates of the last year; Lee and Miller
  2001).

- bootstrap:

  Semi-parametric bootstrap of the parameters.

- n_boot:

  Number of bootstrap fits.

- adjust:

  Re-estimation of `k[t]` for `method = "svd"`: `"dt"`, `"dxt"`, `"e0"`
  or `"none"`.

- choose_period:

  Choose the fitting period (Booth, Maindonald and Smith 2002).

- min_years:

  `choose_period` considers fitting periods of more than `min_years`
  years (`minperiod` of
  [`demography::lca()`](https://pkg.robjhyndman.com/demography/reference/lca.html)).

- rotate:

  Rotate the age pattern of decline `b` towards an ultimate pattern as
  life expectancy at birth rises (Li, Lee and Gerland 2013).

- rotation:

  Settings of the rotation, a named vector: life expectancy at birth at
  which the rotation starts (`e0_start`, default 80) and is complete
  (`e0_end`, default 102), and the exponent of the weights (`power`,
  default 0.5). The defaults are those of Li, Lee and Gerland (2013,
  eq. 2) and MortCast; elements left out keep their defaults.

- family:

  Distribution of the deaths for `method = "poisson"`: `"poisson"`
  (Brouhns, Denuit and Vermunt 2002) or `"nb"` (negative binomial;
  Delwarde, Denuit and Partrat 2007).

## Value

A model specification.

## Details

- `method = "poisson"`: Poisson log-bilinear model (Brouhns, Denuit and
  Vermunt 2002), maximum likelihood by Newton updates;

- `method = "svd"`: original method (Lee and Carter 1992): singular
  value decomposition of centred log rates, with `k[t]` re-estimated in
  a second stage (`adjust`, as in
  [`demography::lca()`](https://pkg.robjhyndman.com/demography/reference/lca.html)):
  so that the fitted deaths match the observed total deaths of each year
  (`"dt"`, Lee and Carter 1992), by Poisson maximum likelihood of the
  deaths by age given `a` and `b` (`"dxt"`, Booth, Maindonald and Smith
  2002), so that the fitted life expectancy at the first age matches the
  observed one (`"e0"`, Lee and Miller 2001), or not at all (`"none"`).

`choose_period = TRUE` fits the model to the last years only, starting
in the year that minimises the ratio of the mean deviance of a fit with
a linear `k[t]` to the mean deviance of the Lee-Carter fit (Booth,
Maindonald and Smith 2002): the longest recent period over which the
decline was close to linear
(`demography::lca(chooseperiod = TRUE, breakmethod = "bms")` gives the
same period). Fitted rates of earlier years are missing.

Compared with
[`demography::lca()`](https://pkg.robjhyndman.com/demography/reference/lca.html),
`k[t]` is centred again after the second stage (a constant shift of `k`
and `a`, with identical fitted rates), and life tables use `a_x = n/2`
at age 0 (as [`life_table()`](life_table.md)), so `adjust = "e0"`
matches a slightly different life expectancy than demography's life
table where age 0 is in the data.

The period index follows a random walk with drift. Draws include the
innovations and the uncertainty of the drift. With `bootstrap = TRUE`
the parameters are re-estimated on `n_boot` Poisson resamples of the
fitted deaths of the fitting period (semi-parametric bootstrap with
fitted deaths of Renshaw and Haberman 2008), which adds the uncertainty
of `a` and `b`. The bootstrap of Brouhns, Denuit and Van Keilegom (2005)
resamples the observed deaths instead:
`bootstrapped(lee_carter(), type = "data")`
([`bootstrapped()`](bootstrapped.md)).

Age groups without deaths in the fitting period get a constant rate (0.5
deaths over the period) and no trend; with `adjust = "e0"` they enter
the life tables at that rate.

**Safeguards of this package for sparse counts.** Where the Poisson
likelihood has no finite maximum (`method = "poisson"` with very few
deaths), the fit falls back to two package-specific estimators
(reference unknown), each with a warning (class `mt_warning_model_fit`;
within [`mt_fit()`](mt_fit.md) it is reported with the stratum): first
to the Poisson fit with a ridge penalty that shrinks `b` towards equal
shares `1 / n_ages`, and if that also fails to the SVD estimates
(`method = "svd"` with `adjust = "dt"`). Their forecasts are not those
of the Poisson Lee-Carter model; wider age groups
([`mt_regroup()`](mt_regroup.md)) or pooling ([`li_lee()`](li_lee.md))
avoid them.

Draws describe the expected rates. Predictive draws
(`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
given the rates and the exposures, the observation model of the Poisson
Lee-Carter model (Brouhns, Denuit and Vermunt 2002); no other
observation noise is added. If intervals are too narrow on real data,
that is a property of the model: see
[`recalibrated()`](recalibrated.md), [`ensemble()`](ensemble.md), or
overdispersed models such as `family = "nb"` below or
`stanmomo("lc", family = "nb")` ([`stanmomo()`](stanmomo.md)).

**Negative binomial version** (`family = "nb"`; Delwarde, Denuit and
Partrat 2007). Following the paper:

- the deaths are negative binomial with mean
  `mu = E exp(a[x] + b[x] k[t])` and variance `mu + phi mu^2`, a Poisson
  distribution mixed over a gamma distributed random effect of each cell
  (Section 2.2). The dispersion `phi` (the paper's `k`) is one parameter
  for all ages and years;

- `a`, `b`, `k` and `phi` are estimated by maximum likelihood (eq. 5):
  uni-dimensional Newton updates of `a`, `b` and `k` (eq. 6-8), each
  given the other parameters, alternated with the estimate of `phi`
  given the fitted deaths;

- `k[t]` follows a random walk with drift (Section 3.5);

- predictive draws (`mt_forecast(predictive = TRUE)`) are negative
  binomial deaths with the estimated dispersion (size `1 / phi`), the
  observation model of Section 2.2.

Package-specific (reference unknown): the Newton updates start from the
Poisson estimates (the paper: `a = 0`, `b = 1`, `k = 0`, or random
values); `phi` is maximised by Newton steps on `log(phi)` with step
halving (the paper: Newton steps on `1 / phi`, eq. 9); the estimate is 0
(the Poisson model) where the derivative of the log-likelihood at
`phi = 0` is not positive; the log-likelihood uses gamma functions, so
deaths need not be whole numbers. The paper states no identification
constraints; `sum(b) = 1` and `sum(k) = 0` of the Poisson model
(Brouhns, Denuit and Vermunt 2002) agree with its Figures 2 and 3. As
for the Poisson model, draws include the innovations of `k` and the
uncertainty of its drift; the paper reports point forecasts only. With
`bootstrap = TRUE` the resampled deaths are negative binomial with the
estimated dispersion, and each draw has the dispersion of its bootstrap
fit (package-specific). The ridge and SVD fallbacks are not used: a fit
that does not converge fails. Like other maximum likelihood estimates of
a variance, `phi` is biased downwards when many parameters are estimated
relative to the number of cells. For France, ages 0-99, the fit to the
Human Mortality Database series in
[`demography::fr.mort`](https://pkg.robjhyndman.com/demography/reference/fr.mort.html)
reproduces the paper's dispersions (Section 3.3), AIC (Section 3.4) and
drifts and standard deviations of `k` (Sections 3.5 and 3.6) within 1%;
the paper's column labelled `sigma^2` holds the standard deviation
(divisor `T - 1`).

`rotate = TRUE` rotates the age pattern of decline as life expectancy
rises (Li, Lee and Gerland 2013; status: published, with the deviations
listed below). In long-lived populations the decline of mortality slows
down at young ages and speeds up at old ages, which a fixed `b[x]` does
not capture over long horizons. For each path and year `t`, the life
expectancy at birth `e0` of the unrotated forecast gives the weight
`w = (0.5 (1 + sin(pi / 2 (2 u - 1))))^power` with
`u = (e0 - e0_start) / (e0_end - e0_start)`; `w` is 0 below `e0_start`
and 1 from `e0_end` (eq. 2 of the paper). The age pattern of that year
is `B[x, t] = (1 - w) b[x] + w bu[x]` (eq. 3). The defaults are those of
the paper: the rotation starts at 80 (about the median life expectancy
of the 20 populations with the highest life expectancy in 2005-2010), is
complete at 102 (the value with which the ratio of the death rates at
age 0 and at ages 15-19 approaches its estimated long-run level of 7.7),
and the exponent is 0.5. The ultimate pattern `bu` has the mean of `b`
over ages 15-64 at all ages below 65, as in the paper. From age 65 it
keeps the shape of `b`, scaled so that the first age group from 65
equals that mean, and it sums to 1 like `b`: the paper only notes that
lowering `b` at young ages raises it at older ages because `b` sums to
1, and this scaling is that of MortCast (Sevcikova et al. 2016). The
index is then solved again, path by path and year by year, so that the
rotated rates have the same life expectancy at birth as the unrotated
ones (eq. 4): rotation changes the age pattern of the forecast and keeps
its life expectancy at birth. Decline thus moves from young ages to the
ages where `bu` exceeds `b`. Where `b` at 65 is already above its mean
over ages 15-64 (as for French women), it moves mostly to ages 15-64.
Fitted rates, and paths with `e0` below `e0_start`, are unchanged.
Rotation needs ages from 0 (life expectancy at birth) with age groups
starting at 15 and 65. With a closed last age group, a warning notes
that life expectancy is then truncated at the highest age. Rotated
forecasts compute life tables for every path and year and are much
slower (about 20 times for single years of age).

Deviations of the rotation from Li, Lee and Gerland (2013):

- the paper forecasts `b` and `k` for both sexes combined (Li and
  Lee 2005) and takes the weights from the life expectancy of both sexes
  combined; `lee_carter()` fits each stratum separately, so the weights
  come from the life expectancy of the stratum (MortCast takes the path
  of life expectancy from a separate model);

- the paper illustrates the method with the median path of life
  expectancy; here every draw is rotated with its own path of life
  expectancy;

- rates change from the jump-off year `T`,
  `log m[x, t] = log m[x, T] + B[x, t] (k*[t] - k[T])`. Eq. (4) writes
  `a[x] + B[x, t] K[t]` with `a` the last observed log rates (as
  `jump_off = "observed"`) or their mean over the years; with
  `jump_off = "fitted"` the fitted rates of year `T` are the starting
  point (package-specific), which keeps the jump-off rates when `B`
  departs from `b`;

- the paper extends the open age group to age 130 with a logistic model
  before fitting (footnote 3); here the age groups of the data are used,
  and life expectancy closes the open interval with `L = l / m`
  ([`life_table()`](life_table.md)).

## References

Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
*Journal of the American Statistical Association* 87(419): 659-671.

Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear
regression approach to the construction of projected lifetables.
*Insurance: Mathematics and Economics* 31(3): 373-393.
[doi:10.1016/S0167-6687(02)00185-3](https://doi.org/10.1016/S0167-6687%2802%2900185-3)

Booth H, Maindonald J, Smith L (2002). Applying Lee-Carter under
conditions of variable mortality decline. *Population Studies* 56(3):
325-336.
[doi:10.1080/00324720215935](https://doi.org/10.1080/00324720215935)

Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
method for forecasting mortality. *Demography* 38(4): 537-549.
[doi:10.1353/dem.2001.0036](https://doi.org/10.1353/dem.2001.0036)

Delwarde A, Denuit M, Partrat C (2007). Negative binomial version of the
Lee-Carter model for mortality forecasting. *Applied Stochastic Models
in Business and Industry* 23(5): 385-401.
[doi:10.1002/asmb.679](https://doi.org/10.1002/asmb.679)

Brouhns N, Denuit M, Van Keilegom I (2005). Bootstrapping the Poisson
log-bilinear model for mortality forecasting. *Scandinavian Actuarial
Journal* 2005(3): 212-224.

Renshaw AE, Haberman S (2008). On simulation-based approaches to risk
measurement in mortality with specific reference to Poisson Lee-Carter
modelling. *Insurance: Mathematics and Economics* 42(2): 797-816.

Li N, Lee R, Gerland P (2013). Extending the Lee-Carter method to model
the rotation of age patterns of mortality decline for long-term
projections. *Demography* 50(6): 2037-2051.
[doi:10.1007/s13524-013-0232-2](https://doi.org/10.1007/s13524-013-0232-2)

Li N, Lee R (2005). Coherent mortality forecasts for a group of
populations: an extension of the Lee-Carter method. *Demography* 42(3):
575-594.
[doi:10.1353/dem.2005.0021](https://doi.org/10.1353/dem.2005.0021)

Sevcikova H, Li N, Kantorova V, Gerland P, Raftery AE (2016).
Age-specific mortality and fertility rates for probabilistic population
projections. In Schoen R (ed.) *Dynamic Demographic Analysis*. Springer,
285-310.
[doi:10.1007/978-3-319-26603-9_15](https://doi.org/10.1007/978-3-319-26603-9_15)

## Examples

``` r
x <- mt_example(deaths = "single")
fit <- mt_fit(lee_carter(), x)
fc <- mt_forecast(fit, h = 10, n_draws = 200, seed = 1)
# \donttest{
# rotation needs ages from 0: French women (HMD, via demography)
if (requireNamespace("demography", quietly = TRUE)) {
  fr <- as_mt_data(demography::fr.mort, series = "female", years = 1970:2006)
  fr <- mt_regroup(fr, age = c(0, 1, seq(5, 100, 5)))
  fc_rot <- mt_forecast(mt_fit(lee_carter(rotate = TRUE), fr), h = 50, n_draws = 100, seed = 1)
  life_expectancy(fc_rot)
}
#>     stratum    sex  label period     mean   q0.025     q0.1     q0.5     q0.9
#>       <int> <char> <char>  <num>    <num>    <num>    <num>    <num>    <num>
#>  1:       1 female female   2007 84.45832 83.97909 84.14176 84.42570 84.84040
#>  2:       1 female female   2008 84.64582 83.98058 84.14838 84.59950 85.21669
#>  3:       1 female female   2009 84.81986 84.07303 84.20940 84.77582 85.50512
#>  4:       1 female female   2010 84.95914 84.04569 84.24819 84.94668 85.64804
#>  5:       1 female female   2011 85.17228 84.19406 84.39411 85.14099 86.06945
#>  6:       1 female female   2012 85.36369 84.29550 84.57521 85.32251 86.25054
#>  7:       1 female female   2013 85.52862 84.25485 84.64479 85.50286 86.57259
#>  8:       1 female female   2014 85.71318 84.34690 84.67043 85.68552 86.81340
#>  9:       1 female female   2015 85.88541 84.46973 84.64564 85.91419 87.05366
#> 10:       1 female female   2016 86.04160 84.41473 84.82095 86.08765 87.20145
#> 11:       1 female female   2017 86.21922 84.48068 84.94156 86.28038 87.39625
#> 12:       1 female female   2018 86.35767 84.40915 85.00913 86.40913 87.58079
#> 13:       1 female female   2019 86.55185 84.48270 85.21895 86.64661 87.86014
#> 14:       1 female female   2020 86.73739 84.59386 85.34588 86.88488 88.21593
#> 15:       1 female female   2021 86.86154 84.42512 85.50065 86.93892 88.28391
#> 16:       1 female female   2022 87.02331 84.72102 85.58790 87.17942 88.46582
#> 17:       1 female female   2023 87.17882 84.52893 85.72703 87.32831 88.81129
#> 18:       1 female female   2024 87.34424 84.46988 85.89195 87.45908 88.97988
#> 19:       1 female female   2025 87.50461 84.85015 86.03837 87.58242 89.08385
#> 20:       1 female female   2026 87.64083 84.90858 86.10474 87.72627 89.32925
#> 21:       1 female female   2027 87.80483 85.23695 86.25216 87.89847 89.50748
#> 22:       1 female female   2028 87.99057 85.60346 86.35903 87.97488 89.43914
#> 23:       1 female female   2029 88.14179 85.67133 86.42451 88.09798 89.69419
#> 24:       1 female female   2030 88.27834 85.71681 86.63274 88.27427 89.95481
#> 25:       1 female female   2031 88.44130 85.77055 86.79347 88.51256 90.18044
#> 26:       1 female female   2032 88.58123 85.83978 86.80131 88.55663 90.17031
#> 27:       1 female female   2033 88.75500 85.89025 86.91354 88.85676 90.58709
#> 28:       1 female female   2034 88.92816 85.86039 87.07892 88.94462 90.73743
#> 29:       1 female female   2035 89.06049 85.90906 87.27126 89.17776 90.87124
#> 30:       1 female female   2036 89.21073 86.11898 87.51358 89.36751 91.02247
#> 31:       1 female female   2037 89.35301 86.21238 87.67122 89.44646 91.30539
#> 32:       1 female female   2038 89.52548 86.20872 87.73986 89.49441 91.54192
#> 33:       1 female female   2039 89.70539 86.27487 87.96464 89.78436 91.90117
#> 34:       1 female female   2040 89.84683 86.38395 88.11563 89.88735 92.10813
#> 35:       1 female female   2041 89.99395 86.85480 88.16030 89.96020 92.17889
#> 36:       1 female female   2042 90.15294 87.11016 88.25911 90.18968 92.35955
#> 37:       1 female female   2043 90.30115 87.20544 88.46409 90.35192 92.40212
#> 38:       1 female female   2044 90.42379 87.30900 88.34859 90.37858 92.47518
#> 39:       1 female female   2045 90.56684 87.39265 88.53738 90.59868 92.69714
#> 40:       1 female female   2046 90.74251 87.68452 88.77240 90.80616 92.93426
#> 41:       1 female female   2047 90.92607 87.92870 88.88210 90.97651 93.05826
#> 42:       1 female female   2048 91.03234 87.83601 88.92326 91.06034 93.36926
#> 43:       1 female female   2049 91.16372 87.79822 89.18376 91.10862 93.48960
#> 44:       1 female female   2050 91.32277 88.06035 89.31103 91.36373 93.77092
#> 45:       1 female female   2051 91.46829 88.03537 89.38917 91.48872 93.95286
#> 46:       1 female female   2052 91.58919 87.91411 89.43010 91.60156 94.11008
#> 47:       1 female female   2053 91.71784 88.32247 89.39964 91.74799 94.15596
#> 48:       1 female female   2054 91.84898 88.25525 89.64994 91.90272 94.28004
#> 49:       1 female female   2055 91.95132 88.26744 89.70716 91.92643 94.40703
#> 50:       1 female female   2056 92.04926 88.41287 89.67235 92.07860 94.61926
#>     stratum    sex  label period     mean   q0.025     q0.1     q0.5     q0.9
#>       <int> <char> <char>  <num>    <num>    <num>    <num>    <num>    <num>
#>       q0.975    at
#>        <num> <num>
#>  1: 84.97797     0
#>  2: 85.37398     0
#>  3: 85.71163     0
#>  4: 85.86549     0
#>  5: 86.23777     0
#>  6: 86.46993     0
#>  7: 86.85593     0
#>  8: 87.11083     0
#>  9: 87.46653     0
#> 10: 87.94733     0
#> 11: 88.06047     0
#> 12: 88.30347     0
#> 13: 88.54289     0
#> 14: 88.87291     0
#> 15: 89.02873     0
#> 16: 89.12558     0
#> 17: 89.42112     0
#> 18: 89.53384     0
#> 19: 89.76067     0
#> 20: 89.93365     0
#> 21: 90.07433     0
#> 22: 90.20059     0
#> 23: 90.47060     0
#> 24: 90.38945     0
#> 25: 90.48946     0
#> 26: 90.92321     0
#> 27: 91.15902     0
#> 28: 91.50720     0
#> 29: 91.64974     0
#> 30: 91.85898     0
#> 31: 92.02428     0
#> 32: 92.12758     0
#> 33: 92.53931     0
#> 34: 92.53557     0
#> 35: 92.65133     0
#> 36: 92.83783     0
#> 37: 92.97059     0
#> 38: 93.15846     0
#> 39: 93.36540     0
#> 40: 93.44503     0
#> 41: 93.76265     0
#> 42: 94.03238     0
#> 43: 94.21914     0
#> 44: 94.28079     0
#> 45: 94.42662     0
#> 46: 94.70687     0
#> 47: 94.90896     0
#> 48: 95.02335     0
#> 49: 95.14042     0
#> 50: 95.19589     0
#>       q0.975    at
#>        <num> <num>
# }
# negative binomial deaths (Delwarde, Denuit and Partrat 2007)
fit_nb <- mt_fit(lee_carter(family = "nb"), x)
mt_parameters(fit_nb)[parameter == "nb_dispersion"]
#>    stratum     parameter  index        value
#>     <char>        <char> <char>        <num>
#> 1:  female nb_dispersion        2.913877e-05
#> 2:    male nb_dispersion        0.000000e+00
```
