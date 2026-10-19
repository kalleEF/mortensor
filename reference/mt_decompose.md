# Decompose differences in mortality by age, cause and subgroup

Splits the difference in a mortality summary between two sets of rates
into contributions of ages, causes of death and population subgroups. It
is a post-processing step for observed data, forecasts, or both:

## Usage

``` r
mt_decompose(
  x,
  y = NULL,
  var = NULL,
  ref = NULL,
  periods = NULL,
  quantity = c("life_expectancy", "prob_dying", "asr", "crude_rate"),
  method = NULL,
  cause = NULL,
  groups = NULL,
  hierarchy = NULL,
  composition = NULL,
  observed = c("data", "fitted", "poisson"),
  pairing = c("independent", "draw"),
  at = NULL,
  to = NULL,
  standard = "esp2013",
  per = 1e+05,
  exposure = NULL,
  n_steps = 20,
  n_draws = 1,
  probs = c(0.025, 0.5, 0.975),
  seed = 1,
  ...
)

# S3 method for class 'mt_decomposition'
summary(
  object,
  by = c("component", "group", "age", "cause"),
  probs = object$probs,
  change = NULL,
  ...
)
```

## Arguments

- x:

  An [mt_data](mt_data.md) object (observed) or an
  [mt_draws](mt_draws.md) object (forecast).

- y:

  Optional second [mt_data](mt_data.md) or [mt_draws](mt_draws.md)
  object with the same age grid and strata.

- var, ref:

  Stratum variable and its reference level for comparisons between
  strata.

- periods:

  Years to compare (see Details). Defaults: the last year for
  comparisons between strata or objects; for changes, the first and last
  observed year of data, or the last observed and the last forecast year
  of draws.

- quantity:

  `"life_expectancy"`, `"prob_dying"`, `"asr"`, `"crude_rate"`, or a
  function `f(m, age)` that takes a matrix of rates `[age, n]` (ages of
  the grid, `age` the [mt_axis](mt_axis.md)) and returns `n` values.

- method:

  `"arriaga"`, `"horiuchi"`, `"stepwise"` or `"kitagawa"` (default:
  `"arriaga"` for life expectancy, `"kitagawa"` for the crude rate,
  `"horiuchi"` otherwise and for pooled composition).

- cause:

  Stratum variable holding causes of death, whose rates add up to all
  causes.

- groups:

  Stratum variable holding population subgroups.

- hierarchy:

  An [`mt_hierarchy()`](mt_hierarchy.md) that declares total and group
  levels of `cause` and `groups` (default: the hierarchy stored in `x`
  or `y`).

- composition:

  `"radix"` (Riffe et al. 2024) or `"pooled"` (specific to this
  package); required with `groups` (see Subgroups).

- observed:

  How observed years held in forecast draws enter: `"data"` (observed
  rates), `"fitted"` (the model's fitted rates) or `"poisson"` (Poisson
  draws of the observed deaths). See Observed and forecast years.

- pairing:

  `"independent"` or `"draw"`: how the draws of `x` and `y` are paired.

- at:

  Starting age (life expectancy at `at`; lower bound of `prob_dying`);
  default the first grid age.

- to:

  Upper age for `"prob_dying"`.

- standard, per:

  Standard population and multiplier for `"asr"` (see
  [`asr()`](derived.md) and Methods); `per` also scales `"crude_rate"`.

- exposure:

  Exposures for forecast years `[age, stratum]` or
  `[age, period, stratum]` (see [`mt_deaths()`](mt_deaths.md)).

- n_steps:

  Number of steps of `"horiuchi"`.

- n_draws:

  For [mt_data](mt_data.md): 1 (the observed rates) or the number of
  Poisson draws.

- probs:

  Quantile levels of the summaries.

- seed:

  Seed of Poisson draws and of the shuffling of `y`.

- ...:

  Passed to the life table (`ax`, `open`, `k`, see
  [`life_table()`](life_table.md)) for life expectancy.

- object:

  An `mt_decomposition` object.

- by:

  Variables of the contributions to keep (`"component"`, `"group"`,
  `"age"`, `"cause"`); the others are summed path by path. `NULL` gives
  the totals.

- change:

  Two periods: the change of each comparison between them (for
  comparisons between strata or objects in several periods), formed path
  by path.

## Value

An `mt_decomposition` object.
[`summary()`](https://rdrr.io/r/base/summary.html) returns contributions
summed (path by path) over the variables not in `by`, with mean and
quantiles; `summary(by = NULL)` returns one row per comparison (its
number in `comparison`) with the values in the two populations, their
difference (mean and quantiles), the jump-off (`jump_off`, an amount in
the unit of the quantity, for example years) where it applies, the sum
of the contributions and the largest difference between the two over
draws (`max_residual`). `summary(change = c(p1, p2))` decomposes the
change between two periods of comparisons between strata or objects (for
example of a gap), path by path.
[`print()`](https://rdrr.io/r/base/print.html) shows the totals.

## Details

- with `var` and `ref`, each level of a stratum variable is compared
  with a reference level (for example the deprivation gap Q5 - Q1 in
  each sex), in each of `periods`;

- with `y`, `y` is compared with `x` for matching strata and periods
  (two forecasts, or observed data and a forecast);

- otherwise the change from the first of `periods` to each later one is
  decomposed for every stratum. For forecast draws, `periods` may mix
  observed years (the training data held in the draws) and forecast
  years.

Forecasts are decomposed draw by draw, so every contribution has an
interval, and sums over ages, causes or subgroups
([`summary()`](https://rdrr.io/r/base/summary.html)) are formed path by
path. Observed data give the decomposition of the observed rates
(`n_draws = 1`) or, with `n_draws > 1`, Poisson draws of the deaths as
in [`mt_observed()`](mt_observed.md). Totals and contributions are means
over draws; [`summary()`](https://rdrr.io/r/base/summary.html) adds
quantiles.

## Observed and forecast years

Forecast draws hold the observed rates of their training years. A change
from an observed year to a forecast year then contains the step from the
observed to the model's fitted rates of that year (the jump-off, see
[`mt_jump_off()`](mt_jump_off.md)) and the sampling noise of the
observed year, besides the forecast change itself. The totals report the
jump-off (`jump_off`): the amount of the difference, in the unit of the
quantity (years for life expectancy), that is due to the observed
instead of the fitted rates of observed years. `observed = "fitted"`
uses the fitted rates of observed years, so that the decomposition
describes the model's change alone; `observed = "poisson"` uses Poisson
draws of the observed deaths, so that their sampling uncertainty enters
the intervals.

Two objects (`y`) are compared draw by draw. Draws of different models
are independent, so by default (`pairing = "independent"`) the draws of
`y` are shuffled before pairing; draws of models with the same seed
would otherwise share their random numbers and give intervals that are
too narrow. Use `pairing = "draw"` when `y` is derived from the draws of
`x` (for example a reconciliation of the same draws).

## Methods

- `"arriaga"` (life expectancy): the symmetric Arriaga decomposition,
  the average of the decompositions from 1 to 2 and (with the sign
  reversed) from 2 to 1, computed on the package's life table (`...`
  passes `ax` and `open` to it). Causes get the age-specific
  contribution in proportion to their rate differences (an implied
  sensitivity, Riffe et al. 2024, eqs. 12-13 of the working paper);
  where the all-cause difference is close to zero, the sensitivity is a
  numerical derivative at the average rates. Exact.

- `"horiuchi"` (any quantity): the line-integral method of Horiuchi,
  Wilmoth and Pletcher (2008) with `n_steps` steps; exact up to a small
  discretisation error, reported in the totals.

- `"stepwise"` (any quantity): stepwise replacement (Andreev, Shkolnikov
  and Begun 2002) from the youngest to the oldest age, averaged over
  both directions of the comparison; with causes, all causes of an age
  are replaced before the next age (DemoDecomp and LEdecomp replace
  cause by cause, which changes the split between causes slightly).
  Exact; results depend on the order. Arriaga's decomposition from 1 to
  2 is stepwise replacement from the oldest to the youngest age, so for
  life expectancy without causes the two symmetric versions are
  identical.

- `"kitagawa"` (crude death rate): the rate effect (age-specific rates
  at the average age structure) and the composition effect (age
  structure at the average rates) of Kitagawa (1955). Exact.

The age-standardised rate (`quantity = "asr"`) is decomposed on the age
groups of the standard population: the rates of each side are first
combined within them with the exposures of that side, as
[`asr()`](derived.md) does (deaths over exposures of the grid cells in a
standard group), so the contributions are by standard age group. The
rate is then linear in the group rates, and `"horiuchi"` and
`"stepwise"` both give the renormalised standard weight times the
difference of the group rate.

Draws whose rates give no valid life table (for example rates far above
any observed value, see the warning of
[`mt_forecast()`](mt_forecast.md)) are left out with a warning, and a
warning reports contributions that do not add up to the difference.

## Subgroups

With `groups`, the populations being compared consist of the levels of a
stratum variable, and the decomposition separates the effect of their
composition. `composition` must be chosen, because the two options
describe different quantities and rest on different grounds:

- `"radix"` (Riffe, Tursun-zade and Trias Llimós 2024, "Arriaga meets
  Kitagawa"): life expectancy (or the probability of dying) of the
  population is the mean of the subgroups' values, weighted by their
  population shares at the starting age (`at`), as when subgroup life
  tables are blended at the radix ("current conditions" in the terms of
  Vaupel 2002, which the authors use). A Kitagawa decomposition of that
  mean gives one composition effect and, for each subgroup, its
  difference times its average share, which is split by age (and cause)
  with `method`. The authors state that this holds when membership is
  fixed from the starting age on and later shares change only through
  mortality: education from age 30 or 35, sex, country of birth. It does
  not hold for areas such as deprivation quintiles or regions, between
  which people move and which are reclassified.

- `"pooled"`: a parametrisation specific to this package (reference
  unknown) for the quantity of the pooled rates, which weight the
  subgroup rates by their population shares at each age (the rates of
  the aggregate strata, as published for whole populations; "current
  rates"). The rates of every subgroup, age and cause and the shares at
  every age are taken as the parameters of the quantity and decomposed
  with the generic line integral of Horiuchi, Wilmoth and
  Pletcher (2008) (`"horiuchi"`, default) or with stepwise replacement
  (`"stepwise"`). The composition effects are reported by age, summed
  over subgroups, whose separate composition terms are not identified.
  Other parametrisations of the same pooled quantity would split the
  difference differently between rates and composition. With many
  subgroups, causes and draws this takes minutes; lower `n_steps` or the
  number of draws for a first look.

Shares come from the exposures: of each observed year, and for forecast
years from `exposure` (default: the last observed year, so composition
does not change in the forecast unless exposures, for example a
population projection, are given).

## Totals among causes and groups

Levels of `cause` and `groups` that are totals of the others (for
example the cause `"all"` or the region `"Germany"`) must not be counted
as parts. They are left out when the hierarchy declares them:
`hierarchy`, or what the object records, the total of
[`mt_group_levels()`](mt_group_levels.md) (`total`, also in fits and
draws of such data) and the hierarchy of
[`mt_aggregate_strata()`](mt_aggregate_strata.md) or
[`reconciled()`](reconciled.md). A level that equals the sum of the
others in the observed deaths without being declared stops the
decomposition.

## References

Arriaga EE (1984). Measuring and explaining the change in life
expectancies. *Demography* 21(1): 83-96.
[doi:10.2307/2061029](https://doi.org/10.2307/2061029)

Kitagawa EM (1955). Components of a difference between two rates.
*Journal of the American Statistical Association* 50(272): 1168-1194.
[doi:10.1080/01621459.1955.10501299](https://doi.org/10.1080/01621459.1955.10501299)

Andreev EM, Shkolnikov VM, Begun AZ (2002). Algorithm for decomposition
of differences between aggregate demographic measures and its
application to life expectancies, healthy life expectancies,
parity-progression ratios and total fertility rates. *Demographic
Research* 7: 499-522.

Horiuchi S, Wilmoth JR, Pletcher SD (2008). A decomposition method based
on a model of continuous change. *Demography* 45(4): 785-801.
[doi:10.1353/dem.0.0033](https://doi.org/10.1353/dem.0.0033)

Riffe T, Tursun-zade R, Trias Llimós S (2024). Arriaga meets Kitagawa.
Life expectancy decomposition with population subgroups. MPIDR Working
Paper WP-2024-029.
[doi:10.4054/MPIDR-WP-2024-029](https://doi.org/10.4054/MPIDR-WP-2024-029)
. Published in *Geographical Analysis* (2026),
[doi:10.1111/gean.70033](https://doi.org/10.1111/gean.70033) .

Vaupel JW (2002). Life expectancy at current rates vs. current
conditions: a reflexion stimulated by Bongaarts and Feeney's "How long
do we live?". *Demographic Research* 7: 365-378.

## Examples

``` r
x <- mt_example(deaths = "single", years = 2000:2015)
# change 2000-2015 in life expectancy at 30, by age, per sex
d <- mt_decompose(x, periods = c(2000, 2015), at = 30)
d
#> <mt_decomposition> life expectancy at 30, method arriaga; 2 comparison(s), 1 draw
#>    comparison    sex   from     to from_value to_value difference
#>         <int> <char> <char> <char>      <num>    <num>      <num>
#> 1:          1 female   2000   2015     53.853   55.443      1.590
#> 2:          2   male   2000   2015     48.585   50.584      1.998
#>    sum_contributions max_residual
#>                <num>        <num>
#> 1:             1.590      4.9e-15
#> 2:             1.998      6.0e-15
summary(d, by = "age")[1:5]
#>       sex   from     to    age        mean      q0.025        q0.5      q0.975
#>    <char> <char> <char> <char>       <num>       <num>       <num>       <num>
#> 1: female   2000   2015     30 0.004438172 0.004438172 0.004438172 0.004438172
#> 2: female   2000   2015     31 0.018028104 0.018028104 0.018028104 0.018028104
#> 3: female   2000   2015     32 0.013227430 0.013227430 0.013227430 0.013227430
#> 4: female   2000   2015     33 0.005466655 0.005466655 0.005466655 0.005466655
#> 5: female   2000   2015     34 0.019223513 0.019223513 0.019223513 0.019223513
# the sex gap in forecast years, with intervals, and its change
fc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 100, seed = 1)
g <- mt_decompose(fc, var = "sex", ref = "male", periods = c(2016, 2025), at = 30)
summary(g, by = NULL)
#>    comparison   from     to period from_value to_value difference
#>         <int> <char> <char>  <num>      <num>    <num>      <num>
#> 1:          1   male female   2016    50.7036 55.53837   4.834771
#> 2:          2   male female   2025    51.8860 56.43574   4.549737
#>    difference_q0.025 difference_q0.5 difference_q0.975 sum_contributions
#>                <num>           <num>             <num>             <num>
#> 1:          4.598310        4.822230          5.116275          4.834771
#> 2:          3.524936        4.490365          5.738382          4.549737
#>    max_residual
#>           <num>
#> 1: 7.105427e-15
#> 2: 6.661338e-15
summary(g, by = "age", change = c(2016, 2025))[1:5]
#>      from     to       period    age          mean       q0.025          q0.5
#>    <char> <char>       <char> <char>         <num>        <num>         <num>
#> 1:   male female 2016 to 2025     30 -0.0024570174 -0.006658112 -0.0028170201
#> 2:   male female 2016 to 2025     31 -0.0022004343 -0.006999818 -0.0025019450
#> 3:   male female 2016 to 2025     32 -0.0052739998 -0.010412948 -0.0057943852
#> 4:   male female 2016 to 2025     33 -0.0035476847 -0.007691005 -0.0039742845
#> 5:   male female 2016 to 2025     34 -0.0001011823 -0.004194642 -0.0000860906
#>         q0.975
#>          <num>
#> 1: 0.002982563
#> 2: 0.004271740
#> 3: 0.001939599
#> 4: 0.001811163
#> 5: 0.005137142
# observed 2015 against the forecast of 2025, with the jump-off (years)
mt_decompose(fc, periods = c(2015, 2025), at = 30)
#> <mt_decomposition> life expectancy at 30, method arriaga; 2 comparison(s), means over 100 draws
#>    comparison    sex   from     to from_value to_value difference
#>         <int> <char> <char> <char>      <num>    <num>      <num>
#> 1:          1 female   2015   2025     55.443   56.436      0.993
#> 2:          2   male   2015   2025     50.584   51.886      1.302
#>    difference_q0.025 difference_q0.5 difference_q0.975 jump_off
#>                <num>           <num>             <num>    <num>
#> 1:             0.442           1.005             1.523   -0.008
#> 2:             0.251           1.337             2.231    0.024
#>    sum_contributions max_residual
#>                <num>        <num>
#> 1:             0.993      4.8e-15
#> 2:             1.302      6.9e-15
#> jump_off: the amount of the difference (in the unit of the quantity) due to observed instead of fitted rates
#> in observed years of the draws (see mt_jump_off()); observed = "fitted" decomposes the model's change alone.
```
