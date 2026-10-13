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
  composition = c("radix", "pooled"),
  at = NULL,
  to = NULL,
  standard = "esp2013",
  per = 1e+05,
  exposure = NULL,
  n_steps = 20,
  n_draws = 1,
  probs = c(0.025, 0.5, 0.975),
  seed = 1
)

# S3 method for class 'mt_decomposition'
summary(
  object,
  by = c("component", "group", "age", "cause"),
  probs = object$probs,
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

  Optional [`mt_hierarchy()`](mt_hierarchy.md): total and group levels
  of `cause` and `groups` are left out (they are sums of the others).

- composition:

  `"radix"` or `"pooled"` (see Subgroups).

- at:

  Starting age (life expectancy at `at`; lower bound of `prob_dying`);
  default the first grid age.

- to:

  Upper age for `"prob_dying"`.

- standard, per:

  Standard population and multiplier for `"asr"` (see
  [`asr()`](derived.md)); `per` also scales `"crude_rate"`.

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

  Seed of the Poisson draws.

- object:

  An `mt_decomposition` object.

- by:

  Variables of the contributions to keep (`"component"`, `"group"`,
  `"age"`, `"cause"`); the others are summed path by path. `NULL` gives
  the totals.

- ...:

  Unused.

## Value

An `mt_decomposition` object.
[`summary()`](https://rdrr.io/r/base/summary.html) returns contributions
summed (path by path) over the variables not in `by`, with mean and
quantiles; `summary(by = NULL)` returns one row per comparison (its
number in `comparison`) with the values in the two populations, their
difference (mean and quantiles), the sum of the contributions and the
largest difference between the two over draws (`max_residual`).
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
in [`mt_observed()`](mt_observed.md). Observed years inside forecast
draws are the same in every draw.

## Methods

- `"arriaga"` (life expectancy): the symmetric Arriaga decomposition,
  the average of the decompositions from 1 to 2 and (with the sign
  reversed) from 2 to 1, computed on the package's life table. Causes
  get the age-specific contribution in proportion to their rate
  differences (an implied sensitivity, Riffe et al. 2024, eqs. 12-13);
  where the all-cause difference is close to zero, the sensitivity is a
  numerical derivative at the average rates. Exact.

- `"horiuchi"` (any quantity): the line-integral method of Horiuchi,
  Wilmoth and Pletcher (2008) with `n_steps` steps; exact up to a small
  discretisation error, reported in the totals.

- `"stepwise"` (any quantity): stepwise replacement (Andreev, Shkolnikov
  and Begun 2002) from the youngest to the oldest age, averaged over
  both directions of the comparison. Exact; results depend on the order.
  Arriaga's decomposition from 1 to 2 is stepwise replacement from the
  oldest to the youngest age, so for life expectancy without causes the
  two symmetric versions give identical results.

- `"kitagawa"` (crude death rate): the rate effect (age-specific rates
  at the average age structure) and the composition effect (age
  structure at the average rates) of Kitagawa (1955). Exact.

## Subgroups

With `groups`, the populations being compared consist of the levels of a
stratum variable (for example deprivation quintiles or education
groups), and the decomposition separates the effect of their
composition:

- `composition = "radix"` (Riffe, Tursun-zade and Trias Llimós 2024,
  "Arriaga meets Kitagawa"): life expectancy (or the probability of
  dying) of the population is the mean of the subgroups' values,
  weighted by their population shares at the starting age (`at`), as
  when subgroup life tables are blended at the radix. A Kitagawa
  decomposition of that mean gives one composition effect and, for each
  subgroup, its difference times its average share, which is split by
  age (and cause) with `method`.

- `composition = "pooled"`: the quantity of the pooled rates, which
  weight the subgroup rates by their population shares at each age (the
  rates of the aggregate strata; Vaupel 2002, "current rates"). The
  rates of every subgroup, age and cause and the shares at every age are
  the parameters of a `"horiuchi"` or `"stepwise"` decomposition; the
  composition effects are reported by age, summed over subgroups, whose
  separate composition terms are not identified.

Shares come from the exposures: of each observed year, and for forecast
years from `exposure` (default: the last observed year, so composition
does not change in the forecast unless exposures are given).

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

Vaupel JW (2002). Life expectancy at current rates vs. current
conditions: a reflexion stimulated by Bongaarts and Feeney's "How long
do we live?". *Demographic Research* 7: 365-378.

## Examples

``` r
x <- mt_example(deaths = "single", years = 2000:2015)
# change 2000-2015 in life expectancy at 30, by age, per sex
d <- mt_decompose(x, periods = c(2000, 2015), at = 30)
d
#> <mt_decomposition> life expectancy at 30, method arriaga; 2 comparison(s), 1 draw(s)
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
# the sex gap in a forecast year, with intervals
fc <- mt_forecast(mt_fit(lee_carter(), x), h = 10, n_draws = 100, seed = 1)
summary(mt_decompose(fc, var = "sex", ref = "male", at = 30), by = NULL)
#>    comparison   from     to period from_value to_value difference
#>         <int> <char> <char>  <num>      <num>    <num>      <num>
#> 1:          1   male female   2025     51.886 56.43574   4.549737
#>    difference_q0.025 difference_q0.5 difference_q0.975 sum_contributions
#>                <num>           <num>             <num>             <num>
#> 1:          3.524936        4.490365          5.738382          4.549737
#>    max_residual
#>           <num>
#> 1: 6.661338e-15
# observed 2015 against the forecast of 2025
mt_decompose(fc, periods = c(2015, 2025), at = 30)
#> <mt_decomposition> life expectancy at 30, method arriaga; 2 comparison(s), 100 draw(s)
#>    comparison    sex   from     to from_value to_value difference
#>         <int> <char> <char> <char>      <num>    <num>      <num>
#> 1:          1 female   2015   2025     55.443   56.436      0.993
#> 2:          2   male   2015   2025     50.584   51.886      1.302
#>    difference_q0.025 difference_q0.5 difference_q0.975 sum_contributions
#>                <num>           <num>             <num>             <num>
#> 1:             0.442           1.005             1.523             0.993
#> 2:             0.251           1.337             2.231             1.302
#>    max_residual
#>           <num>
#> 1:      4.8e-15
#> 2:      6.9e-15
```
