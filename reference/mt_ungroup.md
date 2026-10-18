# Estimate grid-resolution data from grouped observations

Fills the modelling grid (single years of age by default) from deaths
and exposures observed in coarser age groups. The default is the
penalized composite link model (PCLM; Eilers 2007; Rizzi, Gampe and
Eilers 2015) as published: the estimates are the fitted values of the
model, without rescaling.

## Usage

``` r
mt_ungroup(
  data,
  age = NULL,
  dims = 1,
  method = c("pclm", "monotone", "hyman_kannisto", "natural_gompertz", "gam"),
  exposure_method = c("pclm", "monotone"),
  lambda = NULL,
  criterion = c("bic", "aic"),
  kr = 2,
  deg = 3,
  pord = 2,
  lambda_range = c(0.01, 1e+08),
  preserve_totals = FALSE,
  open_split = c("pclm", "stationary"),
  tail = c("kannisto", "gompertz"),
  tail_width = 2
)
```

## Arguments

- data:

  An [mt_data](mt_data.md) object.

- age:

  Target age axis ([mt_axis](mt_axis.md)). By default the grid when it
  was set (with [`mt_set_grid()`](mt_set_grid.md) or
  [`step_grid()`](mt_prep.md)) and is finer than the data, otherwise
  single years of age from the lowest age to the start of the open
  interval.

- dims:

  1 or 2 (for `method = "pclm"`; `"gam"` is always 2D, the other methods
  1D).

- method:

  Method for deaths: `"pclm"` (default), or `"monotone"`, the monotone
  cubic interpolation of Hyman (1983) of the cumulative deaths (counts
  only, 1D, ignores exposures; group totals reproduced exactly).
  Deprecated, without a published source: `"hyman_kannisto"` (monotone
  spline of log group rates with a Kannisto tail), `"natural_gompertz"`
  (natural spline with a log-linear tail) and `"gam"` (2D tensor-product
  GAM of log group rates over age and year with a Kannisto tail, package
  mgcv); they need exposures on the grid and rescale to the
  exposure-weighted group rates. The vignette "Ungrouping age groups:
  methods and diagnostics" compares methods on data with known single
  years; [`mt_plot_ungroup()`](ungroup-diagnostics.md) and
  [`mt_ungroup_check()`](ungroup-diagnostics.md) check a result.

- exposure_method:

  `"pclm"` or `"monotone"` (Hyman 1983).

- lambda:

  Fixed smoothing parameter(s); `NULL` selects them by `criterion`.

- criterion:

  `"bic"` or `"aic"`.

- kr:

  Knot spacing in years (age and period).

- deg, pord:

  Degree of the B-splines and order of the difference penalty.

- lambda_range:

  Search range for the smoothing parameters.

- preserve_totals:

  Rescale the estimates to the observed group totals after the fit
  (default `FALSE`, the published PCLM); see Details.

- open_split:

  Split of the exposure of an open interval of the grid when deaths are
  ungrouped as rates: `"pclm"` (default, the PCLM of the exposures,
  Rizzi, Gampe and Eilers 2015) or `"stationary"` (package-specific
  stationary-population weights); see Details.

- tail, tail_width:

  For the deprecated `method = "gam"`: the parametric curve that takes
  over from the GAM above the start of the last age group (`"kannisto"`,
  logistic, or `"gompertz"`, log-linear), fitted to the GAM rates from
  age 70, and the width (years) of the logistic blend between the two.

## Value

An `mt_data` object with deaths (and exposures) on the grid. Diagnostics
per fit are attached as attribute `"ungroup"`: smoothing parameters,
effective dimension, the largest relative difference between fitted and
observed group totals (before any rescaling), the number of iterations
of the open-interval weights (`open_iterations`, `NA` without an
expanded open interval), and a flag when a smoothing parameter is at the
boundary of the search range.

## Details

- **Exposures** observed coarser than the grid are ungrouped first, as
  counts by the PCLM (`exposure_method = "pclm"`) or by the monotone
  cubic interpolation of Hyman (1983) applied to the cumulative counts
  (`"monotone"`), which reproduces the group totals exactly.

- **Deaths** are ungrouped as **rates** when exposures are available on
  the grid: the exposures enter the composite link as offset (Rizzi,
  Gampe and Eilers 2015), so cohort-size waves stay in the exposures and
  the estimated rate schedule is smooth. Without exposures, deaths are
  ungrouped as counts per year of age (1D only; 2D smoothing of death
  counts cannot represent cohort waves and is refused).

- `dims = 1` fits each stratum and period separately; `dims = 2` fits a
  smooth age-period surface of rates per stratum with two smoothing
  parameters, the PCLM of rates over age and period with exposures of
  Rizzi et al. (2019).

**Group totals.** The PCLM reproduces the total over all groups of a fit
(the penalty leaves the level of the latent curve free), but not each
group total: the gap is reported as `max_rel_diff` in the diagnostics.
With `preserve_totals = TRUE` the estimates within each observed group
are rescaled to the observed total after the fit, so that aggregation
reproduces the data exactly. This rescaling is a post-hoc step of
mortensor without a published source: it creates steps at the group
boundaries, and with `dims = 2` it puts the year-to-year noise of every
group back onto the smooth surface. On England and Wales men (2000-2011,
5-year groups) the median absolute difference of the single-year log
rates between the two settings is about 0.002 with `dims = 1` and 0.008
with `dims = 2`.

**Open age interval.** For counts (exposures, or deaths without
exposures), an open last group is expanded into single years up to the
maximum age `omega` (see [`mt_data()`](mt_data.md)) and the PCLM spreads
its total over them, as with a last interval of width `omega` minus its
lower bound in the ungroup package. For rates, an open interval of the
*grid* (for example 100+) is expanded into single years up to `omega` as
well, and its exposure is split over these years (`open_split`):

- `"pclm"` (default): by the PCLM of the exposures of all grid cells,
  with a group of 0 counts of 15 years beyond `omega`, as Rizzi, Gampe
  and Eilers (2015) ungroup grouped exposures (their maximum age 115 and
  a group 115-130 with 0 counts); the estimates enter the composite link
  of the rates;

- `"stationary"`: with stationary-population weights
  `l(x) = exp(-cumsum(m(x)))` implied by the fitted rates; starting from
  rates that rise by 10% per year of age, the weights and the fit are
  iterated until no weight changes by more than 1e-6 (at most 50
  iterations; the smoothing parameters are selected again once the
  weights have converged). This treatment is specific to mortensor
  (reference unknown); it was the default up to E3.

On French women and men (1997-2006, deaths in 5-year groups and 90+,
grid open at 100+), the rate of 100+ was too high by 0.14 (women) and
0.15 (men) on the log scale with `"stationary"` and by 0.12 and 0.15
with `"pclm"`; the median absolute log error at 90-99 was 0.020 and
0.022 with `"stationary"`, 0.019 and 0.022 with `"pclm"`. On data
simulated with the exposures of a stationary population above 100
([`mt_example()`](mt_example.md)), `"stationary"` was closer (median
absolute log error at 90+ 0.006, against 0.009 with `"pclm"`), because
that generator has the stationary exposures this option assumes.

**Relation to Rizzi, Gampe and Eilers (2015).** As in the paper, the
observed counts are Poisson with means composed from a latent
distribution `exp(B beta)` (times the exposures for rates), estimated by
penalised likelihood with a second-order difference penalty on the
coefficients weighted by `lambda` (their penalised log-likelihood with
`lambda / 2`, by iteratively reweighted least squares), and the
effective dimension is the trace of the hat matrix. The departures:

- basis: the paper uses one coefficient per fine cell (`B` the identity)
  and allows a B-spline basis; here `B` is a cubic B-spline basis, as in
  Eilers (2007) and the ungroup package;

- smoothing parameter: the paper minimises AIC over a grid of
  `log10 lambda`; the default here is BIC, the default of the ungroup
  package, minimised by golden-section search on `log10 lambda` within
  `lambda_range`. `criterion = "aic"` gives the paper's criterion;

- open last interval of counts: the paper sets a maximum age (115) and
  adds a group with 0 counts beyond it (115-130); here the open group
  ends at `omega` without such a group. On French women in 2006 with
  exposures in 5-year groups and 90+, adding a group with 0 counts of 15
  years beyond `omega = 110` changed the median absolute log error of
  the rates at 90-109 from 0.19 (too low) to 0.66 (too high); with
  exposures in single years to 89 and 90+, from 0.57 to 0.14. The split
  of a wide open group from its total is weakly identified either way;

- exposures: as in the paper, they enter the composite link, and grouped
  exposures are ungrouped by the PCLM first; the open interval of the
  grid as described above (`"pclm"` follows the paper);

- iterations: the paper iterates until the coefficients change by less
  than 1e-6; here the iterations stop when the largest change is below
  1e-6 relative to the coefficients or the penalised deviance changes by
  less than 1e-8 relative, with step halving when it rises;

- in 2D (Rizzi et al. 2019) see Defaults below.

When the deaths are grouped more coarsely than the exposures at old ages
(for example deaths 90+ and single-year exposures), the rates within the
last group are extrapolated by the smooth latent curve from the lower
groups, and the comparison with real data is poor at the highest ages.
On French women in 2006 (Human Mortality Database data in the demography
package), with deaths in 5-year groups and 90+ and exposures in single
years to 110+, the median absolute log error of the ungrouped rates at
ages 90-109 against the single-year rates was 0.09 with `dims = 1` (0.02
at ages 30-89). The errors were within 0.1 on the log scale up to age
101, but the rates were too high by 0.3-0.4 at ages 104-109, because the
latent curve continues the rise of the rates below 90 while the observed
rates rise more slowly at the highest ages. With a grid open at 100+
(exposures 100+), the rate of 100+ was too high by 0.14 on the log scale
(median over 1997-2006, see above). With the exposures also grouped in
5-year groups and 90+, the error at 90-109 was 0.19 (rates too low).
With `dims = 2` (1997-2006, ages 30 and over) and smoothing parameters
chosen by BIC, the error at 90-109 was 1.4: BIC chose very small
smoothing parameters (the large counts vary more than Poisson noise) and
the rates inside the wide last group oscillated; with
`lambda = c(100, 10)` fixed it was 0.05. Treat ungrouped rates inside a
wide open group as extrapolations, prefer data with a higher open age,
check the smoothing parameters (attribute `"ungroup"`) and the result
with [`mt_plot_ungroup()`](ungroup-diagnostics.md).

Cells already observed at grid resolution are kept. When deaths are
ungrouped as rates and the grid starts with age 0 as its own interval,
an observed first year of life is kept and left out of the fit, as in
[`mt_smooth()`](mt_smooth.md): infant mortality does not follow the age
pattern of the following ages, and its peak would spill into ages 1-6.

**Defaults.** The basis and penalty follow the defaults of the ungroup
package (version 1.6.3) for one dimension: cubic B-splines (`deg = 3`)
with a knot every `kr = 2` years, a second-order difference penalty
(`pord = 2`) and the smoothing parameter chosen by BIC. The search range
of the smoothing parameter is wider than ungroup's (0.1 to 1e5 in one
dimension, 0.1 to 1000 in two): `lambda_range = c(1e-2, 1e8)`. In two
dimensions ungroup uses fixed smoothing parameters (1, 1) and a knot
every 7 cells by default; here both parameters are chosen by BIC with a
knot every `kr` years and periods.

The methods `"hyman_kannisto"`, `"natural_gompertz"` and `"gam"` are
deprecated: they are heuristics of this package without a published
source. They still work, with a warning of class
`mt_warning_deprecated`.

In backtests use [`step_ungroup()`](mt_prep.md) so that ungrouping, in
particular in 2D, sees only the training window.

## References

Eilers PHC (2007). Ill-posed problems with counts, the composite link
model and penalized likelihood. *Statistical Modelling* 7(3): 239-254.

Hyman JM (1983). Accurate monotonicity preserving cubic interpolation.
*SIAM Journal on Scientific and Statistical Computing* 4(4): 645-654.

Rizzi S, Gampe J, Eilers PHC (2015). Efficient estimation of smooth
distributions from coarsely grouped data. *American Journal of
Epidemiology* 182(2): 138-147.
[doi:10.1093/aje/kwv020](https://doi.org/10.1093/aje/kwv020)

Rizzi S, Halekoh U, Thinggaard M, Engholm G, Christensen N, Johannesen
TB, Lindahl-Jacobsen R (2019). How to estimate mortality trends from
grouped vital statistics. *International Journal of Epidemiology* 48(2):
571-582. [doi:10.1093/ije/dyy183](https://doi.org/10.1093/ije/dyy183)

## See also

[`mt_plot_ungroup()`](ungroup-diagnostics.md) and
[`mt_ungroup_check()`](ungroup-diagnostics.md) for diagnostics.

## Examples

``` r
x <- mt_example(years = 2010:2014)
u <- mt_ungroup(x)
u
#> <mt_data> 2 strata by sex, omega = 110
#>   deaths     single years 30-99, 100+; periods 2010-2014 (width 1); 710 cells
#>   exposure   single years 30-99, 100+; periods 2010-2014 (width 1); 710 cells
#>   grid       <mt_axis: age, 71 intervals> 30 31 32 33 ... 98 99 100+
#>              <mt_axis: period, 5 intervals> 2010 2011 2012 2013 2014
head(attr(u, "ungroup"))
#>    measure stratum period       method lambda_age lambda_period       ed
#>     <char>   <int> <char>       <char>      <num>         <num>    <num>
#> 1:  deaths       1   2010 pclm_rate_1d  18319.323            NA 6.511396
#> 2:  deaths       1   2011 pclm_rate_1d  27689.910            NA 6.082495
#> 3:  deaths       1   2012 pclm_rate_1d   7709.080            NA 7.443061
#> 4:  deaths       1   2013 pclm_rate_1d   2190.986            NA 8.906068
#> 5:  deaths       1   2014 pclm_rate_1d   1862.519            NA 9.109098
#> 6:  deaths       2   2010 pclm_rate_1d  18106.191            NA 6.839731
#>    max_rel_diff at_bound open_iterations open_converged
#>           <num>   <lgcl>           <int>         <lgcl>
#> 1:   0.05830036    FALSE              NA             NA
#> 2:   0.06314078    FALSE              NA             NA
#> 3:   0.04830307    FALSE              NA             NA
#> 4:   0.06514775    FALSE              NA             NA
#> 5:   0.02661996    FALSE              NA             NA
#> 6:   0.04527347    FALSE              NA             NA
# rescaled to the observed group totals (post-hoc step of the package)
ut <- mt_ungroup(x, preserve_totals = TRUE)
```
