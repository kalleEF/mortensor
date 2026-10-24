# Control the jump-off of any model

Forecasts of most models start from their fitted rates in the last
observed year. Where the fit misses the level of the last years (for
example after a recent change of trend), the forecast starts with a step
away from the observations. `anchored()` wraps any model specification
and moves the start of its forecasts to the observed level, per age and
stratum. With the defaults this is the jump-off adjustment of Lee and
Miller (2001): the forecasts start from the observed rates of the last
training year and keep this offset in every forecast year.

## Usage

``` r
anchored(
  spec,
  to = c("observed", "smoothed"),
  years = 1,
  decay = Inf,
  shrink = FALSE
)
```

## Arguments

- spec:

  A model specification.

- to:

  `"observed"` or `"smoothed"`.

- years:

  Number of last training years whose deaths define the offset (default
  1, the last training year; see above for the averages of Stoeldraijer
  et al. 2018).

- decay:

  Years over which the offset fades out (default `Inf`: constant).

- shrink:

  Shrink offsets of ages with few deaths towards the mean of their
  stratum (default `FALSE`).

## Value

A model specification. The fitted object has the offsets in
`$anchor$offset` (`[age, stratum]`, log scale).

## Details

1.  the offset is the log ratio of observed to fitted deaths, summed
    over the last `years` training years (`to = "observed"`: the rates
    before smoothing; `to = "smoothed"`: the rates after preprocessing,
    which differs when the pipeline smooths);

2.  with `shrink = TRUE`, offsets of ages with few deaths are shrunk
    towards the mean offset of their stratum (see below);

3.  the forecast draws are multiplied by `exp(offset * w(h))`, where the
    weight `w(h)` falls linearly from 1 at the first forecast year to 0
    after `decay` years (`decay = Inf` keeps the offset; `decay = 0`
    leaves the model's own jump-off).

The defaults `years = 1, decay = Inf, shrink = FALSE` give, for any
model, the observed jump-off of Lee-Carter models
(`lee_carter(jump_off = "observed")`). The offset multiplies every draw
of a cell by the same factor, so it moves the forecast distribution of
the log rates without changing its width. Models that have no fitted
rates (random walks, time series per age) start from the observed rates
already and are left unchanged; in an [`ensemble()`](ensemble.md) each
draw is anchored with the offset of the component it comes from, which
moves the components relative to each other and can change the width of
the mixture. For [`reconciled()`](reconciled.md) models the bottom
strata are anchored and the aggregates recomputed, so the draws stay
coherent. [`mt_jump_off()`](mt_jump_off.md) reports the step between the
last observed or fitted year and the first forecast year.

**Jump-off options of Stoeldraijer et al. (2018).** The jump-off rates
are the rates of the last year of the fitting period, the jump-off year.
Stoeldraijer et al. compare six options for the Lee-Carter model (their
names in quotes):

- "model values", the fitted rates of the jump-off year (Lee and Carter
  1992): the model itself, without `anchored()` (or `decay = 0`);

- "last observed values" (Lee and Miller 2001): `anchored()` with its
  defaults, for any model; for Lee-Carter also
  `lee_carter(jump_off = "observed")`;

- "average n years observed" (`n` = 2 to 5):
  `lee_carter(jump_off = "observed", jump_off_years = n)`, as the paper
  defines it: the mean of the observed log rates of the last `n` years,
  with the index 0 at the midpoint of these years. `anchored(years = n)`
  is a package-specific variant for any model (reference unknown): its
  offset is the log ratio of the observed to the fitted deaths of these
  years, so the forecast also keeps the change of the fitted rates from
  these years to the last year.

They call a forecast accurate if its out-of-sample errors in historical
data are small (mean absolute error of the forecast life expectancy at
65 one and five years ahead), and robust if it changes little when one
more year of data is added (standard deviation of the changes of the
forecast life expectancy at 65 in 2020 between successive fitting
periods). In eight Western European countries the last observed values
were the most accurate for most countries and an average of observed
years the most robust; the more years averaged, the more robust and the
less accurate the forecast. The model values were the most accurate
where the model fitted the last years well. They recommend moving from
the most accurate jump-off in the first forecast year to the most robust
one by the fifth (Appendix 6, without a formula for the weights). Where
the last observed values are the most accurate and the model values the
most robust, `anchored(decay = 4)` does this: the weight of the observed
offset falls linearly from 1 in the first forecast year to 0 in the
fifth (the linear weights are package-specific). The filtered jump-off
of De Jong and Tickle (2006) is `lee_carter(jump_off = "filtered")`.

**Package-specific options** (reference unknown). Starting from the
observed rates of one year carries their random variation into every
forecast year, which matters where deaths are few. The other settings
are package-specific: averaging the offset over the last `years` years
(see above), fading it out over `decay` years (forecasts more than
`decay` years ahead are then those of the model itself, and during the
fade the offsets add a small apparent trend) and shrinking it
(`shrink = TRUE`). The shrinkage is empirical Bayes within each stratum:
the offset \\r_x\\ of age \\x\\ has the approximate Poisson sampling
variance \\v_x = 1/D_x\\, where \\D_x\\ is the number of deaths in the
last `years` (at least 0.5); the variance \\\tau^2\\ of the true offsets
between the ages of the stratum is estimated with the moment estimator
of DerSimonian and Laird (1986); and each offset becomes \\\hat\mu +
\tau^2 / (\tau^2 + v_x) (r_x - \hat\mu)\\, where \\\hat\mu\\ is the mean
offset weighted by \\1 / (\tau^2 + v_x)\\.

Anchoring helps where the fitted rates of the last year are far from the
observed level; where the model's own jump-off is already close to it,
anchoring changes little. The vignette "Where a forecast starts: the
jump-off" shows the effect in a backtest.

## References

Stoeldraijer L, van Duin C, van Wissen L, Janssen F (2018). Comparing
strategies for matching mortality forecasts to the most recently
observed data: exploring the trade-off between accuracy and robustness.
*Genus* 74: 16.
[doi:10.1186/s41118-018-0041-y](https://doi.org/10.1186/s41118-018-0041-y)

Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
*Journal of the American Statistical Association* 87(419): 659-671.

Lee R, Miller T (2001). Evaluating the performance of the Lee-Carter
method for forecasting mortality. *Demography* 38(4): 537-549.

De Jong P, Tickle L (2006). Extending Lee-Carter mortality forecasting.
*Mathematical Population Studies* 13(1): 1-18.
[doi:10.1080/08898480500452109](https://doi.org/10.1080/08898480500452109)

DerSimonian R, Laird N (1986). Meta-analysis in clinical trials.
*Controlled Clinical Trials* 7(3): 177-188.

## Examples

``` r
x <- mt_example(deaths = "single", years = 1995:2019)
fit <- mt_fit(anchored(lee_carter()), x)
summary(as.vector(fit$anchor$offset))
#>       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#> -0.4028109 -0.0241849  0.0004603 -0.0013737  0.0269568  0.2285328 
fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
mt_jump_off(fc)
#>       sex  label observed   fitted forecast step_observed step_fitted
#>    <char> <char>    <num>    <num>    <num>         <num>       <num>
#> 1: female female 55.95597 55.95828 56.06028    0.10431465   0.1020055
#> 2:   male   male 51.16230 51.17137 51.25410    0.09180217   0.0827333
# package-specific: the mean offset of three years, fading out over ten years
fit3 <- mt_fit(anchored(lee_carter(), years = 3, decay = 10, shrink = TRUE), x)
```
