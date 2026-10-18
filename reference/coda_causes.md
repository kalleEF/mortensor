# Multiple-decrement CoDa model of causes of death (experimental)

`coda_causes()` forecasts death rates by cause with the
multiple-decrement form of the compositional Lee-Carter model, which
Oeppen (2008) proposed and tested on Japanese causes of death. Causes
add up to all causes in every draw. The model is experimental.

## Usage

``` r
coda_causes(
  var = "cause",
  total_level = "all",
  jump_off = c("observed", "fitted"),
  residual = "auto"
)
```

## Arguments

- var:

  Stratum variable with the causes.

- total_level:

  Level of `var` for all causes.

- jump_off:

  `"observed"` (forecasts start from the observed deaths of the last
  year) or `"fitted"` (from the fitted deaths).

- residual:

  Add a residual cause, the part of all-cause mortality not covered by
  the listed causes: `TRUE`, `FALSE`, or `"auto"` (added when the listed
  causes leave more than 0.1% of all-cause mortality unexplained in the
  median cell).

## Value

A model specification.

## Details

Only the abstract of Oeppen (2008) is published. This function
implements the natural multiple-decrement form of [`coda()`](coda.md):

1.  in each group (the strata that differ only in `var`, for example one
    group per sex), the cause-specific rates `m[x, c]` add up to
    all-cause rates `m[x]`, which give life-table deaths `d[x]` with
    radix 1 ([`life_table()`](life_table.md) conventions); the deaths
    are split by cause, `d[x, c] = d[x] m[x, c] / m[x]`
    (multiple-decrement life table);

2.  the composition is the joint distribution of the deaths `d[x, c]`
    over age and cause, plus the survivors at the upper bound of a
    closed age range;

3.  as in [`coda()`](coda.md), the composition is centred by its closed
    geometric mean over the years and transformed by the centred
    log-ratio; the first component of a singular value decomposition
    gives `b[x, c] k[t]`, and `k[t]` follows a random walk with drift;

4.  forecast deaths `d[x, c]` give all-cause rates from the all-cause
    life-table deaths `d[x] = sum_c d[x, c]`, and cause-specific rates
    `m[x, c] = m[x] d[x, c] / d[x]`.

The causes add up to the all-cause rates in every draw, and the strata
of `total_level` get their sum. If the listed causes leave part of the
all-cause deaths unexplained (`residual`), that part is a residual cause
inside the model and counts towards the all-cause rates; it is not a
stratum of the output.

One index `k[t]` drives all ages and causes, so the causes cannot follow
different trends beyond the pattern `b[x, c]`. The extensions of
Kjærgaard, Ergemen, Kallestrup-Lamb, Oeppen and Lindahl-Jacobsen (2019)
are not implemented. In a backtest on
[`mt_example_causes()`](mt_example_causes.md) the model was less
accurate than the deprecated [`cause_shares()`](cause_shares.md);
compare it with [`reconciled()`](reconciled.md) separate models in a
backtest before using it.

Jump-off, the open age interval and zero deaths are treated as in
[`coda()`](coda.md): the all-cause log rate of an open interval is
regressed on `k[t]` (its split by cause comes from the composition), and
zero rates of a cause are replaced by half its smallest positive rate in
the group (as `replace.zeros()` of the R package MortalityForecast).
Draws are expected rates. Predictive draws
(`mt_forecast(predictive = TRUE)`) draw Poisson all-cause deaths given
the all-cause rates and exposures and split them multinomially between
the causes by the cause shares of the draw, so that the causes also add
up to all causes in predictive draws; no other observation noise is
added. If intervals are too narrow on real data, that is a property of
the model (see [`recalibrated()`](recalibrated.md),
[`ensemble()`](ensemble.md), or overdispersed models such as
`stanmomo(family = "nb")`).

## References

Oeppen J (2008). Coherent forecasting of multiple-decrement life tables:
a test using Japanese cause of death data. CODAWORK'08, Universitat de
Girona.

Bergeron-Boucher M-P, Canudas-Romo V, Oeppen J, Vaupel JW (2017).
Coherent forecasts of mortality with compositional data analysis.
*Demographic Research* 37(17): 527-566.
[doi:10.4054/DemRes.2017.37.17](https://doi.org/10.4054/DemRes.2017.37.17)

Kjærgaard S, Ergemen YE, Kallestrup-Lamb M, Oeppen J, Lindahl-Jacobsen R
(2019). Forecasting causes of death by using compositional data
analysis: the case of cancer deaths. *Journal of the Royal Statistical
Society Series C* 68(5): 1351-1370.
[doi:10.1111/rssc.12357](https://doi.org/10.1111/rssc.12357)

## See also

[`reconciled()`](reconciled.md) with a hierarchy of causes,
[`coda()`](coda.md).

## Examples

``` r
x <- mt_example_causes(years = 2005:2019)
fit <- mt_fit(coda_causes(), x)
fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
```
