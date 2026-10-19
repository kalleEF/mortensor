# Multiple-decrement CoDa model of causes of death (experimental)

`coda_causes()` forecasts death rates by cause with the
multiple-decrement form of the compositional Lee-Carter model, which
Oeppen (2008) proposed and tested on Japanese causes of death. Causes
add up to all causes in every draw. The model is experimental. Status:
published, with the deviations listed below; section numbers refer to
Oeppen (2008).

## Usage

``` r
coda_causes(
  var = "cause",
  total_level = "all",
  jump_off = c("observed", "fitted"),
  residual = "auto",
  rank = 1,
  order = c(0, 1, 1),
  drift = TRUE
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

- rank:

  Number of components (1 by default; Oeppen 2008 used 2 for Japan).

- order, drift:

  ARIMA order and drift of each index, as in [`coda()`](coda.md).

## Value

A model specification.

## Details

**Model (Sections 3, 6 and 8.3).** In each group (the strata that differ
only in `var`, for example one group per sex):

1.  the cause-specific rates `m[x, c]` add up to all-cause rates `m[x]`,
    which give life-table deaths `d[x]` with radix 1
    ([`life_table()`](life_table.md) conventions); the deaths are split
    by cause, `d[x, c] = d[x] m[x, c] / m[x]` (relations of the
    multiple-decrement life table, eq. 1 and 7);

2.  the array of deaths by age and cause of each year is unfolded into
    one composition (Section 6), plus the survivors at the upper bound
    of a closed age range;

3.  the compositions are centred by their closed geometric mean over the
    years and transformed by the centred log-ratio; the first `rank`
    components of a singular value decomposition give
    `sum_j b_j[x, c] k_j[t]` (Section 3, steps 1-7);

4.  forecast deaths `d[x, c]` give all-cause rates from the all-cause
    life-table deaths `d[x] = sum_c d[x, c]`, and cause-specific rates
    `m[x, c] = m[x] d[x, c] / d[x]` (eq. 7, Section 8.3).

The causes add up to the all-cause rates in every draw, and the strata
of `total_level` get their sum. If the listed causes leave part of the
all-cause deaths unexplained (`residual`), that part is a residual cause
inside the model and counts towards the all-cause rates; it is not a
stratum of the output (package-specific).

**Indexes and rank.** For Japanese women 1951-1990, Oeppen (2008,
Sections 6 and 7) uses two components and selects ARIMA(0,2,2) models
for the indexes by AICc; Section 3 notes, for the Lee-Carter model, that
one component is usually adequate for advanced economies and that its
index can be forecast by a simple ARIMA model such as the random walk
with drift. These choices belong to one data set; the defaults here are
one component and the index model of [`coda()`](coda.md), ARIMA(0,1,1)
with drift (Bergeron-Boucher et al. 2017, with Oeppen as co-author).
`rank = 2, order = c(0, 2, 2)` gives Oeppen's specification for Japan.
With several components, each index has its own ARIMA model, independent
of the others (the paper does not say how the second index was
forecast).

**Deviations from the paper.** Oeppen (2008) gives point forecasts only:
the jump-off (forecasts start from the observed deaths of the last year,
`jump_off = "observed"`) and the draws (bootstrap of the residuals and
of the index paths, centred on the point forecast) are those of
[`coda()`](coda.md) (Bergeron-Boucher et al. 2017, Appendices B and G).
Oeppen (2008, Section 3) replaces zero counts by lightly smoothing the
deaths of each year with a penalised histogram smoother whose smoothing
parameter is not given; here zero counts are replaced multiplicatively
in the composition of each year (Bergeron-Boucher et al. 2017, Appendix
C; see [`coda()`](coda.md)). How the death rate of the open age interval
follows from the forecast deaths is not stated; here the all-cause log
rate of an open interval is regressed on the indexes and split by cause
through the composition (package-specific, as `open = "regression"` of
[`coda()`](coda.md)).

Draws are expected rates. Predictive draws
(`mt_forecast(predictive = TRUE)`) draw Poisson all-cause deaths given
the all-cause rates and exposures and split them multinomially between
the causes by the cause shares of the draw, so that the causes also add
up to all causes in predictive draws; no other observation noise is
added (package-specific).

The extensions of Kjærgaard, Ergemen, Kallestrup-Lamb, Oeppen and
Lindahl-Jacobsen (2019) are not implemented, nor is the
alpha-transformation that Dong, Shang, Hui and Bruhn (2025) use in place
of the clr for compositions of causes with zero counts. Compare the
model with [`reconciled()`](reconciled.md) separate models in a backtest
before using it.

## References

Oeppen J (2008). Coherent forecasting of multiple-decrement life tables:
a test using Japanese cause of death data. CODAWORK'08, Universitat de
Girona. <https://hdl.handle.net/10256/742>

Bergeron-Boucher M-P, Canudas-Romo V, Oeppen J, Vaupel JW (2017).
Coherent forecasts of mortality with compositional data analysis.
*Demographic Research* 37(17): 527-566.
[doi:10.4054/DemRes.2017.37.17](https://doi.org/10.4054/DemRes.2017.37.17)

Kjærgaard S, Ergemen YE, Kallestrup-Lamb M, Oeppen J, Lindahl-Jacobsen R
(2019). Forecasting causes of death by using compositional data
analysis: the case of cancer deaths. *Journal of the Royal Statistical
Society Series C* 68(5): 1351-1370.
[doi:10.1111/rssc.12357](https://doi.org/10.1111/rssc.12357)

Dong Z, Shang HL, Hui F, Bruhn A (2025). A compositional approach to
modeling cause-specific mortality with zero counts. *Annals of Actuarial
Science* 19(3): 416-441.
[doi:10.1017/S1748499524000320](https://doi.org/10.1017/S1748499524000320)

## See also

[`reconciled()`](reconciled.md) with a hierarchy of causes,
[`coda()`](coda.md).

## Examples

``` r
x <- mt_example_causes(years = 2005:2019)
fit <- mt_fit(coda_causes(), x)
fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
```
