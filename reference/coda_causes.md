# Multiple-decrement CoDa model of causes of death (experimental)

`coda_causes()` forecasts death rates by cause with the
multiple-decrement form of the compositional Lee-Carter model, which
Oeppen (2008) proposed and tested on Japanese causes of death. Causes
add up to all causes in every draw. The options `weights`, `rho` and
`individual` add the extensions of Kjærgaard, Ergemen, Kallestrup-Lamb,
Oeppen and Lindahl-Jacobsen (2019). The model is experimental. Status:
published, with the deviations listed below; section numbers refer to
Oeppen (2008), and to Kjærgaard et al. (2019) for the extensions.

## Usage

``` r
coda_causes(
  var = "cause",
  total_level = "all",
  jump_off = c("observed", "fitted"),
  residual = "auto",
  rank = 1,
  order = c(0, 1, 1),
  drift = TRUE,
  weights = c("equal", "deaths"),
  rho = NULL,
  individual = FALSE
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

- weights:

  Weights of the parts (ages by cause) in the decomposition: `"equal"`
  (Oeppen 2008) or `"deaths"`, their mean life-table deaths (Kjærgaard
  et al. 2019, eq. 10).

- rho:

  `NULL` (equal weights of the years) or the weight of the last year, a
  number between 0 and 1: years weighted by `rho (1 - rho)^(T - t)`
  (Kjærgaard et al. 2019, eq. 11; they use 0.1).

- individual:

  Add one individual component per cause to the joint components
  (Kjærgaard et al. 2019, eq. 12).

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
C, with the package-specific cap `0.5 / z` of `delta` for `z` zeros in a
year; see [`coda()`](coda.md)). How the death rate of the open age
interval follows from the forecast deaths is not stated; here the
all-cause log rate of an open interval is regressed on the indexes and
split by cause through the composition (package-specific, as
`open = "regression"` of [`coda()`](coda.md)).

Draws are expected rates. Predictive draws
(`mt_forecast(predictive = TRUE)`) draw Poisson all-cause deaths given
the all-cause rates and exposures and split them multinomially between
the causes by the cause shares of the draw, so that the causes also add
up to all causes in predictive draws; no other observation noise is
added (package-specific).

**Weights and individual components (Kjærgaard et al. 2019, Section
2.4).** Kjærgaard et al. call the model above CT-CoDA and extend it to
the two-step model 2S-CoDA with two elements, weights of the parts and
years and cause-specific components. They are options here:

- `weights = "deaths"` weights each part (age by cause) by its temporal
  mean of the life-table deaths, `w[x, i] = dbar[x, i] / sum(dbar)` (eq.
  10), so that a large cause counts more than a small one and splitting
  a cause does not double its weight;

- `rho` weights the years `t = 1, ..., T` by `rho (1 - rho)^(T - t)`
  (eq. 11, after Hyndman, Booth and Yasmeen 2013): the last year gets
  the weight `rho`, earlier years less. The paper uses `rho = 0.1`,
  chosen by cross-validation (Section 2.8). Eq. (11) states
  `t = 1, ..., T - 1`; the text gives the last year the weight `rho`,
  which the package follows;

- with weights, the components are the first ones of a generalized
  singular value decomposition (van Loan 1976) with these weights of the
  parts and years: they minimise the sum over parts and years of
  `w[x, i] w[t]` times the squared residual of eq. (9). The centre `a`
  is the unweighted geometric mean, as in CT-CoDA (Section 2.4: 2S-CoDA
  "follows the same first steps"). Loadings have unit length; the index
  of each year is the weighted least-squares fit of that year;

- `individual = TRUE` adds the cause-specific components of eq. (12):
  the residuals of each cause after the joint components (`rank` of
  them; the paper uses one) get the first component of their own
  unweighted SVD, `b_i[x] k_i[t]`, as the paper names it for this step
  (Section 2.4); the weights of eq. (10) and (11) apply to the joint
  components of eq. (9). This follows the joint and individual
  components of Lock et al. (2013) in two steps without iteration, as in
  the paper.

`coda_causes(weights = "deaths", rho = 0.1, individual = TRUE)` is
2S-CoDA; with `jump_off = "fitted"`, forecasts start from the model, as
in eq. (7) of the paper. `rank = 3` is the paper's CT-CoDA (Section
2.8). Every index, joint or individual, has its own ARIMA model of
`order` and `drift`, independent of the others. The paper selects the
order of integration of each index by the augmented Dickey-Fuller test
and the ARMA orders by AIC, with a drift for the rank-one indexes
(Section 2.8); the settings of the test are not given, and the package
uses `order`. Further choices of the package (reference unknown): in the
clr deaths, the individual loadings are centred over all parts, which
leaves the closed composition unchanged; the log rate of an open age
interval is regressed on the joint indexes only, weighted by the weights
of the years; a residual cause gets an individual component, the
survivors of a closed age range none. The paper gives point forecasts
only. Draws refit the same decomposition (with the same weights of the
joint components) to bootstrap samples. With weights of the parts, each
cell takes a residual of its own part: pooled resampling (Appendix B of
Bergeron-Boucher et al. 2017) would move the large residuals of small,
down-weighted parts into the parts that determine the weighted fit, and
on [`mt_example_causes()`](mt_example_causes.md) it about doubled the
innovation standard deviation of the refitted index. With weights of the
years, the year of each residual is drawn with probability proportional
to its weight. [`mt_parameters()`](mt_parameters.md) gives the
individual loadings and indexes as `b_<cause>` and `k_<cause>` with
their ARIMA models and the share of the residual variation of the cause
after the joint components that they explain (`explained_<cause>`), the
weights `w_part` and `w_time`, and the share of the variation of each
cause explained by the joint and the individual components
(`explained_joint`, `explained_individual`; squared Frobenius norms of
the centred clr deaths of the cause, Section 3.3.1 and Fig. 4).

The paper replaces zero counts by half a death; zeros here are replaced
multiplicatively, as above. The paper's third model, VECM-CoDA (Section
2.5: one SVD per cause and a cointegrated vector error correction model
of the cause indexes, with the Johansen test), is not implemented. The
paper's data (France, Institut National d'Etudes Démographiques,
1955-2008; the Netherlands, Statistics Netherlands, 1958-2014; 5-year
age groups from 25) are not available to the package, so its numbers
(Table

1.  are not reproduced. The alpha-transformation that Dong, Shang, Hui
    and Bruhn (2025) use in place of the clr for compositions of causes
    with zero counts is not implemented. Compare the model with
    [`reconciled()`](reconciled.md) separate models in a backtest before
    using it.

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

van Loan CF (1976). Generalizing the singular value decomposition. *SIAM
Journal on Numerical Analysis* 13(1): 76-83.
[doi:10.1137/0713009](https://doi.org/10.1137/0713009)

Hyndman RJ, Booth H, Yasmeen F (2013). Coherent mortality forecasting:
the product-ratio method with functional time series models.
*Demography* 50(1): 261-283.
[doi:10.1007/s13524-012-0145-5](https://doi.org/10.1007/s13524-012-0145-5)

Lock EF, Hoadley KA, Marron JS, Nobel AB (2013). Joint and individual
variation explained (JIVE) for integrated analysis of multiple data
types. *Annals of Applied Statistics* 7(1): 523-542.
[doi:10.1214/12-AOAS597](https://doi.org/10.1214/12-AOAS597)

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
# 2S-CoDA of Kjaergaard et al. (2019): weights and individual components
fit2 <- mt_fit(coda_causes(weights = "deaths", rho = 0.1, individual = TRUE), x)
p <- mt_parameters(fit2)
p[parameter %in% c("explained_joint", "explained_individual")]
#>            stratum            parameter  index      value
#>             <char>               <char> <char>      <num>
#>  1: causes: female      explained_joint    cvd 0.25084751
#>  2: causes: female      explained_joint cancer 0.42160550
#>  3: causes: female      explained_joint  other 0.35442880
#>  4: causes: female explained_individual    cvd 0.16084057
#>  5: causes: female explained_individual cancer 0.09790563
#>  6: causes: female explained_individual  other 0.13299746
#>  7:   causes: male      explained_joint    cvd 0.35871912
#>  8:   causes: male      explained_joint cancer 0.56148267
#>  9:   causes: male      explained_joint  other 0.50032461
#> 10:   causes: male explained_individual    cvd 0.15580499
#> 11:   causes: male explained_individual cancer 0.06960404
#> 12:   causes: male explained_individual  other 0.08919613
```
