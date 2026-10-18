# Compositional data (CoDa) model of life-table deaths

`coda()` is the Lee-Carter model in compositional form (Oeppen 2008) as
published by Bergeron-Boucher, Canudas-Romo, Oeppen and Vaupel (2017).
It models the life-table distribution of deaths by age, which sums to 1.
With `coherent = TRUE` or `common_by`, it is their coherent model for
several populations (CoDa-coherent). Status: published, with the
deviations listed below. Page and equation numbers refer to
Bergeron-Boucher et al. (2017).

## Usage

``` r
coda(
  jump_off = c("observed", "fitted"),
  coherent = FALSE,
  common_by = NULL,
  order = c(0, 1, 1),
  drift = TRUE,
  order_dev = c(1, 1, 0),
  drift_dev = FALSE,
  open = c("auto", "kannisto", "regression")
)
```

## Arguments

- jump_off:

  `"observed"` (forecasts start from the observed deaths of the last
  year) or `"fitted"` (from the fitted deaths).

- coherent:

  Fit the coherent model with a common component per group.

- common_by:

  Stratum variable(s) defining groups with their own common component
  (implies `coherent = TRUE`); `NULL` with `coherent = TRUE` puts all
  strata in one group.

- order:

  ARIMA order `c(p, d, q)` of the index (of the common index in the
  coherent model); `c(0, 1, 1)` by default (Bergeron-Boucher et al.
  2017). `c(0, 1, 0)` with `drift = TRUE` is the random walk with drift
  of Lee and Carter (1992), the default of the R package
  MortalityForecast.

- drift:

  Include a drift (a constant of the differenced index) when `d = 1`;
  ignored otherwise (`d = 0` includes a mean, `d = 2` no constant).

- order_dev, drift_dev:

  ARIMA order and drift of the indexes of the deviations of each stratum
  in the coherent model; `c(1, 1, 0)` without drift by default
  (Bergeron-Boucher et al. 2017).

- open:

  Treatment of the old ages: `"auto"` (the Kannisto extension where the
  age groups allow it, otherwise the regression of the open interval),
  `"kannisto"` or `"regression"`.

## Value

A model specification.

## Details

**Model (Section 2.3, Appendix A).** For one stratum:

1.  the death rates of each year give life-table deaths `d[x, t]` with
    radix 1 at the first age ([`life_table()`](life_table.md)
    conventions);

2.  `a[x]` is the geometric mean of `d[x, t]` over the years, closed to
    1;

3.  the deaths are centred by perturbation (`d[x, t] / a[x]`, closed)
    and transformed by the centred log-ratio (clr), eq. (5);

4.  the first component of a singular value decomposition gives
    `clr(d[x, t] / a[x]) = b[x] k[t]`, eq. (6); `b` and `k` sum to zero,
    and their signs make `k` increase over the fitting period (p. 533);

5.  `k[t]` follows an ARIMA(0,1,1) model with drift, fitted by maximum
    likelihood (`order`, `drift`; p. 536 and Appendix A, p. 555, where
    this order had the best BIC for most of 15 countries);

6.  forecasts are transformed back (inverse clr, perturbation by `a`,
    closure, eq. 7) to life-table deaths and then to death rates.

**Jump-off (Appendix G, p. 566).** With `jump_off = "observed"` (the
default), forecasts are the life-table deaths of the last year `T`
perturbed by the change of the model (eq. 17): `d[x, T + j]` is
proportional to `d[x, T] exp(b[x] (k[T + j] - k[T]))`. With `"fitted"`,
forecasts start from the fitted deaths of year `T`. The point forecast
of the ARIMA(0,1,1) index with drift is shifted so that its first step
equals the drift (eq. 20), which removes the break that the
moving-average term causes in the first forecast year; the package
applies the same correction to any ARIMA(0,1,q) index with drift.

**Coherent model (Section 2.4, Appendix C).** The common component
`B[x] K[t]` is the model above fitted to the average mortality of a
group: the arithmetic mean of the death rates of its strata (equal
weights, p. 557), turned into life-table deaths. Each stratum `i` has
its own centre `a[x, i]` and deviation
`clr(d[x, t, i] / a[x, i] / C[exp(B K)]) = b[x, i] k[t, i]` (eq. 8a),
the first component of the SVD after removing the common component. `K`
follows the model of `order` and `drift` (with eq. 20), each `k[t, i]`
an ARIMA(1,1,0) model without drift (`order_dev`, `drift_dev`; p. 537),
whose point forecast levels off, so that the point forecasts of the
strata do not diverge; their intervals still widen with the horizon.
`order_dev = c(1, 0, 0)` gives an AR(1) process with a mean, which Li
and Lee (2005) suggest for their deviations (as cited on p. 531) and
whose deviations stay bounded in distribution. Groups are all strata, or
each level of `common_by`; if the data carry a hierarchy (see
[`mt_aggregate_strata()`](mt_aggregate_strata.md)), the average uses the
bottom strata only, as [`li_lee()`](li_lee.md) does.

**Draws (Appendix B, pp. 556-557).** Draws carry the uncertainty of the
estimates and of the index forecasts: the residuals of eq. (5) are
resampled cell by cell and added to the fitted component, `b` and `k`
are estimated again (`a` is kept, as in the paper), the ARIMA model is
fitted again to the new index, and future paths of the index are
simulated with innovations resampled from its residuals.
`ceiling(sqrt(n_draws))` refits share the draws (the paper uses 100
refits with 100 paths each). In the coherent model, the common and the
stratum components are simulated independently and added (Appendix B);
the draws of the common component are the same for all strata of a
group. The draws are then centred on the point forecast: each draw is
perturbed by the ratio of the point forecast of the first forecast year
to the median of the draws in that year, part by part. This follows the
authors' code (CoDaForecast.R, supplementary material of the paper); eq.
(19) of the paper shifts the bounds of the intervals of life expectancy
so that their median is the point forecast. Draws are expected rates;
predictive draws (`mt_forecast(predictive = TRUE)`) add Poisson
variation of the deaths given the rates and the exposures, and no other
observation noise. The model weights all parts equally on the clr scale,
so the noise of ages with few deaths enters `k`; its intervals can be
wider than those of [`lee_carter()`](lee_carter.md) (p. 543).

**Old ages (Appendix C, pp. 557-558).** With single-year age groups from
80 (at least to 89) and `open = "auto"` (the default) or `"kannisto"`,
the rates from age 80 are replaced in every year by a Kannisto model
(Thatcher, Kannisto and Vaupel 1998) fitted by Poisson maximum
likelihood to the deaths and exposures of the single ages 80 to 100 of
the grid (the authors' code fits ages 80-100 of the Human Mortality
Database, open interval 110+) and applied up to age 120, where the life
table closes; the single ages 80-120 are parts of the composition. The
rate of an open age interval is the forecast deaths over the forecast
person-years of its parts, and the forecasts of the ages from 80 are
those of the smoothed rates. Other age groups (for example 5-year
groups) do not allow this extension. Then, with `open = "auto"` or
`"regression"`, the deaths of an open interval are the survivors at its
lower bound, and its log death rate gets loadings on the same indexes,
estimated by least squares: `log m[open, t] = a_o + b_o k[t]` (coherent
model: `a_o[i] + B_o K[t] + b_o[i] k[t, i]`, `B_o` from the log of the
average rate of the group); a closed age range adds the survivors at its
upper bound as a part. Both are package-specific (reference unknown).

**Zero deaths (Appendix C, p. 558).** Zero counts of the ages that keep
their own rates (below 80 with the Kannisto extension, all ages
otherwise) are replaced multiplicatively in each year: each zero gets
the share `delta` of the deaths of the year, the other shares are
multiplied by `1 - z delta` (`z` zeros), and `delta` is half the
smallest positive count over all ages and years divided by the deaths of
the year. A warning is given when more than 10% of the cells are zero.
Package safeguards: the imputed share of a year is at most 1/2, and a
year whose deaths at ages from 80 allow no Kannisto fit takes the fit to
all years together (with a warning).

**Deviations from the paper.** Besides the package-specific treatment of
other age groups: the life tables use `a_x = n / 2` at every age (the
authors' code uses 0.06 at age 0); the residuals of the bootstrap are
those of eq. (5), as the paper describes (the authors' code takes them
around the in-sample one-step forecasts of the ARIMA model); the number
of refits follows `n_draws`; the draws are centred as in the authors'
code (see above).

## References

Oeppen J (2008). Coherent forecasting of multiple-decrement life tables:
a test using Japanese cause of death data. CODAWORK'08, Universitat de
Girona. <https://hdl.handle.net/10256/742>

Bergeron-Boucher M-P, Canudas-Romo V, Oeppen J, Vaupel JW (2017).
Coherent forecasts of mortality with compositional data analysis.
*Demographic Research* 37(17): 527-566.
[doi:10.4054/DemRes.2017.37.17](https://doi.org/10.4054/DemRes.2017.37.17)

Thatcher AR, Kannisto V, Vaupel JW (1998). *The Force of Mortality at
Ages 80 to 120*. Odense University Press.

Li N, Lee R (2005). Coherent mortality forecasts for a group of
populations: an extension of the Lee-Carter method. *Demography* 42(3):
575-594.
[doi:10.1353/dem.2005.0021](https://doi.org/10.1353/dem.2005.0021)

Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
*Journal of the American Statistical Association* 87(419): 659-671.

## See also

[`coda_causes()`](coda_causes.md) for causes of death,
[`lee_carter()`](lee_carter.md) for the model of death rates,
[`li_lee()`](li_lee.md) for its coherent version.

## Examples

``` r
x <- mt_example(deaths = "single", years = 1995:2012)
fc <- mt_forecast(mt_fit(coda(), x), h = 5, n_draws = 50, seed = 1)
fc2 <- mt_forecast(mt_fit(coda(coherent = TRUE), x), h = 5, n_draws = 50, seed = 1)
# random walk with drift for the index (Lee and Carter 1992)
fc3 <- mt_forecast(mt_fit(coda(order = c(0, 1, 0)), x), h = 5, n_draws = 50, seed = 1)
```
