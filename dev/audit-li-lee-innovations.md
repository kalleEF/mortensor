# Audit: correlated deviations in li_lee() (coherence plan, item B3)

Question of the owner (decision 4 of `dev/coherence-plan.md`): should the
innovations of the stratum indexes of `li_lee()` be drawn jointly? Li and
Lee (2005) fit one AR(1) per population, and `li_lee()` draws the
innovations independently. Script: `dev/audit-li-lee-innovations.R`
(2026-11-18). The alternative resamples the residual vectors of all strata
of one year together, which keeps their correlation without estimating a
covariance matrix (20 strata, 18 residual years in `ons_data`).

## Data

- `ons_data`: England, sex by deprivation decile (20 strata), 5-year age
  groups, 2001-2019, `li_lee(common_by = "sex")`.
- HMD France (`demography::fr.mort`): women and men, ages 0-99, 1950-2006,
  `li_lee()`.

## Correlation of the residuals of the stratum indexes

Signs aligned so that a positive index raises the rates on average over
age.

| Data | Strata | Years | Mean r | Mean abs(r) | Pairs beyond 2 SE | First eigenvalue share |
|---|---|---|---|---|---|---|
| England deciles | 20 | 18 | 0.05 | 0.19 | 3.7% | 0.17 |
| France sexes | 2 | 56 | 0.49 | 0.49 | 1 of 1 | 0.74 |

In England the deviations of the deciles from the trend of their sex show
no correlation beyond chance (about 5% of pairs would exceed 2 SE). In
France the deviations of women and men from the common trend move
together.

## Forecast spread (2,000 draws; England 10 years, France 20 years ahead)

| Data | Innovations | Mean width of 90% interval of e0 | Mean width of 90% interval of the gap in e0 | Relative width of 90% interval of total deaths |
|---|---|---|---|---|
| England | independent | 3.09 | 0.54 | 0.214 |
| England | joint | 3.09 | 0.53 | 0.215 |
| France | independent | 5.26 | 1.85 | 0.512 |
| France | joint | 5.30 | 1.70 | 0.521 |

## Backtests (500 draws; England origins 2009-2014, h = 5; France origins 1987-1996, h = 10)

| Data | Innovations | Energy e0 joint | Variogram e0 joint | Energy log ASR joint | CRPS gap in e0 | 95% coverage gap in e0 | 95% coverage gap in log ASR | CRPS log rates |
|---|---|---|---|---|---|---|---|---|
| England | independent | 1.88 | 7.83 | 0.205 | 0.416 | 0.30 | 0.14 | 0.0810 |
| England | joint | 1.88 | 7.82 | 0.205 | 0.416 | 0.29 | 0.14 | 0.0809 |
| France | independent | 0.748 | 0.108 | 0.0725 | 1.05 | 0.03 | 0.00 | 0.0870 |
| France | joint | 0.751 | 0.107 | 0.0727 | 1.07 | 0.01 | 0.00 | 0.0869 |

## Decision

No change to `li_lee()`. Joint innovations leave the scores of the joint
distribution and of the gaps unchanged in both data sets; in France they
narrow the gap intervals that are already too narrow.

## Side finding

The intervals of `li_lee()` for gaps between strata are far too narrow in
both backtests (95% coverage 0.14-0.30 in England, 0-0.03 in France),
whatever the innovations. The model fixes `a`, `B`, `b` and the AR(1)
coefficients, so the uncertainty of the stratum indexes covers only their
innovations. This is a property of the model as published; the gap scores
of `mt_backtest(contrast = )` show it, and `recalibrated()` or a model with
parameter uncertainty (for example `apc_bayes(coherent = TRUE)`) address
it.
