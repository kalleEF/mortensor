# Checks of apc_bayes(coherent = TRUE) (coherence plan, A1)

Data: `ons_data`, England, women and men pooled over the deprivation
deciles, ages 30-89 in 5-year groups, 2001-2019. RW1 priors for age,
period and cohort. Run on 2026-11-18.

## Against the branch implementation (`dev/check-coherent-apc.R`)

Reference: `bamp_coherent()` of the bamp branch "coherent_forecasting"
(commit 4363a78), S = 2 path with sampled rho; 20,000 iterations, burn-in
5,000, thinning 10. mortensor: four chains of 5,000 iterations, burn-in
1,000, thinning 5. The branch writes the deviations as +delta and -delta
with sum_j delta_j = 0; mortensor lets them sum to zero over the strata.
The branch's Gamma(1, 0.05) prior on the precision of delta corresponds to
Gamma(1, 0.1) on mortensor's lambda_d, which this comparison uses.

| Quantity | Branch | mortensor |
|---|---|---|
| Fitted rates, relative difference | | median 0.02%, largest 0.28% |
| Mean logit gap between the sexes per year | | equal to three decimals in every year |
| rho, median (90% interval) | 0.71 (0.12-0.98) | 0.68 (0.11-0.98) |
| kappa (age precision), median | 2.32 | 2.32 |
| lambda (period precision), median | 640 | 651 |
| nu (cohort precision), median | 337 | 338 |
| deviation precision on mortensor's scale, median | 94.9 | 96.0 |
| Mean logit gap 10 years ahead, median (90% interval) | -0.44 (-0.85 to -0.06) | -0.45 (-0.91 to 0.00) |
| Rates 10 years ahead, medians | | within 1% |
| Rates 10 years ahead, width of the 90% interval relative to the branch | 1 | 1.03-1.15 |
| Time (seconds) | 136 | 20 |

The remaining difference in the widths comes from the parameterisation:
without the branch's time constraint, the level of the gap carries the
uncertainty of the mean of the deviations. mortensor's chains: Gelman-Rubin
statistics at most 1.001, acceptance 0.95.

## The prior of the deviation precision

The data add about 6e-4 to the rate of the Gamma conditional of the
deviation precision (half the AR(1) quadratic form of the deviations), so
the branch's rate 0.05 sets the posterior. Same model, mortensor's sampler,
3,000 iterations:

| Rate of the prior | Posterior median of lambda_d | Posterior median of rho | Width of the 90% interval of the mean logit gap, 10 years ahead |
|---|---|---|---|
| 0.05 | 189 | 0.69 | 0.62 |
| 0.005 | 1,669 | 0.82 | 0.24 |
| 0.0005 (mortensor's default) | 9,091 | 0.93 | 0.13 |
| 0.00005 | 19,100 | 0.96 | 0.10 |

## Backtest (`dev/coherent-apc-backtest.R`)

Origins 2009-2014, horizon 5 years, 500 draws; gaps of women against men;
results without thresholds. Coherent model: `mcmc_control("fast")`.

| Model | CRPS log rates | 95% coverage log rates | CRPS e30 | CRPS log ASR | Energy e30 joint | Energy log ASR joint | Variogram log ASR joint | CRPS gap e30 | 95% coverage gap e30 | CRPS gap log ASR | 95% coverage gap log ASR |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Lee-Carter (separate) | 0.0383 | 0.686 | 0.205 | 0.0329 | 0.307 | 0.0485 | 0.00030 | 0.0995 | 1.00 | 0.0126 | 1.00 |
| Li-Lee | 0.0378 | 0.746 | 0.192 | 0.0289 | 0.284 | 0.0441 | 0.00125 | 0.0947 | 0.63 | 0.0260 | 0.00 |
| apc_bayes(coherent = TRUE) | 0.0322 | 0.947 | 0.138 | 0.0163 | 0.213 | 0.0250 | 0.00084 | 0.1260 | 1.00 | 0.0150 | 1.00 |

The backtest took 292 seconds, mostly for the coherent model.
