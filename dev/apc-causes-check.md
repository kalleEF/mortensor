# Check of apc_causes() against bamp_multicause() (coherence plan, A2)

Script: `dev/check-apc-causes.R` (2026-11-18). Reference: `bamp_multicause()`
of the bamp branch "coherent_forecasting" (commit 4363a78), 8,000
iterations, burn-in 2,000, thinning 4; mortensor four chains of 5,000,
burn-in 1,000, thinning 5. Data: `ons_data`, England, women, ages 50-89 in
5-year groups, 2001-2019; causes: cancer (all cancers), circulatory
(ischaemic heart disease and stroke), other (the rest of all causes). RW1
priors for age, period and cohort; the branch's Wishart priors.

| Quantity | Branch | mortensor |
|---|---|---|
| Stick-breaking order | other, cancer, circulatory | other, cancer, circulatory |
| Fitted cause shares, last year (8 age groups x 3 causes) | | largest difference 0.0002 |
| Forecast cause shares 10 years ahead, medians | | largest difference 0.0013 (circulatory, 1.4% relative) |
| Correlation of the period trends of the shares (other, cancer), posterior mean | 0.360 | 0.358 |
| Largest Gelman-Rubin statistic | | 1.001 |
| Time (seconds) | 71 | 31 |
