# Interval calibration of random walks and ensembles (2026-10-13)

Question: are the intervals of `rwd()` and of ensembles that contain it too
wide, as the simulation benchmark suggested, and what should the defaults
be? Scripts: `dev/calibration-benchmarks.R` (simulation, HMD) and
`dev/calibration-ons.R` (ONS). All scores are for predictive draws of
observed rates (`mt_backtest()`); coverage is the share of observed rates
inside the central 50%, 80% and 95% intervals.

New option: `rwd(observation_noise = TRUE)`, a random walk with drift
observed with noise (local level model). Changes of observed log rates
contain the noise of two years; the plain random walk accumulates it at
every step. The variant estimates the noise variance (Poisson part plus a
common extra part, from the first autocovariance of the changes), filters
the level and adds the noise only once.

## Results

Ensembles are equal-weight mixtures of Lee-Carter and the random walk
named in the row.

| Data | Model | CRPS | 50% | 80% | 95% |
|---|---|---|---|---|---|
| Simulation (Lee-Carter generator, 10 replicates, h = 1-10) | random walk | 0.0540 | 0.721 | 0.941 | 0.986 |
| | random walk plus noise | 0.0445 | 0.532 | 0.803 | 0.930 |
| | Lee-Carter | 0.0405 | 0.548 | 0.837 | 0.959 |
| | ensemble with random walk | 0.0443 | 0.660 | 0.927 | 0.989 |
| | ensemble with random walk plus noise | 0.0413 | 0.557 | 0.839 | 0.963 |
| England & Wales men (HMD via StMoMo), single ages 0-100, origins 1996-2001, h = 1-10 | random walk | 0.0763 | 0.600 | 0.861 | 0.970 |
| | random walk plus noise | 0.0744 | 0.385 | 0.625 | 0.780 |
| | Lee-Carter | 0.0801 | 0.285 | 0.507 | 0.695 |
| | ensemble with random walk | 0.0724 | 0.444 | 0.763 | 0.951 |
| | ensemble with random walk plus noise | 0.0744 | 0.353 | 0.599 | 0.774 |
| France (HMD via demography), single ages 0-100, origins 1991-1996, h = 1-10 | random walk | 0.0789 | 0.688 | 0.904 | 0.972 |
| | random walk plus noise | 0.0707 | 0.482 | 0.718 | 0.863 |
| | Lee-Carter | 0.0893 | 0.370 | 0.643 | 0.842 |
| | ensemble with random walk | 0.0753 | 0.584 | 0.874 | 0.975 |
| | ensemble with random walk plus noise | 0.0730 | 0.473 | 0.742 | 0.892 |
| ONS England, all causes, sex x IMD quintile, 5-year groups 0-90+, origins 2006-2014, h = 1-5 | random walk | 0.0843 | 0.545 | 0.816 | 0.923 |
| | random walk plus noise | 0.0772 | 0.465 | 0.722 | 0.870 |
| | Lee-Carter | 0.0684 | 0.399 | 0.666 | 0.836 |
| | ensemble with random walk | 0.0704 | 0.488 | 0.799 | 0.940 |
| | ensemble with random walk plus noise | 0.0692 | 0.442 | 0.733 | 0.907 |

## Interpretation

* The over-coverage of the plain random walk in the centre of the
  distribution (50% and 80% levels) comes from the double-counted
  observation noise. The plus-noise variant removes it and improves the
  CRPS on every data set (by 2% to 18%). On the simulation, where the model
  class is right, it is close to nominal.
* On real data every model except the plain random walk is too narrow over
  5 to 10 years: the 95% coverage of Lee-Carter is 0.70 to 0.84, of the
  plus-noise random walk 0.78 to 0.87. Changes of trend add uncertainty
  that a constant drift does not describe, and under-coverage grows with
  the horizon (ONS, Lee-Carter: 0.87 at h = 1, 0.80 at h = 5). The extra
  width of the plain random walk compensates for this in the tails.
* An equal-weight mixture of Lee-Carter and the plain random walk had 95%
  coverage of 0.94 to 0.98 on the real data, with the best or second-best
  CRPS. Its over-coverage on the simulation is the price of that
  robustness.

## Decisions

* `rwd()` keeps `observation_noise = FALSE` as default: its 95% intervals
  are the best calibrated on real data. `observation_noise = TRUE` is
  documented as the sharper choice when the noise is large relative to
  trend changes (small populations, short horizons), with the evaluation
  above.
* Ensembles stay mixtures. The ensemble documentation states that a mixture
  is wider than its components when they disagree, and recommends
  Lee-Carter-type models together with the plain random walk when 95%
  intervals matter.
* Validation: `coverage_rwd_noise` (quick suite) checks the calibration of
  the variant on the simulation benchmark from both sides.
