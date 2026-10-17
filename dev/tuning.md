# Tuning and recalibration (2026-10-17)

Question: does choosing settings by forecast accuracy on the last years of
the training data (the hold-out workflow often used with
`demography::fdm()`) improve forecasts out of sample, and what does
`recalibrated()` add? Script: `dev/tuning-benchmark.R` (about 25 minutes),
results in `dev/tuning-benchmark.rds`.

Design: rolling-origin backtests, h = 10. England & Wales men (HMD via
StMoMo), single ages 30-100, origins 1996-2001; France men (HMD via
demography), single ages 0-95, 1950-2006, origins 1986, 1991, 1996. Every
tuning and calibration uses only the training data of its outer origin
(`tuned()` and `recalibrated()` inside `mt_backtest()`).

* fdm: `hyndman_ullah()` with 6 components, smoothed over age (default) or
  unsmoothed; tuned over the number of components (1, 2, 3, 6) and the
  knots of the smoother (10, 30) by a hold-out of the last 10 training
  years (RMSE of the median log rates), or by 3 rolling inner origins with
  h = 5 (CRPS, one-standard-error rule); tuned and calibrated.
* Lee-Carter: all years, calibrated, and the fitting window (last 15, 20,
  30 or 60 years, `step_window()`) tuned by the 10-year hold-out.

## Log death rates

| Data | Model | CRPS | RMSE | 80% | 95% |
|---|---|---|---|---|---|
| E&W | Lee-Carter | 0.0687 | 0.116 | 0.546 | 0.746 |
| | Lee-Carter, calibrated | 0.0657 | 0.116 | 0.664 | 0.877 |
| | Lee-Carter, window tuned | 0.0573 | 0.0996 | 0.528 | 0.733 |
| | fdm, smoothed | 0.0641 | 0.105 | 0.459 | 0.630 |
| | fdm, unsmoothed | 0.0644 | 0.108 | 0.525 | 0.699 |
| | fdm tuned (hold-out) | 0.0657 | 0.106 | 0.455 | 0.622 |
| | fdm tuned (rolling, one-SE) | 0.0662 | 0.107 | 0.463 | 0.637 |
| | fdm tuned (hold-out), calibrated | 0.0613 | 0.106 | 0.731 | 0.898 |
| France | Lee-Carter | 0.109 | 0.202 | 0.473 | 0.690 |
| | Lee-Carter, calibrated | 0.107 | 0.202 | 0.673 | 0.817 |
| | Lee-Carter, window tuned | 0.106 | 0.190 | 0.477 | 0.659 |
| | fdm, smoothed | 0.0960 | 0.169 | 0.477 | 0.653 |
| | fdm, unsmoothed | 0.0882 | 0.159 | 0.534 | 0.715 |
| | fdm tuned (hold-out) | 0.104 | 0.183 | 0.484 | 0.675 |
| | fdm tuned (rolling, one-SE) | 0.0964 | 0.170 | 0.493 | 0.664 |
| | fdm tuned (hold-out), calibrated | 0.0969 | 0.183 | 0.747 | 0.909 |

## Life expectancy (at 30 for E&W, at birth for France; years)

| Data | Model | CRPS | RMSE | 80% | 95% |
|---|---|---|---|---|---|
| E&W | Lee-Carter | 0.473 | 0.829 | 0.550 | 0.917 |
| | Lee-Carter, calibrated | 0.442 | 0.831 | 0.917 | 1.000 |
| | Lee-Carter, window tuned | 0.287 | 0.532 | 0.783 | 1.000 |
| | fdm, smoothed | 0.585 | 0.812 | 0.083 | 0.250 |
| | fdm tuned (hold-out) | 0.607 | 0.831 | 0.083 | 0.267 |
| | fdm tuned (hold-out), calibrated | 0.495 | 0.870 | 0.500 | 0.867 |
| France | Lee-Carter | 0.575 | 0.908 | 0.167 | 0.667 |
| | Lee-Carter, window tuned | 0.473 | 0.774 | 0.333 | 0.733 |
| | fdm, smoothed | 0.709 | 0.961 | 0.100 | 0.167 |
| | fdm tuned (hold-out) | 0.791 | 1.080 | 0.133 | 0.167 |
| | fdm tuned (hold-out), calibrated | 0.703 | 1.160 | 0.533 | 0.800 |

Settings chosen by the hold-out: fdm 6 or 3 components alternating between
consecutive origins with 10 knots (E&W), 1 or 3 components with 30 knots
(France); Lee-Carter window 15 years at every E&W origin and 2 of 3 French
origins.

## Interpretation

* Tuning the fitting window of Lee-Carter by a hold-out of the last ten
  years helped on both data sets, most for life expectancy (E&W 0.29
  against 0.47 years). The choice was stable across origins: the window
  carries real information about the pace of recent decline (Booth,
  Maindonald and Smith 2002).
* Tuning fdm over components and smoothing by the same hold-out did not
  help and was worse on France (CRPS 0.104 against 0.096). The choice
  alternated between origins one year apart: with one 10-year window,
  the differences between settings were mostly noise of that window.
  Rolling inner origins with the one-standard-error rule were close to the
  default. Tuning pays where the setting changes the forecast
  substantively and consistently (the window), less where settings differ
  by second-order details of the age pattern.
* Smoothing over age (the Hyndman-Ullah method) did not improve fdm here
  and narrowed its intervals. The unsmoothed variant counted the
  observation noise twice (in the resampled model errors and in the
  predictive draws); its extra width improved coverage by accident, as for
  the plain random walk (`dev/calibration.md`).
* All models were too narrow over ten years, fdm the most for life
  expectancy (95% coverage 0.17 to 0.27). Calibration on the last ten
  training years brought the rates to 0.90-0.91 and life expectancy to
  0.80-0.87, and gave the best fdm CRPS on E&W.
* With six and three origins of one population each, the differences are
  indicative; the validation check `tuning_ew` guards the E&W results.
