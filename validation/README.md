# Validation pipeline

A fixed set of accuracy, calibration and runtime benchmarks run against
thresholds, so that every iteration of the package is checked against the
same evidence. It complements the unit tests (`tests/testthat`), which check
behaviour; the validation checks statistical quality on Human Mortality
Database data (through `StMoMo::EWMaleData` and `demography::fr.mort`) and on
simulated data with known truth.

| File | Content |
|---|---|
| `checks.R` | The checks; each returns metrics as a data frame (`metric`, `value`) |
| `thresholds.csv` | Pass criteria per check and metric (`<=` or `>=`) with a note on the value when the threshold was set |
| `run_validation.R` | Runner: selects a suite, runs the checks, compares with thresholds, writes `results.csv` and `report.md`, exits with status 1 on a breach |
| `history.csv` | Results recorded with `--record` (one row per metric and run), committed with the change that produced them |

## Running

```sh
# quick suite on the working tree (about a minute)
Rscript validation/run_validation.R --load-all
# full suite on the installed package, recorded in history.csv
R CMD INSTALL . && Rscript validation/run_validation.R --suite full --record
# check that a breach makes the pipeline fail
Rscript validation/run_validation.R --load-all --selftest
```

The `unit_tests` check of the full suite runs the test files of the working
tree against the installed package, so install the current tree first and
do not change it during the run.

GitHub Actions (`.github/workflows/validation.yaml`) runs the quick suite on
every push and pull request and the full suite every Monday and on demand
(workflow dispatch). The report is shown in the job summary and uploaded as
an artifact.

## Changing thresholds

Thresholds encode the accepted quality. Tighten them when a change improves
a metric; loosen them only with a written reason in the `note` column and in
the pull request.
