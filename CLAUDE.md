# mortensor: notes for Claude Code sessions

Read first: `dev/progress.md` (status per work package), the work package row
in `dev/implementation-plan.md`, and `dev/revamp-plan.md` for the design.

## Layout

| Path | Content |
|---|---|
| `R/data-*.R` | Axes, `mt_data` container, regrouping, validation, exposures, adapters, simulated example |
| `R/harmonise-*.R` | PCLM engine and `mt_ungroup()`, smoothing and preprocessing steps |
| `R/model-*.R` | Model specifications, fitting, forecast draws |
| `R/eval-*.R` | Backtests and scoring |
| `R/derive-*.R` | Life tables and derived quantities |
| `R/validation.R` | Benchmarks shared by slow tests and `validation/` |
| `validation/` | Validation pipeline (thresholds, runner, history); see its README |
| `dev/` | Plans, review, progress log, figures (not built into the package) |

## Conventions

- R >= 4.3.0. Imports: data.table, Matrix, methods, parallel, splines, stats,
  utils. Everything else is in Suggests and guarded with `need_pkg()` or
  `requireNamespace()`; tests use `skip_if_not_installed()`.
- Exported verbs carry the `mt_` prefix; model constructors are nouns
  (`lee_carter()`, `rwd()`, `tucker3()`, `apc_bayes()`). A test checks that
  exports do not clash with base, stats, utils, data.table, dplyr, forecast,
  generics, demography or StMoMo.
- Errors and warnings: `mt_abort(msg, class = "<topic>")` and
  `mt_warn(msg, class = "<topic>")` give classes `mt_error_<topic>` and
  `mt_warning_<topic>`; tests check classes, not messages. Messages only via
  `mt_inform()` (silent unless `options(mortensor.verbose = TRUE)`).
- data.table inside the package: column names used in non-standard
  evaluation go into `utils::globalVariables()` in `R/mortensor-package.R`.
  Never name a local variable like a column used in the same `[` call.
- Randomness: functions that simulate take `seed = NULL` and use
  `with_seed()`; parallel code uses `parallel` with L'Ecuyer streams.
- Code adapted from other packages carries a header with source, version,
  licence and authors.
- No `print()`/`cat()` for messages outside print methods; plots return
  ggplot objects; no files written outside `tempdir()`.

## Commands

```sh
Rscript -e 'roxygen2::roxygenise()'                 # docs and NAMESPACE
Rscript -e 'devtools::test()'                       # unit tests
MORTENSOR_SLOW_TESTS=true Rscript -e 'devtools::test()'   # with benchmarks
Rscript validation/run_validation.R --load-all      # quick validation suite
R CMD build . && R CMD check --as-cran mortensor_*.tar.gz
```

## Workflow

- One work package per commit or pull request; update `dev/progress.md` and
  `NEWS.md` in the same change.
- New statistical code gets a validation check or a slow test with a
  threshold, not only unit tests.
- Data from other packages (HMD extracts in StMoMo and demography) is loaded
  in tests and vignettes, never copied into the repository. The owner's
  German data stays outside the repository.
