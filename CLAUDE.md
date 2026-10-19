# mortensor: notes for Claude Code sessions

Read first: `dev/progress.md` (status per work package), the work package row
in `dev/implementation-plan.md`, and `dev/revamp-plan.md` for the design.

## Layout

| Path | Content |
|---|---|
| `R/data-*.R` | Axes, `mt_data` container, regrouping, validation, exposures, adapters (HMD, ONS, GENESIS), simulated examples, `ons_data` documentation |
| `data/`, `data-raw/` | `ons_data`, `ons_cause_groups` (`data-raw/ons_data.R`); `de_deaths`, `de_population`, `de_cause_groups` (`data-raw/de_data.R`) |
| `R/harmonise-*.R` | PCLM engine and `mt_ungroup()`, smoothing and preprocessing steps |
| `R/model-*.R` | Model specifications, fitting, forecast draws |
| `R/hierarchy.R` | Hierarchies of strata, aggregation, reconciliation |
| `R/eval-*.R` | Backtests and scoring |
| `R/derive-*.R` | Life tables and derived quantities |
| `dev/` | Plans, review, progress log, calibration study, package comparison script (not built into the package) |

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
- Plot snapshots (vdiffr) are made with ggplot2 >= 4.0; regenerate them
  after an intended change or a ggplot2 release (`testthat::snapshot_review()`).
- Precomputed vignette results: `vignettes/comparison-results.rds` comes
  from `dev/compare-packages.R` (about one hour); the bamp vignette from
  `dev/precompute-vignettes.R`.
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
R CMD INSTALL . && Rscript -e 'lintr::lint_package()'   # lint
R CMD build . && R CMD check --as-cran mortensor_*.tar.gz
```

## Workflow

- CI budget: pushes to development branches start no workflows. Run the
  unit tests and lint locally before each commit. Pull requests (when
  ready for review) and pushes to main run the Ubuntu checks; Windows and
  R-devel run monthly; macOS only by hand (`workflow_dispatch` with
  `platforms = all`).

- One work package per commit or pull request; update `dev/progress.md` and
  `NEWS.md` in the same change.
- Documentation is part of every feature. After a feature works (tests,
  lint), review and update in the same change: roxygen docs
  (`roxygenise()`), README, NEWS, the vignettes that use or should use it
  (re-knit precomputed ones with `dev/precompute-vignettes.R`), the function
  map (`vignettes/articles/function-map.Rmd`), the FAQ, and `_pkgdown.yml`
  (reference index and articles). Then build the pkgdown site locally
  (`pkgdown::build_site()`, offline: `options(pkgdown.internet = FALSE)`)
  and check that every exported function is indexed and every page
  renders. The site is deployed from `main` by `.github/workflows/pkgdown.yaml`.
- The package is public and labelled as experimental, vibe-coded software
  (README notice and lifecycle badge, DESCRIPTION, package help, startup
  message). Keep these notices.
- Correctness first: run the unit tests; new methods get tests that check
  them against the published method or a reference implementation
  (published numbers, the authors' code, an established package), and
  exact identities where they hold. Forecast performance (coverage, CRPS,
  accuracy) depends on the data and is checked by the user with
  backtests; the repository has no performance benchmarks with
  thresholds.
- Data from other packages (HMD extracts in StMoMo and demography) is loaded
  in tests and vignettes, never copied into the repository. The owner's
  German data stays outside the repository. Exceptions, with the owner's
  agreement: `ons_data` (ONS, Open Government Licence) and the German
  regional data `de_deaths`/`de_population` (Destatis, Datenlizenz
  Deutschland by-2-0), built by the scripts in `data-raw/`.
