# Progress log

Status per work package of `dev/implementation-plan.md`. Branch for this
work: `claude/vigilant-albattani-2ssw1x` (plays the role of `dev-0.2`).

| WP | Status | Commit | Notes |
|---|---|---|---|
| P0.1, P0.2, P0.5 | done | 41058ba | 0.1.1: generics and `%\|\|%` fixed, GPL-3, baseline in `dev/baseline/` |
| P0.3 | done | this branch | Workflows: R-CMD-check (R 4.3, oldrel, release, devel; macOS, Windows), test-coverage, lint, validation, pkgdown |
| P0.4 | dropped | | 0.1.x vignettes not revived; the clean break replaces them |
| P0.6 | done | this branch | `CLAUDE.md`, this file |
| Q0.1 | done | this branch | `validation/` pipeline with quick and full suites |
| D1.1 | done | this branch | Skeleton, conditions, checks, name-clash test |
| D1.2 | done | this branch | `mt_axis`, `parse_age_labels()` |
| D1.3 | partial | this branch | Strata table; hierarchy sums checked in `mt_validate()`; reconciliation later (P5.4) |
| D1.4, D1.5 | done | this branch | `mt_data()`, `mt_add_*()`, grid, `mt_array()`, `mt_window()`, `mt_granularity()`, `mt_validate()` |
| D1.6 | done | this branch | `mt_regroup()`, `mt_operator()` |
| D1.7 | done | this branch | `mt_as_exposure()` (mid-year, January 1) |
| D1.8 | done | this branch | Native PCLM (1D, 2D), rate-based with stationary weights in the open interval |
| D1.9 | done | this branch | `mt_ungroup()`; E&W benchmark 0.015 / 0.024 (1D), 0.015 / 0.023 (2D) |
| D1.12 | partial | this branch | `as_mt_data()` for data frames, StMoMoData, demogdata; `mt_example()` |

## Decisions taken during implementation

- 2026-10-12: Exported verbs use the `mt_` prefix; `mt_obs()` returns an
  observation table (the class name in the plan is not used as a class).
- 2026-10-12: The open age interval is expanded into single years up to
  omega inside the PCLM; for rates the exposure of the open interval is split
  with stationary-population weights from the current fit (three outer
  iterations). This replaces the artificial bin of `ungroup`.
- 2026-10-12: Default smoothing selection BIC on log10 lambda in [-2, 8],
  knots every 2 years, cubic B-splines, second-order differences
  (alternatives tested on E&W: third-order differences and wider knots gave
  the same accuracy within 0.003).
- 2026-10-12: 2D ungrouping takes about 90 s for 71 ages x 32 years; it is
  an option, 1D per year is the default.

## Open questions

See `dev/implementation-plan.md` Section 17.
