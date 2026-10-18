#' Estimate grid-resolution data from grouped observations
#'
#' Fills the modelling grid (single years of age by default) from deaths and
#' exposures observed in coarser age groups. The default is the penalized
#' composite link model (PCLM; Eilers 2007; Rizzi, Gampe and Eilers 2015) as
#' published: the estimates are the fitted values of the model, without
#' rescaling.
#'
#' * **Exposures** observed coarser than the grid are ungrouped first, as
#'   counts by the PCLM (`exposure_method = "pclm"`) or by the monotone cubic
#'   interpolation of Hyman (1983) applied to the cumulative counts
#'   (`"monotone"`), which reproduces the group totals exactly.
#' * **Deaths** are ungrouped as **rates** when exposures are available on
#'   the grid: the exposures enter the composite link as offset (Rizzi, Gampe
#'   and Eilers 2015), so cohort-size waves stay in the exposures and the
#'   estimated rate schedule is smooth. Without exposures, deaths are
#'   ungrouped as counts per year of age (1D only; 2D smoothing of death
#'   counts cannot represent cohort waves and is refused).
#' * `dims = 1` fits each stratum and period separately; `dims = 2` fits a
#'   smooth age-period surface of rates per stratum with two smoothing
#'   parameters, the PCLM of rates over age and period with exposures of
#'   Rizzi et al. (2019).
#'
#' **Group totals.** The PCLM reproduces the total over all groups of a fit
#' (the penalty leaves the level of the latent curve free), but not each
#' group total: the gap is reported as `max_rel_diff` in the diagnostics.
#' With `preserve_totals = TRUE` the estimates within each observed group
#' are rescaled to the observed total after the fit, so that aggregation
#' reproduces the data exactly. This rescaling is a post-hoc step of
#' mortensor without a published source: it creates steps at the
#' group boundaries, and with `dims = 2` it puts the year-to-year noise of
#' every group back onto the smooth surface. On England and Wales men
#' (2000-2011, 5-year groups) the median absolute difference of the
#' single-year log rates between the two settings is about 0.002 with
#' `dims = 1` and 0.008 with `dims = 2`.
#'
#' **Open age interval.** For counts (exposures, or deaths without
#' exposures), an open last group is expanded into single years up to the
#' maximum age `omega` (see [mt_data()]) and the PCLM spreads its total over
#' them, as with a last interval of width `omega` minus its lower bound in
#' the ungroup package. For rates, an open interval of the *grid* (for
#' example 100+) is expanded into single years up to `omega` as well, and
#' its exposure is split over these years (`open_split`):
#'
#' * `"pclm"` (default): by the PCLM of the exposures of all grid cells,
#'   with a group of 0 counts of 15 years beyond `omega`, as Rizzi, Gampe
#'   and Eilers (2015) ungroup grouped exposures (their maximum age 115 and
#'   a group 115-130 with 0 counts); the estimates enter the composite link
#'   of the rates;
#' * `"stationary"`: with stationary-population weights `l(x) =
#'   exp(-cumsum(m(x)))` implied by the fitted rates; starting from rates
#'   that rise by 10% per year of age, the weights and the fit are iterated
#'   until no weight changes by more than 1e-6 (at most 50 iterations; the
#'   smoothing parameters are selected again once the weights have
#'   converged). This treatment is specific to mortensor (reference
#'   unknown); it was the default up to E3.
#'
#' On French women and men (1997-2006, deaths in 5-year groups and 90+,
#' grid open at 100+), the rate of 100+ was too high by 0.14 (women) and
#' 0.15 (men) on the log scale with `"stationary"` and by 0.12 and 0.15 with
#' `"pclm"`; the median absolute log error at 90-99 was 0.020 and 0.022
#' with `"stationary"`, 0.019 and 0.022 with `"pclm"`. On data simulated
#' with the exposures of a stationary population above 100 ([mt_example()]),
#' `"stationary"` was closer (median absolute log error at 90+ 0.006,
#' against 0.009 with `"pclm"`), because that generator has the stationary
#' exposures this option assumes.
#'
#' **Relation to Rizzi, Gampe and Eilers (2015).** As in the paper, the
#' observed counts are Poisson with means composed from a latent
#' distribution `exp(B beta)` (times the exposures for rates), estimated by
#' penalised likelihood with a second-order difference penalty on the
#' coefficients weighted by `lambda` (their penalised log-likelihood with
#' `lambda / 2`, by iteratively reweighted least squares), and the effective
#' dimension is the trace of the hat matrix. The departures:
#'
#' * basis: the paper uses one coefficient per fine cell (`B` the identity)
#'   and allows a B-spline basis; here `B` is a cubic B-spline basis, as in
#'   Eilers (2007) and the ungroup package;
#' * smoothing parameter: the paper minimises AIC over a grid of `log10
#'   lambda`; the default here is BIC, the default of the ungroup package,
#'   minimised by golden-section search on `log10 lambda` within
#'   `lambda_range`. `criterion = "aic"` gives the paper's criterion;
#' * open last interval of counts: the paper sets a maximum age (115) and
#'   adds a group with 0 counts beyond it (115-130); here the open group
#'   ends at `omega` without such a group. On French women in 2006 with
#'   exposures in 5-year groups and 90+, adding a group with 0 counts of 15
#'   years beyond `omega = 110` changed the median absolute log error of the
#'   rates at 90-109 from 0.19 (too low) to 0.66 (too high); with exposures
#'   in single years to 89 and 90+, from 0.57 to 0.14. The split of a wide
#'   open group from its total is weakly identified either way;
#' * exposures: as in the paper, they enter the composite link, and grouped
#'   exposures are ungrouped by the PCLM first; the open interval of the
#'   grid as described above (`"pclm"` follows the paper);
#' * iterations: the paper iterates until the coefficients change by less
#'   than 1e-6; here the iterations stop when the largest change is below
#'   1e-6 relative to the coefficients or the penalised deviance changes by
#'   less than 1e-8 relative, with step halving when it rises;
#' * in 2D (Rizzi et al. 2019) see Defaults below.
#'
#' When the deaths are grouped more coarsely than the exposures at old ages
#' (for example deaths 90+ and single-year exposures), the rates within the
#' last group are extrapolated by the smooth latent curve from the lower
#' groups, and the comparison with real data is poor at the highest ages.
#' On French women in 2006 (Human Mortality Database data in the
#' demography package), with deaths in 5-year groups and 90+ and exposures
#' in single years to 110+, the median absolute log error of the ungrouped
#' rates at ages 90-109 against the single-year rates was 0.09 with
#' `dims = 1` (0.02 at ages 30-89). The errors were within 0.1 on the log
#' scale up to age 101, but the rates were too high by 0.3-0.4 at ages
#' 104-109, because the latent curve continues the rise of the rates below
#' 90 while the observed rates rise more slowly at the highest ages. With a grid
#' open at 100+ (exposures 100+), the rate of 100+ was too high by 0.14 on
#' the log scale (median over 1997-2006, see above). With the exposures also
#' grouped in 5-year groups and 90+, the error at 90-109 was 0.19 (rates too
#' low). With `dims = 2` (1997-2006, ages 30 and over)
#' and smoothing parameters chosen by BIC, the error at 90-109 was 1.4: BIC
#' chose very small smoothing parameters (the large counts vary more than
#' Poisson noise) and the rates inside the wide last group oscillated; with
#' `lambda = c(100, 10)` fixed it was 0.05. Treat ungrouped rates inside a
#' wide open group as extrapolations, prefer data with a higher open age,
#' check the smoothing parameters (attribute `"ungroup"`) and the result
#' with [mt_plot_ungroup()].
#'
#' Cells already observed at grid resolution are kept. When deaths are
#' ungrouped as rates and the grid starts with age 0 as its own interval, an
#' observed first year of life is kept and left out of the fit, as in
#' [mt_smooth()]: infant mortality does not follow the age pattern of the
#' following ages, and its peak would spill into ages 1-6.
#'
#' **Defaults.** The basis and penalty follow the defaults of the ungroup
#' package (version 1.6.3) for one dimension: cubic B-splines (`deg = 3`)
#' with a knot every `kr = 2` years, a second-order difference penalty
#' (`pord = 2`) and the smoothing parameter chosen by BIC. The search range
#' of the smoothing parameter is wider than ungroup's (0.1 to 1e5 in one
#' dimension, 0.1 to 1000 in two): `lambda_range = c(1e-2, 1e8)`. In two
#' dimensions ungroup uses fixed smoothing parameters (1, 1) and a knot every
#' 7 cells by default; here both parameters are chosen by BIC with a knot
#' every `kr` years and periods.
#'
#' The methods `"hyman_kannisto"`, `"natural_gompertz"` and `"gam"` are
#' deprecated: they are heuristics of this package without a published
#' source. They still work, with a warning of class
#' `mt_warning_deprecated`.
#'
#' In backtests use `step_ungroup()` so that ungrouping, in particular in 2D,
#' sees only the training window.
#'
#' @seealso [mt_plot_ungroup()] and [mt_ungroup_check()] for diagnostics.
#'
#' @param data An [mt_data] object.
#' @param age Target age axis ([mt_axis]). By default the grid when it was
#'   set (with [mt_set_grid()] or `step_grid()`) and is finer than the data,
#'   otherwise single years of age from the lowest age to the start of the
#'   open interval.
#' @param dims 1 or 2 (for `method = "pclm"`; `"gam"` is always 2D, the
#'   other methods 1D).
#' @param method Method for deaths: `"pclm"` (default), or `"monotone"`, the
#'   monotone cubic interpolation of Hyman (1983) of the cumulative deaths
#'   (counts only, 1D, ignores exposures; group totals reproduced exactly).
#'   Deprecated, without a published source: `"hyman_kannisto"` (monotone
#'   spline of log group rates with a Kannisto tail), `"natural_gompertz"`
#'   (natural spline with a log-linear tail) and `"gam"` (2D tensor-product
#'   GAM of log group rates over age and year with a Kannisto tail, package
#'   mgcv); they need exposures on the grid and rescale to the
#'   exposure-weighted group rates. The vignette "Ungrouping age groups:
#'   methods and diagnostics" compares methods on data with known single
#'   years; [mt_plot_ungroup()] and [mt_ungroup_check()] check a result.
#' @param exposure_method `"pclm"` or `"monotone"` (Hyman 1983).
#' @param lambda Fixed smoothing parameter(s); `NULL` selects them by
#'   `criterion`.
#' @param criterion `"bic"` or `"aic"`.
#' @param kr Knot spacing in years (age and period).
#' @param deg,pord Degree of the B-splines and order of the difference
#'   penalty.
#' @param lambda_range Search range for the smoothing parameters.
#' @param preserve_totals Rescale the estimates to the observed group totals
#'   after the fit (default `FALSE`, the published PCLM); see Details.
#' @param open_split Split of the exposure of an open interval of the grid
#'   when deaths are ungrouped as rates: `"pclm"` (default, the PCLM of the
#'   exposures, Rizzi, Gampe and Eilers 2015) or `"stationary"`
#'   (package-specific stationary-population weights); see Details.
#' @param tail,tail_width For the deprecated `method = "gam"`: the parametric
#'   curve that takes over from the GAM above the start of the last age group
#'   (`"kannisto"`, logistic, or `"gompertz"`, log-linear), fitted to the GAM
#'   rates from age 70, and the width (years) of the logistic blend between
#'   the two.
#' @return An `mt_data` object with deaths (and exposures) on the grid.
#'   Diagnostics per fit are attached as attribute `"ungroup"`: smoothing
#'   parameters, effective dimension, the largest relative difference between
#'   fitted and observed group totals (before any rescaling), the number of
#'   iterations of the open-interval weights (`open_iterations`, `NA` without
#'   an expanded open interval), and a flag when a smoothing parameter is at
#'   the boundary of the search range.
#' @references
#' Eilers PHC (2007). Ill-posed problems with counts, the composite link model
#' and penalized likelihood. *Statistical Modelling* 7(3): 239-254.
#'
#' Hyman JM (1983). Accurate monotonicity preserving cubic interpolation.
#' *SIAM Journal on Scientific and Statistical Computing* 4(4): 645-654.
#'
#' Rizzi S, Gampe J, Eilers PHC (2015). Efficient estimation of smooth
#' distributions from coarsely grouped data. *American Journal of
#' Epidemiology* 182(2): 138-147. \doi{10.1093/aje/kwv020}
#'
#' Rizzi S, Halekoh U, Thinggaard M, Engholm G, Christensen N, Johannesen TB,
#' Lindahl-Jacobsen R (2019). How to estimate mortality trends from grouped
#' vital statistics. *International Journal of Epidemiology* 48(2): 571-582.
#' \doi{10.1093/ije/dyy183}
#' @examples
#' x <- mt_example(years = 2010:2014)
#' u <- mt_ungroup(x)
#' u
#' head(attr(u, "ungroup"))
#' # rescaled to the observed group totals (post-hoc step of the package)
#' ut <- mt_ungroup(x, preserve_totals = TRUE)
#' @export
mt_ungroup <- function(data, age = NULL, dims = 1,
                       method = c("pclm", "monotone", "hyman_kannisto", "natural_gompertz", "gam"),
                       exposure_method = c("pclm", "monotone"),
                       lambda = NULL,
                       criterion = c("bic", "aic"), kr = 2, deg = 3, pord = 2,
                       lambda_range = c(1e-2, 1e8), preserve_totals = FALSE, open_split = c("pclm", "stationary"),
                       tail = c("kannisto", "gompertz"), tail_width = 2) {
  check_class(data, "mt_data", "data")
  method <- match.arg(method)
  open_split <- match.arg(open_split)
  tail <- match.arg(tail)
  check_number(tail_width, "tail_width", lower = 0.1)
  exposure_method <- match.arg(exposure_method)
  criterion <- match.arg(criterion)
  check_number(dims, "dims", lower = 1, upper = 2, integerish = TRUE)
  if (!is.logical(preserve_totals) || length(preserve_totals) != 1L || is.na(preserve_totals)) {
    mt_abort("`preserve_totals` must be TRUE or FALSE.", class = "argument")
  }
  if (method %in% c("hyman_kannisto", "natural_gompertz", "gam")) {
    warn_deprecated_ungroup(sprintf("mt_ungroup(method = \"%s\")", method))
  }
  if (dims == 2 && method %in% c("monotone", "hyman_kannisto", "natural_gompertz")) {
    mt_warn(sprintf("dims = 2 is ignored by method = \"%s\", which works per stratum and period.", method),
            class = "ungroup_dims")
  }
  g <- mt_grid(data)
  if (!is.null(age)) {
    check_class(age, "mt_axis", "age")
    g$age <- age
  } else if (is.null(data$grid) || !any(vapply(data$obs, needs_ungroup, logical(1), g = g))) {
    g$age <- single_year_grid(data)
  }
  ctl <- list(lambda = lambda, criterion = criterion, kr = kr, deg = deg, pord = pord,
              lambda_range = lambda_range, preserve_totals = preserve_totals, infant = TRUE, open_split = open_split)
  diag_rows <- list()
  if (!is.null(data$obs$exposure) && needs_ungroup(data$obs$exposure, g)) {
    res <- ungroup_counts_1d(data$obs$exposure, g, data$omega, ctl, exposure_method, "exposure")
    data$obs$exposure <- res$obs
    diag_rows <- c(diag_rows, list(res$diag))
  }
  if (!is.null(data$obs$deaths) && needs_ungroup(data$obs$deaths, g)) {
    E <- data$obs$exposure
    if (method == "monotone") {
      res <- ungroup_counts_1d(data$obs$deaths, g, data$omega, ctl, "monotone", "deaths")
    } else if (method %in% c("hyman_kannisto", "natural_gompertz", "gam")) {
      if (is.null(E) || needs_ungroup(E, g)) {
        mt_abort(sprintf("method = \"%s\" needs exposures on the grid.", method), class = "ungroup")
      }
      res <- list(obs = ungroup_rates_spline(data$obs$deaths, E, g, data$omega, method, tail, tail_width),
                  diag = data.table(measure = "deaths", stratum = NA_integer_, period = "all", method = method,
                                    at_bound = FALSE))
    } else if (!is.null(E) && !needs_ungroup(E, g)) {
      res <- ungroup_rates(data$obs$deaths, E, g, data$omega, ctl, dims)
    } else {
      if (dims == 2) {
        mt_abort(paste0("2D ungrouping of death counts without exposures is not supported: ",
                        "smoothing counts over time cannot represent cohort-size waves. ",
                        "Add exposures on the grid or use dims = 1."), class = "ungroup")
      }
      mt_warn(paste0("Deaths are ungrouped as counts without exposures; the split of the open ",
                     "interval is weakly identified."), class = "ungroup_counts")
      res <- ungroup_counts_1d(data$obs$deaths, g, data$omega, ctl, "pclm", "deaths")
    }
    data$obs$deaths <- res$obs
    diag_rows <- c(diag_rows, list(res$diag))
  }
  dg <- rbindlist(diag_rows, fill = TRUE)
  if (nrow(dg) && any(dg$at_bound %in% TRUE)) {
    mt_warn(sprintf("The smoothing parameter is at the boundary of the search range in %d of %d fits.",
                    sum(dg$at_bound %in% TRUE), nrow(dg)), class = "lambda_bound")
  }
  if (nrow(dg) && "open_converged" %in% names(dg) && any(dg$open_converged %in% FALSE)) {
    mt_warn(sprintf(paste("The weights of the open age interval did not converge in %d of %d fits (50 iterations);",
                          "the rates in the open interval are approximate."),
                    sum(dg$open_converged %in% FALSE), nrow(dg)), class = "ungroup_open")
  }
  data$grid <- g
  attr(data, "ungroup") <- dg
  data
}

# Warning for the deprecated ungrouping and disaggregation heuristics
warn_deprecated_ungroup <- function(what, alternative = NULL) {
  alternative <- alternative %||% paste("the penalized composite link model (method = \"pclm\") or the monotone cubic",
                                        "interpolation of cumulative counts of Hyman (1983) (method = \"monotone\")")
  mt_warn(sprintf("%s is deprecated: it is a heuristic of this package without a published source. Use %s.",
                  what, alternative), class = "deprecated", call = NULL)
}

needs_ungroup <- function(obs, g) {
  idx <- obs_grid_index(obs, g)
  any(!is.na(idx$a1) & (idx$a1 != idx$a2 | idx$p1 != idx$p2))
}

# Counts (exposures, or deaths without exposures), 1D per stratum and period
ungroup_counts_1d <- function(obs, g, omega, ctl, method, measure) {
  idx <- obs_grid_index(obs, g, measure)
  obs <- cbind(obs, idx)[!is.na(a1)]
  if (any(obs$p1 != obs$p2)) {
    mt_abort(sprintf("The %s table has periods coarser than the grid; 1D ungrouping needs grid periods.", measure),
             class = "ungroup")
  }
  fa <- fine_axis(g$age, omega)
  x <- (fa$lower + fa$upper) / 2
  bas <- pclm_basis(x, min(fa$lower), max(fa$upper), ctl$kr, ctl$deg, ctl$pord)
  width <- fa$upper - fa$lower
  open_cells <- is.infinite(obs$age_hi) & obs$a1 != obs$a2
  if (any(open_cells)) {
    mt_warn(sprintf(paste0("The open %s interval from age %s is split into %d grid cells from its total alone; ",
                           "the age pattern inside it is extrapolated and can be badly wrong. Prefer data with ",
                           "a higher open age or a grid whose open interval starts at %s."),
                    measure, min(obs$age_lo[open_cells]), max(obs$a2 - obs$a1 + 1), min(obs$age_lo[open_cells])),
            class = "ungroup_open_counts", call = NULL)
  }
  out <- list()
  dg <- list()
  for (key in split(seq_len(nrow(obs)), list(obs$stratum, obs$p1), drop = TRUE)) {
    o <- obs[key]
    hi <- pmin(o$age_hi, omega)
    if (all(o$a1 == o$a2)) {
      out[[length(out) + 1L]] <- o[, list(stratum, a = a1, p = p1, value, provenance)]
      next
    }
    C <- link_matrix(o$age_lo, hi, fa$lower, fa$upper)
    cov <- colSums(C) > 0
    if (method == "monotone") {
      est <- monotone_split(o$age_lo, hi, o$value, fa$lower[cov], fa$upper[cov])
      fit <- list(lambda = NA_real_, ed = NA_real_, at_bound = FALSE, mu = o$value)
    } else {
      fit <- pclm_fit(o$value, C[, cov, drop = FALSE], width[cov], bas$B[cov, , drop = FALSE], list(bas$P),
                      ctl$lambda, ctl$criterion, ctl$lambda_range)
      est <- width[cov] * fit$gamma
    }
    rel <- max(abs(fit$mu / pmax(o$value, 1e-12) - 1)[o$value > 0], 0)
    if (ctl$preserve_totals) est <- rescale_to_groups(est, C[, cov, drop = FALSE], o$value)
    gi <- fa$grid[cov]
    v <- tapply(est, gi, sum)
    exact <- o[a1 == a2]
    res <- data.table(stratum = o$stratum[1], a = as.integer(names(v)), p = o$p1[1], value = as.numeric(v),
                      provenance = method)
    res[exact, on = c(a = "a1"), `:=`(value = i.value, provenance = i.provenance)]
    out[[length(out) + 1L]] <- res
    dg[[length(dg) + 1L]] <- data.table(measure = measure, stratum = o$stratum[1], period = g$period$labels[o$p1[1]],
                                        method = method, lambda_age = fit$lambda[1], lambda_period = NA_real_,
                                        ed = fit$ed,
                                        max_rel_diff = rel, at_bound = isTRUE(fit$at_bound))
  }
  list(obs = grid_table(rbindlist(out), g), diag = rbindlist(dg))
}

# Deaths as rates with grid exposures as offset. With ctl$infant, age 0
# observed as its own grid cell is kept and left out of the fit (fine axis and
# basis from age 1), as in mt_smooth()
ungroup_rates <- function(D, E, g, omega, ctl, dims) {
  idx <- obs_grid_index(D, g, "deaths")
  D <- cbind(D, idx)[!is.na(a1)]
  Earr <- obs_to_matrix(E, g)
  fa <- fine_axis(g$age, omega)
  xa <- (fa$lower + fa$upper) / 2
  n_open <- sum(fa$grid == length(g$age) & fa$open)
  bas <- pclm_basis(xa, min(fa$lower), max(fa$upper), ctl$kr, ctl$deg, ctl$pord)
  sep0 <- isTRUE(ctl$infant) && g$age$lower[1] == 0 && g$age$upper[1] == 1 && length(g$age) > 3
  if (sep0) {
    k1 <- fa$grid > 1L
    fa1 <- list(lower = fa$lower[k1], upper = fa$upper[k1], grid = fa$grid[k1], open = fa$open)
    bas1 <- pclm_basis(xa[k1], min(fa1$lower), max(fa1$upper), ctl$kr, ctl$deg, ctl$pord)
  }
  # rows of age 0 set aside: only when age 0 is never part of a wider cell
  infant <- function(o) {
    inf <- sep0 & o$a1 == 1L & o$a2 == 1L & o$p1 == o$p2
    if (all(inf) || any(o$a1 == 1L & !inf)) inf[] <- FALSE
    inf
  }
  out <- list()
  dg <- list()
  for (s in unique(D$stratum)) {
    Ds <- D[stratum == s]
    if (dims == 1) {
      if (any(Ds$p1 != Ds$p2)) {
        mt_abort("Deaths have periods coarser than the grid; use dims = 2.", class = "ungroup")
      }
      for (p in sort(unique(Ds$p1))) {
        o <- Ds[p1 == p]
        e_grid <- Earr[, p, s]
        inf <- infant(o)
        if (any(inf)) {
          r <- ungroup_rate_block(o[!inf], list(e_grid), fa1, bas1$B, list(bas1$P), ctl, omega, n_open)
        } else {
          r <- ungroup_rate_block(o, list(e_grid), fa, bas$B, list(bas$P), ctl, omega, n_open)
        }
        out[[length(out) + 1L]] <- data.table(stratum = s, a = c(o$a1[inf], r$a), p = p,
                                              value = c(o$value[inf], r$deaths),
                                              provenance = c(o$provenance[inf], r$provenance))
        dg[[length(dg) + 1L]] <- data.table(measure = "deaths", stratum = s, period = g$period$labels[p],
                                            method = paste0(ctl$method_label %||% "pclm_rate", "_1d"),
                                            lambda_age = r$fit$lambda[1],
                                            lambda_period = NA_real_, ed = r$fit$ed,
                                            max_rel_diff = r$rel, at_bound = isTRUE(r$fit$at_bound),
                                            open_iterations = r$open_iterations, open_converged = r$open_converged)
      }
    } else {
      ps <- seq(min(Ds$p1), max(Ds$p2))
      xp <- (g$period$lower[ps] + g$period$upper[ps]) / 2
      bp <- pclm_basis(xp, g$period$lower[min(ps)], g$period$upper[max(ps)], ctl$kr, ctl$deg, ctl$pord)
      inf <- infant(Ds)
      fb <- if (any(inf)) list(fa = fa1, bas = bas1) else list(fa = fa, bas = bas)
      B2 <- kronecker(bp$B, fb$bas$B)
      pen <- list(kronecker(diag(ncol(bp$B)), fb$bas$P), kronecker(bp$P, diag(ncol(fb$bas$B))))
      r <- ungroup_rate_block(Ds[!inf], lapply(ps, function(p) Earr[, p, s]), fb$fa, B2, pen, ctl, omega, n_open,
                              ps)
      out[[length(out) + 1L]] <- data.table(stratum = s, a = c(Ds$a1[inf], r$a), p = c(Ds$p1[inf], r$p),
                                            value = c(Ds$value[inf], r$deaths),
                                            provenance = c(Ds$provenance[inf], r$provenance))
      dg[[length(dg) + 1L]] <- data.table(measure = "deaths", stratum = s, period = "all",
                                          method = paste0(ctl$method_label %||% "pclm_rate", "_2d"),
                                          lambda_age = r$fit$lambda[1],
                                          lambda_period = r$fit$lambda[2],
                                          ed = r$fit$ed, max_rel_diff = r$rel, at_bound = isTRUE(r$fit$at_bound),
                                          open_iterations = r$open_iterations, open_converged = r$open_converged)
    }
  }
  list(obs = grid_table(rbindlist(out), g), diag = rbindlist(dg, fill = TRUE))
}

# One stratum: a block of periods (one for 1D). e_list: grid exposures per period.
ungroup_rate_block <- function(o, e_list, fa, B, pen, ctl, omega, n_open, ps = NULL) {
  ps <- ps %||% o$p1[1]
  nf <- length(fa$lower)
  ng <- max(fa$grid)
  # composite link: rows = observed cells; columns = fine cells x periods
  rows <- lapply(seq_len(nrow(o)), function(i) {
    cols_a <- which(fa$lower >= o$age_lo[i] - 1e-9 & fa$upper <= pmin(o$age_hi[i], omega) + 1e-9)
    cols_p <- match(seq(o$p1[i], o$p2[i]), ps)
    as.vector(outer(cols_a, (cols_p - 1L) * nf, `+`))
  })
  C <- Matrix::sparseMatrix(i = rep(seq_along(rows), lengths(rows)), j = unlist(rows), x = 1,
                            dims = c(nrow(o), nf * length(ps)))
  e_grid <- do.call(cbind, e_list)
  if (anyNA(e_grid)) mt_abort("Exposures are missing for grid cells needed by the ungrouping.", class = "ungroup")
  open_idx <- which(fa$grid == ng & fa$open)
  w <- matrix(1, nf, length(ps))
  stationary <- length(open_idx) && identical(ctl$open_split, "stationary")
  if (length(open_idx)) {
    if (stationary) {
      # starting weights: rates rising by 10% per year of age in the open interval
      w[open_idx, ] <- open_weights(sum(o$value) / sum(e_grid) * exp(0.1 * seq_along(open_idx)))
    } else {
      # the exposure of the open interval split by the PCLM of the exposures
      # of the grid cells (Rizzi, Gampe and Eilers 2015)
      for (j in seq_along(ps)) w[open_idx, j] <- open_weights_pclm(e_grid[, j], fa, ctl)
    }
  }
  fit_w <- function(w, lambda, beta) {
    off <- as.vector(e_grid[fa$grid, , drop = FALSE] * w)
    cov <- Matrix::colSums(C) > 0 & off > 0
    fit <- pclm_fit(o$value, C[, cov, drop = FALSE], off[cov], B[cov, , drop = FALSE], pen, lambda,
                    ctl$criterion, ctl$lambda_range, beta = beta)
    list(fit = fit, off = off, cov = cov)
  }
  r <- fit_w(w, ctl$lambda, NULL)
  open_it <- NA_integer_
  open_conv <- NA
  if (stationary) {
    # the exposure of the open interval is split with the stationary weights
    # of the fitted rates: weights and fit are iterated to a fixed point, with
    # the smoothing parameters fixed; they are then selected again with the
    # converged weights, until the selection no longer changes
    open_it <- 0L
    open_conv <- TRUE
    select <- is.null(ctl$lambda)
    at_bound <- isTRUE(r$fit$at_bound)
    for (round in seq_len(if (select) 3L else 1L)) {
      lam <- r$fit$lambda
      repeat {
        rate_m <- matrix(exp(as.numeric(B %*% r$fit$beta)), nf, length(ps))
        w_new <- w
        for (j in seq_along(ps)) w_new[open_idx, j] <- open_weights(rate_m[open_idx, j])
        if (max(abs(w_new - w)) < 1e-6) break
        if (open_it >= 50L) {
          open_conv <- FALSE
          break
        }
        open_it <- open_it + 1L
        w <- w_new
        r <- fit_w(w, lam, r$fit$beta)
      }
      if (!select || !open_conv || round == 3L) break
      r2 <- fit_w(w, NULL, r$fit$beta)
      at_bound <- isTRUE(r2$fit$at_bound)
      if (all(abs(log10(r2$fit$lambda / lam)) < 0.01)) break
      r <- r2
    }
    r$fit$at_bound <- at_bound
  }
  fit <- r$fit
  off <- r$off
  cov <- r$cov
  deaths_f <- rep(0, length(off))
  deaths_f[cov] <- off[cov] * fit$gamma
  rel <- max(abs(fit$mu / pmax(o$value, 1e-12) - 1)[o$value > 0], 0)
  if (ctl$preserve_totals) deaths_f[cov] <- rescale_to_groups(deaths_f[cov], C[, cov, drop = FALSE], o$value)
  dm <- matrix(deaths_f, nf, length(ps))
  dgrid <- rowsum(dm, fa$grid)
  covered <- rowsum(matrix(as.numeric(Matrix::colSums(C) > 0), nf, length(ps)), fa$grid) > 0
  # grid cells of the fine axis (from age 1 when age 0 is set aside)
  ga <- sort(unique(fa$grid))
  a <- rep(ga, length(ps))
  p <- rep(ps, each = length(ga))
  keep <- as.vector(covered)
  prov <- rep(ctl$provenance %||% "pclm", length(a))
  vals <- as.vector(dgrid)
  if (!isFALSE(ctl$keep_observed)) {
    exact <- o[a1 == a2 & p1 == p2]
    m <- match(paste(exact$a1, exact$p1), paste(a, p))
    vals[m] <- exact$value
    prov[m] <- "observed"
  }
  list(a = a[keep], p = p[keep], deaths = vals[keep], provenance = prov[keep], fit = fit, rel = rel,
       open_iterations = open_it, open_converged = open_conv)
}

# Weights of the single years of the open grid interval from the PCLM of
# the exposures of all grid cells of one period (exposures as counts, offset
# the widths of the fine cells, the open cell expanded to omega), as Rizzi,
# Gampe and Eilers (2015) ungroup grouped exposures before they enter the
# composite link of the rates
open_weights_pclm <- function(e, fa, ctl) {
  open_idx <- which(fa$grid == max(fa$grid) & fa$open)
  cells <- sort(unique(fa$grid))
  if (!(e[max(cells)] > 0)) return(rep(1 / length(open_idx), length(open_idx)))
  # single years from omega on, in a group with 0 counts
  zb <- zero_bin(max(fa$upper))
  lower <- c(fa$lower, zb$lower)
  upper <- c(fa$upper, zb$upper)
  grid <- c(fa$grid, rep(max(cells) + 1L, length(zb$lower)))
  C <- 1 * outer(c(cells, max(cells) + 1L), grid, "==")
  width <- upper - lower
  bas <- pclm_basis((lower + upper) / 2, min(lower), max(upper), ctl$kr, ctl$deg, ctl$pord)
  fit <- pclm_fit(c(e[cells], 0), C, width, bas$B, list(bas$P), ctl$lambda, ctl$criterion, ctl$lambda_range)
  est <- (width * fit$gamma)[open_idx]
  est / sum(est)
}

# Group with 0 counts beyond the maximum age omega that closes an open last
# interval of counts, as in Rizzi, Gampe and Eilers (2015): they set the
# maximum age to 115 and add a group 115-130 with 0 counts. Its single years
# from omega to omega + width
zero_bin <- function(omega, width = 15) {
  lo <- seq(omega, omega + width - 1)
  list(lower = lo, upper = lo + 1)
}

rescale_to_groups <- function(est, C, y) {
  C <- as.matrix(C)
  fitted <- as.numeric(C %*% est)
  f <- ifelse(fitted > 0, y / fitted, 0)
  grp <- apply(C, 2, function(col) if (any(col > 0)) which(col > 0)[1] else NA_integer_)
  est * ifelse(is.na(grp), 1, f[grp])
}

monotone_split <- function(lo, hi, value, fine_lo, fine_hi) {
  ord <- order(lo)
  lo <- lo[ord]
  hi <- hi[ord]
  value <- value[ord]
  if (any(abs(hi[-length(hi)] - lo[-1]) > 1e-9)) {
    mt_abort("Monotone ungrouping needs contiguous groups.", class = "ungroup")
  }
  xs <- c(lo, hi[length(hi)])
  cs <- c(0, cumsum(value))
  f <- stats::splinefun(xs, cs, method = "hyman")
  pmax(diff(f(c(fine_lo, fine_hi[length(fine_hi)]))), 0)
}

obs_to_matrix <- function(obs, g) {
  idx <- obs_grid_index(obs, g)
  ok <- !is.na(idx$a1) & idx$a1 == idx$a2 & idx$p1 == idx$p2
  arr <- array(NA_real_, c(length(g$age), length(g$period), max(obs$stratum)))
  arr[cbind(idx$a1[ok], idx$p1[ok], obs$stratum[ok])] <- obs$value[ok]
  arr
}

grid_table <- function(x, g) {
  out <- x[, list(stratum, age_lo = g$age$lower[a], age_hi = g$age$upper[a], period_lo = g$period$lower[p],
                  period_hi = g$period$upper[p], value, provenance)]
  setorderv(out, c("stratum", "period_lo", "age_lo"))
  out[]
}
