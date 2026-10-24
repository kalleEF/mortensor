#' Hierarchies of strata
#'
#' Describes how strata add up, for example women and men to both sexes,
#' regions to the country and deprivation quintiles to all quintiles. Several
#' variables give a crossed (grouped) structure: with `sex`, `region` and
#' `quintile` every combination of a level or a total is a node, from the
#' bottom series (one sex, one region, one quintile) to the national total
#' for both sexes. Intermediate groups of levels, such as federal states
#' grouped into East and West, are given with `groups`.
#'
#' Deaths add up over every hierarchy variable. Exposures add up over
#' variables that divide the population (sex, regions, deprivation
#' quintiles: `additive = "deaths"`, the default) and are shared by the
#' levels of variables that do not (causes of death: `additive = "rates"`,
#' whose rates add up to the all-cause rate). Variables of both kinds can be
#' crossed, for example regions and causes: deaths are summed over both,
#' exposures over regions only, and forecasts are reconciled on the deaths
#' scale with these exposures.
#'
#' @param ... One argument per aggregated stratum variable: the label of its
#'   total level (for example `region = "Germany"`), or a list with `total`,
#'   optionally `groups` (a named list of intermediate groups, each a vector
#'   of bottom levels) and `additive`.
#' @param additive Default additivity of the variables: `"deaths"` or
#'   `"rates"`.
#' @return An object of class `mt_hierarchy`.
#' @seealso [mt_aggregate_strata()], [mt_validate()], [mt_reconcile()],
#'   [reconciled()].
#' @examples
#' h <- mt_hierarchy(sex = "both", region = list(total = "DE", groups = list(East = c("BE", "SN"))),
#'                   quintile = "all")
#' h
#' # regions crossed with causes of death, which share the exposure
#' mt_hierarchy(region = "DE", cause = list(total = "all", additive = "rates"))
#' @export
mt_hierarchy <- function(..., additive = c("deaths", "rates")) {
  if (!missing(additive) && length(additive) != 1) {
    mt_abort(paste("`additive` must be one value; give the additivity of single variables in their list, for",
                   "example cause = list(total = \"all\", additive = \"rates\")."), class = "hierarchy")
  }
  additive <- match.arg(additive)
  specs <- list(...)
  if (!length(specs) || is.null(names(specs)) || any(!nzchar(names(specs)))) {
    mt_abort("mt_hierarchy() needs named arguments, one per stratum variable.", class = "hierarchy")
  }
  vars <- lapply(specs, function(s) {
    if (is.character(s) && length(s) == 1) return(list(total = s, groups = list(), additive = additive))
    if (is.list(s) && !is.null(s$total)) {
      a <- s$additive %||% additive
      if (!a %in% c("deaths", "rates") || length(a) != 1) {
        mt_abort("`additive` of a variable must be \"deaths\" or \"rates\".", class = "hierarchy")
      }
      g <- s$groups %||% list()
      gn <- names(g)
      bad <- !is.list(g) || (length(g) && (is.null(gn) || any(!nzchar(gn)) || anyDuplicated(gn) > 0 ||
                                             any(vapply(g, function(m) !is.atomic(m) || !length(m) || anyNA(m), TRUE))))
      if (bad) mt_abort("`groups` must be a named list of non-empty vectors of bottom levels.", class = "hierarchy")
      if (any(as.character(unlist(g)) %in% c(s$total, gn)) || s$total %in% gn) {
        mt_abort("Groups hold bottom levels: not the total, not other groups, and no group named like the total.",
                 class = "hierarchy")
      }
      return(list(total = s$total, groups = g, additive = a))
    }
    mt_abort("Each hierarchy variable needs a total label or list(total, groups, additive).", class = "hierarchy")
  })
  add <- vapply(vars, `[[`, "", "additive")
  # deaths add up over every variable; exposures add up over the variables
  # with additive = "deaths" and are shared by the levels of the others.
  # Hierarchies of shared-exposure variables only are reconciled on rates.
  structure(list(vars = vars, additive = if (all(add == "rates")) "rates" else "deaths"), class = "mt_hierarchy")
}

#' @export
print.mt_hierarchy <- function(x, ...) {
  cat("<mt_hierarchy>\n")
  for (v in names(x$vars)) {
    g <- x$vars[[v]]$groups
    a <- x$vars[[v]]$additive %||% x$additive
    cat(sprintf("  %s: total \"%s\"%s; %s\n", v, x$vars[[v]]$total,
                if (length(g)) sprintf("; groups %s", paste(names(g), collapse = ", ")) else "",
                if (a == "deaths") "deaths and exposures add up" else "deaths add up, exposures are shared"))
  }
  invisible(x)
}

as_hierarchy <- function(h) {
  if (inherits(h, "mt_hierarchy")) return(h)
  if (is.list(h) && !is.null(h$var) && !is.null(h$total)) {
    args <- stats::setNames(list(h$total), h$var)
    return(do.call(mt_hierarchy, c(args, list(additive = h$additive %||% "deaths"))))
  }
  mt_abort("`hierarchy` must be an mt_hierarchy() object.", class = "hierarchy")
}

# Membership of nodes in terms of bottom strata: S [node, bottom]. With
# exposure = TRUE, the membership for exposures: variables whose levels share
# one exposure (additive = "rates") contribute one representative level
membership <- function(nodes, bottom, h, exposure = FALSE) {
  vars <- names(h$vars)
  others <- setdiff(intersect(names(nodes), names(bottom)), c("stratum", "label", vars))
  members <- function(v, level) {
    hv <- h$vars[[v]]
    if (identical(level, hv$total)) return(NULL)
    if (level %in% names(hv$groups)) return(hv$groups[[level]])
    level
  }
  shared <- vars[vapply(vars, function(v) identical(h$vars[[v]]$additive %||% h$additive, "rates"), logical(1))]
  S <- matrix(0, nrow(nodes), nrow(bottom))
  for (i in seq_len(nrow(nodes))) {
    ok <- rep(TRUE, nrow(bottom))
    for (v in vars) {
      m <- members(v, nodes[[v]][i])
      if (exposure && v %in% shared) {
        lev <- unique(bottom[[v]])
        m <- (if (is.null(m)) lev else intersect(lev, m))[1]
      }
      if (!is.null(m)) ok <- ok & bottom[[v]] %in% m
    }
    for (o in others) ok <- ok & bottom[[o]] == nodes[[o]][i]
    S[i, ok] <- 1
  }
  S
}

# Exposures of the aggregate strata in the data against those the hierarchy
# implies (sums of the parts, or the shared exposure of additive = "rates")
check_hierarchy_exposure <- function(data, hs) {
  if (is.null(data$obs$exposure) || all(hs$is_bottom)) return(invisible())
  E <- tryCatch(mt_array(data, "exposure"), error = function(e) NULL)
  if (is.null(E)) return(invisible())
  d <- dim(E)
  Eb <- matrix(E[, , hs$bottom, drop = FALSE], d[1] * d[2])
  agg <- which(!hs$is_bottom)
  implied <- Eb %*% t(hs$Sx[agg, , drop = FALSE])
  observed <- matrix(E[, , agg, drop = FALSE], d[1] * d[2])
  ok <- is.finite(implied) & is.finite(observed) & implied > 0
  if (!any(ok)) return(invisible())
  dev <- stats::median(abs(observed[ok] / implied[ok] - 1))
  if (dev > 0.01) {
    mt_warn(sprintf(paste("The exposures of the aggregate strata differ from those the hierarchy implies by %.0f%%",
                          "(median): exposures add up over variables with additive = \"deaths\" and are shared",
                          "over variables with additive = \"rates\". For causes of death use",
                          "list(total = ..., additive = \"rates\") in mt_hierarchy()."), 100 * dev),
            class = "hierarchy_exposure")
  }
  invisible()
}

# Structure of the strata present: S [stratum row, bottom stratum], positions of the bottom rows
hierarchy_structure <- function(strata, h) {
  h <- as_hierarchy(h)
  vars <- names(h$vars)
  miss <- setdiff(vars, names(strata))
  if (length(miss)) mt_abort(sprintf("Hierarchy variable(s) not in the strata: %s.", paste(miss, collapse = ", ")),
                             class = "hierarchy")
  agg_levels <- lapply(vars, function(v) c(h$vars[[v]]$total, names(h$vars[[v]]$groups)))
  names(agg_levels) <- vars
  is_bottom <- Reduce(`&`, lapply(vars, function(v) !strata[[v]] %in% agg_levels[[v]]))
  if (!any(is_bottom)) mt_abort("No bottom-level strata found.", class = "hierarchy")
  S <- membership(strata, strata[is_bottom], h)
  if (any(rowSums(S) == 0)) mt_abort("Some aggregate strata have no bottom strata.", class = "hierarchy")
  Sx <- membership(strata, strata[is_bottom], h, exposure = TRUE)
  list(S = S, Sx = Sx, bottom = which(is_bottom), is_bottom = is_bottom, hierarchy = h)
}

# Members of the groups of a hierarchy are bottom levels of the data: a group
# without any is an error (its strata would have no data), members that are
# not levels (typing errors, or levels removed from the data) are reported
check_groups <- function(bottom, h) {
  for (v in names(h$vars)) {
    lev <- unique(as.character(bottom[[v]]))
    for (g in names(h$vars[[v]]$groups)) {
      m <- unique(as.character(h$vars[[v]]$groups[[g]]))
      miss <- setdiff(m, lev)
      if (length(miss) == length(m)) {
        mt_abort(sprintf("Group \"%s\" of `%s` has no level in the data (%s); the levels are %s.", g, v,
                         paste(miss, collapse = ", "), paste(utils::head(lev, 10), collapse = ", ")),
                 class = "hierarchy")
      }
      if (length(miss)) {
        mt_warn(sprintf("Group \"%s\" of `%s`: %s not in the data; the group holds %s.", g, v,
                        paste(miss, collapse = ", "), paste(intersect(m, lev), collapse = ", ")), class = "hierarchy")
      }
    }
  }
  invisible()
}

hierarchy_nodes <- function(bottom, h, by) {
  vars <- names(h$vars)
  others <- setdiff(by, vars)
  level_sets <- lapply(vars, function(v) c(unique(bottom[[v]]), names(h$vars[[v]]$groups), h$vars[[v]]$total))
  names(level_sets) <- vars
  other_sets <- lapply(others, function(o) unique(bottom[[o]]))
  names(other_sets) <- others
  nodes <- do.call(CJ, c(level_sets, other_sets, list(sorted = FALSE)))
  setcolorder(nodes, by)
  nodes[, label := do.call(paste, c(lapply(by, function(b) nodes[[b]]), sep = "."))]
  nodes[]
}

#' Add aggregate strata to data
#'
#' Sums deaths, exposures and populations of the bottom strata into every
#' node of the hierarchy (all combinations of levels, groups and totals of
#' the hierarchy variables). Cells are aggregated where all member strata
#' observe the same age and period interval; other cells are dropped with a
#' warning. Aggregate strata that are already in the data are kept as
#' observed unless `replace = TRUE`; use [mt_validate()] with the hierarchy
#' to check that they agree with the sums. Members of `groups` that are not
#' levels of the data give a warning; a group without any is an error.
#'
#' @param data An [mt_data] object with (at least) the bottom strata.
#' @param hierarchy An [mt_hierarchy()].
#' @param replace Replace aggregate strata that are already in the data.
#' @return An `mt_data` object with all nodes as strata; the hierarchy is
#'   stored in `$hierarchy`.
#' @examples
#' x <- mt_example_regions(years = 2010:2014)
#' xa <- mt_aggregate_strata(x, mt_hierarchy(sex = "both", region = "all", quintile = "all"))
#' nrow(mt_strata(xa))
#' @export
mt_aggregate_strata <- function(data, hierarchy, replace = FALSE) {
  check_class(data, "mt_data", "data")
  h <- as_hierarchy(hierarchy)
  st <- data$strata
  by <- data$by
  hs <- hierarchy_structure(st, h)
  bottom <- st[hs$is_bottom]
  check_groups(bottom, h)
  tl <- find_total_levels(strata_subset(data, hs$is_bottom), intersect(names(h$vars), by))
  if (nrow(tl)) {
    mt_warn(sprintf(paste("Level \"%s\" of `%s` equals the sum of the other levels in %.0f%% of cells; as a bottom",
                          "level it is counted twice in the aggregates. Declare it as the total in mt_hierarchy()",
                          "or remove it with mt_subset()."), tl$level[1], tl$var[1], 100 * tl$share[1]),
            class = "aggregate_total")
  }
  nodes <- hierarchy_nodes(bottom, h, by)
  key_of <- function(d) do.call(paste, c(lapply(by, function(v) d[[v]]), sep = "\r"))
  keep <- if (replace) st[hs$is_bottom] else st
  new_nodes <- nodes[!key_of(nodes) %in% key_of(keep)]
  S_new <- membership(new_nodes, bottom, h)
  Sx_new <- membership(new_nodes, bottom, h, exposure = TRUE)
  final <- rbind(keep[, c(by, "label"), with = FALSE], new_nodes)
  final[, stratum := seq_len(.N)]
  setcolorder(final, c("stratum", by, "label"))
  old_to_new <- stats::setNames(match(key_of(keep), key_of(final)), keep$stratum)
  dropped <- 0L
  for (m in names(data$obs)) {
    o <- data$obs[[m]]
    o <- o[o$stratum %in% keep$stratum]
    ob <- o[o$stratum %in% bottom$stratum]
    add <- vector("list", nrow(new_nodes))
    Sm <- if (m == "deaths") S_new else Sx_new
    for (i in seq_len(nrow(new_nodes))) {
      mem <- bottom$stratum[Sm[i, ] > 0]
      agg <- ob[ob$stratum %in% mem, list(value = sum(value), n = .N),
                by = c("age_lo", "age_hi", "period_lo", "period_hi")]
      complete <- agg$n == length(mem)
      dropped <- dropped + sum(!complete)
      agg <- agg[complete]
      if (nrow(agg)) {
        agg[, `:=`(stratum = nrow(keep) + i, provenance = "aggregated", n = NULL)]
        add[[i]] <- agg
      }
    }
    o <- copy(o)
    o[, stratum := as.integer(old_to_new[as.character(stratum)])]
    o <- rbind(o, rbindlist(add), fill = TRUE)
    setcolorder(o, c("stratum", "age_lo", "age_hi", "period_lo", "period_hi", "value", "provenance"))
    setorderv(o, c("stratum", "period_lo", "age_lo"))
    data$obs[[m]] <- o
  }
  if (dropped > 0) {
    mt_warn(sprintf("%d aggregate cell(s) dropped because not all member strata observe them.", dropped),
            class = "aggregate_incomplete")
  }
  data$strata <- final
  data$hierarchy <- h
  data
}

#' Reconcile forecast draws over a hierarchy
#'
#' Makes forecasts of nested or crossed strata add up, draw by draw. For each
#' age, year and draw, the base forecasts `y` of all strata in the hierarchy
#' (bottom series and aggregates) are projected onto the coherent subspace,
#' `y - W U (U' W U)^-1 U' y`, where `U' y = 0` are the aggregation
#' constraints (each aggregate equals the sum of its bottom series). With
#' `additive = "deaths"` in the hierarchy, the projection works on expected
#' deaths (rates times exposures; exposures of aggregates are the sums of
#' their bottom series), and reconciled rates are reconciled deaths divided
#' by exposures. Methods choose `W`:
#'
#' * `"ols"`: identity (Hyndman, Ahmed, Athanasopoulos and Shang 2011);
#' * `"structural"`: diagonal with the number of bottom series in each node
#'   (structural scaling, Athanasopoulos, Hyndman, Kourentzes and
#'   Petropoulos 2017);
#' * `"wls"`: diagonal with the variance of each stratum's base draws
#'   (weighted least squares, Hyndman, Lee and Wang 2016; see below);
#' * `"mint"`: covariance of the base draws across strata, shrunk towards
#'   its diagonal with the intensity of Schäfer and Strimmer (2005) by
#'   default (minimum-trace reconciliation, Wickramasuriya, Athanasopoulos
#'   and Hyndman 2019; see below);
#' * `"bottom_up"`: aggregates are the sums of the bottom forecasts.
#'
#' **`W` of `"wls"` and `"mint"`.** In the papers, `W` holds the variances
#' (Hyndman, Lee and Wang 2016) or the shrunk covariance matrix
#' (Wickramasuriya et al. 2019) of the in-sample one-step-ahead residuals of
#' the base models. Here `W` is the variance (`"wls"`) or the
#' Schäfer-Strimmer shrunk covariance (`"mint"`) of the base forecast draws
#' across draws, computed separately in each age and forecast year cell: the
#' variance of the model's own `h`-step forecast of that cell (of expected
#' deaths, or of rates when the levels share exposures; without the Poisson
#' variation of the deaths unless predictive draws are reconciled). It
#' changes with the horizon and follows the uncertainty the model reports,
#' not its past errors: a stratum whose model is too narrow is changed too
#' little. For strata fitted by separate models the base draws are
#' independent across strata, so the shrunk covariance is close to its
#' diagonal and `"mint"` gives nearly the result of `"wls"`; it differs
#' where a joint model (for example [li_lee()]) simulates correlated
#' strata, or after the Schaake shuffle. [reconciled()] with `covariance =
#' "errors"` estimates `W` from the errors of past forecasts instead, the
#' input of the papers.
#'
#' **Non-negative forecasts.** Projections can make bottom values negative
#' where the base forecasts disagree strongly. For the draws where this
#' happens, the bottom values `b` are those of the non-negative
#' reconciliation of Wickramasuriya, Turlach and Hyndman (2020): they
#' minimise the same weighted distance `(y - S b)' W^-1 (y - S b)` subject to
#' `b >= 0` (their Eq. 3), and the aggregates are their sums. The quadratic
#' program is solved exactly by their block principal pivoting algorithm
#' (their Algorithm 3.1: all infeasible values exchanged while their number
#' falls, at most 3 times without a fall, then one value at a time),
#' started from the unconstrained solution, with the projection above for
#' each set of bottom values held at zero. Infeasibilities below 1e-10 times
#' the largest absolute base forecast of the draw count as zero (the paper
#' uses 1e-12 in absolute terms). Draws without negative values keep the
#' unconstrained solution, which is then the optimum. A warning reports when
#' more than 5% of the draws have bottom values at zero, because the
#' distributions of the affected strata (often small causes) are then
#' truncated at zero. Earlier versions of mortensor held negative values at
#' zero and projected the rest again, a heuristic without a guaranteed
#' optimum; it was removed.
#'
#' **Fixed nodes.** `fixed` keeps the forecasts of some nodes, as in the
#' reconciliation with immutable forecasts of Zhang, Kang, Panagiotelis and
#' Li (2023). The package gives the fixed nodes a near-zero variance in `W`
#' (1e-10 times the largest variance); with `"mint"` the other nodes get
#' their covariance conditional on the fixed nodes. The projection is then
#' the weighted least-squares solution with the fixed nodes as equality
#' constraints, up to a relative change of the fixed nodes of about 1e-9,
#' and the other strata absorb the whole gap. With negative bottom values,
#' the same weights enter the non-negative quadratic program, so the fixed
#' nodes are kept as equality constraints there as well, as long as
#' non-negative parts can add up to them. When they cannot (a fixed node
#' below zero, or a node fixed together with parts whose forecasts do not
#' add up to it), the fixed nodes move, with a warning.
#'
#' **Proportional.** `"proportional"` keeps the `fixed` nodes and scales the
#' bottom strata below each of them by one factor per draw, so that they add
#' up to the fixed node; nodes between them are recomputed as sums. This is
#' a single-level top-down method (Gross and Sohl 1990) with forecast
#' proportions (Athanasopoulos, Ahmed and Hyndman 2009): the proportions are
#' those of the bottom forecasts below the fixed node in the same draw,
#' while Athanasopoulos et al. (2009) derive them level by level down the
#' hierarchy (the same for two levels). It never sets a part to zero. Fixed
#' nodes must not share parts.
#'
#' **Dependence between strata.** Strata forecast by separate models have
#' independent draws. With `dependence = "none"` (default), the draws are
#' reconciled as they are, so sums of separately fitted parts have narrower
#' intervals than they would if the parts moved together. Joint models (for
#' example [li_lee()], [coherent_fdm()] and [coda()] with `coherent = TRUE`)
#' carry the dependence between strata in the model itself; their draws are
#' never reordered (`"schaake"` and `"empirical"` then give a warning and are
#' ignored). Three methods give separately fitted strata a dependence before
#' the reconciliation, all by reordering their draws. The Schaake shuffle
#' takes the dependence from the observed changes of the training years,
#' either in every age, forecast year and stratum (`dependence =
#' "schaake"`) or for the whole path of each stratum (`dependence =
#' "empirical"`), both below. Ensemble copula coupling (Schefzik,
#' Thorarinsdottir and Gneiting 2013) takes it from the draws of a joint
#' model ([mt_couple()] for draws before `mt_reconcile()`, and [coupled()]
#' or `reconciled(dependence = "ecc")` for models).
#'
#' **The Schaake shuffle** (`dependence = "schaake"`; Clark, Gangopadhyay,
#' Hay, Rajagopalan and Wilby 2004, Section 3b). Status: published, adapted
#' to mortality forecasts as listed below. In the paper, the ensemble
#' members of each station and variable are reordered so that their ranks
#' follow those of observations from as many dates of the historical record,
#' the same dates for all stations and variables (their eqs. 5-12); the
#' later lead times take the following days of the same dates, which keeps
#' the persistence over lead times. Here a margin is an age, a forecast year
#' and a stratum. The template value of a date, a training year `tau`, for
#' the margin of forecast year `h` (counted from the last training year) is
#' the observed change of the log death rate from `tau` to `tau + h`: one
#' observed trajectory per date, over consecutive years. The dates are the
#' `n` training years whose next `H` years are observed (`H` the number of
#' forecast years, so `n` is the number of training years minus `H`), from
#' the training data of the draws only. In every margin the draws receive
#' the sorted values of that margin in the rank order of the template values
#' of their dates. Every margin keeps exactly its values, and the draws take
#' the rank correlations of the observed changes between strata, between
#' ages and between forecast years. Ties among template values are ranked
#' at random, as in the paper (Section 4c); margins with a missing template
#' value keep the order of their draws, as the paper's code skips vectors
#' with missing values (appendix), with a warning. Schefzik et al. (2013,
#' Section 3.4) describe the shuffle as the empirical copula of the
#' historical record applied to samples of the margins; ensemble copula
#' coupling takes the copula from the draws of a model instead.
#'
#' Package-specific (reference unknown): (1) the template holds changes of
#' log rates from the start year of each date, where the paper takes the
#' observed values themselves (from dates within seven days of the
#' forecast date), so that the ranks follow the variation of past changes
#' around the trend of mortality; (2) deaths of zero count as 0.5
#' in the template; (3) the paper has as many dates as ensemble members.
#' With more draws than dates, the draws are taken in a random order in
#' blocks of `n`, and each block is shuffled with all `n` dates; a last
#' block with fewer draws uses a random subset of distinct dates, as do
#' fewer draws than dates (the paper selects its dates at random). The
#' draws of a block take the ranks of the template among themselves, and
#' the same dates recur in every block.
#'
#' The shuffle assumes that the dependence of the training years holds in
#' the forecast years (Clark et al. 2004, Section 5). It sets the
#' dependence between the ages and years of a stratum as well, so measures
#' over ages, such as the life expectancy of a single stratum, change:
#' Lee-Carter draws move all ages of a stratum together, and after the
#' shuffle they move together as much as the observed changes did, which
#' narrows the intervals of life expectancy where these are less correlated
#' between ages. The observed changes carry the Poisson variation of the
#' deaths, which is independent between strata and ages; where deaths are
#' few, it weakens the dependence the template gives. Templates are short:
#' with 19 training years and 5 forecast years there are 14 dates, with 16
#' forecast years 3.
#'
#' **The Schaake shuffle of stratum paths** (`dependence = "empirical"`).
#' Status: the Schaake shuffle of Clark et al. (2004, Section 3b) with the
#' stratum as margin; the margin and the summary that ranks the paths are
#' package-specific (see below). The path of a stratum in one draw, all its
#' ages and forecast years, is one ensemble member and is moved as a whole.
#' The members of each stratum are ranked by a summary of their path, and
#' they are reordered so that their ranks follow those of the same summary
#' in the observed data of `n` distinct training years (the dates), the
#' same dates for all strata, as in the paper. Every stratum keeps exactly
#' its set of paths: its forecast distribution and the dependence between
#' its ages and years stay those of the model, and with `method =
#' "bottom_up"` so do all scores of the bottom strata ([reconciled()] keeps
#' this for predictive draws as well). The pairing
#' of the paths of different strata changes, so that strata whose summaries
#' moved together in the same training years move together in the
#' forecast. In the terms of Schefzik et al. (2013, Section 3.4), the
#' summaries of the strata take the empirical copula of their observed
#' changes. Ben Taieb, Taylor and Hyndman (2021) reorder samples of the
#' bottom series of a hierarchy in the same way, with the empirical copula
#' of in-sample residuals, and sum them bottom-up; Panagiotelis,
#' Gamakumara, Athanasopoulos and Hyndman (2023, Sections 1 and 7.3)
#' discuss that method and note that it needs as many draws as training
#' observations.
#'
#' The summary is the log age-standardised death rate `z`, with one standard
#' population for all strata: the mean exposure of each age over the strata
#' and training years of the draws, so that ages enter with their weight in
#' the population. The template value of a stratum and training year `tau`
#' is the observed annual change `z[tau] - z[tau - 1]`, so `n` is the number
#' of training years minus one. A draw is ranked by the mean annual change of
#' its summary over the forecast years, `(z[H] - z[0]) / H`, where `z[H]` is
#' the summary in the last forecast year and `z[0]` that of the observed
#' rates of the last training year. Both ranks so refer to annual changes of
#' the summary; `z[0]` is the same for all draws of a stratum, so the ranks
#' of the draws are those of `z[H]`. Draws and dates are matched in blocks
#' as for `"schaake"` (item 3 above); ties among template values are ranked
#' at random (Clark et al. 2004, Section 4c); deaths of zero count as 0.5.
#' With fewer than three training years (one observed change) the draws keep
#' their order, with a warning. A stratum with a missing template value
#' keeps the order of its draws, with a warning.
#'
#' Package-specific (reference unknown): the stratum as margin, the summary
#' and its standard population, the ranking of the draws by their mean
#' annual change, and the blocks with more draws than dates.
#'
#' The three methods differ in what they reorder. `"schaake"` reorders every
#' age, forecast year and stratum separately: it replaces the dependence
#' between the ages and years of a stratum with that of the observed changes
#' as well, which changes measures over ages and years of a single stratum
#' (life expectancy, cumulative changes). `"empirical"` leaves each stratum
#' to its model and gives only the strata a dependence, from one summary per
#' stratum and year, which carries less Poisson variation than the rates of
#' single ages. `"ecc"` takes the dependence of all cells from the draws of a
#' joint model, in place of observed years. Each assumes that the
#' dependence it takes holds in the forecast years. The template of
#' `"empirical"` holds the changes of single years, and the draws are ranked
#' by their changes over all forecast years: the strata of one draw take the
#' ranks that their annual changes had in one training year, so they move
#' together as far as their annual changes did. The observed changes carry
#' the Poisson variation of the deaths, which is independent between strata
#' and weakens the dependence where deaths are few.
#'
#' In the package's backtests of bottom-up Lee-Carter forecasts, on data
#' simulated with correlated period shocks, on the deprivation quintiles of
#' England ([ons_data]) and on the German Laender, `"empirical"` widened the
#' intervals of the totals towards their nominal coverage and lowered their
#' CRPS, with the scores of every stratum unchanged; on
#' [mt_example_regions()], whose separate forecasts already covered the
#' totals, it made the intervals of the totals too wide. The vignettes
#' "Coherent forecasts: totals and their parts" and "Dependence between
#' groups" show such backtests.
#'
#' The one-variable interface with `var` and `total` is still accepted. It
#' reconciles deaths when the exposures of the parts add up to that of the
#' total (sexes, regions) and rates when each part has the exposure of the
#' total (causes of death) or there are no exposures; give a hierarchy with
#' `additive` to choose.
#'
#' @param draws An [mt_draws] object with the strata of the hierarchy.
#' @param hierarchy An [mt_hierarchy()] (or, for one variable, `var` and
#'   `total`).
#' @param method `"wls"` (default), `"mint"`, `"structural"`, `"ols"`,
#'   `"bottom_up"` or `"proportional"` (with `fixed`). The vignette
#'   "Coherent forecasts: totals and their parts" compares them.
#'   Score-optimal reconciliation (`"energy_score"`) is learnt from past
#'   forecasts and is available in [reconciled()] only.
#' @param shrink Shrinkage intensity for `"mint"` (number between 0 and 1),
#'   or `"auto"`.
#' @param exposure Exposures `[age, period, stratum]` or `[age, stratum]` for
#'   the forecast years (default: last training year). Only bottom strata are
#'   used, and they are enough; aggregates are summed.
#' @param dependence `"none"` (default: keep the draws as they are),
#'   `"schaake"` (the Schaake shuffle of every age, year and stratum of the
#'   draws of separately fitted strata before the reconciliation) or
#'   `"empirical"` (the Schaake shuffle of the whole paths of each stratum);
#'   see Details. Draws of joint models are not reordered: both then give a
#'   warning and are ignored. Both need the training data of the draws (as
#'   [mt_forecast()] returns them); `"schaake"` needs at least two more
#'   training years than forecast years, `"empirical"` at least three
#'   training years.
#' @param seed Random seed for the Schaake shuffle (order of the draws,
#'   choice of dates, ties).
#' @param fixed Strata whose forecasts are kept unchanged (labels, or
#'   `"total"` for the grand total of the hierarchy), for example a trusted
#'   all-cause forecast that the causes must add up to. The projection gives
#'   these nodes near-zero variance, so the other strata absorb the whole gap;
#'   `$reconciliation$max_relative_change_fixed` reports how closely they
#'   were kept; they move, with a warning, when their own forecasts do not
#'   add up (a node fixed together with all its parts) or when every bottom
#'   stratum below them would have to be negative. Needs a projection
#'   method or `"proportional"`.
#' @param var,total,parts One-variable form: stratum variable and its total
#'   level (`parts` is ignored: all other levels are parts).
#' @return The reconciled `mt_draws`, with diagnostics in `$reconciliation`
#'   (method, largest relative gap before reconciliation, share of bottom
#'   values held at zero, mean shrinkage intensity). After the Schaake
#'   shuffle, `$schaake` holds the number of dates (`dates`), their start
#'   years (`years`), the start year of the date of each draw (`draw_year`),
#'   the block of each draw (`block`), the number of blocks (`blocks`), the
#'   number of shuffled margins (`margins`) and the share of margins with
#'   tied template values (`tied`); `$component` is then dropped, since the
#'   values of one draw come from different base draws. After the shuffle of
#'   stratum paths, `$empirical` holds the same items, with the years of the
#'   observed changes as `years` and the number of reordered strata as
#'   `margins`, and also the standard population (`standard`) and the
#'   template (`template`, the observed annual changes of the log
#'   age-standardised rate `[stratum, year]`); `$component` then gives the
#'   component of each stratum and draw (ensembles).
#' @references
#' Hyndman RJ, Ahmed RA, Athanasopoulos G, Shang HL (2011). Optimal
#' combination forecasts for hierarchical time series. *Computational
#' Statistics & Data Analysis* 55(9): 2579-2589.
#'
#' Hyndman RJ, Lee AJ, Wang E (2016). Fast computation of reconciled
#' forecasts for hierarchical and grouped time series. *Computational
#' Statistics & Data Analysis* 97.
#'
#' Athanasopoulos G, Hyndman RJ, Kourentzes N, Petropoulos F (2017).
#' Forecasting with temporal hierarchies. *European Journal of Operational
#' Research* 262(1): 60-74.
#'
#' Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast
#' reconciliation for hierarchical and grouped time series through trace
#' minimization. *Journal of the American Statistical Association*
#' 114(526): 804-819. \doi{10.1080/01621459.2018.1448825}
#'
#' Wickramasuriya SL, Turlach BA, Hyndman RJ (2020). Optimal non-negative
#' forecast reconciliation. *Statistics and Computing* 30: 1167-1182.
#' \doi{10.1007/s11222-020-09930-0}
#'
#' Zhang B, Kang Y, Panagiotelis A, Li F (2023). Optimal reconciliation
#' with immutable forecasts. *European Journal of Operational Research*
#' 308(2): 650-660. \doi{10.1016/j.ejor.2022.11.035}
#'
#' Gross CW, Sohl JE (1990). Disaggregation methods to expedite product line
#' forecasting. *Journal of Forecasting* 9(3): 233-254.
#'
#' Athanasopoulos G, Ahmed RA, Hyndman RJ (2009). Hierarchical forecasts for
#' Australian domestic tourism. *International Journal of Forecasting*
#' 25: 146-166.
#'
#' Clark M, Gangopadhyay S, Hay L, Rajagopalan B, Wilby R (2004). The Schaake
#' shuffle: a method for reconstructing space-time variability in forecasted
#' precipitation and temperature fields. *Journal of Hydrometeorology* 5(1):
#' 243-262.
#'
#' Schefzik R, Thorarinsdottir TL, Gneiting T (2013). Uncertainty
#' quantification in complex simulation models using ensemble copula
#' coupling. *Statistical Science* 28(4): 616-640. \doi{10.1214/13-STS443}
#'
#' Ben Taieb S, Taylor JW, Hyndman RJ (2021). Hierarchical probabilistic
#' forecasting of electricity demand with smart meter data. *Journal of the
#' American Statistical Association* 116(533): 27-43.
#' \doi{10.1080/01621459.2020.1736081}
#'
#' Panagiotelis A, Gamakumara P, Athanasopoulos G, Hyndman RJ (2023).
#' Probabilistic forecast reconciliation: properties, evaluation and score
#' optimisation. *European Journal of Operational Research* 306(2):
#' 693-706. \doi{10.1016/j.ejor.2022.07.040}
#'
#' Schäfer J, Strimmer K (2005). A shrinkage approach to large-scale
#' covariance matrix estimation and implications for functional genomics.
#' *Statistical Applications in Genetics and Molecular Biology* 4(1):
#' Article 32. \doi{10.2202/1544-6115.1175}
#' @examples
#' x <- mt_example_regions(years = 2005:2019, regions = 3)
#' h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
#' xa <- mt_aggregate_strata(x, h)
#' fc <- mt_forecast(mt_fit(lee_carter(), xa), h = 3, n_draws = 50, seed = 1)
#' rc <- mt_reconcile(fc, h, method = "wls")
#' rc$reconciliation
#' # the strata take the dependence of the observed changes of 12 training years
#' rs <- mt_reconcile(fc, h, method = "wls", dependence = "schaake", seed = 1)
#' rs$schaake$years
#' # whole paths of each stratum, ranked by the annual changes of 14 training years
#' re <- mt_reconcile(fc, h, method = "bottom_up", dependence = "empirical", seed = 1)
#' re$empirical$template[1:3, 1:5]
#' @export
mt_reconcile <- function(draws, hierarchy = NULL,
                         method = c("wls", "mint", "structural", "ols", "bottom_up", "proportional"),
                         shrink = "auto", exposure = NULL, dependence = c("none", "schaake", "empirical"),
                         seed = NULL, fixed = NULL, var = NULL, total = NULL, parts = NULL) {
  check_class(draws, "mt_draws", "draws")
  if (identical(method, "energy_score")) {
    mt_abort(paste("`method = \"energy_score\"` learns the reconciliation from past forecasts and observations,",
                   "which forecast draws do not hold: use reconciled(spec, hierarchy, method = \"energy_score\"),",
                   "which makes them inside the training data."), class = "argument")
  }
  method <- match.arg(method)
  dependence <- match.arg(dependence)
  dependence <- joint_dependence(dependence, draws$spec)
  if (is.character(hierarchy) && length(hierarchy) == 1 && is.null(var)) {
    var <- hierarchy
    hierarchy <- NULL
  }
  if (is.null(hierarchy)) {
    if (is.null(var) || is.null(total)) mt_abort("Give a hierarchy or `var` and `total`.", class = "argument")
    if (!var %in% names(draws$strata)) mt_abort("`var` is not a stratum variable.", class = "argument")
    if (!total %in% draws$strata[[var]]) mt_abort("No total with parts found for the hierarchy.", class = "argument")
    hierarchy <- do.call(mt_hierarchy, c(stats::setNames(list(total), var),
                                         list(additive = one_var_additive(draws, var, total, exposure))))
  }
  h <- as_hierarchy(hierarchy)
  hs <- tryCatch(hierarchy_structure(draws$strata, h), mt_error_hierarchy = function(e) {
    mt_abort(conditionMessage(e), class = "argument")
  })
  fixed_idx <- fixed_nodes(fixed, draws$strata, h, hs, method)
  arr <- draws$rates
  if (dependence != "none") {
    dep <- with_seed(seed, dependence_reorder(arr, dependence, train_of(draws)))
    warn_schaake(dep$info)
    arr <- dep$arr
    if (dependence == "schaake") {
      # the values of one draw come from different base draws
      draws$component <- NULL
      draws$nb_size <- NULL
      draws$schaake <- dep$info
    } else {
      # each path of a stratum comes from one base draw
      draws$component <- permute_component(draws$component, dep$perm)
      draws$nb_size <- permute_rows(draws$nb_size, dep$perm)
      draws$empirical <- dep$info
    }
  }
  rec <- reconcile_array(arr, hs, h$additive, method, shrink, exposure %||% draws$exposure_last, dim(arr),
                         fixed = fixed_idx)
  draws$rates <- rec$rates
  draws$reconciliation <- rec$info
  draws
}

# Whether the draws of a specification carry the dependence between strata:
# joint models, wrappers around them, and ensembles whose components all do
# (each draw of an ensemble comes from one component for all strata)
spec_is_joint <- function(spec) {
  if (!inherits(spec, "mt_spec")) return(FALSE)
  if (isTRUE(spec$joint)) return(TRUE)
  # coupled draws carry the dependence of their template
  if (inherits(spec, "mt_spec_coupled")) return(TRUE)
  comps <- spec$params$components
  if (inherits(spec, "mt_spec_ensemble") && length(comps)) return(all(vapply(comps, spec_is_joint, logical(1))))
  if (inherits(spec$params$base, "mt_spec")) return(spec_is_joint(spec$params$base))
  FALSE
}

# The dependence used for draws of `spec` (or of any of a list of specs, such
# as the specs fitted below wrappers): the paths of joint models are never
# reordered, since their draws carry the dependence of the model
joint_dependence <- function(dependence, spec) {
  specs <- if (inherits(spec, "mt_spec")) list(spec) else spec
  if (dependence %in% c("empirical", "schaake") && any(vapply(specs, spec_is_joint, logical(1)))) {
    mt_warn(sprintf(paste("The draws of joint models carry the dependence between strata and are not reordered:",
                          "`dependence = \"%s\"` is ignored."), dependence), class = "reconcile")
    return("none")
  }
  dependence
}

# Additivity of the one-variable form of mt_reconcile(), from the exposures:
# "deaths" when the exposures of the parts add up to that of their total
# (sexes, regions), "rates" when each part has the exposure of the total
# (causes of death) or there are no exposures
one_var_additive <- function(draws, var, total, exposure) {
  ex <- exposure %||% draws$exposure_last
  if (is.null(ex)) return("rates")
  E <- tryCatch(suppressWarnings(expand_exposure(ex, dim(draws$rates), dimnames(draws$rates))),
                error = function(e) NULL)
  # exposures of the parts only: they add up
  if (is.null(E)) return("deaths")
  st <- draws$strata
  others <- setdiff(names(st), c("stratum", "label", var))
  key <- if (length(others)) do.call(paste, c(unname(as.list(st)[others]), sep = "\r")) else rep("", nrow(st))
  ratio <- vapply(which(st[[var]] == total), function(i) {
    j <- which(key == key[i] & st[[var]] != total)
    if (length(j)) sum(E[, , j]) / sum(E[, , i]) else NA_real_
  }, 0)
  if (all(is.na(ratio)) || stats::median(ratio, na.rm = TRUE) < 1.5) "deaths" else "rates"
}

# Indices of the strata whose forecasts are kept in reconciliation: labels,
# or "total" for the grand total of the hierarchy
fixed_nodes <- function(fixed, strata, h, hs, method) {
  if (method == "proportional" && !length(fixed)) {
    mt_abort("`method = \"proportional\"` scales the parts of fixed nodes: give `fixed`.", class = "argument")
  }
  if (is.null(fixed) || !length(fixed)) return(integer())
  if (method == "bottom_up") {
    mt_abort("`fixed` needs a projection method (\"wls\", \"mint\", \"structural\", \"ols\") or \"proportional\".",
             class = "argument")
  }
  idx <- integer()
  if ("total" %in% fixed) {
    tot <- Reduce(`&`, lapply(names(h$vars), function(v) strata[[v]] == h$vars[[v]]$total))
    idx <- which(tot)
    if (!length(idx)) mt_abort("The strata hold no grand total of the hierarchy.", class = "argument")
  }
  lab <- setdiff(fixed, "total")
  bad <- setdiff(lab, strata$label)
  if (length(bad)) mt_abort(sprintf("Unknown stratum label(s) in `fixed`: %s.", paste(bad, collapse = ", ")),
                            class = "argument")
  idx <- union(idx, match(lab, strata$label))
  if (all(hs$is_bottom[idx]) && length(idx) == sum(hs$is_bottom)) {
    mt_abort("`fixed` cannot hold every bottom stratum.", class = "argument")
  }
  if (method == "proportional") {
    agg <- idx[!hs$is_bottom[idx]]
    if (!length(agg)) mt_abort("`method = \"proportional\"` needs an aggregate among `fixed`.", class = "argument")
    if (any(colSums(hs$S[agg, , drop = FALSE] > 0) > 1)) {
      mt_abort(paste("`method = \"proportional\"` needs fixed nodes without common parts (for example all",
                     "causes of every region, but not also the national total)."), class = "argument")
    }
  }
  idx
}

# Core: arr [age, period, stratum, draw] of rates. `W_cell`: optional
# function(age, period) giving the weights of "wls" (a vector) or "mint"
# (list(W, lambda)) from base forecast errors, or NULL for the draws
reconcile_array <- function(arr, hs, additive, method, shrink, exposure, d, fixed = integer(), W_cell = NULL) {
  S <- hs$S
  b <- hs$bottom
  is_agg <- !hs$is_bottom
  agg_rows <- which(is_agg)
  nb <- length(b)
  E <- expand_exposure(exposure, d, dimnames(arr), hs)
  if (additive == "deaths") {
    Sx <- hs$Sx %||% S
    for (t in seq_len(d[2])) E[, t, ] <- matrix(E[, t, b], d[1]) %*% t(Sx)
    Y <- arr * as.vector(E)
  } else {
    Y <- arr
  }
  Ut <- NULL
  if (length(agg_rows)) {
    Ut <- matrix(0, length(agg_rows), d[3])
    Ut[cbind(seq_along(agg_rows), agg_rows)] <- 1
    Ut[, b] <- -S[agg_rows, , drop = FALSE]
  }
  gap_before <- 0
  fixed_change <- 0
  zero_share <- 0
  draw_zero <- logical(d[4])
  shrink_used <- numeric(0)
  nnls_failed <- 0L
  out <- Y
  for (a in seq_len(d[1])) {
    for (t in seq_len(d[2])) {
      Yc <- t(matrix(Y[a, t, , ], nrow = d[3]))
      if (length(agg_rows)) {
        tot_b <- Yc[, b, drop = FALSE] %*% t(S[agg_rows, , drop = FALSE])
        gap_before <- max(gap_before, max(abs(Yc[, agg_rows] - tot_b) / pmax(abs(tot_b), 1e-12)))
      }
      Yb <- if (method == "bottom_up" || !length(agg_rows)) {
        Yc[, b, drop = FALSE]
      } else if (method == "proportional") {
        # parts of each fixed node scaled in proportion to their forecasts
        Ypb <- Yc[, b, drop = FALSE]
        for (fi in fixed[!hs$is_bottom[fixed]]) {
          mem <- which(S[fi, ] > 0)
          tot <- Ypb[, mem, drop = FALSE] %*% S[fi, mem]
          Ypb[, mem] <- Ypb[, mem, drop = FALSE] * as.vector(Yc[, fi] / pmax(tot, 1e-300))
        }
        Ypb
      } else {
        # diagonal weights as vectors, MinT as a matrix; fixed nodes get a
        # near-zero variance
        We <- if (!is.null(W_cell)) W_cell(a, t)
        W <- if (!is.null(We)) {
          if (method == "mint") {
            shrink_used <- c(shrink_used, We$lambda)
            We$W
          } else {
            We
          }
        } else {
          switch(method,
                 ols = rep(1, d[3]),
                 structural = rowSums(S),
                 wls = pmax(apply(Yc, 2, stats::var), 1e-20),
                 mint = {
                   sh <- shrink_cov(Yc, shrink)
                   shrink_used <- c(shrink_used, sh$lambda)
                   sh$W
                 })
        }
        W <- keep_weights(W, fixed)
        Yp <- project_rows(Yc, W, Ut)
        # non-negative reconciliation (Wickramasuriya, Turlach and Hyndman
        # 2020) for the draws with negative bottom values
        neg <- which(rowSums(Yp[, b, drop = FALSE] < 0) > 0)
        if (length(neg)) {
          nn <- nnls_reconcile(Yc[neg, , drop = FALSE], W, Ut, b, Yp[neg, , drop = FALSE])
          Yp[neg, ] <- nn$X
          zero_share <- zero_share + sum(nn$held)
          draw_zero[neg] <- draw_zero[neg] | rowSums(nn$held) > 0
          nnls_failed <- nnls_failed + sum(!nn$converged)
        }
        Yp[, b, drop = FALSE]
      }
      zero_share <- zero_share + sum(Yb < 0)
      Yb <- pmax(Yb, 0)
      out[a, t, , ] <- t(Yb %*% t(S))
      if (length(fixed)) {
        y0 <- Y[a, t, fixed, ]
        fixed_change <- max(fixed_change, max(abs(out[a, t, fixed, ] - y0) / pmax(abs(y0), 1e-12)))
      }
    }
  }
  zero_share <- zero_share / (d[1] * d[2] * d[4] * nb)
  if (nnls_failed > 0) {
    mt_warn(sprintf(paste("The non-negative reconciliation did not reach its optimum in %d draw(s) and cell(s);",
                          "their negative bottom values were set to zero."), nnls_failed), class = "reconcile")
  }
  if (length(fixed) && fixed_change > 1e-6) {
    mt_warn(sprintf(paste("Fixed nodes moved by up to %.1f%%: they cannot all be kept when their forecasts do not",
                          "add up (a node fixed together with all its parts) or when the parts below them would",
                          "have to be negative."), 100 * fixed_change), class = "reconcile")
  }
  if (mean(draw_zero) > 0.05) {
    mt_warn(sprintf(paste("In %.0f%% of the draws some bottom forecasts were set to zero to keep the aggregates%s;",
                          "their distributions are truncated. Check the base forecasts (mt_coherence()) or use",
                          "method = \"proportional\" with fixed nodes."), 100 * mean(draw_zero),
                    if (length(fixed)) " and fixed nodes" else ""), class = "reconcile")
  }
  rates <- if (additive == "deaths") out / as.vector(E) else out
  rates <- pmax(rates, 1e-12)
  list(rates = rates, exposure = E,
       info = list(method = method, max_relative_gap_before = gap_before, share_set_to_zero = zero_share,
                   mean_shrinkage = if (length(shrink_used)) mean(shrink_used) else NA_real_,
                   share_draws_with_zeros = mean(draw_zero),
                   max_relative_change_fixed = if (length(fixed)) fixed_change else NA_real_))
}

# Weights with (almost) zero variance for the nodes in `keep`, so that a
# projection leaves them unchanged: their variance becomes 1e-10 times the
# largest variance. In a matrix, the other nodes get their covariance
# conditional on the kept nodes and no covariance with them: the inverse of
# this block is the block of the inverse of W, so that the projection is the
# weighted least-squares solution with the kept nodes as equality
# constraints (up to the near-zero variance). Diagonal weights need no
# conditioning
keep_weights <- function(W, keep = integer()) {
  if (!length(keep)) return(W)
  tiny <- 1e-10 * max(if (is.matrix(W)) diag(W) else W)
  if (!is.matrix(W)) {
    W[keep] <- tiny
    return(W)
  }
  o <- setdiff(seq_len(ncol(W)), keep)
  Wk <- W[keep, keep, drop = FALSE]
  Wok <- W[o, keep, drop = FALSE]
  Wc <- W[o, o, drop = FALSE] - Wok %*% tryCatch(solve(Wk, t(Wok)), error = function(e) MASS_ginv(Wk) %*% t(Wok))
  out <- matrix(0, ncol(W), ncol(W))
  out[o, o] <- (Wc + t(Wc)) / 2
  out[cbind(keep, keep)] <- tiny
  out
}

# Projection of the rows of Y (draws x strata) onto the subspace where the
# linear constraints Ct x = 0 hold, with weights W (a vector of diagonal
# weights or a matrix): x = y - W C (C' W C)^-1 C' y. With the aggregation
# constraints Ut this is the reconciliation of Wickramasuriya, Athanasopoulos
# and Hyndman (2019) in the form of their Lemma 1 (Eq. 2 of Wickramasuriya,
# Turlach and Hyndman 2020). With `multipliers = TRUE` the Lagrange
# multipliers (C' W C)^-1 C' y [constraint, row] are returned as well
project_rows <- function(Y, W, Ct, multipliers = FALSE) {
  WC <- if (is.matrix(W)) W %*% t(Ct) else W * t(Ct)
  M <- Ct %*% WC
  rhs <- Ct %*% t(Y)
  lam <- tryCatch(solve(M, rhs), error = function(e) MASS_ginv(M) %*% rhs)
  X <- Y - t(WC %*% lam)
  if (multipliers) list(X = X, lambda = lam) else X
}

# Non-negative reconciliation of Wickramasuriya, Turlach and Hyndman (2020):
# the bottom values b of each row y of Y minimise (y - S b)' W^-1 (y - S b)
# subject to b >= 0 (their Eq. 3). Solved by their block principal pivoting
# algorithm (Algorithm 3.1, Judice and Pires): the bottom values are split
# into a free set F and a set G held at zero; for a split, the solution with
# b_G = 0 (Eq. 9) and the multipliers g_G of the held values (Eq. 10) come
# from the projection with the constraints U' x = 0 and x_G = 0 (the
# multiplier of x_j = 0 is -g_j). Infeasible values (b_j < 0 in F, g_j < 0
# in G) are exchanged all at once while their number falls, up to p = 3
# times without a fall, and then one at a time (the last index, Murty's
# rule), which ends after finitely many steps. Infeasibilities smaller than
# 1e-10 times the largest absolute base forecast of the row (for g, in the
# units of y: g_j times the weight of j) count as zero. The start is the
# unconstrained solution `X0` (F = all bottom values), as in the paper for
# diagonal weights. Rows with the same split are solved together. Returns
# the solutions of all nodes `X`, the held bottom values `held` [row,
# bottom] and whether each row reached the optimum (`converged`)
nnls_reconcile <- function(Y, W, Ut, b, X0, p_bar = 3L, max_iter = NULL) {
  nr <- nrow(Y)
  nb <- length(b)
  max_iter <- max_iter %||% (10L * nb + 100L)
  wb <- if (is.matrix(W)) diag(W)[b] else W[b]
  tol <- 1e-10 * pmax(apply(abs(Y), 1, max), 1e-300)
  X <- X0
  xb <- X0[, b, drop = FALSE]
  g <- matrix(0, nr, nb)
  held <- matrix(FALSE, nr, nb)
  ninf <- rep(nb + 1L, nr)
  p <- rep(p_bar, nr)
  converged <- rep(FALSE, nr)
  for (it in seq_len(max_iter + 1L)) {
    inf <- (!held & xb < -tol) | (held & g * rep(wb, each = nr) < -tol)
    n_inf <- rowSums(inf)
    converged <- n_inf == 0L
    active <- which(!converged)
    if (!length(active) || it > max_iter) break
    for (i in active) {
      if (n_inf[i] < ninf[i]) {
        ninf[i] <- n_inf[i]
        p[i] <- p_bar
      } else if (p[i] >= 1L) {
        p[i] <- p[i] - 1L
      } else {
        # backup rule: exchange only the last infeasible variable
        r <- max(which(inf[i, ]))
        inf[i, ] <- FALSE
        inf[i, r] <- TRUE
      }
    }
    # exchange: infeasible free values are held at zero, infeasible held ones freed
    held[active, ] <- xor(held[active, , drop = FALSE], inf[active, , drop = FALSE])
    key <- do.call(paste0, as.data.frame(held[active, , drop = FALSE] + 0L))
    for (k in unique(key)) {
      rr <- active[key == k]
      hg <- which(held[rr[1], ])
      Ct <- Ut
      if (length(hg)) {
        Eg <- matrix(0, length(hg), ncol(Ut))
        Eg[cbind(seq_along(hg), b[hg])] <- 1
        Ct <- rbind(Ut, Eg)
      }
      pr <- project_rows(Y[rr, , drop = FALSE], W, Ct, multipliers = TRUE)
      X[rr, ] <- pr$X
      X[rr, b[hg]] <- 0
      xb[rr, ] <- X[rr, b, drop = FALSE]
      g[rr, ] <- 0
      if (length(hg)) g[rr, hg] <- -t(pr$lambda[nrow(Ut) + seq_along(hg), , drop = FALSE])
    }
  }
  if (any(!converged)) {
    # safeguard: negative bottom values of rows without the optimum set to
    # zero (the caller computes the aggregates as their sums)
    bad <- which(!converged)
    xb_bad <- pmax(xb[bad, , drop = FALSE], 0)
    X[bad, b] <- xb_bad
    held[bad, ] <- xb_bad == 0
  }
  list(X = X, held = held, converged = converged, iterations = it - 1L)
}

# The Schaake shuffle (Clark et al. 2004, Section 3b; the empirical copula
# reading of Schefzik et al. 2013, Section 3.4). `arr` [age, period,
# stratum, draw]: base draws, whose period h is h years after the last
# training year; `y` [age, training period, stratum]: observed log rates of
# the training years (schaake_log_rates()); `periods`: the training years.
# The template of date i (start year tau_i, the i-th of the n = T - H
# training years whose next H years are observed) holds, for each margin
# (age, horizon h, stratum), the change of the log rate from tau_i to
# tau_i + h: consecutive years for the horizons, as the paper's consecutive
# days for the lead times. For each margin, the draws receive the sorted
# base values in the rank order of the template values of their dates
# (schaake_ranks()). Returns the reordered draws, the base draw each value
# comes from (`draw`, as dim(arr)) and `info`: the n start years, the start
# year of each draw's date, the block of each draw
schaake_reorder <- function(arr, y, periods) {
  d <- dim(arr)
  M <- d[4]
  H <- d[2]
  n <- dim(y)[2] - H
  if (n < 2) {
    mt_abort(sprintf(paste("The Schaake shuffle needs at least %d training years for a forecast of %d years: the",
                           "template takes the changes over 1 to %d years from at least two distinct years."),
                     H + 2, H, H), class = "reconcile")
  }
  y <- match_training(arr, y)
  # template: one row per margin (age, horizon, stratum), one column per date
  L <- prod(d[1:3])
  tpl <- array(NA_real_, c(d[1:3], n))
  for (h in seq_len(H)) {
    tpl[, h, , ] <- aperm(y[, h + seq_len(n), , drop = FALSE] - y[, seq_len(n), , drop = FALSE], c(1, 3, 2))
  }
  Tm <- matrix(tpl, L, n)
  B <- matrix(arr, L, M)
  sr <- schaake_ranks(B, Tm)
  out <- array(B[cbind(rep(seq_len(L), M), as.vector(sr$src))], d, dimnames = dimnames(arr))
  miss <- array(sr$na_margin & !sr$base_na, d[1:3])
  info <- list(method = "Schaake shuffle", dates = n, years = periods[seq_len(n)], draw_year = periods[sr$date],
               block = sr$block, blocks = sr$blocks, margins = sum(sr$ok), tied = sr$tied,
               uncoupled = dimnames(arr)[[3]][apply(miss, 3, any)])
  list(arr = out, draw = array(sr$src, d), info = info)
}

# The reordering of the Schaake shuffle for K margins (Clark et al. 2004,
# eqs. 5-12). `V` [margin, draw]: the values by which the members of each
# margin are ranked; `Tm` [margin, date]: the template values of the n
# dates. The draws are taken in a random order in blocks of n, each block
# shuffled with the n distinct dates; with fewer than n draws in a block
# (also M < n), its dates are a random subset of distinct dates
# (package-specific for M > n: the paper has as many dates as members). In
# each margin and block, the member with the k-th smallest value goes to the
# draw whose date has the k-th smallest template value among the dates of
# the block (ties ranked at random, Section 4c). Margins with a missing
# template or base value keep the order of the base draws (Clark et al.
# 2004, appendix). Returns the base draw of each reordered member (`src`
# [margin, draw]), the date and the block of each draw, the margins
# reordered (`ok`) and the share of them with tied template values
schaake_ranks <- function(V, Tm) {
  K <- nrow(V)
  M <- ncol(V)
  n <- ncol(Tm)
  Tm[!is.finite(Tm)] <- NA
  na_margin <- rowSums(is.na(Tm)) > 0
  base_na <- rowSums(is.na(V)) > 0
  ok <- !na_margin & !base_na
  # the j-th draw of a random order is in block blk[j] and has date dates[j]
  perm <- sample.int(M)
  blk <- as.integer(ceiling(seq_len(M) / n))
  dates <- integer(M)
  for (b in unique(blk)) {
    j <- which(blk == b)
    dates[j] <- if (length(j) == n) seq_len(n) else sample.int(n, length(j))
  }
  src <- matrix(rep(seq_len(M), each = K), K)
  tied <- 0
  if (any(ok)) {
    Lo <- sum(ok)
    Vp <- V[ok, perm, drop = FALSE]
    Tp <- Tm[ok, dates, drop = FALSE]
    # one group per margin and block: in each group, the member of rank k
    # among the base values goes to the draw whose template value has rank k
    # among the template values (ties ranked at random)
    g <- (rep(blk, each = Lo) - 1) * Lo + rep(seq_len(Lo), M)
    ob <- order(g, Vp)
    ot <- order(g, Tp, stats::runif(Lo * M))
    Sn <- matrix(0L, Lo, M)
    Sn[ot] <- rep(perm, each = Lo)[ob]
    src[ok, perm] <- Sn
    # share of margins with tied template values among the n dates
    To <- Tm[ok, , drop = FALSE]
    st <- matrix(To[order(rep(seq_len(Lo), n), To)], Lo, n, byrow = TRUE)
    tied <- mean(rowSums(st[, -1, drop = FALSE] == st[, -n, drop = FALSE]) > 0)
  }
  date <- block <- integer(M)
  date[perm] <- dates
  block[perm] <- blk
  list(src = src, date = date, block = block, blocks = max(blk), ok = ok, na_margin = na_margin,
       base_na = base_na, tied = tied)
}

# The Schaake shuffle of stratum paths (dependence = "empirical"): the
# shuffle of schaake_ranks() with the stratum as margin, each draw's path of
# a stratum (all ages and forecast years) one member, moved as a whole.
# Members are ranked by the mean annual change of the log age-standardised
# rate z over the forecast years, (z_H - z_0) / H, with z_0 that of the
# observed rates of the last training year (the same for all draws, so the
# ranks are those of z_H); the template value of a stratum and training year
# tau is the observed change z_tau - z_(tau - 1), so the n = T - 1 dates are
# the training years after the first. The standard population of z is the
# mean exposure by age over the strata and training years (`E` [age,
# training period, stratum]). `y`: observed log rates of the training years
# (schaake_log_rates()). With fewer than 3 training years the draws keep
# their order, with a warning. Returns the reordered draws, the base draw of
# each stratum path (`perm` [stratum, draw]) and `info`
empirical_reorder <- function(arr, y, E, periods) {
  d <- dim(arr)
  H <- d[2]
  y <- match_training(arr, y)
  E <- match_training(arr, E)
  nt <- dim(y)[2]
  n <- nt - 1L
  if (n < 2) {
    mt_warn(sprintf(paste("`dependence = \"empirical\"` needs at least 3 training years (2 observed annual changes);",
                          "the draws have %d and keep their order."), nt), class = "reconcile", call = NULL)
    return(list(arr = arr, perm = NULL, info = NULL))
  }
  w <- apply(E, 1, mean, na.rm = TRUE)
  w <- w / sum(w)
  # log age-standardised rates [period, stratum] of the training years and
  # [stratum, draw] of the draws in the last forecast year
  z <- matrix(log(colSums(matrix(exp(y), d[1]) * w)), nt)
  Tm <- t(z[-1, , drop = FALSE] - z[-nt, , drop = FALSE])
  zh <- matrix(log(colSums(matrix(arr[, H, , , drop = FALSE], d[1]) * w)), d[3])
  V <- (zh - z[nt, ]) / H
  V[!is.finite(V)] <- NA
  sr <- schaake_ranks(V, Tm)
  out <- arr
  for (s in which(sr$ok)) out[, , s, ] <- arr[, , s, sr$src[s, ], drop = FALSE]
  years <- periods[-1]
  labels <- dimnames(arr)[[3]] %||% as.character(seq_len(d[3]))
  dimnames(Tm) <- list(stratum = labels, year = as.character(years))
  info <- list(method = "Schaake shuffle of stratum paths", dates = n, years = years, draw_year = years[sr$date],
               block = sr$block, blocks = sr$blocks, margins = sum(sr$ok), tied = sr$tied,
               uncoupled = labels[sr$na_margin & !sr$base_na],
               standard = stats::setNames(w, dimnames(arr)[[1]]), template = Tm)
  list(arr = out, perm = sr$src, info = info)
}

# Training data `A` [age, training period, stratum] matched to the ages and
# strata of the draws `arr` by their labels
match_training <- function(arr, A) {
  dn <- dimnames(arr)
  dy <- dimnames(A)
  if (!is.null(dn[[3]]) && !is.null(dy[[3]])) {
    ia <- if (!is.null(dn[[1]]) && !is.null(dy[[1]])) match(dn[[1]], dy[[1]]) else seq_len(dim(A)[1])
    is <- match(dn[[3]], dy[[3]])
    if (anyNA(ia) || anyNA(is)) {
      mt_abort("The training data of the draws lack some of their ages or strata.", class = "reconcile")
    }
    A <- A[ia, , is, drop = FALSE]
  }
  if (!identical(dim(A)[c(1, 3)], dim(arr)[c(1, 3)])) {
    mt_abort("The training data of the draws do not match their ages and strata.", class = "reconcile")
  }
  A
}

# Observed log rates [age, period, stratum] of the training years for the
# templates of the Schaake shuffle: deaths over exposures, with deaths of 0
# set to 0.5 (package-specific) so that cells without deaths have a value
schaake_log_rates <- function(observed, E, D = NULL, dependence = "schaake") {
  Dobs <- if (!is.null(observed)) observed * E else D
  if (is.null(Dobs) || is.null(E)) {
    mt_abort(sprintf("`dependence = \"%s\"` needs the observed deaths and exposures of the training years.",
                     dependence), class = "reconcile")
  }
  y <- log(pmax(Dobs, 0.5) / E)
  y[!is.finite(y)] <- NA
  y
}

# Draws `arr` [age, period, stratum, draw] of separately fitted strata given
# the dependence of `dependence` ("schaake", "empirical" or "none"), with the
# template from the training data `train` of their fit or draws
# (train_of()). Returns the reordered draws, the base draw of each value
# (`draw`, as dim(arr); NULL without reordering), the base draw of each
# stratum path (`perm` [stratum, draw], "empirical" only) and the report of
# the reordering (`info`)
dependence_reorder <- function(arr, dependence, train) {
  none <- list(arr = arr, draw = NULL, perm = NULL, info = NULL)
  if (identical(dependence, "schaake")) {
    y <- schaake_log_rates(train$observed, train$train_E, train$train_D)
    sr <- schaake_reorder(arr, y, train$periods %||% seq_len(dim(y)[2]))
    return(list(arr = sr$arr, draw = sr$draw, perm = NULL, info = sr$info))
  }
  if (identical(dependence, "empirical")) {
    y <- schaake_log_rates(train$observed, train$train_E, train$train_D, "empirical")
    er <- empirical_reorder(arr, y, train$train_E, train$periods %||% seq_len(dim(y)[2]))
    if (is.null(er$perm)) return(none)
    d <- dim(arr)
    return(list(arr = er$arr, draw = array(rep(as.vector(er$perm), each = d[1] * d[2]), d), perm = er$perm,
                info = er$info))
  }
  none
}

# Training data of a fit or of forecast draws for dependence_reorder()
train_of <- function(x) {
  list(observed = x$observed, train_E = x$train_E, train_D = x$train_D,
       periods = if (inherits(x, "mt_draws")) x$observed_periods else x$periods)
}

warn_schaake <- function(info) {
  if (!length(info$uncoupled)) return(invisible())
  mt_warn(sprintf(paste("The template of the %s has missing observed values for %d stratum/strata (%s):",
                        "their draws keep the order of the base model where the template is missing."),
                  info$method, length(info$uncoupled), paste(utils::head(info$uncoupled, 5), collapse = ", ")),
          class = "reconcile", call = NULL)
}

# Component of each draw (ensembles) after the paths of each stratum were
# permuted (`perm` [stratum, draw], from empirical_reorder()): a matrix
# [stratum, draw]
permute_component <- function(component, perm) {
  if (is.null(component) || is.null(perm)) return(component)
  ns <- nrow(perm)
  cm <- if (is.matrix(component)) component else matrix(component, ns, length(component), byrow = TRUE)
  permute_rows(cm, perm)
}

# Values per stratum and draw [stratum, draw] (components, sizes of the
# negative binomial distribution of the deaths; sizes [age, stratum, draw]
# with a dispersion by age) after the paths of each stratum were permuted
permute_rows <- function(x, perm) {
  if (is.null(x) || is.null(perm)) return(x)
  if (length(dim(x)) == 3L) {
    for (a in seq_len(dim(x)[1])) x[a, , ] <- permute_rows(matrix(x[a, , ], dim(x)[2]), perm)
    return(x)
  }
  matrix(x[cbind(rep(seq_len(nrow(perm)), ncol(perm)), as.vector(perm))], nrow(perm))
}

# Exposures [age, period, stratum] for the forecast: from an mt_data object
# or an array, matched to the forecast by the names of ages, periods and
# strata (dn: dimnames of the forecast rates); a single period is used for
# every forecast period. With the structure `hs` of a hierarchy, exposures
# of the bottom strata only are enough: those of the aggregates are their
# sums (one shared exposure for variables with additive = "rates")
expand_exposure <- function(exposure, d, dn = NULL, hs = NULL) {
  if (!is.null(hs)) {
    Eb <- bottom_exposure(exposure, d, dn, hs)
    if (!is.null(Eb)) {
      E <- array(NA_real_, d[1:3], dimnames = dn[1:3])
      Sx <- hs$Sx %||% hs$S
      for (t in seq_len(d[2])) E[, t, ] <- matrix(Eb[, t, ], d[1]) %*% t(Sx)
      return(E)
    }
  }
  E <- exposure
  if (inherits(E, "mt_data")) {
    if (is.null(E$obs$exposure)) {
      mt_abort("The `exposure` data have no exposures; convert populations with mt_as_exposure().",
               class = "argument")
    }
    E <- mt_array(E, "exposure")
  }
  if (is.data.frame(E)) {
    mt_abort(paste("Give `exposure` as an mt_data object (for example mt_as_exposure(mt_data(population = df,",
                   "by = ...))) or as an array from mt_array()."), class = "argument")
  }
  if (length(dim(E)) == 2) {
    E <- array(E, c(dim(E)[1], 1, dim(E)[2]), dimnames = if (!is.null(dimnames(E))) {
      list(dimnames(E)[[1]], NULL, dimnames(E)[[2]])
    })
  }
  if (length(dim(E)) != 3) mt_abort("`exposure` must be [age, stratum] or [age, period, stratum].", class = "argument")
  en <- dimnames(E)
  if (!is.null(dn) && !is.null(en) && !is.null(en[[1]]) && !is.null(en[[3]])) {
    ia <- match(dn[[1]], en[[1]])
    is <- match(dn[[3]], en[[3]])
    if (anyNA(ia) || anyNA(is)) {
      miss <- c(dn[[1]][is.na(ia)], dn[[3]][is.na(is)])
      mt_abort(sprintf("`exposure` lacks ages or strata of the forecast: %s.", paste(utils::head(miss, 6),
                                                                                       collapse = ", ")),
               class = "argument")
    }
    ip <- if (dim(E)[2] == 1) rep(1L, d[2]) else match(dn[[2]], en[[2]])
    if (anyNA(ip)) {
      mt_abort(sprintf("`exposure` lacks forecast years: %s.", paste(utils::head(dn[[2]][is.na(ip)], 6),
                                                                    collapse = ", ")), class = "argument")
    }
    E <- E[ia, ip, is, drop = FALSE]
    dimnames(E) <- dn[1:3]
    return(E)
  }
  if (is.null(en) && !is.null(dn)) {
    mt_warn(paste("The `exposure` array has no dimnames and is matched to the forecast by position; name its",
                  "dimensions (age, period and stratum labels, as mt_array() returns) to match by name."),
            class = "exposure_order")
  }
  if (dim(E)[2] == 1) E <- E[, rep(1, d[2]), , drop = FALSE]
  if (!all(dim(E) == d[1:3])) mt_abort("`exposure` does not match the draws.", class = "argument")
  E
}

# Exposures [age, period, bottom stratum] when `exposure` holds the bottom
# strata of the hierarchy but not all strata of the forecast (NULL otherwise)
bottom_exposure <- function(exposure, d, dn, hs) {
  E <- exposure
  if (inherits(E, "mt_data")) {
    if (is.null(E$obs$exposure)) return(NULL)
    E <- mt_array(E, "exposure")
  }
  k <- length(dim(E))
  if (!k %in% 2:3) return(NULL)
  b <- hs$bottom
  if (length(b) == d[3]) return(NULL)
  en <- dimnames(E)[[k]]
  if (is.null(en) || is.null(dn[[3]])) {
    # unnamed: one column per bottom stratum
    if (dim(E)[k] != length(b)) return(NULL)
  } else {
    if (all(dn[[3]] %in% en) || !all(dn[[3]][b] %in% en)) return(NULL)
    sel <- match(dn[[3]][b], en)
    E <- if (k == 2) E[, sel, drop = FALSE] else E[, , sel, drop = FALSE]
  }
  expand_exposure(E, c(d[1:2], length(b)), if (!is.null(dn)) c(dn[1:2], list(dn[[3]][b])))
}

# Shrinkage covariance towards the diagonal: the sample correlations shrunk
# towards the identity with the intensity of Schaefer and Strimmer (2005),
# lambda = sum_{i != j} Var(r_ij) / sum_{i != j} r_ij^2, as in MinT
# (Wickramasuriya et al. 2019; 2020, Section 5). `center = TRUE` (draws):
# the sample covariance of the centred rows with divisor n - 1. `center =
# FALSE` (base forecast errors, assumed to have mean zero): (1 / n) sum e e',
# the sample estimator of the error covariance in MinT, with r_ij the mean
# of the products of the standardised errors and Var(r_ij) = sum_k (w_kij -
# mean w_ij)^2 / (n (n - 1))
shrink_cov <- function(X, shrink = "auto", center = TRUE) {
  n <- nrow(X)
  V <- if (center) stats::cov(X) else crossprod(X) / n
  s <- sqrt(pmax(diag(V), 1e-30))
  if (is.numeric(shrink)) {
    lambda <- shrink
  } else {
    Z <- sweep(if (center) sweep(X, 2, colMeans(X)) else X, 2, s, "/")
    R <- crossprod(Z) / (if (center) n - 1 else n)
    Z2 <- Z^2
    v_r <- if (center) {
      (crossprod(Z2) - crossprod(Z)^2 / n) * n / (n - 1)^3
    } else {
      (crossprod(Z2) - crossprod(Z)^2 / n) / (n * (n - 1))
    }
    off <- row(R) != col(R)
    lambda <- sum(v_r[off]) / max(sum(R[off]^2), 1e-30)
    lambda <- min(max(lambda, 0), 1)
  }
  W <- (1 - lambda) * V + lambda * diag(diag(V), ncol(V))
  diag(W) <- pmax(diag(W), 1e-20)
  list(W = W, lambda = lambda)
}

#' Models with reconciled forecasts
#'
#' Wraps a model specification so that its forecasts are reconciled over a
#' hierarchy every time it is forecast, also inside [mt_backtest()]. If the
#' training data contain only the bottom strata, the aggregate strata are
#' added with [mt_aggregate_strata()] before fitting, so that the forecasts
#' cover every node; exposures for the forecast (`exposure` of
#' [mt_forecast()]) can then be given for the bottom strata only, and those
#' of the aggregates are their sums. Predictive draws (used by backtests)
#' are coherent as well: given the reconciled rates, the Poisson variation
#' of the deaths of the bottom strata is drawn independently per cell and
#' draw (none for models of observed rates such as [rwd()]; for [ensemble()]
#' and [stratified()] models as for the component or group each draw comes
#' from; negative binomial for `stanmomo(family = "nb")`), and the
#' aggregates are the sums of the bottom strata.
#'
#' **Bottom-up.** With `method = "bottom_up"`, the base model is fitted to
#' the bottom strata only, and the aggregates are the sums of their draws
#' (of deaths, or of rates for variables with `additive = "rates"`). Wrappers
#' inside the base model ([ensemble()] weights, [tuned()], [recalibrated()])
#' then run their inner backtests on the bottom strata. The fit describes the
#' aggregates by sums as well: their deaths and exposures in the training
#' years, and their fitted rates (fitted deaths of the bottom strata over the
#' exposures). Levels of [stratified()] models that are aggregates of the
#' hierarchy, such as `stratified("cause", all = ...)`, have no strata in
#' this fit and are skipped. Joint models that pool the bottom strata and
#' give each aggregate parameters of its own ([li_lee()], [cae()],
#' [coherent_fdm()], `coda(coherent = TRUE)`) estimate the same parameters of
#' the bottom strata as when fitted to all strata; their draws of the bottom
#' strata have the same distribution and differ by the random numbers drawn.
#' The other methods fit the base model to all strata, and so does
#' bottom-up for models whose fit of the bottom strata uses the aggregates:
#' [tucker3()] and [coda_multiway()] decompose all strata together, [saint()]
#' takes its trend from reference strata, and [coda_causes()] needs the
#' all-cause stratum. Observation noise of predictive draws is drawn before
#' the reordering of `dependence`, for the base draws, and moves with their
#' values; reordering and noise are independent, so this gives the
#' distribution of noise drawn after the reordering, and with the same seed
#' every bottom stratum keeps exactly the predictive draws of `dependence =
#' "none"`.
#'
#' Whether coherence costs accuracy, and which method and base model work
#' best, depends on the data: compare them in a backtest at every level of
#' the hierarchy (`mt_score_table(..., hierarchy = )`), as in the vignette
#' "Coherent forecasts: totals and their parts".
#'
#' **Score-optimal reconciliation.** `method = "energy_score"` follows
#' Panagiotelis et al. (2023). The reconciled draws are
#' \eqn{\tilde y = S (d + G \hat y)}, where \eqn{\hat y} is a base draw of
#' all nodes, \eqn{S} sums the bottom strata to the nodes, and \eqn{d} and
#' \eqn{G} minimise the total energy score of past forecasts. Each time the
#' model is fitted (also at every origin of an outer [mt_backtest()]), the
#' base model is refitted at `n_origins` origins in the last `years` of the
#' training data, as in [recalibrated()], and forecasts the next `horizon`
#' years. Every origin, horizon and age with observed values for all strata
#' is one training forecast. \eqn{d} and \eqn{G} are estimated separately for
#' each age, from the training forecasts of all origins and horizons, on the
#' scale of the other methods (deaths, or rates when the levels share
#' exposures). The estimation is stochastic gradient descent (Adam) started
#' from OLS reconciliation, as in the package ProbReco. The score is that
#' of the predictive distribution the forecast gives: reconciled draws plus
#' the observation noise of the bottom strata. In the forecast, bottom
#' values below zero are set to zero and the aggregates are their sums.
#' The estimates and a table of the training are in the fit (see Value).
#'
#' As in the paper (Section 6, Algorithm 1), \eqn{d} and \eqn{G} are free
#' (no constraint such as \eqn{G S = I}), the objective is the total energy
#' score with exponent 1 estimated from pairs of base draws (their Eq. 5),
#' and Adam starts from OLS reconciliation with their settings (`eta =
#' 0.001`, `beta1 = 0.9`, `beta2 = 0.999`, `epsilon = 1e-8`, their Section
#' 7). The gradient is written out analytically; ProbReco differentiates the
#' same estimate automatically; with the same draws both follow the same
#' optimisation path. The departures from the paper:
#'
#' * training forecasts: Algorithm 1 refits the base models on rolling
#'   windows of the length of the final fit at consecutive origins (250
#'   one-step forecasts in their simulations), and the objective (their Eq.
#'   4) is that of one horizon. Here the inner models use all training years
#'   up to each of `n_origins` origins (expanding windows, shorter than the
#'   final fit), and the forecasts of horizons 1 to `horizon` enter one
#'   objective, so that one \eqn{d} and \eqn{G} serve all horizons;
#' * one \eqn{d} and \eqn{G} per age: every age is a hierarchy of its own;
#' * the score includes the observation noise of the bottom strata, so it
#'   is the score of the predictive distribution of the deaths or rates;
#' * `draws_per_step` is 50 by default (the paper uses 250 draws per score
#'   in its simulations), and the draws come from a fixed set of inner draws
#'   (see below);
#' * stopping: the paper stops when the change in all gradients is less
#'   than 10% of `eta`; ProbReco and mortensor stop when no parameter moves
#'   by more than `tol` (1e-4, 10% of the default `eta`) in a step, or after
#'   `max_iter` steps (500, the default of ProbReco);
#' * the rescaling before the optimisation and the bottom values set to
#'   zero in the forecast are package-specific (reference unknown). The
#'   paper does not treat negative values; the non-negative quadratic
#'   program of the projection methods (see [mt_reconcile()]) does not apply
#'   to the map \eqn{S (d + G \hat y)}.
#'
#' Two details differ from ProbReco. (1) Each step estimates the energy
#' score from two samples of `draws_per_step` base draws of every training
#' forecast, which ProbReco simulates independently. Here they are taken
#' from the `n_draws` inner draws without replacement when there are at
#' least `2 * draws_per_step` of them, and otherwise from two disjoint
#' random halves, so that a draw is never paired with itself (such pairs
#' would contribute a distance of zero to the score). (2) Before the
#' optimisation, the bottom values are divided by their mean observed value
#' and every node, as an input, by its own mean observed value; the
#' estimates are transformed back afterwards. This leaves the optimum
#' unchanged, but Adam moves every parameter by steps of about `eta`, so the
#' rescaling changes the path of the optimisation, and since it stops after
#' `max_iter` steps or when the steps fall below `tol`, also the estimates,
#' relative to ProbReco, which optimises on the original scale. The
#' rescaling is package-specific; it makes the steps comparable across ages
#' and between small and large nodes.
#'
#' The method has many parameters and few training forecasts. \eqn{G} has
#' one parameter per bottom stratum and node at each age (20 x 54 for 2
#' sexes x 2 regions x 5 quintiles with all aggregates); the defaults give
#' about 20 training forecasts per age, which overlap. The only
#' regularisation is the limited number of steps from the OLS start
#' (`control$max_iter`, `control$eta`). \eqn{d} and \eqn{G} then fit the
#' training forecasts much better than OLS, mostly through level corrections
#' of single strata. These help in new years only if the errors of the
#' training years recur: they do not after a change of trend, and the inner
#' models are fitted to shorter series than the final model. In the package's
#' backtests on simulated regional data and on the ONS data by sex and
#' deprivation quintile, the method was less accurate than `"wls"` at every
#' level of the hierarchy, also with fewer and smaller steps. The
#' applications in the paper had longer training windows. Compare the method
#' with `"wls"` or `"mint"` in a backtest before using it.
#'
#' The training refits the base model at `n_origins` origins; the
#' optimisation takes about 3 seconds per age for 54 strata with the
#' defaults. `control$max_iter = 0` gives OLS reconciliation at the cost of
#' the inner forecasts.
#'
#' **Weights from past errors.** `covariance = "errors"` with `method =
#' "wls"` or `"mint"`. MinT weights the base forecasts by the covariance
#' matrix of their `h`-step forecast errors (Wickramasuriya, Turlach and
#' Hyndman 2020, Section 2.1, after Wickramasuriya et al. 2019); Hyndman et
#' al. (2011, Section 3) note that reconciled intervals need this matrix.
#' With `covariance = "errors"` it is estimated from errors of past
#' forecasts. Each time the model is fitted (also at every origin of an
#' outer [mt_backtest()]), the base model is refitted at `n_origins` inner
#' origins in the last `years` of the training data, as for
#' `"energy_score"`, and forecasts the next `horizon` years (cut at the last
#' training year) with `control$n_draws` draws and the observed exposures.
#' The error of a stratum is its observed value minus the mean of its
#' draws, on the scale of the reconciliation (deaths, or rates when the
#' levels share exposures). For each age and inner horizon `h`, `"mint"`
#' takes the shrinkage estimator of MinT: the sample covariance of the `n`
#' error vectors `e`, \eqn{\hat W = n^{-1} \sum e e'} (errors of mean zero),
#' shrunk towards its diagonal, \eqn{\lambda \, diag(\hat W) + (1 - \lambda)
#' \hat W}, with \eqn{\lambda = \sum_{i \ne j} \widehat{Var}(r_{ij}) /
#' \sum_{i \ne j} r_{ij}^2} from the correlations \eqn{r_{ij}} of the errors
#' (Schäfer and Strimmer 2005; Wickramasuriya et al. 2019; Wickramasuriya et
#' al. 2020, Section 5: the estimator "shrinks the sample correlation matrix
#' towards an identity matrix"). Panagiotelis et al. (2023, Table 1) give
#' the same form and write \eqn{\lambda} with the covariances. `"wls"` takes
#' the diagonal, the mean squared errors. `shrink` fixes \eqn{\lambda};
#' `shrink = 1` gives the result of `"wls"`. The draws are then projected as
#' with the default weights, as Panagiotelis et al. (2023, Section 7.3,
#' MinTShr) do for probabilistic forecasts. The errors are reported in the
#' fit (see Value).
#'
#' Deviations: the papers estimate the weights from the in-sample
#' one-step-ahead residuals of the base models and take them as
#' proportional for all horizons (Wickramasuriya et al. 2019). Here the
#' `h`-step errors of inner forecasts estimate the weights of each horizon,
#' one set per age (every age is a hierarchy of its own). Package-specific
#' (reference unknown): forecast years beyond `horizon` take the weights of
#' the longest inner horizon (a projection does not change when its weights
#' are multiplied by a constant, so this is the papers' proportionality
#' beyond `horizon`); an inner horizon with fewer than 2 error vectors with
#' observed values for every stratum takes the weights of the next shorter
#' one; ages without any take the weights from the draws, with a warning.
#' Each age and horizon has at most `n_origins` errors, from overlapping
#' forecasts; with many strata the shrinkage intensity is then often high,
#' and `"mint"` near `"wls"` with the error variances.
#'
#' @param spec Base model specification (fitted to all strata, or to the
#'   bottom strata with `method = "bottom_up"`, see Details).
#' @param hierarchy An [mt_hierarchy()]. When the data hold the aggregate
#'   strata, their exposures are compared with those the hierarchy implies,
#'   with a warning (class `mt_warning_hierarchy_exposure`) when they differ,
#'   for example causes of death declared without `additive = "rates"`.
#' @param method A method of [mt_reconcile()], or `"energy_score"`
#'   (score-optimal reconciliation, see Details).
#' @param shrink See [mt_reconcile()]; with `covariance = "errors"` the
#'   shrinkage intensity for the errors.
#' @param dependence `"none"` (default: the draws are reconciled as they
#'   are), `"ecc"`, `"schaake"`, `"empirical"` or `"auto"` (the same as
#'   `"none"`, kept for compatibility). With `"ecc"` the draws of the base
#'   model take the dependence of the draws of `template` before
#'   reconciliation (ensemble copula coupling, Schefzik, Thorarinsdottir and
#'   Gneiting 2013): the same as `reconciled(coupled(spec, template), ...)`,
#'   see [coupled()]. With `"schaake"` the draws of every age, forecast year
#'   and stratum take the ranks of the observed changes of distinct training
#'   years (the Schaake shuffle of Clark et al. 2004), which also replaces the
#'   dependence between the ages and years of each stratum. With
#'   `"empirical"` the whole paths of each stratum take the ranks of the
#'   observed annual changes of its log age-standardised rate in distinct
#'   training years (the Schaake shuffle with the stratum as margin): every
#'   stratum keeps its paths, and only their pairing across strata changes.
#'   [mt_reconcile()] describes both. Their templates come from the training
#'   data of each fit, also inside [mt_backtest()]. Draws of joint models
#'   such as [li_lee()] and [tucker3()], which carry the dependence between
#'   strata, are never reordered by `"schaake"` or `"empirical"` (they then
#'   give a warning and are ignored). The observation noise of predictive
#'   draws is independent across cells either way.
#' @param covariance `"draws"` (default: the weights of `"wls"` and `"mint"`
#'   from the base draws of each age and forecast year, see
#'   [mt_reconcile()]) or `"errors"` (from the errors of inner forecasts in
#'   the training data, see Details).
#' @param template With `dependence = "ecc"`: the model specification whose
#'   draws give the dependence, usually a joint model such as [cae()] or
#'   [li_lee()].
#' @param fixed Strata whose forecasts are kept (see [mt_reconcile()]); not
#'   with `"energy_score"`.
#' @param years,n_origins,horizon Training of `method = "energy_score"` and
#'   of `covariance = "errors"`: the length of the training window (the last
#'   `years` observed years), the number of inner origins in it (as in
#'   [recalibrated()]) and the number of years each inner forecast covers
#'   (horizons 1 to `horizon`, cut at the last training year).
#' @param control Settings of the training for `method = "energy_score"`
#'   (`n_draws` also for `covariance = "errors"`): a
#'   list with any of `n_draws` (draws per inner forecast, default 100),
#'   `draws_per_step` (draws per training forecast in each step, 50),
#'   `max_iter` (number of steps, 500), `eta` (learning rate, 0.001),
#'   `beta1` (0.9), `beta2` (0.999), `tol` (1e-4) and `epsilon` (1e-8), as in
#'   `ProbReco::scoreopt.control()`. The optimisation stops early when no
#'   parameter moves by more than `tol` in a step, so `tol` must be below
#'   `eta` (0 runs all steps). `max_iter = 0` gives OLS reconciliation.
#' @param seed Random seed of the training for `method = "energy_score"`
#'   and `covariance = "errors"` (inner forecasts, observation noise and the
#'   draws of each step), so that fits are reproducible. The Schaake
#'   shuffles of `dependence` use the `seed` of [mt_forecast()].
#'
#' @return A model specification. With `method = "energy_score"` its fit
#'   holds `$reconcile$scoreopt`: `d` `[bottom stratum, age]`, `G` `[bottom
#'   stratum, stratum, age]` on the scale of the reconciliation, the inner
#'   `origins`, the `control` settings and a `table` per age with the number
#'   of training forecasts, the scale of the bottom strata, the number of
#'   steps, and the mean in-sample energy score with OLS (`energy_ols`) and
#'   with the estimated \eqn{d} and \eqn{G} (`energy_opt`). With
#'   `covariance = "errors"` its fit holds `$reconcile$errors`: the errors
#'   `[age, horizon, stratum, origin]`, the number of error vectors with
#'   observed values for every stratum per age and horizon (`n`) and the
#'   inner `origins`. With `method = "bottom_up"`, `$reconcile$bottom_only`
#'   tells whether the base model was fitted to the bottom strata only. Forecasts
#'   with `dependence = "schaake"` report the shuffle in `$schaake`, those
#'   with `dependence = "empirical"` in `$empirical` (see [mt_reconcile()]).
#' @references
#' Panagiotelis A, Gamakumara P, Athanasopoulos G, Hyndman RJ (2023).
#' Probabilistic forecast reconciliation: properties, evaluation and score
#' optimisation. *European Journal of Operational Research* 306(2):
#' 693-706. \doi{10.1016/j.ejor.2022.07.040}
#'
#' Schefzik R, Thorarinsdottir TL, Gneiting T (2013). Uncertainty
#' quantification in complex simulation models using ensemble copula
#' coupling. *Statistical Science* 28(4): 616-640. \doi{10.1214/13-STS443}
#'
#' Clark M, Gangopadhyay S, Hay L, Rajagopalan B, Wilby R (2004). The Schaake
#' shuffle: a method for reconstructing space-time variability in forecasted
#' precipitation and temperature fields. *Journal of Hydrometeorology* 5(1):
#' 243-262.
#'
#' Hyndman RJ, Ahmed RA, Athanasopoulos G, Shang HL (2011). Optimal
#' combination forecasts for hierarchical time series. *Computational
#' Statistics & Data Analysis* 55(9): 2579-2589.
#'
#' Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast
#' reconciliation for hierarchical and grouped time series through trace
#' minimization. *Journal of the American Statistical Association*
#' 114(526): 804-819. \doi{10.1080/01621459.2018.1448825}
#'
#' Wickramasuriya SL, Turlach BA, Hyndman RJ (2020). Optimal non-negative
#' forecast reconciliation. *Statistics and Computing* 30: 1167-1182.
#' \doi{10.1007/s11222-020-09930-0}
#'
#' Schäfer J, Strimmer K (2005). A shrinkage approach to large-scale
#' covariance matrix estimation and implications for functional genomics.
#' *Statistical Applications in Genetics and Molecular Biology* 4(1):
#' Article 32. \doi{10.2202/1544-6115.1175}
#' @examples
#' x <- mt_example_regions(years = 2000:2019, regions = 3)
#' h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
#' fit <- mt_fit(reconciled(lee_carter(), h, method = "wls"), mt_aggregate_strata(x, h))
#' fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
#' \donttest{
#' # the Schaake shuffle of the separate Lee-Carter forecasts, then MinT with
#' # the covariance of the errors of forecasts from 4 inner origins
#' sm <- reconciled(lee_carter(), h, method = "mint", dependence = "schaake", covariance = "errors",
#'                  years = 8, n_origins = 4, horizon = 3)
#' fit_sm <- mt_fit(sm, mt_aggregate_strata(x, h))
#' fit_sm$reconcile$errors$n[1:3, ]
#' fc_sm <- mt_forecast(fit_sm, h = 5, n_draws = 50, seed = 1)
#' fc_sm$schaake$dates
#' # score-optimal reconciliation, trained on forecasts of 2012-2019 from 4 origins
#' es <- reconciled(lee_carter(), h, method = "energy_score", years = 8, n_origins = 4, horizon = 3,
#'                  control = list(max_iter = 200))
#' fit_es <- mt_fit(es, mt_aggregate_strata(x, h))
#' fit_es$reconcile$scoreopt$table
#' fc_es <- mt_forecast(fit_es, h = 5, n_draws = 50, seed = 1)
#' }
#' @export
reconciled <- function(spec, hierarchy,
                       method = c("wls", "mint", "structural", "ols", "bottom_up", "proportional", "energy_score"),
                       shrink = "auto", dependence = c("none", "ecc", "schaake", "empirical", "auto"), fixed = NULL,
                       years = 10, n_origins = 5, horizon = 5, control = list(), seed = 1, template = NULL,
                       covariance = c("draws", "errors")) {
  check_class(spec, "mt_spec", "spec")
  method <- match.arg(method)
  dependence <- match.arg(dependence)
  covariance <- match.arg(covariance)
  if (dependence == "auto") dependence <- "none"
  if (covariance == "errors" && !method %in% c("wls", "mint")) {
    mt_abort("`covariance = \"errors\"` sets the weights of `method = \"wls\"` or `\"mint\"`.", class = "argument")
  }
  if (dependence == "ecc") {
    if (is.null(template)) mt_abort("`dependence = \"ecc\"` needs a `template` model.", class = "argument")
    spec <- coupled(spec, template)
  } else if (!is.null(template)) {
    mt_abort("`template` is used with `dependence = \"ecc\"` only.", class = "argument")
  }
  dependence <- joint_dependence(dependence, spec)
  if (length(fixed) && method == "bottom_up") {
    mt_abort("`fixed` needs a projection method (\"wls\", \"mint\", \"structural\", \"ols\") or \"proportional\".",
             class = "argument")
  }
  if (method == "proportional" && !length(fixed)) {
    mt_abort("`method = \"proportional\"` scales the parts of fixed nodes: give `fixed`.", class = "argument")
  }
  params <- list(base = spec, hierarchy = as_hierarchy(hierarchy), method = method, shrink = shrink,
                 dependence = dependence, fixed = fixed, covariance = covariance)
  min_periods <- spec$min_periods
  if (method == "energy_score" || covariance == "errors") {
    if (length(fixed) && method == "energy_score") {
      mt_abort(paste("`fixed` is not available with `method = \"energy_score\"`, which learns the map of all nodes;",
                     "use \"wls\", \"mint\" or \"proportional\" to keep fixed nodes."), class = "argument")
    }
    check_number(years, "years", lower = 2, integerish = TRUE)
    check_number(n_origins, "n_origins", lower = 1, integerish = TRUE)
    check_number(horizon, "horizon", lower = 1, integerish = TRUE)
    if (!is.null(seed)) check_number(seed, "seed", integerish = TRUE)
    if (covariance == "errors" && n_origins < 2) {
      mt_abort("`covariance = \"errors\"` needs at least 2 inner origins (`n_origins`).", class = "argument")
    }
    params <- c(params, list(years = years, n_origins = min(n_origins, years - 1), horizon = horizon,
                             control = scoreopt_control(control), seed = seed))
    min_periods <- spec$min_periods + years
  }
  new_mt_spec("reconciled", params, min_periods = min_periods,
              label = sprintf("%s, reconciled (%s%s)", spec$label, method,
                              if (covariance == "errors") ", error covariance" else ""))
}

#' @export
mt_fit.mt_spec_reconciled <- function(spec, data, prep = NULL, cores = 1, ...) {
  p <- spec$params
  hs <- hierarchy_structure(data$strata, p$hierarchy)
  if (all(hs$is_bottom)) data <- mt_aggregate_strata(data, p$hierarchy) else check_hierarchy_exposure(data, hs)
  hs <- hierarchy_structure(data$strata, p$hierarchy)
  # bottom-up uses the forecasts of the bottom strata only: the base model is
  # fitted to them, unless its fit of the bottom strata uses the aggregates
  bottom_only <- identical(p$method, "bottom_up") && !length(p$fixed) && fits_bottom_alone(p$base)
  if (bottom_only) {
    bottom <- strata_subset(data, hs$bottom)
    bottom$hierarchy <- p$hierarchy
    fit <- mt_fit(p$base, bottom, prep = prep, cores = cores, ...)
    if (!identical(fit$strata$label, data$strata$label[hs$bottom])) {
      mt_abort("The base model was fitted to other strata than the bottom strata of the hierarchy.",
               class = "reconcile")
    }
  } else {
    fit <- mt_fit(p$base, data, prep = prep, cores = cores, ...)
    hs <- hierarchy_structure(fit$strata, p$hierarchy)
  }
  # the fitted specs can be joint where the spec could not tell (tuned())
  p$dependence <- joint_dependence(p$dependence, c(list(fit$spec), fit$inner_specs))
  strata <- if (bottom_only) data$strata else fit$strata
  rc <- list(hierarchy = p$hierarchy, method = p$method, shrink = p$shrink, dependence = p$dependence,
             structure = hs, fixed = fixed_nodes(p$fixed, strata, p$hierarchy, hs, p$method),
             covariance = p$covariance %||% "draws", bottom_only = bottom_only)
  # score-optimal reconciliation: d and G per age from inner forecasts
  if (identical(p$method, "energy_score")) rc$scoreopt <- scoreopt_train(p, data, prep, fit, hs, cores)
  # weights of WLS and MinT from the errors of inner forecasts
  if (identical(rc$covariance, "errors")) rc$errors <- errors_train(p, data, prep, fit, hs, cores)
  fields <- c(list(reconcile = rc), if (bottom_only) node_fields(fit, strata, hs))
  wrap_fit(fit, spec, "mt_fit_reconciled", fields)
}

# Whether the fit of the bottom strata of a hierarchy does not depend on the
# aggregates: models fitted per stratum, and joint models that pool the
# bottom strata only and give each aggregate parameters of its own
# (li_lee(), cae(), coherent_fdm(), coda(coherent = TRUE)). The fits of
# tucker3() and coda_multiway() decompose all strata together, saint()
# takes its trend from reference strata (often aggregates), and
# coda_causes() needs the all-cause stratum: they are fitted to all strata.
# Wrappers, ensembles and stratified models decide by the models inside
# them; the candidates of tuned() are created from its grid (FALSE when that
# fails)
fits_bottom_alone <- function(spec) {
  if (!inherits(spec, "mt_spec")) return(TRUE)
  if (spec$model %in% c("tucker3", "coda_multiway", "saint", "coda_causes")) return(FALSE)
  # relative_model(): reference strata of the data (often aggregates), or a
  # reference population of its own
  if (spec$model == "relative") return(inherits(spec$params$reference, "mt_data"))
  p <- spec$params
  inner <- switch(spec$model,
                  ensemble = p$components,
                  stratified = c(p$specs, if (!is.null(p$default)) list(p$default)),
                  coupled = list(p$base, p$template),
                  tuned = tryCatch(lapply(seq_len(nrow(p$grid)), function(i) {
                    tune_candidate(p$spec_fun, as.list(p$grid[i, , drop = FALSE]), NULL)$spec
                  }), error = function(e) NULL),
                  if (inherits(p$base, "mt_spec")) list(p$base) else list())
  if (is.null(inner)) return(FALSE)
  all(vapply(inner, fits_bottom_alone, logical(1)))
}

# Fields of a fit to the bottom strata of a hierarchy (structure `hs`) for
# all its nodes (`strata`): deaths of each node summed over its bottom
# strata, exposures summed or shared (additive = "rates"), observed rates
# their ratio. The fits of the base model per stratum (fits, components,
# group_fits) refer to the bottom strata and are left out; the base fit
# keeps them (unwrap_fit(), bottom_base_fit())
node_fields <- function(fit, strata, hs) {
  agg <- which(!hs$is_bottom)
  up <- function(A, M) {
    if (is.null(A)) return(NULL)
    d <- dim(A)
    dn <- dimnames(A)
    out <- array(NA_real_, c(d[1:2], nrow(strata)),
                 dimnames = stats::setNames(list(dn[[1]], dn[[2]], strata$label), names(dn)))
    out[, , hs$bottom] <- A
    if (length(agg)) {
      for (t in seq_len(d[2])) out[, t, agg] <- matrix(A[, t, ], d[1]) %*% t(M[agg, , drop = FALSE])
    }
    out
  }
  E <- up(fit$train_E, hs$Sx %||% hs$S)
  obs <- up(fit$raw_D %||% fit$train_D, hs$S) / E
  obs[, , hs$bottom] <- fit$observed
  list(strata = strata, observed = obs, train_D = up(fit$train_D, hs$S), raw_D = up(fit$raw_D, hs$S),
       train_E = E, exposure_last = E[, dim(E)[2], , drop = FALSE], fits = NULL, components = NULL,
       group_fits = NULL, groups = NULL)
}

# The fit of the base model to the bottom strata of a reconciled model with
# `method = "bottom_up"` (`fit$reconcile$bottom_only`), also below wrappers
# around it (anchored(), recalibrated(), coupled(), bootstrapped()); NULL for
# other fits
bottom_base_fit <- function(fit) {
  k <- match("mt_fit_reconciled", class(fit))
  if (is.na(k) || !isTRUE(fit$reconcile$bottom_only)) return(NULL)
  for (cl in class(fit)[seq_len(k - 1)]) fit <- unwrap_fit(fit, cl)
  unwrap_fit(fit, "mt_fit_reconciled")
}

#' @export
fitted.mt_fit_reconciled <- function(object, ...) {
  base <- bottom_base_fit(object)
  if (is.null(base)) return(NextMethod())
  hs <- object$reconcile$structure
  fb <- stats::fitted(base)
  # fitted deaths of the bottom strata summed to the nodes
  out <- node_fields(list(train_D = fb * base$train_E, train_E = base$train_E, observed = fb), object$strata,
                     hs)$observed
  out[, , hs$bottom] <- fb
  out
}

#' @export
mt_forecast.mt_fit_reconciled <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                                          exposure = NULL, ...) {
  # the base fit without the hierarchy (unwrap_fit() restores its own fields),
  # so that its draws are those of the base model, not already summed
  base_fit <- unwrap_fit(fit, "mt_fit_reconciled")
  rc <- fit$reconcile
  hs <- rc$structure
  bottom_only <- isTRUE(rc$bottom_only)
  out <- with_seed(seed, {
    fc <- mt_forecast(base_fit, h = if (is.null(to)) h else to - max(fit$periods), n_draws = n_draws)
    dn <- dimnames(fc$rates)
    dn[[3]] <- fit$strata$label
    d <- c(dim(fc$rates)[1:2], nrow(fit$strata), n_draws)
    E <- expand_exposure(exposure %||% fit$exposure_last, d, dn, hs)
    res <- if (bottom_only) {
      bottom_up_forecast(fc, base_fit, rc$dependence, hs, rc$hierarchy$additive, E, predictive)
    } else {
      dep <- dependence_reorder(fc$rates, rc$dependence, train_of(base_fit))
      rec <- if (identical(rc$method, "energy_score")) {
        scoreopt_apply(dep$arr, hs, rc$hierarchy$additive, rc$scoreopt, E, d)
      } else {
        W_cell <- if (!is.null(rc$errors)) function(a, t) errors_weights(rc$errors, a, t, rc$method, rc$shrink)
        reconcile_array(dep$arr, hs, rc$hierarchy$additive, rc$method, rc$shrink, E, d,
                        fixed = rc$fixed %||% integer(), W_cell = W_cell)
      }
      rates <- rec$rates
      if (predictive) {
        # observation noise of the bottom strata given the reconciled rates,
        # from the component or group of the base draw each value comes from,
        # independent across cells; the aggregates are their sums
        rates <- add_observation_noise(rates, base_fit, fc, hs$bottom, rec$exposure, draw = dep$draw)
        rates <- aggregate_bottom(rates, rec$exposure, hs, rc$hierarchy$additive)
      }
      list(rates = rates, info = rec$info, dep = dep)
    }
    res$fc <- fc
    res
  })
  dep <- out$dep
  warn_schaake(dep$info)
  res <- out$fc
  if (bottom_only) {
    # the draws describe all nodes, with the training data of the fit
    full <- new_mt_draws(out$rates, fit$age, res$periods, fit$strata, model = fit$spec$label, origin = res$origin,
                         fit = fit)
    for (n in c("strata", "observed", "fitted", "train_D", "train_E", "exposure_last")) res[n] <- list(full[[n]])
  }
  res$rates <- out$rates
  dimnames(res$rates) <- dn
  res$model <- fit$spec$label
  res$predictive <- predictive
  res$reconciliation <- out$info
  res$hierarchy <- rc$hierarchy
  # the component (ensembles) and the size of the negative binomial
  # distribution of each draw: per stratum and draw after the shuffle of
  # stratum paths, none after the Schaake shuffle (the values of one draw
  # come from different base draws)
  comp <- res$component
  size <- res$nb_size
  if (identical(dep$info$method, "Schaake shuffle")) {
    comp <- NULL
    size <- NULL
    res$schaake <- dep$info
  }
  if (!is.null(dep$perm)) {
    comp <- permute_component(comp, dep$perm)
    size <- permute_rows(size, dep$perm)
    res$empirical <- dep$info
  }
  res$component <- if (bottom_only) node_rows(comp, hs) else comp
  res$nb_size <- if (bottom_only) node_rows(size, hs) else size
  res
}

# Bottom-up forecast of all nodes (structure `hs`, exposures `E` of all
# nodes) from the draws `fc` of a base model fitted to the bottom strata.
# Observation noise of predictive draws is drawn first, independently per
# cell given the base rates, and every value keeps its noise when the draws
# are reordered for the dependence between the strata (dependence_reorder(),
# which ranks the draws without noise): reordering and noise commute in
# distribution, and so every bottom stratum keeps exactly the draws it has
# with dependence = "none". The aggregates are the sums of the bottom strata
bottom_up_forecast <- function(fc, base_fit, dependence, hs, additive, E, predictive) {
  arr <- fc$rates
  d <- dim(arr)
  neg <- !is.na(arr) & arr < 0
  arr <- pmax(arr, 0)
  vals <- if (predictive) add_observation_noise(arr, base_fit, fc, seq_len(d[3]), E[, , hs$bottom, drop = FALSE])
  dep <- dependence_reorder(arr, dependence, train_of(base_fit))
  vals <- if (is.null(vals)) dep$arr else if (is.null(dep$draw)) vals else take_draws(vals, dep$draw)
  full <- array(NA_real_, c(d[1:2], nrow(hs$S), d[4]))
  full[, , hs$bottom, ] <- vals
  rates <- aggregate_bottom(full, E, hs, additive)
  if (!predictive) rates <- pmax(rates, 1e-12)
  info <- list(method = "bottom_up", max_relative_gap_before = NA_real_, share_set_to_zero = mean(neg),
               mean_shrinkage = NA_real_, share_draws_with_zeros = mean(apply(neg, 4, any)),
               max_relative_change_fixed = NA_real_)
  list(rates = rates, info = info, dep = dep)
}

# Draws `A` [age, period, stratum, draw] reordered value by value: the value
# of cell (a, t, s) in draw j is that of draw `draw[a, t, s, j]`
take_draws <- function(A, draw) {
  d <- dim(A)
  L <- prod(d[1:3])
  array(A[rep(seq_len(L), d[4]) + (as.vector(draw) - 1) * L], d, dimnames = dimnames(A))
}

# Rows [bottom stratum, draw] of a matrix for all nodes of a hierarchy
# (aggregates NA; [age, stratum, draw] for sizes with a dispersion by age);
# vectors (one value per draw) and NULL are kept
node_rows <- function(x, hs) {
  if (length(dim(x)) == 3L) {
    out <- array(NA_real_, c(dim(x)[1], nrow(hs$S), dim(x)[3]))
    out[, hs$bottom, ] <- x
    return(out)
  }
  if (!is.matrix(x)) return(x)
  out <- matrix(NA, nrow(hs$S), ncol(x))
  storage.mode(out) <- storage.mode(x)
  out[hs$bottom, ] <- x
  out
}

# Base forecast errors for the weights of WLS and MinT (`covariance =
# "errors"`; Wickramasuriya et al. 2020, Section 2.1: the covariance of the
# h-step base forecast errors). The base
# model is refitted at the inner origins in the last `years` training years
# (as for "energy_score"), forecasts the next `horizon` years (cut at the
# last training year) with the observed exposures, and the error of each
# node is its observed value minus the mean of its draws, on the scale of
# the reconciliation (deaths, or rates when the levels share exposures).
# Returns the errors [age, horizon, node, origin], the number of complete
# error vectors (all nodes observed) per age and horizon, and the origins
errors_train <- function(p, data, prep, fit, hs, cores = 1) {
  ctl <- scoreopt_control(p$control)
  g <- mt_grid(data)
  last <- max(g$period$lower)
  first_ok <- min(g$period$lower) + p$base$min_periods - 1
  origins <- calibration_origins(last, p$years, p$n_origins)
  origins <- origins[origins >= first_ok]
  if (length(origins) < 2) {
    mt_abort(sprintf(paste("`covariance = \"errors\"` needs at least 2 inner origins with %d training years each:",
                           "give a shorter `years` or more training years."), p$base$min_periods),
             class = "model_data")
  }
  p0 <- p
  # the mean of the draws of each cell does not depend on their order
  p0$dependence <- "none"
  pieces <- with_seed(p$seed, par_lapply(origins, function(o) {
    tryCatch(scoreopt_origin(p0, data, o, last, prep, ctl$n_draws, hs, fit$strata, p$hierarchy$additive,
                             noise = FALSE), error = function(e) e)
  }, cores, streams = TRUE))
  bad <- vapply(pieces, inherits, logical(1), what = "error")
  if (any(bad)) {
    mt_warn(sprintf("Inner origin(s) %s of the error covariance failed: %s", paste(origins[bad], collapse = ", "),
                    conditionMessage(pieces[[which(bad)[1]]])), class = "reconcile")
  }
  keep <- !bad & !vapply(pieces, is.null, logical(1))
  pieces <- pieces[keep]
  if (length(pieces) < 2) mt_abort("Fewer than 2 inner origins of the error covariance succeeded.", class = "reconcile")
  K <- p$horizon
  err <- array(NA_real_, c(length(fit$age), K, nrow(fit$strata), length(pieces)),
               dimnames = list(age = fit$age$labels, horizon = seq_len(K), stratum = fit$strata$label,
                               origin = origins[keep]))
  for (i in seq_along(pieces)) {
    pc <- pieces[[i]]
    k <- dim(pc$Y)[2]
    err[, seq_len(k), , i] <- pc$Y - rowMeans(pc$X, dims = 3)
  }
  err[!is.finite(err)] <- NA
  n <- rowSums(apply(!is.na(err), c(1, 2, 4), all), dims = 2)
  dimnames(n) <- dimnames(err)[1:2]
  none <- rowSums(n >= 2) == 0
  if (any(none)) {
    mt_warn(sprintf(paste("Age(s) %s have fewer than 2 inner forecasts with observed values for every stratum: their",
                          "weights come from the base draws."), paste(utils::head(fit$age$labels[none], 6),
                                                                      collapse = ", ")), class = "reconcile")
  }
  list(errors = err, n = n, origins = origins[keep])
}

# Weights of forecast year `t` at age `a` from the inner errors: those of
# the longest inner horizon up to `t` with at least 2 complete error vectors
# (forecast years beyond the inner horizons use the last one). "wls": the
# mean squared errors; "mint": shrink_cov() of the errors, not centred.
# NULL when the age has no such horizon
errors_weights <- function(errs, a, t, method, shrink = "auto") {
  K <- dim(errs$errors)[2]
  for (k in rev(seq_len(min(t, K)))) {
    if (errs$n[a, k] < 2) next
    X <- t(matrix(errs$errors[a, k, , ], dim(errs$errors)[3]))
    X <- X[rowSums(is.na(X)) == 0, , drop = FALSE]
    if (method == "wls") return(pmax(colMeans(X^2), 1e-20))
    return(shrink_cov(X, shrink, center = FALSE))
  }
  NULL
}

#' Aggregate forecast draws over strata or ages
#'
#' Computes draws for coarser groups from draws on the grid, path by path:
#' over strata along a hierarchy (bottom-up: deaths are summed and divided
#' by summed exposures) and over age groups (exposure-weighted mean rates).
#'
#' @param draws An [mt_draws] object.
#' @param hierarchy Optional [mt_hierarchy()]; aggregate strata are added
#'   bottom-up from the bottom strata.
#' @param age Optional coarser [mt_axis] (cells must nest).
#' @param exposure Exposures `[age, period, stratum]` or `[age, stratum]`
#'   (default: last training year).
#' @return An `mt_draws` object.
#' @examples
#' x <- mt_example(deaths = "single", years = 2000:2019)
#' fc <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 20, seed = 1)
#' fa <- mt_aggregate_draws(fc, hierarchy = mt_hierarchy(sex = "both"),
#'                          age = age_groups(c(30, seq(40, 100, 10))))
#' fa
#' @export
mt_aggregate_draws <- function(draws, hierarchy = NULL, age = NULL, exposure = NULL) {
  check_class(draws, "mt_draws", "draws")
  d <- dim(draws$rates)
  st <- draws$strata
  h <- if (!is.null(hierarchy)) as_hierarchy(hierarchy)
  hs0 <- if (!is.null(h)) hierarchy_structure(st, h)
  E <- expand_exposure(exposure %||% draws$exposure_last, d, dimnames(draws$rates), hs0)
  X <- draws$rates * as.vector(E)
  if (!is.null(hierarchy)) {
    by <- setdiff(names(st), c("stratum", "label"))
    check_groups(st[hs0$is_bottom], h)
    nodes <- hierarchy_nodes(st[hs0$is_bottom], h, by)
    nodes[, stratum := seq_len(.N)]
    setcolorder(nodes, c("stratum", by, "label"))
    hs <- hierarchy_structure(nodes, h)
    bpos <- match(nodes$label[hs$bottom], st$label)
    Xb <- X[, , bpos, , drop = FALSE]
    Eb <- E[, , bpos, drop = FALSE]
    Xn <- array(NA_real_, c(d[1], d[2], nrow(nodes), d[4]))
    En <- array(NA_real_, c(d[1], d[2], nrow(nodes)))
    for (t in seq_len(d[2])) {
      En[, t, ] <- matrix(Eb[, t, ], d[1]) %*% t(hs$Sx %||% hs$S)
      for (j in seq_len(d[4])) Xn[, t, , j] <- matrix(Xb[, t, , j], d[1]) %*% t(hs$S)
    }
    X <- Xn
    E <- En
    st <- nodes
  }
  ages <- draws$age
  if (!is.null(age)) {
    check_class(age, "mt_axis", "age")
    ia <- axis_locate(age, ages$lower, ages$upper)
    if (anyNA(ia)) mt_abort("The grid ages do not nest in `age`.", class = "argument")
    A <- matrix(0, length(age), length(ages))
    A[cbind(ia, seq_along(ia))] <- 1
    keep <- rowSums(A) > 0
    A <- A[keep, , drop = FALSE]
    X <- apply(X, 2:4, function(v) as.numeric(A %*% v))
    E <- apply(E, 2:3, function(v) as.numeric(A %*% v))
    if (length(dim(X)) == 3) X <- array(X, c(1, dim(X)))
    if (length(dim(E)) == 2) E <- array(E, c(1, dim(E)))
    ages <- mt_axis(age$lower[keep], age$upper[keep], age$labels[keep], kind = "age")
  }
  rates <- X / as.vector(E)
  out <- draws
  dimnames(rates) <- list(age = ages$labels, period = as.character(draws$periods), stratum = st$label, draw = NULL)
  out$rates <- rates
  out$age <- ages
  out$strata <- st
  out$exposure_last <- E[, 1, , drop = FALSE]
  out$observed <- NULL
  out$fitted <- NULL
  # components per stratum and draw (reordered ensembles) do not describe sums of strata
  if (is.matrix(out$component)) out$component <- NULL
  out
}

#' Check the coherence of forecast draws
#'
#' For every aggregate stratum of a hierarchy, the ratio of the sum of its
#' bottom strata to the aggregate itself, on the deaths scale (rates times
#' the exposures of the last observed year, or `exposure`) or on the rate
#' scale for variables whose levels share an exposure (causes of death).
#' Coherent forecasts give ratios of 1 in every draw; separately forecast
#' strata show how far they disagree.
#'
#' @param draws An [mt_draws] object with the strata of the hierarchy.
#' @param hierarchy An [mt_hierarchy()].
#' @param exposure Optional exposures (see [mt_deaths()]).
#' @return `mt_coherence()`: a `data.table` with the aggregate strata,
#'   `period`, and the median and the 2.5% and 97.5% quantiles of the ratio
#'   over ages and draws. `mt_plot_coherence()`: a ggplot object with the
#'   distribution of the median ratios of the aggregates per year and the
#'   range of the ratios of single paths.
#' @examples
#' x <- mt_example_regions(years = 2005:2019, regions = 2)
#' h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
#' fc <- mt_forecast(mt_fit(lee_carter(), mt_aggregate_strata(x, h)), h = 3, n_draws = 30, seed = 1)
#' head(mt_coherence(fc, h))
#' head(mt_coherence(mt_reconcile(fc, h, method = "wls"), h))
#' @export
mt_coherence <- function(draws, hierarchy, exposure = NULL) {
  check_class(draws, "mt_draws", "draws")
  h <- as_hierarchy(hierarchy)
  hs <- hierarchy_structure(draws$strata, h)
  d <- dim(draws$rates)
  E <- expand_exposure(exposure %||% draws$exposure_last, d, dimnames(draws$rates), hs)
  by_deaths <- identical(h$additive, "deaths")
  V <- if (by_deaths) draws$rates * as.vector(E) else draws$rates
  agg <- which(!hs$is_bottom)
  if (!length(agg)) {
    mt_abort(paste("The draws hold no aggregate strata of the hierarchy, only its bottom strata, so there is no",
                   "sum to check. Forecast the aggregates as well: fit the data aggregated with",
                   "mt_aggregate_strata(), or fit reconciled(), whose draws include the aggregates."),
             class = "hierarchy")
  }
  b <- hs$bottom
  rows <- lapply(agg, function(i) {
    w <- hs$S[i, ]
    part <- apply(sweep(V[, , b, , drop = FALSE], 3, w, "*"), c(1, 2, 4), sum)
    # [age, period, draw] also with one period or one draw
    r <- part / array(V[, , i, ], dim(part))
    q <- t(apply(matrix(aperm(r, c(1, 3, 2)), ncol = d[2]), 2, stats::quantile, c(0.025, 0.5, 0.975),
                 na.rm = TRUE))
    data.table(stratum = i, period = draws$periods, q0.025 = q[, 1], q0.5 = q[, 2], q0.975 = q[, 3])
  })
  out <- merge(draws$strata, rbindlist(rows), by = "stratum", sort = FALSE)
  out[, stratum := NULL]
  out[]
}

#' @rdname mt_coherence
#' @export
mt_plot_coherence <- function(draws, hierarchy, exposure = NULL) {
  need_ggplot()
  tab <- mt_coherence(draws, hierarchy, exposure)
  rng <- tab[, list(lo = min(q0.025), hi = max(q0.975)), by = period]
  ggplot2::ggplot(tab, ggplot2::aes(x = factor(period), y = q0.5)) +
    ggplot2::geom_linerange(data = rng, ggplot2::aes(y = NULL, ymin = lo, ymax = hi), colour = "grey60",
                            linewidth = 2, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "#D7191C") +
    ggplot2::geom_boxplot(fill = "#66C2A5", alpha = 0.6, outlier.size = 0.6) +
    ggplot2::labs(x = NULL, y = "Sum of parts / aggregate", title = "Coherence of the forecast",
                  subtitle = sprintf(paste("%s; boxes: medians over ages and draws per aggregate stratum; grey: range",
                                           "of the 2.5%%-97.5%% quantiles of single paths; 1 = coherent"),
                                     draws$model)) +
    theme_mortensor()
}
