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
      return(list(total = s$total, groups = s$groups %||% list(), additive = a))
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
#' to check that they agree with the sums.
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
#' * `"mint"`: covariance of the base draws across strata, shrunk towards
#'   its diagonal (Schaefer-Strimmer intensity by default). This follows the
#'   minimum-trace reconciliation of Wickramasuriya et al. (2019), with the
#'   covariance of the forecast draws in place of in-sample residuals;
#' * `"wls"`: diagonal with the variance of each stratum's draws;
#' * `"structural"`: diagonal with the number of bottom series in each node;
#' * `"ols"`: identity;
#' * `"bottom_up"`: aggregates are the sums of the bottom forecasts.
#'
#' Projections can make bottom values negative where the base forecasts
#' disagree strongly. Such values are held at zero and the remaining bottom
#' strata of the draw are projected again, until none is negative (a simple
#' form of the non-negative reconciliation of Wickramasuriya et al. 2020);
#' the aggregates are recomputed from the bottom series, so the result is
#' coherent and non-negative, and fixed nodes are kept. A warning reports
#' when more than 5% of the draws needed this, because the distributions of
#' the affected strata (often small causes) are then truncated at zero.
#'
#' `"proportional"` keeps the `fixed` nodes and scales the bottom strata
#' below each of them in proportion to their own forecasts, draw by draw
#' (top-down allocation with forecast proportions, Athanasopoulos, Ahmed
#' and Hyndman 2009). It never sets a part to zero and keeps the shape of
#' the parts' distributions; it is the usual choice when causes of death
#' must add up to a trusted all-cause forecast. Fixed nodes must not share
#' parts.
#'
#' Draws of strata forecast by separate models are independent, so their sums
#' would be too narrow. With `dependence = "empirical"` (default), the paths
#' of each stratum are first reordered so that their ranks follow the
#' co-movement of the strata in the training data (an empirical copula, the
#' "Schaake shuffle" of Clark et al. 2004): every draw is matched to a
#' training year, and the stratum whose average year-to-year change of log
#' rates was relatively large in that year receives a relatively high path.
#' The marginal distribution of each stratum is unchanged.
#'
#' The one-variable interface with `var` and `total` is still accepted.
#'
#' @param draws An [mt_draws] object with the strata of the hierarchy.
#' @param hierarchy An [mt_hierarchy()] (or, for one variable, `var` and
#'   `total`).
#' @param method `"wls"` (default), `"mint"`, `"structural"`, `"ols"`,
#'   `"bottom_up"` or `"proportional"` (with `fixed`). The vignette
#'   "Coherent forecasts: totals and their parts" compares them.
#' @param shrink Shrinkage intensity for `"mint"` (number between 0 and 1),
#'   or `"auto"`.
#' @param exposure Exposures `[age, period, stratum]` or `[age, stratum]` for
#'   the forecast years (default: last training year). Only bottom strata are
#'   used; aggregates are summed.
#' @param dependence `"empirical"` or `"none"` (keep the draws as they are).
#' @param seed Random seed for the reordering.
#' @param fixed Strata whose forecasts are kept unchanged (labels, or
#'   `"total"` for the grand total of the hierarchy), for example a trusted
#'   all-cause forecast that the causes must add up to. The projection gives
#'   these nodes zero variance, so the other strata absorb the whole gap;
#'   `$reconciliation$max_relative_change_fixed` reports how closely they
#'   were kept; they move, with a warning, when their own forecasts do not
#'   add up (a node fixed together with all its parts) or when every bottom
#'   stratum below them would have to be negative. Needs a projection
#'   method or `"proportional"`.
#' @param var,total,parts One-variable form: stratum variable and its total
#'   level (`parts` is ignored: all other levels are parts).
#' @return The reconciled `mt_draws`, with diagnostics in `$reconciliation`
#'   (method, largest relative gap before reconciliation, share of bottom
#'   values held at zero, mean shrinkage intensity).
#' @references
#' Wickramasuriya SL, Athanasopoulos G, Hyndman RJ (2019). Optimal forecast
#' reconciliation for hierarchical and grouped time series through trace
#' minimization. *Journal of the American Statistical Association*
#' 114(526): 804-819. \doi{10.1080/01621459.2018.1448825}
#'
#' Athanasopoulos G, Ahmed RA, Hyndman RJ (2009). Hierarchical forecasts for
#' Australian domestic tourism. *International Journal of Forecasting*
#' 25: 146-166.
#'
#' Wickramasuriya SL, Turlach BA, Hyndman RJ (2020). Optimal non-negative
#' forecast reconciliation. *Statistics and Computing* 30: 1167-1182.
#' \doi{10.1007/s11222-020-09930-0}
#'
#' Clark M, Gangopadhyay S, Hay L, Rajagopalan B, Wilby R (2004). The Schaake
#' shuffle: a method for reconstructing space-time variability in forecasted
#' precipitation and temperature fields. *Journal of Hydrometeorology* 5(1):
#' 243-262.
#'
#' Schaefer J, Strimmer K (2005). A shrinkage approach to large-scale
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
#' @export
mt_reconcile <- function(draws, hierarchy = NULL,
                         method = c("wls", "mint", "structural", "ols", "bottom_up", "proportional"),
                         shrink = "auto", exposure = NULL, dependence = c("empirical", "none"), seed = NULL,
                         fixed = NULL, var = NULL, total = NULL, parts = NULL) {
  check_class(draws, "mt_draws", "draws")
  method <- match.arg(method)
  dependence <- match.arg(dependence)
  if (is.character(hierarchy) && length(hierarchy) == 1 && is.null(var)) {
    var <- hierarchy
    hierarchy <- NULL
  }
  if (is.null(hierarchy)) {
    if (is.null(var) || is.null(total)) mt_abort("Give a hierarchy or `var` and `total`.", class = "argument")
    if (!var %in% names(draws$strata)) mt_abort("`var` is not a stratum variable.", class = "argument")
    if (!total %in% draws$strata[[var]]) mt_abort("No total with parts found for the hierarchy.", class = "argument")
    hierarchy <- do.call(mt_hierarchy, c(stats::setNames(list(total), var),
                                         list(additive = if (is.null(exposure)) "rates" else "deaths")))
  }
  h <- as_hierarchy(hierarchy)
  hs <- tryCatch(hierarchy_structure(draws$strata, h), mt_error_hierarchy = function(e) {
    mt_abort(conditionMessage(e), class = "argument")
  })
  fixed_idx <- fixed_nodes(fixed, draws$strata, h, hs, method)
  arr <- draws$rates
  if (dependence == "empirical" && !is.null(draws$train_D)) {
    Dc <- if (!is.null(draws$observed)) draws$observed * draws$train_E else draws$train_D
    arr <- with_seed(seed, copula_reorder(arr, Dc, draws$train_E)$arr)
  }
  rec <- reconcile_array(arr, hs, h$additive, method, shrink, exposure %||% draws$exposure_last, dim(arr),
                         fixed = fixed_idx)
  draws$rates <- rec$rates
  draws$reconciliation <- rec$info
  draws
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

# Core: arr [age, period, stratum, draw] of rates
reconcile_array <- function(arr, hs, additive, method, shrink, exposure, d, fixed = integer()) {
  S <- hs$S
  b <- hs$bottom
  is_agg <- !hs$is_bottom
  agg_rows <- which(is_agg)
  nb <- length(b)
  E <- expand_exposure(exposure, d, dimnames(arr))
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
        # diagonal weights as vectors, MinT as a matrix
        W <- switch(method,
                    ols = rep(1, d[3]),
                    structural = rowSums(S),
                    wls = pmax(apply(Yc, 2, stats::var), 1e-20),
                    mint = {
                      sh <- shrink_cov(Yc, shrink)
                      shrink_used <- c(shrink_used, sh$lambda)
                      sh$W
                    })
        WU <- if (is.matrix(W)) W %*% t(Ut)
        Yp <- project_rows(Yc, W, Ut, fixed, WU)[, b, drop = FALSE]
        # non-negative: bottom values below zero are held at zero and the rest
        # projected again, so that fixed nodes and aggregates are kept; draws
        # with the same held strata are projected together
        held <- matrix(FALSE, nrow(Yp), nb)
        active <- which(rowSums(Yp < 0) > 0)
        for (it in seq_len(nb)) {
          if (!length(active)) break
          held[active, ] <- held[active, , drop = FALSE] | Yp[active, , drop = FALSE] < 0
          active <- active[rowSums(held[active, , drop = FALSE]) < nb]
          key <- do.call(paste0, as.data.frame(held[active, , drop = FALSE] + 0L))
          for (k in unique(key)) {
            rr <- active[key == k]
            hb <- b[held[rr[1], ]]
            Yk <- Yc[rr, , drop = FALSE]
            Yk[, hb] <- 0
            Yp[rr, ] <- project_rows(Yk, W, Ut, union(fixed, hb), WU)[, b, drop = FALSE]
          }
          active <- active[rowSums(Yp[active, , drop = FALSE] < 0) > 0]
        }
        zero_share <- zero_share + sum(held)
        draw_zero <- draw_zero | rowSums(held) > 0
        Yp
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

# Projection of the rows of Y (draws x strata) onto the coherent subspace with
# weights W; nodes in `keep` get (almost) zero variance, so the projection
# leaves them unchanged
project_rows <- function(Y, W, Ut, keep = integer(), WU = NULL) {
  Utt <- t(Ut)
  if (!is.matrix(W)) {
    # diagonal weights given as a vector
    if (length(keep)) W[keep] <- 1e-10 * max(W)
    WU <- W * Utt
  } else {
    # W with the rows and columns of `keep` set to zero and a tiny diagonal,
    # updated from the full product W %*% t(Ut)
    WU <- WU %||% (W %*% Utt)
    if (length(keep)) {
      WU <- WU - W[, keep, drop = FALSE] %*% Utt[keep, , drop = FALSE]
      WU[keep, ] <- 1e-10 * max(diag(W)) * Utt[keep, , drop = FALSE]
    }
  }
  M <- Ut %*% WU
  rhs <- Ut %*% t(Y)
  lam <- tryCatch(solve(M, rhs), error = function(e) MASS_ginv(M) %*% rhs)
  Y - t(WU %*% lam)
}

# Reorder the predictive noise of each bottom stratum, cell by cell, along the
# training-year ranks used for the paths (keeps co-movement of observation noise)
noise_reorder <- function(deaths, D, E, tau) {
  d <- dim(deaths)
  y <- log(pmax(D, 0.5) / E)
  nt <- dim(y)[2]
  ch <- apply(y[, -1, , drop = FALSE] - y[, -nt, , drop = FALSE], c(2, 3), mean)
  out <- deaths
  for (s in seq_len(d[3])) {
    r <- rank(ch[tau, s], ties.method = "random")
    for (a in seq_len(d[1])) {
      for (t in seq_len(d[2])) out[a, t, s, ] <- sort(deaths[a, t, s, ])[r]
    }
  }
  out
}

# Empirical copula: reorder whole paths of each stratum so that their ranks follow
# the strata's average one-year changes of log rates in sampled training years
copula_reorder <- function(arr, D, E) {
  d <- dim(arr)
  y <- log(pmax(D, 0.5) / E)
  nt <- dim(y)[2]
  if (nt < 3) return(list(arr = arr, years = NULL))
  ch <- apply(y[, -1, , drop = FALSE] - y[, -nt, , drop = FALSE], c(2, 3), mean)
  tau <- sample.int(nt - 1, d[4], replace = TRUE)
  out <- arr
  for (s in seq_len(d[3])) {
    v <- apply(log(pmax(arr[, , s, , drop = FALSE], 1e-300)), 4, mean)
    if (anyNA(v)) next
    paths <- order(v)
    r <- rank(ch[tau, s], ties.method = "random")
    out[, , s, ] <- arr[, , s, paths[r], drop = FALSE]
  }
  list(arr = out, years = tau)
}

# Exposures [age, period, stratum] for the forecast: from an mt_data object
# or an array, matched to the forecast by the names of ages, periods and
# strata (dn: dimnames of the forecast rates); a single period is used for
# every forecast period
expand_exposure <- function(exposure, d, dn = NULL) {
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

# Shrinkage covariance towards the diagonal (Schaefer and Strimmer 2005, target D)
shrink_cov <- function(X, shrink = "auto") {
  n <- nrow(X)
  V <- stats::cov(X)
  s <- sqrt(pmax(diag(V), 1e-30))
  if (is.numeric(shrink)) {
    lambda <- shrink
  } else {
    Z <- sweep(sweep(X, 2, colMeans(X)), 2, s, "/")
    R <- crossprod(Z) / (n - 1)
    Z2 <- Z^2
    v_r <- (crossprod(Z2) - crossprod(Z)^2 / n) * n / (n - 1)^3
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
#' added with [mt_aggregate_strata()] before fitting, so that the base model
#' forecasts every node. Predictive draws (used by backtests) are coherent
#' as well: Poisson and extra-Poisson variation are drawn for the bottom
#' strata only and the aggregates are their sums.
#'
#' Whether coherence costs accuracy, and which method and base model work
#' best, depends on the data: compare them in a backtest at every level of
#' the hierarchy (`mt_score_table(..., hierarchy = )`), as in the vignette
#' "Coherent forecasts: totals and their parts".
#'
#' @param spec Base model specification (fitted to all strata).
#' @param hierarchy An [mt_hierarchy()].
#' @param method,shrink See [mt_reconcile()].
#' @param dependence `"auto"` (default: `"none"` for models fitted jointly
#'   across strata, whose draws already carry their dependence, such as
#'   [li_lee()] and [tucker3()]; `"empirical"` for models fitted per stratum),
#'   `"empirical"` or `"none"`. With `"empirical"` the same training years
#'   order the paths before reconciliation and the predictive noise of the
#'   bottom strata, so aggregates of predictive draws keep the observed
#'   co-movement.
#' @param fixed Strata whose forecasts are kept (see [mt_reconcile()]).
#'
#' @return A model specification.
#' @examples
#' x <- mt_example_regions(years = 2000:2019, regions = 3)
#' h <- mt_hierarchy(sex = "both", region = "all", quintile = "all")
#' fit <- mt_fit(reconciled(lee_carter(), h, method = "wls"), mt_aggregate_strata(x, h))
#' fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
#' @export
reconciled <- function(spec, hierarchy, method = c("wls", "mint", "structural", "ols", "bottom_up", "proportional"),
                          shrink = "auto", dependence = c("auto", "empirical", "none"), fixed = NULL) {
  check_class(spec, "mt_spec", "spec")
  method <- match.arg(method)
  dependence <- match.arg(dependence)
  if (dependence == "auto") dependence <- if (isTRUE(spec$joint)) "none" else "empirical"
  if (length(fixed) && method == "bottom_up") {
    mt_abort("`fixed` needs a projection method (\"wls\", \"mint\", \"structural\", \"ols\") or \"proportional\".",
             class = "argument")
  }
  if (method == "proportional" && !length(fixed)) {
    mt_abort("`method = \"proportional\"` scales the parts of fixed nodes: give `fixed`.", class = "argument")
  }
  sp <- new_mt_spec("reconciled", list(base = spec, hierarchy = as_hierarchy(hierarchy), method = method,
                                       shrink = shrink, dependence = dependence, fixed = fixed),
                    min_periods = spec$min_periods, label = sprintf("%s, reconciled (%s)", spec$label, method))
  sp
}

#' @export
mt_fit.mt_spec_reconciled <- function(spec, data, prep = NULL, cores = 1, ...) {
  p <- spec$params
  hs <- hierarchy_structure(data$strata, p$hierarchy)
  if (all(hs$is_bottom)) data <- mt_aggregate_strata(data, p$hierarchy)
  fit <- mt_fit(p$base, data, prep = prep, cores = cores, ...)
  hs <- hierarchy_structure(fit$strata, p$hierarchy)
  fit$reconcile <- list(hierarchy = p$hierarchy, method = p$method, shrink = p$shrink, dependence = p$dependence,
                        structure = hs, fixed = fixed_nodes(p$fixed, fit$strata, p$hierarchy, hs, p$method))
  wrap_fit(fit, spec, "mt_fit_reconciled")
}

#' @export
mt_forecast.mt_fit_reconciled <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                                          exposure = NULL, ...) {
  base_fit <- unwrap_fit(fit, "mt_fit_reconciled")
  rc <- fit$reconcile
  out <- with_seed(seed, {
    fc <- mt_forecast(base_fit, h = if (is.null(to)) h else to - max(fit$periods), n_draws = n_draws)
    d <- dim(fc$rates)
    E <- expand_exposure(exposure %||% fit$exposure_last, d, dimnames(fc$rates))
    arr <- fc$rates
    tau <- NULL
    if (identical(rc$dependence, "empirical")) {
      cr <- copula_reorder(arr, fit$raw_D %||% fit$train_D, fit$train_E)
      arr <- cr$arr
      tau <- cr$years
    }
    rec <- reconcile_array(arr, rc$structure, rc$hierarchy$additive, rc$method, rc$shrink, E, d,
                           fixed = rc$fixed %||% integer())
    rates <- rec$rates
    if (predictive) {
      hs <- rc$structure
      b <- hs$bottom
      Eb <- rec$exposure[, , b, drop = FALSE]
      mb <- rates[, , b, , drop = FALSE]
      if (isTRUE(base_fit$spec$smooth)) {
        sd <- array(fit$noise_sd[, b, drop = FALSE][, rep(seq_along(b), each = d[2])], c(d[1], d[2], length(b)))
        mb <- mb * exp(stats::rnorm(length(mb), 0, as.vector(sd)))
        mu <- mb * as.vector(Eb)
        deaths_b <- array(stats::rpois(length(mu), mu), dim(mu))
      } else {
        deaths_b <- mb * as.vector(Eb)
      }
      if (!is.null(tau)) deaths_b <- noise_reorder(deaths_b, (fit$raw_D %||% fit$train_D)[, , b, drop = FALSE],
                                                   fit$train_E[, , b, drop = FALSE], tau)
      # deaths add up over strata (and are divided by the summed exposures), or
      # rates add up (causes of death sharing one exposure)
      by_deaths <- rc$hierarchy$additive == "deaths"
      vb <- if (by_deaths) deaths_b else deaths_b / as.vector(Eb)
      full <- array(NA_real_, d)
      for (j in seq_len(d[4])) {
        for (t in seq_len(d[2])) full[, t, , j] <- matrix(vb[, t, , j], d[1]) %*% t(hs$S)
      }
      rates <- if (by_deaths) full / as.vector(rec$exposure) else full
    }
    list(fc = fc, rates = rates, info = rec$info)
  })
  res <- out$fc
  res$rates <- out$rates
  dimnames(res$rates) <- dimnames(out$fc$rates)
  res$model <- fit$spec$label
  res$predictive <- predictive
  res$reconciliation <- out$info
  res$hierarchy <- rc$hierarchy
  res
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
  E <- expand_exposure(exposure %||% draws$exposure_last, d, dimnames(draws$rates))
  st <- draws$strata
  X <- draws$rates * as.vector(E)
  if (!is.null(hierarchy)) {
    h <- as_hierarchy(hierarchy)
    hs0 <- hierarchy_structure(st, h)
    by <- setdiff(names(st), c("stratum", "label"))
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
  E <- expand_exposure(exposure %||% draws$exposure_last, d, dimnames(draws$rates))
  by_deaths <- identical(h$additive, "deaths")
  V <- if (by_deaths) draws$rates * as.vector(E) else draws$rates
  agg <- which(!hs$is_bottom)
  b <- hs$bottom
  rows <- lapply(agg, function(i) {
    w <- hs$S[i, ]
    part <- apply(sweep(V[, , b, , drop = FALSE], 3, w, "*"), c(1, 2, 4), sum)
    r <- part / V[, , i, ]
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
