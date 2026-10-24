#' Relative model of subpopulations and a reference population
#'
#' Joint model of subpopulations (for example deprivation quintiles) and the
#' population they belong to (for example England): the relative model of
#' Villegas and Haberman (2014, Section 2.4). The reference population has a
#' Lee-Carter model, with a cohort effect if `cohort = TRUE`; each
#' subpopulation's rates are the reference rates times an age pattern of
#' level differences and an index of improvement differences. Section and
#' equation numbers refer to that paper.
#'
#' @section Model (Section 2.4, eqs. 16-18):
#' The reference rates are `log mu'[x, t] = alpha'[x] + beta'[x] kappa'[t] +
#' gamma'[t - x]` (eq. 16), the age-period-cohort model of Renshaw and
#' Haberman (2006) with cohort effects `gamma'`; with `cohort = FALSE` the
#' cohort effects are left out and the reference has the Lee-Carter model
#' (the paper's "LC reference"). The rates of subpopulation `g` are
#' `log mu[x, t, g] = log mubar[x, t] + alpha[x, g] + beta[x] kappa[t, g]`
#' (eq. 17). `exp(alpha[x, g])` is the average ratio of the mortality of `g`
#' to the reference (level differences); `kappa[t, g]` is the index of the
#' improvements of `g` relative to the reference (a falling index: faster
#' improvement), with the age pattern `beta[x]` common to all
#' subpopulations. With `beta[x] > 0`, subpopulations keep their order as
#' long as their indexes keep it (Section 2.4).
#'
#' `mubar[x, t]` is the geometric mean of the reference rates over the
#' single ages of the subpopulations' age group (eq. 18). The reference can
#' cover more ages and years than the subpopulations (Section 2.4): the
#' paper used England and Wales at ages 10-99 in 1961-2009 as the reference
#' for English deprivation quintiles at ages 50-84 in 1981-2007.
#'
#' @section Reference:
#' `reference` is either
#'
#' * strata of the data: stratum labels or a named list of values of stratum
#'   variables, as in [saint()], for example `list(quintile = "England")`
#'   after [mt_aggregate_strata()]. The reference strata of a group of `by`
#'   are pooled (sum of deaths and exposures) into one reference population.
#'   A single reference stratum of a group is forecast by the reference
#'   model; all other strata of the group are subpopulations. Reference and
#'   subpopulations then share ages and years, and `mubar` is the reference
#'   rate of the same age group;
#' * an [mt_data] object with the reference population, which can cover more
#'   ages and years (Section 2.4). With `by`, its strata are matched to the
#'   groups of `by` by the values of these variables; the strata that match
#'   a group are pooled. Without `by`, all its strata are pooled into one
#'   reference. Every age group of the subpopulations must be a union of age
#'   groups of the reference, and every year of the subpopulations a year of
#'   the reference. All strata of the data are subpopulations. Preprocessing
#'   (`mt_fit(prep = )`) applies to the data, not to the reference.
#'
#' Package-specific (reference unknown): where the reference has age groups
#' wider than one year inside an age group of the subpopulations, `mubar` is
#' their geometric mean weighted by the widths (eq. 18 for single ages);
#' years of the reference after the last year of the data are not used, so
#' that backtests use no later data (the paper's reference ends two years
#' after its subpopulation data, Section 3).
#'
#' @section Estimation (Sections 2.5 and 3.1):
#' Deaths of the reference and of the subpopulations are independent Poisson
#' counts, and all parameters maximise the total log-likelihood (eq. 25), so
#' the subpopulations also inform the reference parameters. The paper uses
#' "straightforward extensions" of the Newton-Raphson procedure of Brouhns,
#' Denuit and Vermunt (2002). Here each iteration updates `alpha'`,
#' `kappa'`, `beta'`, `gamma'`, `alpha`, `kappa` and `beta` in turn by
#' one-dimensional Newton steps of the total log-likelihood, started from
#' the Poisson Lee-Carter fit of the reference and a singular value
#' decomposition of the subpopulations' log ratios to it. `alpha'` takes its
#' steps from the reference deaths only: `alpha[x, g]` absorbs any change of
#' `alpha'` in the subpopulation part, so its score there vanishes at the
#' maximum and the steps reach the same maximum (package-specific).
#' Package-specific (reference unknown): with cohort effects these
#' iterations converge slowly, as known for this model (Section 2.5), and
#' their fixed point is found with Anderson acceleration (Anderson 1965;
#' Walker and Ni 2011), which extrapolates from the last ten iterations.
#' The iterations stop when one iteration changes no fitted log rate by
#' more than `1e-9` (at most 5,000 iterations; a warning of class
#' `mt_warning_relative_model` otherwise). The fit reproduces
#' `StMoMo::rh(approxConst = TRUE)` for a reference without
#' subpopulations.
#'
#' Constraints: `sum(beta') = 1`, `sum(kappa') = 0`, `sum(gamma') = 0` over
#' the cohorts (eqs. 19-21), `sum(beta) = 1` and `sum(kappa[, g]) = 0` for
#' every `g` (eqs. 22-23). With `cohort = TRUE`, the cohort effects also
#' have no linear trend, `sum((c - cbar) gamma'[c]) = 0` (eq. 26), the
#' approximate constraint of Hunt and Villegas (2015). It is imposed after
#' every iteration by their approximately invariant transformation, as in
#' `StMoMo::rh(approxConst = TRUE)` (Villegas, Kaishev and Millossovich 2018),
#' the package of the first author: `alpha'`, `beta'`, `kappa'` and `gamma'`
#' change so that the cohort effects lose their linear trend and the fitted
#' rates stay the same where `kappa'` is linear. The transformation is
#' skipped where `kappa'` has no trend (package-specific safeguard).
#'
#' [mt_parameters()] reports, per part, the log-likelihood (with the
#' `log d!` terms of eqs. 24-25), the effective number of parameters and the
#' BIC of Section 3.1, `L - 0.5 nu log N` with `N` the number of cells
#' (higher is better). The subpopulation part has `nu = k m + k + n m - 1 -
#' m` parameters for `k` age groups, `n` years and `m` subpopulations,
#' which gives the 171 of Table 1; the Lee-Carter reference has `2 k' + n' -
#' 2` (227 in Table 1), the age-period-cohort reference `2 k' + n' + (k' +
#' n' - 1) - 4`, counting every cohort (363 for the dimensions of Table 1,
#' which reports 353).
#'
#' @section Forecasts (Sections 2.4 and 3.3):
#' Forecasts start from the fitted rates of the last year (eq. 32). The
#' reference index `kappa'` follows a random walk with drift (Section 2.4);
#' `trend_from` estimates it from the years from `trend_from` only, as the
#' paper did with the years after 1980, where the index was close to linear
#' (Section 3.3, following Haberman and Renshaw 2009). The subpopulation
#' indexes follow, by `deviation`:
#'
#' * `"mrwd"`: a multivariate random walk with drift with correlated normal
#'   innovations (eq. 4), the paper's model (Sections 2.4 and 3.3);
#' * `"rwd"`: independent random walks with drift (eq. 3). Section 2.1 notes
#'   that each subpopulation has the same forecast distribution under both
#'   processes, and that differences between subpopulations get wider
#'   intervals from independent walks; Section 3.1 compares both in a
#'   backtest (for the independent Lee-Carter model). With the same seed,
#'   the draws of the first subpopulation of a group are identical for both;
#' * `"constant"`: the indexes stay at their values of the last year, so
#'   all subpopulations improve as the reference (Scenario 2 of Section 3.4,
#'   eq. 33: level differences only).
#'
#' The drifts are the mean year-to-year changes and the innovation
#' (co)variances those of the changes around them, with divisor the number
#' of changes minus one (as `StMoMo::mrwd()`). Draws include the innovations
#' only: parameters, drifts and covariances are fixed, as in the paper,
#' whose intervals "ignore any provision for parameter uncertainty"
#' (Section 3.1). The reference index and the subpopulation indexes are
#' independent, separate time series models as in the paper. Each group of
#' `by` (for example each sex) is a separate model, as the paper fits women
#' and men separately. Predictive draws (`mt_forecast(predictive = TRUE)`)
#' add Poisson variation of the deaths given the rates.
#'
#' The drifting indexes let gaps between subpopulations widen or narrow
#' without limit; [li_lee()], [coherent_fdm()] and [saint()] keep the
#' deviations from a common trend stationary.
#'
#' Cohort effects: forecasts need no cohort beyond the fitted ones while
#' the horizon is at most the distance between the youngest age of the
#' subpopulations and the youngest age of the reference (the paper's case,
#' Section 2.4). Later cohorts follow an ARIMA(1,1,0) model with drift of
#' the fitted cohort effects, drawn with its innovations (package-specific:
#' the paper does not extrapolate cohorts; this is the default of
#' `StMoMo::forecast.fitStMoMo()`).
#'
#' @section Data:
#' The reference needs single-year consecutive ages without an open
#' interval for `cohort = TRUE` (the paper uses single ages, Section 3),
#' positive exposures, and deaths at every age (and in every cohort).
#' Subpopulations need deaths in every age group and every year.
#'
#' @param reference The reference population: stratum labels, a named list
#'   of values of stratum variables (for example `list(quintile =
#'   "England")`), or an [mt_data] object (see Details).
#' @param by Stratum variable(s) with a model each, for example `"sex"`.
#'   `NULL`: one model for all strata.
#' @param cohort Whether the reference model has cohort effects (eq. 16;
#'   the paper's "APC reference"). `FALSE`: Lee-Carter reference.
#' @param deviation Forecast of the subpopulation indexes: `"mrwd"`,
#'   `"rwd"` or `"constant"` (see Details).
#' @param trend_from `NULL` (all years), or the first year of the reference
#'   index used for its random walk with drift.
#' @return A model specification (joint over strata). The fit
#'   (`fit$fits$models`, one per group) holds the reference parameters
#'   (`ref`: `alpha`, `beta`, `kappa`, `gamma` and the time series
#'   parameters), the subpopulation parameters (`alpha`, `beta`, `kappa`,
#'   `drift`, `sigma`) and the log-likelihoods, effective numbers of
#'   parameters and BIC of both parts.
#' @references
#' Villegas AM, Haberman S (2014). On the modeling and forecasting of
#' socioeconomic mortality differentials: an application to deprivation and
#' mortality in England. *North American Actuarial Journal* 18(1): 168-193.
#' \doi{10.1080/10920277.2013.866034}
#'
#' Renshaw AE, Haberman S (2006). A cohort-based extension to the Lee-Carter
#' model for mortality reduction factors. *Insurance: Mathematics and
#' Economics* 38(3): 556-570. \doi{10.1016/j.insmatheco.2005.12.001}
#'
#' Hunt A, Villegas AM (2015). Robustness and convergence in the Lee-Carter
#' model with cohort effects. *Insurance: Mathematics and Economics* 64:
#' 186-202. \doi{10.1016/j.insmatheco.2015.05.004}
#'
#' Brouhns N, Denuit M, Vermunt JK (2002). A Poisson log-bilinear regression
#' approach to the construction of projected lifetables. *Insurance:
#' Mathematics and Economics* 31(3): 373-393.
#' \doi{10.1016/S0167-6687(02)00185-3}
#'
#' Haberman S, Renshaw A (2009). On age-period-cohort parametric mortality
#' rate projections. *Insurance: Mathematics and Economics* 45(2): 255-270.
#' \doi{10.1016/j.insmatheco.2009.07.006}
#'
#' Villegas AM, Kaishev VK, Millossovich P (2018). StMoMo: an R package for
#' stochastic mortality modeling. *Journal of Statistical Software* 84(3):
#' 1-38. \doi{10.18637/jss.v084.i03}
#'
#' Anderson DG (1965). Iterative procedures for nonlinear integral
#' equations. *Journal of the ACM* 12(4): 547-560. \doi{10.1145/321296.321305}
#'
#' Walker HF, Ni P (2011). Anderson acceleration for fixed-point iterations.
#' *SIAM Journal on Numerical Analysis* 49(4): 1715-1735.
#' \doi{10.1137/10078356X}
#' @seealso [saint()], another model of populations relative to a
#'   reference; [cae()], whose populations also keep drifting indexes.
#' @examples
#' x <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile"),
#'              age = "age_group", value = c(deaths = "all_causes", population = "population"))
#' x <- mt_as_exposure(x, type = "midyear")
#' x <- mt_group_levels(x, "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)),
#'                      name = "quintile")
#' x <- mt_subset(x, ages = c(50, 85))
#' xa <- mt_aggregate_strata(x, mt_hierarchy(quintile = list(total = "England")))
#' fit <- mt_fit(relative_model(reference = list(quintile = "England"), by = "sex"), xa)
#' p <- mt_parameters(fit)
#' p[parameter == "bic"]
#' fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
#' @export
relative_model <- function(reference, by = NULL, cohort = FALSE, deviation = c("mrwd", "rwd", "constant"),
                           trend_from = NULL) {
  deviation <- match.arg(deviation)
  check_flag(cohort, "cohort")
  external <- inherits(reference, "mt_data")
  ok_ref <- external || (is.character(reference) && length(reference) > 0 && !anyNA(reference)) ||
    (is.list(reference) && length(reference) > 0 && !is.null(names(reference)) && all(nzchar(names(reference))))
  if (!ok_ref) {
    mt_abort(paste("`reference` must name the reference strata (stratum labels, or a named list of values of",
                   "stratum variables, for example list(quintile = \"England\")) or be an mt_data object with",
                   "the reference population."), class = "argument")
  }
  if (!is.null(by) && (!is.character(by) || !length(by) || anyNA(by))) {
    mt_abort("`by` must be NULL or names of stratum variables, for example \"sex\".", class = "argument")
  }
  if (!is.null(trend_from)) check_number(trend_from, "trend_from", integerish = TRUE)
  new_mt_spec("relative", list(reference = reference, by = by, cohort = cohort, deviation = deviation,
                               trend_from = trend_from),
              joint = TRUE, min_periods = 5,
              label = paste0("Relative model (", if (cohort) "age-period-cohort" else "Lee-Carter", " reference",
                             if (!is.null(by)) sprintf(", by %s", paste(by, collapse = ", ")),
                             switch(deviation, mrwd = "", rwd = ", independent random walks",
                                    constant = ", constant deviations"),
                             if (!is.null(trend_from)) sprintf(", reference trend from %d", as.integer(trend_from)),
                             ")"))
}

#' @export
fit_model.mt_spec_relative <- function(spec, D, E, ages, periods, strata = NULL, hierarchy = NULL, ...) {
  p <- spec$params
  ns <- dim(D)[3]
  if (is.null(strata) || nrow(strata) != ns) {
    mt_abort("relative_model() needs the strata of the data.", class = "argument")
  }
  inputs <- if (inherits(p$reference, "mt_data")) {
    relative_external(p, strata, ages, periods)
  } else {
    relative_internal(p, D, E, ages, periods, strata)
  }
  models <- lapply(seq_along(inputs$groups), function(g) {
    inp <- inputs$ref[[g]]
    tgt <- which(inputs$group == g & inputs$target)
    if (p$cohort) relative_check_single(inp$lower, inp$upper)
    relative_check_deaths(inp$D, D[, , tgt, drop = FALSE], p$cohort, inp$labels, ages$labels, strata$label[tgt])
    m <- relative_fit_group(inp$D, inp$E, D[, , tgt, drop = FALSE], E[, , tgt, drop = FALSE], inp$W, inp$jj,
                            p$cohort)
    m$targets <- tgt
    m$ref <- c(m$ref, list(labels = inp$labels, lower = inp$lower, years = inp$periods))
    m$ref <- c(m$ref, relative_reference_series(m$ref, p$trend_from, p$cohort))
    m$sub <- c(m$sub, relative_index_series(m$sub$kappa))
    # the youngest reference age whose rates the forecasts need: with a
    # reference stratum in the forecasts, all reference ages
    used <- if (any(inputs$group == g & inputs$sole)) 1L else min(which(colSums(inp$W) > 0))
    m$ref$first_used <- used
    m
  })
  if (!all(vapply(models, `[[`, logical(1), "converged"))) {
    mt_warn("Relative model: the Newton-Raphson iterations did not converge within 5,000 iterations.",
            class = "relative_model")
  }
  fitted <- array(NA_real_, dim(D))
  for (g in seq_along(models)) {
    m <- models[[g]]
    lmb <- m$W %*% m$ref$eta[, m$jj, drop = FALSE]
    for (j in seq_along(m$targets)) {
      fitted[, , m$targets[j]] <- exp(lmb + m$sub$alpha[, j] + outer(m$sub$beta, m$sub$kappa[, j]))
    }
    for (s in which(inputs$group == g & inputs$sole)) fitted[, , s] <- exp(lmb)
  }
  list(models = models, group = inputs$group, groups = inputs$groups, target = inputs$target, sole = inputs$sole,
       deviation = p$deviation, cohort = p$cohort, fitted = fitted)
}

#' @export
simulate_model.mt_spec_relative <- function(spec, fit, h, n_draws, ...) {
  na <- nrow(fit$fitted)
  ns <- length(fit$group)
  out <- array(NA_real_, c(na, h, ns, n_draws))
  for (g in seq_along(fit$models)) {
    m <- fit$models[[g]]
    r <- m$ref
    na_r <- length(r$alpha)
    nt_r <- length(r$kappa)
    # random numbers: reference index, then cohort effects, then
    # subpopulation indexes (the same order for every `deviation`)
    kp <- r$kappa[nt_r] + outer(r$drift * seq_len(h), rep(1, n_draws)) +
      r$sd * relative_cumsum(matrix(stats::rnorm(h * n_draws), h, n_draws))
    gam <- NULL
    if (fit$cohort) gam <- relative_cohort_paths(r, h, n_draws)
    kt <- relative_target_paths(m$sub, h, n_draws, fit$deviation)
    sole <- which(fit$group == g & fit$sole)
    for (j in seq_len(n_draws)) {
      er <- r$alpha + outer(r$beta, kp[, j])
      if (!is.null(gam)) {
        # cohort index of reference age i in forecast year k: n_c + 1 + k - i
        ci <- outer(-seq_len(na_r), seq_len(h), `+`) + length(r$gamma) + 1L
        er <- er + matrix(c(r$gamma, gam[, j])[ci], na_r, h)
      }
      lmb <- m$W %*% er
      for (i in seq_along(m$targets)) {
        out[, , m$targets[i], j] <- exp(lmb + m$sub$alpha[, i] + m$sub$beta %o% kt[i, , j])
      }
      for (s in sole) out[, , s, j] <- exp(lmb)
    }
  }
  out
}

# Inputs ------------------------------------------------------------------------

# Reference strata of the data: groups of `by`, pooled reference strata of
# each group (as saint()), on the grid of the data
relative_internal <- function(p, D, E, ages, periods, strata) {
  st <- saint_strata(list(reference = p$reference, by = p$by, ecm = NULL), strata, dim(D)[3])
  na <- dim(D)[1]
  nt <- dim(D)[2]
  ref <- lapply(seq_along(st$groups), function(g) {
    use <- st$group == g & st$reference
    list(D = matrix(apply(D[, , use, drop = FALSE], 1:2, sum), na, nt),
         E = matrix(apply(E[, , use, drop = FALSE], 1:2, sum), na, nt),
         W = diag(na), jj = seq_len(nt), labels = ages$labels, lower = ages$lower, upper = ages$upper,
         periods = periods)
  })
  list(groups = st$groups, group = st$group, target = !st$sole, sole = st$sole, ref = ref)
}

# Reference population of an mt_data object: strata matched to the groups of
# `by` and pooled, years up to the last year of the data, and the weights of
# its age groups in each age group of the data (eq. 18)
relative_external <- function(p, strata, ages, periods) {
  rd <- p$reference
  md <- model_inputs(rd, list(min_periods = 2, label = "The reference population of relative_model()"))
  rp <- md$periods
  keep <- rp <= max(periods)
  if (!all(periods %in% rp[keep])) {
    mt_abort(sprintf(paste("The reference population must cover every year of the data (Villegas and Haberman 2014,",
                           "Section 2.4); the data have years %d-%d, the reference %d-%d."), min(periods),
                     max(periods), min(rp), max(rp)), class = "model_data", call = NULL)
  }
  rs <- rd$strata
  ns <- nrow(strata)
  if (is.null(p$by)) {
    group <- rep(1L, ns)
    groups <- "all"
    rgroup <- rep(1L, nrow(rs))
  } else {
    bad <- setdiff(p$by, setdiff(names(strata), c("stratum", "label")))
    bad_r <- setdiff(p$by, setdiff(names(rs), c("stratum", "label")))
    if (length(bad) || length(bad_r)) {
      mt_abort(sprintf("`by` names stratum variable(s) %s that the %s does not have.",
                       paste(unique(c(bad, bad_r)), collapse = ", "), if (length(bad)) "data" else "reference"),
               class = "argument", call = NULL)
    }
    key_of <- function(s) do.call(paste, c(lapply(p$by, function(v) as.character(s[[v]])), sep = "."))
    key <- key_of(strata)
    groups <- unique(key)
    group <- match(key, groups)
    rgroup <- match(key_of(rs), groups)
    miss <- groups[!seq_along(groups) %in% rgroup]
    if (length(miss)) {
      mt_abort(sprintf("The reference population has no stratum for group(s) %s of `by`.",
                       paste(miss, collapse = ", ")), class = "argument", call = NULL)
    }
  }
  W <- relative_age_weights(ages, md$age)
  jj <- match(periods, rp[keep])
  ref <- lapply(seq_along(groups), function(g) {
    use <- which(rgroup == g)
    list(D = matrix(apply(md$D[, keep, use, drop = FALSE], 1:2, sum), length(md$age)),
         E = matrix(apply(md$E[, keep, use, drop = FALSE], 1:2, sum), length(md$age)),
         W = W, jj = jj, labels = md$age$labels, lower = md$age$lower, upper = md$age$upper,
         periods = rp[keep])
  })
  list(groups = groups, group = group, target = rep(TRUE, ns), sole = rep(FALSE, ns), ref = ref)
}

# Weights of the reference age groups in each age group of the data: the
# widths of the reference groups inside it over its width (1 / n for single
# ages, eq. 18); an open age group of the data must be an open age group of
# the reference
relative_age_weights <- function(ages, rages) {
  lo <- ages$lower
  up <- ages$upper
  rlo <- rages$lower
  rup <- rages$upper
  W <- matrix(0, length(lo), length(rlo))
  for (x in seq_along(lo)) {
    if (!is.finite(up[x])) {
      i <- which(abs(rlo - lo[x]) < 1e-9 & !is.finite(rup))
    } else {
      i <- which(rlo > lo[x] - 1e-9 & rup < up[x] + 1e-9)
    }
    ok <- length(i) > 0 && (!is.finite(up[x]) ||
                              (abs(min(rlo[i]) - lo[x]) < 1e-9 && abs(max(rup[i]) - up[x]) < 1e-9 &&
                                 abs(sum(rup[i] - rlo[i]) - (up[x] - lo[x])) < 1e-9))
    if (!ok) {
      mt_abort(sprintf(paste("Age group %s of the data is not a union of age groups of the reference population",
                             "(ages %s to %s); restrict the ages of the data with mt_subset(ages = )."),
                       ages$labels[x], rages$labels[1], rages$labels[length(rlo)]), class = "model_data",
               call = NULL)
    }
    W[x, i] <- if (is.finite(up[x])) (rup[i] - rlo[i]) / (up[x] - lo[x]) else 1
  }
  W
}

# Cohort effects need consecutive single ages without an open interval
relative_check_single <- function(lower, upper) {
  if (any(!is.finite(upper)) || any(abs(upper - lower - 1) > 1e-9) ||
      (length(lower) > 1 && any(abs(diff(lower) - 1) > 1e-9))) {
    mt_abort(paste("relative_model(cohort = TRUE) needs a reference with consecutive single-year ages and no open",
                   "interval: the cohort effects follow the year of birth (Villegas and Haberman 2014 use single",
                   "ages). Drop the open interval with mt_subset(ages = ), or use a reference in single years",
                   "(an mt_data object as `reference`)."), class = "model_data", call = NULL)
  }
  invisible()
}

# Deaths in every age (and cohort) of the reference and in every age group and
# year of each subpopulation: otherwise a parameter has no finite estimate
relative_check_deaths <- function(Dr, D, cohort, rlabels, labels, targets) {
  if (any(rowSums(Dr) <= 0)) {
    mt_abort(sprintf("The reference population has no deaths at age(s) %s; restrict its ages.",
                     paste(utils::head(rlabels[rowSums(Dr) <= 0], 5), collapse = ", ")), class = "model_data",
             call = NULL)
  }
  if (cohort) {
    cid <- col(Dr) - row(Dr) + nrow(Dr)
    if (any(rowsum(as.vector(Dr), as.vector(cid)) <= 0)) {
      mt_abort("A cohort of the reference population has no deaths; restrict its ages or years.",
               class = "model_data", call = NULL)
    }
  }
  if (dim(D)[3] == 0) return(invisible())
  by_age <- apply(D, c(1, 3), sum)
  by_year <- apply(D, c(2, 3), sum)
  if (any(by_age <= 0) || any(by_year <= 0)) {
    s <- which(colSums(by_age <= 0) > 0 | colSums(by_year <= 0) > 0)[1]
    mt_abort(sprintf(paste("Subpopulation %s has an age group or a year without deaths; use wider age groups",
                           "(mt_regroup()) or fewer ages (mt_subset(ages = ))."), targets[s]),
             class = "model_data", call = NULL)
  }
  invisible()
}

# Estimation ----------------------------------------------------------------

# Joint Poisson maximum likelihood of one group (Villegas and Haberman 2014,
# eqs. 16-18 and 25): reference deaths Dr, exposures Er [ref age, ref year];
# subpopulations D, E [age, year, subpopulation]; W [age, ref age] the
# weights of eq. 18; jj the reference years of the data's years. One
# iteration (relative_step()): Newton steps of alpha' (reference deaths),
# kappa', beta' and gamma' (both parts), alpha, kappa and beta
# (subpopulations), then the constraints (19)-(23) and, with cohort effects,
# the approximately invariant transformation of Hunt and Villegas (2015) for
# (26), as in StMoMo:::fit.rh() with approxConst = TRUE. The fixed point of
# these iterations is found with Anderson acceleration (Anderson 1965;
# Walker and Ni 2011, Algorithm AA with the last `memory` differences); an
# extrapolation with non-finite rates restarts it. The iterations stop when
# one iteration changes no linear predictor by more than `tol`.
relative_fit_group <- function(Dr, Er, D, E, W, jj, cohort, tol = 1e-9, maxit = 5000, accelerate = TRUE,
                               memory = 10) {
  na_r <- nrow(Dr)
  nt_r <- ncol(Dr)
  na <- dim(D)[1]
  nt <- dim(D)[2]
  ns <- dim(D)[3]
  nc <- na_r + nt_r - 1L
  pr <- list(Dr = Dr, Er = Er, D = D, E = E, W = W, W2 = W^2, jj = jj, cohort = cohort, na_r = na_r, nt_r = nt_r,
             na = na, nt = nt, ns = ns, nc = nc, cid = col(Dr) - row(Dr) + na_r,
             # cell (i, j) of a reference matrix in row (cohort) j - i + na_r and
             # column i of a matrix [cohort, age], for sums over cohorts
             cpos = as.vector(col(Dr) - row(Dr) + na_r + (row(Dr) - 1) * nc),
             tbar = seq_len(nt_r) - (nt_r + 1) / 2, xbar = seq_len(na_r) - (na_r + 1) / 2,
             cbar = seq_len(nc) - (nc + 1) / 2)
  len <- c(a_r = na_r, b_r = na_r, k_r = nt_r, g_r = if (cohort) nc else 0L, alpha = na * ns, beta = na,
           kappa = nt * ns)
  pr$pos <- split(seq_len(sum(len)), rep(factor(names(len), names(len)), len))
  x <- relative_pack(relative_start(pr))
  lp_x <- relative_lp(x, pr)
  dF <- dG <- f_old <- g_old <- x_safe <- NULL
  converged <- FALSE
  for (it in seq_len(maxit)) {
    g <- relative_step(x, pr)
    lp_g <- relative_lp(g, pr)
    if (!all(is.finite(lp_g))) {
      # an extrapolation that leads to non-finite rates: restart from the
      # iteration it replaced
      if (is.null(x_safe)) mt_abort("The relative model diverged.", class = "model_fit", call = NULL)
      x <- x_safe
      lp_x <- relative_lp(x, pr)
      dF <- dG <- f_old <- g_old <- x_safe <- NULL
      next
    }
    if (max(abs(lp_g - lp_x)) < tol) {
      x <- g
      converged <- TRUE
      break
    }
    # Anderson acceleration with the last `memory` differences of the steps
    f <- g - x
    if (accelerate && !is.null(f_old)) {
      dF <- cbind(dF, f - f_old)
      dG <- cbind(dG, g - g_old)
      if (ncol(dF) > memory) {
        dF <- dF[, -1, drop = FALSE]
        dG <- dG[, -1, drop = FALSE]
      }
    }
    f_old <- f
    g_old <- g
    x <- g
    lp_x <- lp_g
    x_safe <- NULL
    if (!is.null(dF)) {
      gam <- qr.coef(qr(dF), f)
      gam[is.na(gam)] <- 0
      cand <- as.vector(g - dG %*% gam)
      lp_c <- relative_lp(cand, pr)
      if (all(is.finite(lp_c))) {
        x_safe <- g
        x <- cand
        lp_x <- lp_c
      } else {
        dF <- dG <- NULL
      }
    }
  }
  v0 <- x
  th <- relative_unpack(v0, pr)
  # exact constraints at the end (rates unchanged): the transformation of the
  # last iteration changes sum(beta') and sum(kappa')
  th <- relative_normalise(th, pr, transform = FALSE)
  er <- relative_eta_ref(th, pr)
  ll_r <- sum(Dr * (er + log(Er)) - Er * exp(er) - lgamma(Dr + 1))
  ll_s <- if (ns > 0) {
    es <- relative_eta_sub(er, th, pr)
    sum(D * (es + log(E)) - E * exp(es) - lgamma(D + 1))
  } else {
    0
  }
  npar <- relative_npar(na_r, nt_r, cohort, na, nt, ns)
  npar_r <- npar[["ref"]]
  npar_s <- npar[["sub"]]
  n_r <- length(Dr)
  n_s <- length(D)
  list(ref = list(alpha = th$a_r, beta = th$b_r, kappa = th$k_r, gamma = if (cohort) th$g_r, eta = er,
                  fitted = exp(er)),
       sub = list(alpha = th$alpha, beta = th$beta, kappa = th$kappa),
       W = W, jj = jj, loglik_ref = ll_r, loglik_sub = ll_s, npar_ref = npar_r, npar_sub = npar_s, n_ref = n_r,
       n_sub = n_s, bic_ref = ll_r - 0.5 * npar_r * log(n_r),
       bic_sub = if (ns > 0) ll_s - 0.5 * npar_s * log(n_s) else NA_real_, iterations = it, converged = converged)
}

# Effective numbers of parameters (Section 3.1, footnote 3): the reference
# part less the constraints (19)-(21) and (26), the subpopulation part less
# (22) and the m constraints (23)
relative_npar <- function(na_r, nt_r, cohort, na, nt, ns) {
  c(ref = 2 * na_r + nt_r - 2 + if (cohort) na_r + nt_r - 1 - 2 else 0, sub = na * ns + na + nt * ns - 1 - ns)
}

# Parameters as one vector and back
relative_pack <- function(th) c(th$a_r, th$b_r, th$k_r, th$g_r, th$alpha, th$beta, th$kappa)

relative_unpack <- function(v, pr) {
  g <- function(nm) v[pr$pos[[nm]]]
  list(a_r = g("a_r"), b_r = g("b_r"), k_r = g("k_r"), g_r = if (pr$cohort) g("g_r"),
       alpha = matrix(g("alpha"), pr$na, pr$ns), beta = g("beta"), kappa = matrix(g("kappa"), pr$nt, pr$ns))
}

# Log rates of the reference [ref age, ref year] and of the subpopulations
# [age, year, subpopulation] (eqs. 16-18)
relative_eta_ref <- function(th, pr) {
  e <- th$a_r + outer(th$b_r, th$k_r)
  if (pr$cohort) e <- e + th$g_r[pr$cid]
  e
}

relative_eta_sub <- function(er, th, pr) {
  lmb <- pr$W %*% er[, pr$jj, drop = FALSE]
  out <- array(0, c(pr$na, pr$nt, pr$ns))
  for (s in seq_len(pr$ns)) out[, , s] <- lmb + th$alpha[, s] + outer(th$beta, th$kappa[, s])
  out
}

# Linear predictors of all cells (convergence criterion)
relative_lp <- function(v, pr) {
  th <- relative_unpack(v, pr)
  er <- relative_eta_ref(th, pr)
  c(er, if (pr$ns > 0) relative_eta_sub(er, th, pr))
}

# Starting values: the Poisson Lee-Carter fit of the reference (Brouhns et
# al. 2002), no cohort effects, and for the subpopulations the mean log
# ratio to the reference and the first singular vectors of the centred log
# ratios (sum(beta) = 1)
relative_start <- function(pr) {
  st <- tryCatch(lc_poisson(pr$Dr, pr$Er), mt_error_model_fit = function(e) lc_svd(pr$Dr, pr$Er))
  th <- list(a_r = st$a, b_r = st$b, k_r = st$k, g_r = if (pr$cohort) numeric(pr$nc),
             alpha = matrix(0, pr$na, pr$ns), beta = rep(1 / pr$na, pr$na), kappa = matrix(0, pr$nt, pr$ns))
  if (pr$ns > 0) {
    er <- relative_eta_ref(th, pr)
    y <- log_rates(pr$D, pr$E) - as.vector(pr$W %*% er[, pr$jj, drop = FALSE])
    th$alpha <- matrix(apply(y, c(1, 3), mean), pr$na, pr$ns)
    sv <- svd(matrix(sweep(y, c(1, 3), th$alpha), pr$na), nu = 1, nv = 1)
    if (abs(sum(sv$u)) > 1e-6) {
      th$beta <- sv$u[, 1] / sum(sv$u)
      th$kappa <- matrix(sv$d[1] * sv$v[, 1] * sum(sv$u), pr$nt, pr$ns)
    }
  }
  th
}

# One iteration of one-dimensional Newton steps of the total log-likelihood,
# followed by the constraints
relative_step <- function(v, pr) {
  th <- relative_unpack(v, pr)
  Dr <- pr$Dr
  Er <- pr$Er
  D <- pr$D
  E <- pr$E
  W <- pr$W
  jj <- pr$jj
  na <- pr$na
  nt <- pr$nt
  ns <- pr$ns
  # deaths of the subpopulations minus fitted deaths, and fitted deaths,
  # summed over subpopulations [age, year]
  sub_sums <- function(er) {
    if (ns == 0) return(list(R = matrix(0, na, nt), H = matrix(0, na, nt)))
    Dh <- E * exp(relative_eta_sub(er, th, pr))
    list(R = rowSums(D - Dh, dims = 2), H = rowSums(Dh, dims = 2))
  }
  # alpha': reference deaths (alpha[x, g] absorbs alpha' in the other part,
  # whose score of alpha' vanishes where that of alpha does)
  Dh_r <- Er * exp(relative_eta_ref(th, pr))
  th$a_r <- th$a_r + rowSums(Dr - Dh_r) / rowSums(Dh_r)
  # kappa'
  er <- relative_eta_ref(th, pr)
  Dh_r <- Er * exp(er)
  ss <- sub_sums(er)
  wb <- as.vector(W %*% th$b_r)
  sk <- colSums((Dr - Dh_r) * th$b_r)
  ik <- colSums(Dh_r * th$b_r^2)
  sk[jj] <- sk[jj] + colSums(ss$R * wb)
  ik[jj] <- ik[jj] + colSums(ss$H * wb^2)
  th$k_r <- th$k_r + sk / ik
  # beta'
  er <- relative_eta_ref(th, pr)
  Dh_r <- Er * exp(er)
  ss <- sub_sums(er)
  kk <- matrix(th$k_r, pr$na_r, pr$nt_r, byrow = TRUE)
  sb <- rowSums((Dr - Dh_r) * kk) + as.vector(crossprod(W, ss$R) %*% th$k_r[jj])
  ib <- rowSums(Dh_r * kk^2) + as.vector(crossprod(pr$W2, ss$H) %*% th$k_r[jj]^2)
  if (max(abs(th$k_r)) > 1e-10) th$b_r <- th$b_r + sb / ib
  # gamma'
  if (pr$cohort) {
    er <- relative_eta_ref(th, pr)
    Dh_r <- Er * exp(er)
    ss <- sub_sums(er)
    Rr <- Dr - Dh_r
    Hr <- Dh_r
    Rr[, jj] <- Rr[, jj] + crossprod(W, ss$R)
    Hr[, jj] <- Hr[, jj] + crossprod(pr$W2, ss$H)
    th$g_r <- th$g_r + relative_by_cohort(Rr, pr) / relative_by_cohort(Hr, pr)
  }
  # alpha, kappa, beta of the subpopulations
  if (ns > 0) {
    er <- relative_eta_ref(th, pr)
    Dh <- E * exp(relative_eta_sub(er, th, pr))
    th$alpha <- th$alpha + matrix(apply(D - Dh, c(1, 3), sum) / apply(Dh, c(1, 3), sum), na, ns)
    Dh <- E * exp(relative_eta_sub(er, th, pr))
    for (s in seq_len(ns)) {
      R <- matrix(D[, , s] - Dh[, , s], na, nt)
      th$kappa[, s] <- th$kappa[, s] + colSums(R * th$beta) / colSums(matrix(Dh[, , s], na, nt) * th$beta^2)
    }
    Dh <- E * exp(relative_eta_sub(er, th, pr))
    num <- 0
    den <- 0
    for (s in seq_len(ns)) {
      ks <- matrix(th$kappa[, s], na, nt, byrow = TRUE)
      num <- num + rowSums(matrix(D[, , s] - Dh[, , s], na, nt) * ks)
      den <- den + rowSums(matrix(Dh[, , s], na, nt) * ks^2)
    }
    # beta is not identified while every kappa is zero (package-specific
    # safeguard: no step then)
    if (max(abs(th$kappa)) > 1e-10) th$beta <- th$beta + num / den
  }
  relative_pack(relative_normalise(th, pr, transform = TRUE))
}

# Sums of a reference matrix [ref age, ref year] over the cells of each cohort
relative_by_cohort <- function(M, pr) {
  Z <- numeric(pr$nc * pr$na_r)
  Z[pr$cpos] <- M
  rowSums(matrix(Z, pr$nc))
}

# Constraints (19)-(23), which leave the rates unchanged, and with
# `transform` the approximately invariant transformation of Hunt and Villegas
# (2015) for (26), as in StMoMo; it is skipped where kappa' has no trend
relative_normalise <- function(th, pr, transform) {
  if (pr$ns > 0) {
    cs <- colMeans(th$kappa)
    th$alpha <- th$alpha + outer(th$beta, cs)
    th$kappa <- sweep(th$kappa, 2, cs)
    sb <- sum(th$beta)
    th$beta <- th$beta / sb
    th$kappa <- th$kappa * sb
  }
  c1 <- mean(th$k_r)
  th$a_r <- th$a_r + th$b_r * c1
  th$k_r <- th$k_r - c1
  s2 <- sum(th$b_r)
  th$b_r <- th$b_r / s2
  th$k_r <- th$k_r * s2
  if (pr$cohort) {
    c3 <- mean(th$g_r)
    th$a_r <- th$a_r + c3
    th$g_r <- th$g_r - c3
    if (transform) {
      K <- sum(pr$tbar * th$k_r) / sum(pr$tbar^2)
      e <- -sum(pr$cbar * th$g_r) / sum(pr$cbar^2)
      if (abs(K) > 1e-10 && abs(K - e) > 1e-10) {
        th$a_r <- th$a_r + e * pr$xbar
        th$b_r <- K / (K - e) * th$b_r - e / (K - e)
        th$k_r <- (K - e) / K * th$k_r
        th$g_r <- th$g_r + e * pr$cbar
      }
    }
  }
  th
}

# Time series ---------------------------------------------------------------

# Random walk with drift of the reference index (from `trend_from`), and the
# ARIMA(1,1,0) model with drift of the cohort effects (package-specific, the
# default of StMoMo), used for cohorts after the last fitted one
relative_reference_series <- function(ref, trend_from, cohort) {
  years <- ref$years
  use <- if (is.null(trend_from)) rep(TRUE, length(years)) else years >= trend_from
  if (sum(use) < 3) {
    mt_abort(sprintf("`trend_from` (%d) leaves fewer than three years of the reference index (years %d-%d).",
                     as.integer(trend_from), min(years), max(years)), class = "argument", call = NULL)
  }
  dk <- diff(ref$kappa[use])
  drift <- mean(dk)
  out <- list(drift = drift, sd = sqrt(sum((dk - drift)^2) / max(length(dk) - 1, 1)),
              trend_years = range(years[use]))
  if (cohort) {
    dg <- diff(ref$gamma)
    a <- tryCatch(stats::arima(dg, order = c(1, 0, 0), include.mean = TRUE), error = function(e) NULL)
    out$cohort_arima <- if (!is.null(a)) {
      list(phi = unname(a$coef[["ar1"]]), mu = unname(a$coef[["intercept"]]), sd = sqrt(a$sigma2),
           last = dg[length(dg)])
    }
  }
  out
}

# Multivariate random walk with drift of the subpopulation indexes (eq. 4):
# drift and covariance of the year-to-year changes (divisor: number of
# changes minus one)
relative_index_series <- function(kappa) {
  ns <- ncol(kappa)
  if (ns == 0) return(list(drift = numeric(), sigma = matrix(0, 0, 0), root = matrix(0, 0, 0), last = numeric()))
  dk <- kappa[-1, , drop = FALSE] - kappa[-nrow(kappa), , drop = FALSE]
  drift <- colMeans(dk)
  R <- sweep(dk, 2, drift)
  sigma <- crossprod(R) / max(nrow(dk) - 1, 1)
  list(drift = drift, sigma = sigma, root = mvn_root(sigma), last = kappa[nrow(kappa), ])
}

# Cumulative sums of the columns of a matrix
relative_cumsum <- function(Z) {
  matrix(apply(Z, 2, cumsum), nrow(Z))
}

# Paths of the subpopulation indexes [subpopulation, horizon, draw]: "mrwd"
# correlated normal innovations (lower triangular root, so that the first
# subpopulation has the same draws as with "rwd"), "rwd" independent ones,
# "constant" the last values
relative_target_paths <- function(sub, h, n_draws, deviation) {
  ns <- length(sub$drift)
  if (ns == 0) return(array(0, c(0, h, n_draws)))
  out <- array(sub$last, c(ns, h, n_draws))
  if (deviation == "constant") return(out)
  root <- if (deviation == "mrwd") sub$root else diag(sqrt(diag(sub$sigma)), ns)
  U <- 1 * upper.tri(diag(h), diag = TRUE)
  for (j in seq_len(n_draws)) {
    eps <- root %*% matrix(stats::rnorm(ns * h), ns, h)
    out[, , j] <- sub$last + (eps + sub$drift) %*% U
  }
  out
}

# Cohort effects of the h cohorts after the last fitted one [cohort, draw]:
# ARIMA(1,1,0) with drift, drawn only where the forecasts need them
relative_cohort_paths <- function(r, h, n_draws) {
  if (h < r$first_used) return(matrix(0, h, n_draws))
  a <- r$cohort_arima
  if (is.null(a)) {
    mt_abort(paste("The forecast needs cohorts after the last fitted one, and the ARIMA(1,1,0) model of the cohort",
                   "effects could not be estimated; shorten the horizon."), class = "model_fit", call = NULL)
  }
  out <- matrix(0, h, n_draws)
  last <- r$gamma[length(r$gamma)]
  for (j in seq_len(n_draws)) {
    d <- a$last
    g <- last
    for (k in seq_len(h)) {
      d <- a$mu + a$phi * (d - a$mu) + a$sd * stats::rnorm(1)
      g <- g + d
      out[k, j] <- g
    }
  }
  out
}

# Inspection ----------------------------------------------------------------

relative_parameters <- function(fit) {
  f <- fit$fits
  lab <- fit$strata$label
  ages <- fit$age$labels
  yrs <- as.character(fit$periods)
  row <- function(stratum, parameter, index, value) {
    data.table(stratum = stratum, parameter = parameter, index = as.character(index), value = as.numeric(value))
  }
  rbindlist(lapply(seq_along(f$models), function(g) {
    m <- f$models[[g]]
    r <- m$ref
    ryrs <- as.character(r$years)
    rl <- paste("reference", f$groups[g])
    sl <- paste("subpopulations", f$groups[g])
    ref <- rbind(row(rl, "alpha", r$labels, r$alpha), row(rl, "beta", r$labels, r$beta),
                 row(rl, "kappa", ryrs, r$kappa),
                 if (!is.null(r$gamma)) {
                   row(rl, "gamma", as.numeric(ryrs[1]) - r$lower[length(r$lower)] + seq_along(r$gamma) - 1,
                       r$gamma)
                 },
                 row(rl, c("drift", "sd"), "", c(r$drift, r$sd)),
                 if (!is.null(r$cohort_arima)) {
                   row(rl, c("ar_gamma", "drift_gamma", "sd_gamma"), "",
                       c(r$cohort_arima$phi, r$cohort_arima$mu, r$cohort_arima$sd))
                 },
                 row(rl, c("loglik", "npar", "nobs", "bic"), "", c(m$loglik_ref, m$npar_ref, m$n_ref, m$bic_ref)))
    tl <- lab[m$targets]
    sub <- if (length(m$targets)) {
      s <- m$sub
      sdv <- sqrt(diag(s$sigma))
      cr <- s$sigma / outer(sdv, sdv)
      rbind(row(sl, "beta", ages, s$beta),
            row(sl, c("loglik", "npar", "nobs", "bic"), "", c(m$loglik_sub, m$npar_sub, m$n_sub, m$bic_sub)),
            rbindlist(lapply(seq_along(m$targets), function(i) {
              rbind(row(tl[i], "alpha", ages, s$alpha[, i]), row(tl[i], "kappa", yrs, s$kappa[, i]),
                    row(tl[i], c("drift", "sd"), "", c(s$drift[i], sdv[i])), row(tl[i], "cor", tl, cr[i, ]))
            })))
    }
    rbind(ref, sub)
  }))
}

relative_assumptions <- function(spec) {
  p <- spec$params
  paste(
    "Relative model (Villegas and Haberman 2014): log m(x,t) of the reference =",
    if (isTRUE(p$cohort)) "alpha'(x) + beta'(x) kappa'(t) + gamma'(t-x) (Renshaw and Haberman 2006, cohort effects" else
      "alpha'(x) + beta'(x) kappa'(t) (Lee-Carter);",
    if (isTRUE(p$cohort)) "without linear trend as in Hunt and Villegas 2015);",
    "log m(x,t,g) of each subpopulation = log of the geometric mean of the reference rates over the age group",
    "+ alpha(x,g) + beta(x) kappa(t,g). All parameters by joint Poisson maximum likelihood of the reference and",
    "the subpopulations. kappa' follows a random walk with drift",
    if (!is.null(p$trend_from)) sprintf("estimated from %d;", as.integer(p$trend_from)) else ";",
    switch(p$deviation, mrwd = "kappa(t,g) a multivariate random walk with drift (correlated innovations),",
           rwd = "kappa(t,g) independent random walks with drift,",
           constant = "kappa(t,g) stays at its last value (level differences only),"),
    if (p$deviation != "constant") "so gaps between subpopulations can widen or narrow without limit." else
      "so subpopulations improve as the reference.",
    if (isTRUE(p$cohort)) "Cohorts after the last fitted one follow an ARIMA(1,1,0) model with drift.",
    "Draws include the innovations; parameters, drifts and covariances are fixed. Predictive draws add Poisson",
    "variation of the deaths given the rates. Jump-off: fitted rates.")
}
