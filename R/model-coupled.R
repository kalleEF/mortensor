#' Give separate forecasts the dependence of a joint model
#'
#' `coupled()` wraps a model (the base model) and a joint model (the
#' template): its forecasts keep the margins of the base model and take the
#' dependence between ages, years and strata from the template, by ensemble
#' copula coupling (ECC; Schefzik, Thorarinsdottir and Gneiting 2013).
#' `mt_couple()` does the same with two sets of forecast draws. Status:
#' published (the ECC-R variant), with the package-specific choices listed
#' below; section and equation numbers refer to the paper.
#'
#' **Method (Sections 4.1 and 4.2).** Both models are fitted to the same data
#' and forecast with the same number of draws `M`. A margin is one cell of
#' the forecast: an age, a forecast year and a stratum. For each margin, the
#' `M` values of the base model are sorted, and draw `m` receives the value
#' whose rank among them equals the rank of draw `m` of the template among
#' the template's values (eq. 4.3). Each margin keeps exactly the values of
#' the base model. The draws take the multivariate rank structure (the
#' empirical copula) of the template, and so its Spearman rank correlation
#' between any two margins. The margins of the paper are weather variables
#' at locations and lead times; here they are the cells of the grid, so the
#' template sets the dependence between strata, between ages and between
#' years. The template's draws must be exchangeable (Section 4); draws of a
#' fitted model are.
#'
#' **Variant (Section 4.3).** The base model's draws are a random sample
#' from its forecast distribution of each margin, so the method is ECC-R,
#' with the sample taken by the base model's simulation. The paper
#' recommends ECC-Q, which uses the equidistant quantiles `F^-1(m / (M +
#' 1))` of each margin and did better in its case study (Section 5.3); it
#' needs the quantile function of each margin, which models that forecast by
#' simulation do not give.
#'
#' **Ties.** Template values that are tied within a margin get their ranks
#' at random, as the paper suggests (Section 4.1, footnote 6), with the
#' random number generator of the forecast (`seed`). Tied base values need
#' nothing, since only their sorted values are used. Draws of expected rates
#' seldom have ties; predictive draws often have them, since deaths are
#' counts. Margins without template draws (a stratum where the template
#' failed) keep the order of the base model, with a warning.
#'
#' **Package-specific choices.**
#'
#' * `coupled()` couples the draws of the models before observation noise:
#'   for models of expected rates (such as [lee_carter()] and [cae()]) the
#'   expected rates. Predictive draws (`mt_forecast(predictive = TRUE)`)
#'   then add the base model's Poisson (or negative binomial) deaths given
#'   the coupled rates, independent across cells, as for every model of the
#'   package. Without noise the margins are those of the base model exactly;
#'   with noise they are the same in distribution. The paper applies the
#'   ranks to the forecast of the observed quantity itself.
#' * Around a reconciled base model (`coupled(reconciled(...), ...)`), the
#'   bottom strata are coupled and the aggregates recomputed, so the draws
#'   stay coherent. Coupling before reconciling,
#'   `reconciled(spec, hierarchy, dependence = "ecc", template = )`, is the
#'   order that gives the totals the template's dependence between the
#'   parts.
#'
#' **What changes.** Coupling changes no single cell: the forecast
#' distribution of each cell is that of the base model, so scores of single
#' cells, such as their CRPS or the coverage of their intervals, change only
#' by simulation noise. It changes quantities computed from several cells:
#' totals over strata, gaps between strata, life expectancy (over ages) and
#' cumulative changes (over years). The template should be a joint model
#' ([li_lee()], [cae()], [coherent_fdm()], [tucker3()], `coda(coherent =
#' TRUE)`); a template fitted per stratum has independent strata and gives a
#' warning (class `mt_warning_coupled`).
#'
#' The Schaake shuffle (Clark et al. 2004) reorders draws in the same way by
#' the ranks of past observations; ECC takes the ranks from a model's
#' forecast of the same years (Schefzik et al. 2013, Section 3.4). The
#' package-specific `reconciled(dependence = "empirical")` is related to the
#' Schaake shuffle.
#'
#' `mt_couple()` couples two sets of draws with the same ages, years, strata
#' and number of draws. Both should describe the same quantity (both
#' predictive or both not). Couple before reconciling: draws reconciled
#' before give a warning, because each stratum is reordered separately and
#' their totals no longer add up.
#'
#' @param spec The base model specification (margins).
#' @param template For `coupled()`, a model specification whose draws give
#'   the dependence, usually a joint model of the same strata; for
#'   `mt_couple()`, its forecast draws ([mt_draws]).
#' @param draws Forecast draws ([mt_draws]) that give the margins.
#' @param seed Random seed for the ranks of tied template values (`NULL`
#'   leaves the random number generator as it is).
#' @return `coupled()`: a model specification; its fit holds the fit of the
#'   template in `$coupling$template`, and its forecasts report the coupling
#'   in `$coupling` (template, number of margins and share of margins with
#'   tied template values). `mt_couple()`: the draws with the values of each
#'   margin reordered, and the coupling in `$coupling`.
#' @references
#' Schefzik R, Thorarinsdottir TL, Gneiting T (2013). Uncertainty
#' quantification in complex simulation models using ensemble copula
#' coupling. *Statistical Science* 28(4): 616-640. \doi{10.1214/13-STS443}
#'
#' Clark M, Gangopadhyay S, Hay L, Rajagopalan B, Wilby R (2004). The Schaake
#' shuffle: a method for reconstructing space-time variability in forecasted
#' precipitation and temperature fields. *Journal of Hydrometeorology* 5(1):
#' 243-262.
#' @seealso [cae()], [li_lee()]; [reconciled()] with `dependence = "ecc"`.
#' @examples
#' x <- mt_subset(mt_example_regions(years = 2000:2019, regions = 2), ages = c(60, Inf))
#' fit <- mt_fit(coupled(lee_carter(), template = cae()), x)
#' fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
#' fc$coupling
#' # the same with draws
#' fc_lc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 100, seed = 1)
#' fc_cae <- mt_forecast(mt_fit(cae(), x), h = 5, n_draws = 100, seed = 2)
#' fc2 <- mt_couple(fc_lc, fc_cae, seed = 3)
#' @export
coupled <- function(spec, template) {
  check_class(spec, "mt_spec", "spec")
  check_class(template, "mt_spec", "template")
  if (!spec_is_joint(template)) {
    mt_warn(paste("The template's draws are independent between strata (it is not a joint model): the coupled",
                  "forecasts then have independent strata. Use a joint model such as cae() or li_lee()."),
            class = "coupled")
  }
  new_mt_spec("coupled", list(base = spec, template = template), joint = spec$joint,
              min_periods = max(spec$min_periods, template$min_periods), smooth = spec$smooth,
              label = sprintf("%s, coupled to %s (ECC)", spec$label, template$label))
}

#' @export
mt_fit.mt_spec_coupled <- function(spec, data, prep = NULL, cores = 1, ...) {
  p <- spec$params
  fit <- mt_fit(p$base, data, prep = prep, cores = cores, ...)
  tpl <- mt_fit(p$template, data, prep = prep, cores = cores, ...)
  if (!setequal(tpl$strata$label, fit$strata$label) || !identical(tpl$age$labels, fit$age$labels) ||
      !identical(tpl$periods, fit$periods)) {
    mt_abort(paste("The template is fitted to other strata, ages or years than the base model. Wrap both in the",
                   "same way, for example reconciled(coupled(spec, template), hierarchy)."), class = "coupled")
  }
  wrap_fit(fit, spec, "mt_fit_coupled", list(coupling = list(template = tpl)))
}

#' @export
mt_forecast.mt_fit_coupled <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                                       exposure = NULL, ...) {
  base_fit <- unwrap_fit(fit, "mt_fit_coupled")
  tpl_fit <- fit$coupling$template
  if (!is.null(to)) {
    if (!missing(h)) mt_warn("Both `h` and `to` are given; `to` is used.", class = "argument")
    h <- to - max(fit$periods)
  }
  out <- with_seed(seed, {
    # the base model first, so that its draws are those of mt_forecast(base, seed = seed)
    fc <- suppressWarnings(mt_forecast(base_fit, h = h, n_draws = n_draws, exposure = exposure))
    ft <- suppressWarnings(mt_forecast(tpl_fit, h = h, n_draws = n_draws))
    targets <- calibration_targets(fit)
    idx <- match(fc$strata$label, ft$strata$label)
    ec <- ecc_reorder(fc$rates, ft$rates[, , idx, , drop = FALSE], targets)
    arr <- ec$rates
    hs <- fit$reconcile$structure
    E <- expand_exposure(exposure %||% fit$exposure_last, dim(arr), dimnames(arr), hs)
    # observation noise given the coupled rates, from the base draw each value
    # comes from (ensembles: its component)
    if (predictive) arr <- add_observation_noise(arr, base_fit, fc, targets, E, draw = ec$draw)
    if (!is.null(hs)) arr <- aggregate_bottom(arr, E, hs, fit$reconcile$hierarchy$additive)
    list(fc = fc, rates = arr, info = c(list(template = ft$model), ec$info))
  })
  warn_uncoupled(out$info$uncoupled)
  res <- out$fc
  res$rates <- out$rates
  dimnames(res$rates) <- dimnames(out$fc$rates)
  check_draws_plausible(res$rates, fit$age, fit$strata, fit$spec$label)
  res$model <- fit$spec$label
  res$spec <- fit$spec
  res$predictive <- predictive
  # after coupling, the values of one draw come from different draws of the base model
  res$component <- NULL
  res$nb_size <- NULL
  res$coupling <- out$info[c("template", "method", "margins", "tied")]
  res
}

#' @rdname coupled
#' @export
mt_couple <- function(draws, template, seed = NULL) {
  check_class(draws, "mt_draws", "draws")
  check_class(template, "mt_draws", "template")
  d <- dim(draws$rates)
  idx <- match(draws$strata$label, template$strata$label)
  if (anyNA(idx) || !identical(draws$age$labels, template$age$labels) ||
      !isTRUE(all.equal(as.numeric(draws$periods), as.numeric(template$periods))) ||
      dim(template$rates)[4] != d[4]) {
    mt_abort("`draws` and `template` need the same ages, years, strata and number of draws.", class = "coupled")
  }
  if (!identical(isTRUE(draws$predictive), isTRUE(template$predictive))) {
    mt_warn(paste("One of `draws` and `template` is predictive and the other is not: couple draws that describe the",
                  "same quantity."), class = "coupled")
  }
  if (!is.null(draws$reconciliation)) {
    mt_warn(paste("`draws` were reconciled: coupling reorders each stratum separately, so their totals no longer",
                  "add up. Couple before reconciling."), class = "coupled")
  }
  ec <- with_seed(seed, ecc_reorder(draws$rates, template$rates[, , idx, , drop = FALSE]))
  warn_uncoupled(ec$info$uncoupled)
  draws$rates[] <- ec$rates
  draws$component <- NULL
  draws$nb_size <- NULL
  draws$coupling <- c(list(template = template$model), ec$info[c("method", "margins", "tied")])
  draws
}

# Ensemble copula coupling (Schefzik, Thorarinsdottir and Gneiting 2013,
# eq. 4.3): for each margin (age, year, stratum among `targets`), the sorted
# base values are given to the draws in the order of the template's values,
# whose ties are ranked at random (footnote 6). Margins with missing values
# in either array keep the base values in their order. Returns the rates,
# the base draw each value comes from (`draw`, same dimensions) and a
# summary (`info`: number of coupled margins, share with tied template
# values, strata the template has no draws for)
ecc_reorder <- function(arr, tpl, targets = seq_len(dim(arr)[3])) {
  d <- dim(arr)
  M <- d[4]
  draw <- array(rep(seq_len(M), each = prod(d[1:3])), d)
  info <- list(method = "ECC-R", margins = 0L, tied = 0, uncoupled = character())
  if (M < 2 || !length(targets)) return(list(rates = arr, draw = draw, info = info))
  B <- matrix(arr[, , targets, , drop = FALSE], ncol = M)
  Tm <- matrix(tpl[, , targets, , drop = FALSE], ncol = M)
  ok <- rowSums(is.na(B)) == 0 & rowSums(is.na(Tm)) == 0
  miss_t <- array(rowSums(is.na(Tm)) > 0 & rowSums(is.na(B)) == 0, c(d[1:2], length(targets)))
  info$uncoupled <- dimnames(arr)[[3]][targets][apply(miss_t, 3, any)]
  L <- sum(ok)
  if (!L) return(list(rates = arr, draw = draw, info = info))
  Bo <- B[ok, , drop = FALSE]
  To <- Tm[ok, , drop = FALSE]
  # o[l, k]: the base draw with the k-th smallest value; r[l, m]: rank of template draw m
  o <- t(apply(Bo, 1, order))
  r <- t(apply(To, 1, rank, ties.method = "random"))
  li <- rep(seq_len(L), M)
  src <- matrix(o[cbind(li, as.vector(r))], L)
  B[ok, ] <- Bo[cbind(li, as.vector(src))]
  Dr <- matrix(draw[, , targets, , drop = FALSE], ncol = M)
  Dr[ok, ] <- src
  arr[, , targets, ] <- B
  draw[, , targets, ] <- Dr
  info$margins <- L
  info$tied <- mean(apply(To, 1, anyDuplicated) > 0)
  list(rates = arr, draw = draw, info = info)
}

warn_uncoupled <- function(labels) {
  if (!length(labels)) return(invisible())
  mt_warn(sprintf(paste("The template has no draws for %d stratum/strata (%s): their draws keep the order of the",
                        "base model."), length(labels), paste(utils::head(labels, 5), collapse = ", ")),
          class = "coupled", call = NULL)
}
