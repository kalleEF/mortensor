#' Bayesian age-period-cohort model (bamp)
#'
#' Wrapper around the \pkg{bamp} package (Schmid and Held 2007): binomial
#' likelihood for deaths given exposures, logit link, random walk priors of
#' first or second order for age, period and cohort effects, optional
#' overdispersion, MCMC with four chains. Forecasts extrapolate period and
#' cohort effects as random walks for every posterior draw
#' (`bamp::predict_apc()`), so the draws include parameter uncertainty.
#'
#' The wrapper makes bamp behave like every other model:
#'
#' * **Resolution.** With `resolution = "observed"` (default) the model is
#'   fitted where the deaths were observed: 5-year age groups with annual
#'   periods use `periods_per_agegroup = 5` (Knorr-Held and Rainer 2001).
#'   Exposures are aggregated to the death cells. Age groups must have equal
#'   widths; an open last group is treated as one more group of that width.
#' * **Grid.** Draws for age groups are mapped to the modelling grid with
#'   age-specific shape factors (single-year rate relative to its group rate,
#'   exposure-weighted to 1 within each group) estimated from the
#'   preprocessed training data of the last `shape_years` years. Use
#'   `step_ungroup()` in `prep` when deaths are grouped.
#' * **Convergence.** `bamp::checkConvergence()` is called explicitly; a
#'   stratum that has not converged is refitted with twice as many iterations
#'   up to `max_refits` times, then reported as failed in the diagnostics.
#' * **Integers.** bamp needs whole numbers; deaths and exposures are rounded
#'   and the largest relative rounding change is recorded.
#'
#' Heterogeneity variants (`"rw1+het"`, `"rw2+het"`) are not offered because
#' `bamp::predict_apc()` does not forecast them in bamp 3.0.0.
#'
#' @section Coherent model (`coherent = TRUE`):
#' Fits the strata of each group of `common_by` (all strata without it) in
#' one model, so that their forecasts do not diverge:
#'
#' `logit p[s, i, j] = mu + a[s] + theta[s, i] + phi[j] + psi[k] + d[s, j]`
#'
#' for stratum `s`, age group `i`, year `j` and cohort `k`, with binomial
#' deaths given the exposures. The period effect `phi` and the cohort
#' effect `psi` are shared and carry the common trend; their random-walk
#' priors (orders `period`, `cohort`) are scaled to unit generalised
#' variance (Sørbye and Rue 2014) and forecast as random walks. Each stratum
#' has its own age effect `theta[s, ]` (random walk of order `age`, one
#' precision for all strata) and level `a[s]`. The deviations `d[s, ]` of
#' the strata from the shared period effect are stationary AR(1) series with
#' a common coefficient `rho` (uniform prior on 0 to 1, sampled) and
#' innovation precision; they sum to zero over the strata in every year, or,
#' with `reference`, the deviation of the reference stratum is zero. Their
#' forecasts revert towards zero, so the gaps between strata stay bounded.
#' `cohort_deviation = TRUE` adds deviations `c[s, k]` of the same form
#' along the cohorts. The precisions have Gamma priors with shape 1 and
#' rates 0.5 (age) and 5e-4 (period, cohort, deviations). For the deviations
#' the branch uses the rate 0.05, which dominates typical data: the
#' posterior of their precision and the width of the forecast gaps are then
#' set by the prior.
#'
#' **Sources.** The model is the coherent model of the development branch
#' "coherent_forecasting" of bamp (`bamp_coherent()`, by Chris
#' Kypridemos), rewritten in this package because that branch is not
#' released; when a bamp release contains it, `apc_bayes(coherent = TRUE)`
#' will call bamp and the native sampler will be deprecated. Its published
#' anchor is the two-population Bayesian APC model of Cairns et al. (2011):
#' a common period trend and a mean-reverting AR(1) spread with level zero,
#' population-specific age effects, a single-stage Bayesian fit, and two
#' parameterisations, the spread between equal populations (their Section
#' 4.5, the default here) and the spread from a dominant population (Section
#' 4.4, `reference`). It differs from that paper, and these parts are
#' package-specific (reference unknown): the binomial-logit likelihood and
#' the random-walk smoothing priors of bamp (Schmid and Held 2007) in place
#' of the Poisson likelihood with a random walk with drift and AR(2) cohort
#' processes; one cohort effect shared by all strata (the alternative
#' cohort hypothesis of their Section 4.3), with optional AR(1) cohort
#' deviations in place of an AR(2) cohort spread; independent innovations
#' of the common trend and the deviations (the paper: correlated); the
#' priors; and more than two strata.
#'
#' **Sampler.** Four chains of a Metropolis-Hastings sampler: all effects are
#' proposed jointly from a Gaussian approximation of their conditional
#' posterior (one Newton step from the current state, with sparse Cholesky
#' factors; Rue and Held 2005, Section 4.4), the precisions are drawn from
#' their Gamma conditionals, `rho` by a random-walk Metropolis step on the
#' logit scale whose step size adapts during the burn-in only. Convergence
#' is checked with the Gelman-Rubin statistic of the hyperparameters, the
#' last period effect and the level of every stratum in the last year
#' (`rhat` below 1.1); a group that has not converged is refitted with
#' twice as many iterations up to `max_refits` times. Overdispersion is not
#' available in the coherent model.
#'
#' Draws are posterior draws of the rates. Predictive draws
#' (`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths
#' given the rates of each draw and the exposures, as for the other models
#' of expected rates; no other observation noise is added. bamp's binomial
#' likelihood and its overdispersion term (`overdispersion = TRUE`) are not
#' carried into the predictive draws: the overdispersion only enters the
#' posterior of the rates. If intervals are too narrow on real data, that is
#' a property of the model (see [recalibrated()], [ensemble()], or
#' overdispersed models such as `stanmomo(family = "nb")`).
#'
#' @param age,period,cohort Prior for each effect: `"rw1"` or `"rw2"`. The
#'   defaults use first-order random walks for period and cohort: with
#'   second-order random walks the extrapolated trends can make forecast
#'   intervals span orders of magnitude. First-order random walks carry the
#'   last period effect forward without a trend, so forecasts are biased when
#'   mortality keeps declining. The vignette "Bayesian age-period-cohort
#'   models with bamp" compares both.
#' @param overdispersion Include cell-level overdispersion.
#' @param mcmc MCMC settings from [mcmc_control()].
#' @param resolution `"observed"` or `"grid"`.
#' @param shape_years Number of recent years used for the shape factors.
#' @param max_refits Refits with longer chains if not converged.
#' @param coherent Fit the strata of each group jointly with shared period
#'   and cohort effects and mean-reverting stratum deviations (see the
#'   section "Coherent model"); needs at least two strata per group.
#' @param common_by Stratum variable(s) defining groups with their own
#'   common trend (`coherent = TRUE`), for example `"sex"` for the
#'   deprivation groups of each sex. `NULL`: one group of all strata.
#' @param reference With `coherent = TRUE`: `NULL` (the deviations sum to
#'   zero over the strata of a group) or a named list of values of stratum
#'   variables that selects one reference stratum in every group, for
#'   example `list(quintile = "England")`, whose deviation is zero.
#' @param cohort_deviation With `coherent = TRUE`: add AR(1) deviations of
#'   the strata along the cohorts.
#' @param seed Seed of the native sampler (`coherent = TRUE`); `NULL` draws
#'   one from R's random number generator, so that `set.seed()` makes fits
#'   reproducible.
#' @return A model specification.
#' @references
#' Knorr-Held L, Rainer E (2001). Prognosis of lung cancer mortality in
#' West Germany: a case study in Bayesian prediction. *Biostatistics* 2(1):
#' 109-129.
#'
#' Schmid VJ, Held L (2007). Bayesian age-period-cohort modeling and
#' prediction: BAMP. *Journal of Statistical Software* 21(8): 1-15.
#' \doi{10.18637/jss.v021.i08}
#'
#' Cairns AJG, Blake D, Dowd K, Coughlan GD, Khalaf-Allah M (2011). Bayesian
#' stochastic mortality modelling for two populations. *ASTIN Bulletin*
#' 41(1): 29-59. \doi{10.2143/AST.41.1.2084385}
#'
#' Gelman A, Rubin DB (1992). Inference from iterative simulation using
#' multiple sequences. *Statistical Science* 7(4): 457-472.
#' \doi{10.1214/ss/1177011136}
#'
#' Rue H, Held L (2005). *Gaussian Markov Random Fields: Theory and
#' Applications*. Chapman and Hall/CRC, Boca Raton.
#'
#' Sørbye SH, Rue H (2014). Scaling intrinsic Gaussian Markov random field
#' priors in spatial modelling. *Spatial Statistics* 8: 39-51.
#' \doi{10.1016/j.spasta.2013.06.004}
#' @examples
#' \donttest{
#' if (requireNamespace("bamp", quietly = TRUE)) {
#'   x <- mt_example(years = 2005:2019)
#'   fit <- mt_fit(apc_bayes(mcmc = mcmc_control("test")), x, prep = mt_prep(step_ungroup()))
#'   fc <- mt_forecast(fit, h = 5, n_draws = 200, seed = 1)
#' }
#' }
#' # coherent model of the two sexes (native sampler, no bamp needed)
#' \donttest{
#' x <- mt_example(years = 2005:2019)
#' fit <- mt_fit(apc_bayes(coherent = TRUE, mcmc = mcmc_control("test"), seed = 1), x,
#'               prep = mt_prep(step_ungroup()))
#' }
#' @export
apc_bayes <- function(age = "rw2", period = "rw1", cohort = "rw1", overdispersion = FALSE,
                      mcmc = mcmc_control(), resolution = c("observed", "grid"), shape_years = 5,
                      max_refits = 1, coherent = FALSE, common_by = NULL, reference = NULL,
                      cohort_deviation = FALSE, seed = NULL) {
  check_flag(coherent, "coherent")
  check_flag(cohort_deviation, "cohort_deviation")
  if (!coherent && (!is.null(common_by) || !is.null(reference) || cohort_deviation)) {
    mt_abort("`common_by`, `reference` and `cohort_deviation` need `coherent = TRUE`.", class = "argument")
  }
  if (coherent && overdispersion) {
    mt_abort("Overdispersion is not available in the coherent model.", class = "argument")
  }
  if (!is.null(reference) && (!is.list(reference) || is.null(names(reference)))) {
    mt_abort("`reference` must be a named list of values of stratum variables.", class = "argument")
  }
  if (!coherent) need_pkg("bamp", "for apc_bayes()")
  for (nm in c("age", "period", "cohort")) {
    v <- get(nm)
    if (!v %in% c("rw1", "rw2")) {
      mt_abort(sprintf("`%s` must be \"rw1\" or \"rw2\" (heterogeneity variants cannot be forecast with bamp).", nm),
               class = "argument")
    }
  }
  resolution <- match.arg(resolution)
  if (coherent) {
    return(new_mt_spec("apc_coherent", list(age = age, period = period, cohort = cohort, mcmc = mcmc,
                                            resolution = resolution, shape_years = shape_years,
                                            max_refits = max_refits, common_by = common_by, reference = reference,
                                            cohort_deviation = cohort_deviation, seed = seed),
                       joint = TRUE, min_periods = 5, label = "Bayesian APC, coherent"))
  }
  new_mt_spec("apc_bayes", list(age = age, period = period, cohort = cohort, overdispersion = overdispersion,
                                mcmc = mcmc, resolution = resolution, shape_years = shape_years,
                                max_refits = max_refits),
              min_periods = 5, label = "Bayesian APC (bamp)")
}

#' MCMC settings for apc_bayes()
#'
#' Presets: `"fast"` (5,000 iterations, burn-in 1,000, thinning 5; suitable
#' for backtests), `"final"` (bamp's automatic choice of chain length), and
#' `"test"` (1,000 iterations, for examples and unit tests only).
#'
#' @param preset `"fast"`, `"final"` or `"test"`.
#' @param iterations,burn_in,thin Override the preset.
#' @param tuning Tuning iterations.
#' @param parallel Run the four chains in parallel. bamp runs them in a
#'   socket cluster without seeded random number streams, so results differ
#'   slightly between runs even when a seed is set.
#' @return A list.
#' @examples
#' mcmc_control("fast")
#' @export
mcmc_control <- function(preset = c("fast", "final", "test"), iterations = NULL, burn_in = NULL, thin = NULL,
                         tuning = NULL, parallel = TRUE) {
  preset <- match.arg(preset)
  base <- switch(preset,
                 fast = list(number_of_iterations = 5000, burn_in = 1000, step = 5, tuning = 500),
                 final = list(number_of_iterations = "auto", burn_in = "auto", step = "auto", tuning = 500),
                 test = list(number_of_iterations = 1000, burn_in = 200, step = 2, tuning = 100))
  if (!is.null(iterations)) base$number_of_iterations <- iterations
  if (!is.null(burn_in)) base$burn_in <- burn_in
  if (!is.null(thin)) base$step <- thin
  if (!is.null(tuning)) base$tuning <- tuning
  list(options = base, parallel = parallel, preset = preset)
}

#' @export
mt_fit.mt_spec_apc_bayes <- function(spec, data, prep = NULL, cores = 1, ...) {
  check_class(data, "mt_data", "data")
  need_pkg("bamp", "for apc_bayes()")
  pr <- mt_prep_apply_raw(prep, data)
  d <- pr$data
  md <- model_inputs(d, spec)
  D_obs <- raw_deaths(pr$raw, md) %||% md$D
  p <- spec$params
  src <- if (p$resolution == "observed") data else d
  cells <- apc_cells(src, md$periods)
  shape <- apc_shape(md, cells$age, p$shape_years)
  st <- d$strata
  fits <- vector("list", nrow(st))
  diag <- vector("list", nrow(st))
  for (s in seq_len(nrow(st))) {
    t0 <- proc.time()[["elapsed"]]
    res <- tryCatch(apc_fit_one(spec, cells$D[, , s], cells$E[, , s], cells$ppa), error = function(e) e)
    secs <- proc.time()[["elapsed"]] - t0
    if (inherits(res, "error")) {
      diag[[s]] <- data.table(stratum = st$label[s], status = "failed", message = conditionMessage(res),
                              seconds = secs)
      next
    }
    res$shape <- shape$r[, s]
    res$group <- shape$group
    res$fitted <- res$fitted_groups[shape$group, , drop = FALSE] * shape$r[, s]
    fits[[s]] <- res
    diag[[s]] <- data.table(stratum = st$label[s], status = "ok", message = res$message, seconds = secs)
  }
  diagnostics <- rbindlist(diag)
  if (all(diagnostics$status == "failed")) {
    mt_abort(sprintf("bamp failed in every stratum: %s", diagnostics$message[1]), class = "model_fit")
  }
  structure(list(spec = spec, prep = prep, fits = fits, diagnostics = diagnostics, age = md$age,
                 periods = md$periods, strata = st, omega = d$omega, observed = D_obs / md$E,
                 train_D = md$D, raw_D = raw_deaths(pr$raw, md), train_E = md$E,
                 exposure_last = md$E[, length(md$periods), , drop = FALSE], level_totals = data$level_totals,
                 apc_age = cells$age),
            class = "mt_fit")
}

# Deaths and exposures at the death resolution: arrays [age group, period, stratum]
apc_cells <- function(x, periods) {
  D <- x$obs$deaths
  D <- D[D$period_lo >= min(periods) & D$period_hi <= max(periods) + 1]
  if (any(D$period_hi - D$period_lo != 1)) mt_abort("apc_bayes() needs annual periods.", class = "model_data")
  ag <- unique(D[, c("age_lo", "age_hi")])
  ag <- ag[order(age_lo)]
  if (anyDuplicated(ag$age_lo)) {
    mt_abort("Age groups of the deaths change over time; regroup them first with mt_regroup().", class = "model_data")
  }
  axis <- mt_axis(ag$age_lo, ag$age_hi, kind = "age")
  w <- unique((ag$age_hi - ag$age_lo)[is.finite(ag$age_hi)])
  if (length(w) != 1L) {
    mt_abort("bamp needs age groups of equal width; regroup the deaths with mt_regroup().", class = "model_data")
  }
  g <- list(age = axis, period = single_years(min(periods), max(periods)))
  xx <- x
  xx$obs <- list(deaths = D)
  xx <- mt_set_grid(mt_window(xx, min(periods), max(periods)), age = axis, period = g$period)
  Dm <- mt_array(xx, "deaths")
  Ex <- x$obs$exposure
  ia <- axis_locate(axis, Ex$age_lo, Ex$age_hi)
  ip <- axis_locate(g$period, Ex$period_lo, Ex$period_hi)
  inside <- Ex$age_lo >= min(axis$lower) & Ex$age_hi <= max(axis$upper) & !is.na(ip)
  if (any(inside & is.na(ia))) mt_abort("Exposures are coarser than the deaths.", class = "model_data")
  ok <- inside & !is.na(ia)
  Em <- array(0, dim(Dm))
  agg <- data.table(a = ia[ok], p = ip[ok], s = Ex$stratum[ok], v = Ex$value[ok])
  agg <- agg[, list(v = sum(v)), by = c("a", "p", "s")]
  Em[cbind(agg$a, agg$p, agg$s)] <- agg$v
  if (anyNA(Dm) || any(Em <= 0)) mt_abort("Deaths or exposures are incomplete at the death resolution.",
                                          class = "model_data")
  list(D = Dm, E = Em, age = axis, ppa = w)
}

# Shape factors: grid rate relative to its group rate, exposure-weighted to 1 within groups
apc_shape <- function(md, group_axis, years) {
  ga <- md$age
  group <- findInterval(ga$lower, group_axis$lower)
  if (any(group == 0) || any(ga$upper > group_axis$upper[pmax(group, 1)] + 1e-9)) {
    mt_abort("The grid does not nest in the age groups of the deaths.", class = "model_data")
  }
  nt <- length(md$periods)
  use <- max(1, nt - years + 1):nt
  ns <- dim(md$D)[3]
  r <- matrix(1, length(group), ns)
  for (s in seq_len(ns)) {
    D <- md$D[, use, s, drop = FALSE]
    E <- md$E[, use, s, drop = FALSE]
    m <- rowSums(D) / rowSums(E)
    mg <- tapply(rowSums(D), group, sum) / tapply(rowSums(E), group, sum)
    rr <- pmax(m, 1e-12) / mg[as.character(group)]
    El <- md$E[, nt, s]
    norm <- tapply(El * rr, group, sum) / tapply(El, group, sum)
    r[, s] <- rr / norm[as.character(group)]
  }
  list(r = r, group = group)
}

# Gelman-Rubin potential scale reduction factors (R-hat, point estimates) of
# the identified quantities, as reported by bamp::checkConvergence()
apc_psrf <- function(m) {
  out <- tryCatch(utils::capture.output(bamp::checkConvergence(m, info = TRUE, level = 1)),
                  error = function(e) character())
  i <- grep("Identified quantities", out)
  j <- grep("Raw per-effect", out)
  if (!length(i)) return(NULL)
  lines <- out[(i + 1):(if (length(j)) j - 1 else length(out))]
  m <- regmatches(lines, regexec("^\\s+(\\S+)\\s+([0-9.]+)", lines))
  m <- m[lengths(m) == 3]
  if (!length(m)) return(NULL)
  stats::setNames(as.numeric(vapply(m, `[`, "", 3)), vapply(m, `[`, "", 2))
}

apc_fit_one <- function(spec, D, E, ppa) {
  p <- spec$params
  cases <- round(t(D))
  pop <- round(t(E))
  rounding <- max(abs(pop - t(E)) / pmax(t(E), 1))
  opts <- p$mcmc$options
  msg <- character()
  for (attempt in 0:p$max_refits) {
    utils::capture.output(m <- suppressWarnings(bamp::bamp(
      cases, pop, age = p$age, period = p$period, cohort = p$cohort, overdisp = p$overdispersion,
      periods_per_agegroup = ppa, mcmc.options = opts, parallel = p$mcmc$parallel, verbose = FALSE)))
    conv <- isTRUE(tryCatch(bamp::checkConvergence(m), error = function(e) FALSE))
    if (conv) break
    if (is.numeric(opts$number_of_iterations)) {
      opts$number_of_iterations <- 2 * opts$number_of_iterations
      opts$burn_in <- 2 * opts$burn_in
    }
    msg <- c(msg, sprintf("not converged (attempt %d)", attempt + 1))
  }
  if (!conv) {
    mt_abort(sprintf("bamp chains did not converge after %d refit(s); use longer chains (mcmc_control()).",
                     p$max_refits), class = "convergence")
  }
  p0 <- bamp::predict_apc(m, periods = 0)
  fitted_groups <- t(apply(p0$samples$pr, c(1, 2), mean))
  list(model = m, fitted_groups = fitted_groups, rounding = rounding, psrf = apc_psrf(m),
       message = paste(c(msg, if (rounding > 0.001) sprintf("exposures rounded (max change %.2g)", rounding)),
                       collapse = "; "))
}

#' @export
simulate_model.mt_spec_apc_bayes <- function(spec, fit, h, n_draws, ...) {
  pr <- bamp::predict_apc(fit$model, periods = h)$samples$pr
  np <- dim(pr)[1]
  fut <- pr[(np - h + 1):np, , , drop = FALSE]
  k <- dim(fut)[3]
  pick <- if (n_draws <= k) sample.int(k, n_draws) else sample.int(k, n_draws, replace = TRUE)
  grp <- fut[, , pick, drop = FALSE]
  out <- array(NA_real_, c(length(fit$group), h, n_draws))
  for (j in seq_len(n_draws)) {
    out[, , j] <- t(grp[, fit$group, j, drop = FALSE][, , 1]) * fit$shape
  }
  out
}

# Age, period and cohort effects of a bamp fit: posterior medians and 90%
# intervals on the logit scale (bamp's identification); see mt_parameters()
apc_effects <- function(fit) {
  check_class(fit, "mt_fit", "fit")
  if (!inherits(fit$spec, "mt_spec_apc_bayes")) mt_abort("Not an apc_bayes() fit.", class = "argument")
  rbindlist(lapply(seq_along(fit$fits), function(s) {
    f <- fit$fits[[s]]
    if (is.null(f)) return(NULL)
    smp <- f$model$samples
    rbindlist(lapply(c("age", "period", "cohort"), function(e) {
      se <- smp[[e]]
      M <- if (is.list(se) && !is.data.frame(se)) do.call(rbind, lapply(se, as.matrix)) else as.matrix(se)
      q <- t(apply(M, 2, stats::quantile, c(0.05, 0.5, 0.95)))
      data.table(stratum = fit$strata$label[s], effect = e, index = seq_len(nrow(q)), q0.05 = q[, 1],
                 q0.5 = q[, 2], q0.95 = q[, 3], convention = attr(f$model$effects, "gauge_convention") %||% "")
    }))
  }))
}

#' Choose the random-walk orders of a Bayesian APC model by DIC
#'
#' Fits [apc_bayes()] with every combination of first- and second-order
#' random walks for age, period and cohort (and, optionally, with and
#' without overdispersion) and ranks the variants by the deviance
#' information criterion (Spiegelhalter et al. 2002) that bamp reports,
#' summed over strata. DIC
#' measures in-sample fit; to choose by forecast accuracy, compare the
#' candidates with [mt_backtest()] instead. Selection on all data followed by
#' a backtest of the selected model reuses the test years once.
#'
#' @param data An [mt_data] object.
#' @param prep Optional [mt_prep()] pipeline.
#' @param orders Random-walk orders to try.
#' @param overdispersion Values of `overdispersion` to try.
#' @param mcmc MCMC settings, see [mcmc_control()].
#' @param ... Further arguments of [apc_bayes()].
#' @return A list of class `mt_apc_selection` with `best` (the chosen
#'   specification) and `table` (DIC, effective number of parameters `pd`
#'   and deviance at the posterior mean per candidate, sorted by DIC; DIC =
#'   deviance at the mean + 2 pd).
#' @references
#' Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A (2002). Bayesian
#' measures of model complexity and fit. *Journal of the Royal Statistical
#' Society B* 64(4): 583-639.
#' @examples
#' \donttest{
#' if (requireNamespace("bamp", quietly = TRUE)) {
#'   x <- mt_example(years = 2005:2019)
#'   sel <- apc_select(x, prep = mt_prep(step_ungroup()), orders = "rw1", mcmc = mcmc_control("test"))
#'   sel$table
#' }
#' }
#' @export
apc_select <- function(data, prep = NULL, orders = c("rw1", "rw2"), overdispersion = FALSE,
                       mcmc = mcmc_control("fast"), ...) {
  need_pkg("bamp", "for apc_select()")
  check_class(data, "mt_data", "data")
  grid <- expand.grid(age = orders, period = orders, cohort = orders, overdispersion = overdispersion,
                      stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    spec <- apc_bayes(age = g$age, period = g$period, cohort = g$cohort, overdispersion = g$overdispersion,
                      mcmc = mcmc, ...)
    fit <- tryCatch(mt_fit(spec, data, prep = prep), error = function(e) NULL)
    dev <- if (is.null(fit)) NULL else lapply(Filter(Negate(is.null), fit$fits), function(f) f$model$deviance)
    ok <- !is.null(fit) && length(dev) == nrow(fit$strata) && all(vapply(dev, function(d) !is.null(d$DIC), TRUE))
    data.table(age = g$age, period = g$period, cohort = g$cohort, overdispersion = g$overdispersion,
               dic = if (ok) sum(vapply(dev, function(d) d$DIC, 0)) else NA_real_,
               pd = if (ok) sum(vapply(dev, function(d) d$pD, 0)) else NA_real_,
               deviance_at_mean = if (ok) sum(vapply(dev, function(d) d$deviance.mean, 0)) else NA_real_,
               converged = ok, spec = list(spec))
  })
  tab <- rbindlist(rows)
  if (!any(tab$converged)) mt_abort("No candidate could be fitted in every stratum.", class = "model_fit")
  setorderv(tab, "dic", na.last = TRUE)
  structure(list(best = tab$spec[[1]], table = tab[, !"spec"]), class = "mt_apc_selection")
}

#' @export
print.mt_apc_selection <- function(x, ...) {
  b <- x$best$params
  cat(sprintf("<mt_apc_selection> best: age %s, period %s, cohort %s%s\n", b$age, b$period, b$cohort,
              if (isTRUE(b$overdispersion)) ", overdispersion" else ""))
  print(x$table, ...)
  invisible(x)
}
