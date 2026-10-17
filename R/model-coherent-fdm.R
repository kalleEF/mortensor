# Product-ratio functional demographic model.
# Adapted from demography 2.0.1 (coherentfdm() and forecast.fdmpr()),
# GPL (>= 3), by Rob J Hyndman, with contributions from Heather Booth,
# Leonie Tickle and John Maindonald. Changes: the product is the geometric
# mean of the bottom strata of each group (aggregates of a hierarchy get a
# ratio but do not enter the product), the rates of each stratum are
# smoothed over age with smooth.demogdata() before the decomposition, and
# forecasts are simulated paths (product and ratios simulated independently,
# rescaled path by path so that the ratios of each group multiply to one).

#' Coherent functional demographic model (product-ratio, demography)
#'
#' The product-ratio method of Hyndman, Booth and Yasmeen (2013), with the
#' functional data models of the demography package: for a group of strata
#' (for example women and men, or the deprivation quintiles of one sex)
#'
#' * the **product**, the geometric mean of the strata's rates, is modelled
#'   with a functional demographic model ([hyndman_ullah()]) whose scores are
#'   forecast with ARIMA models, so it carries the common trend;
#' * the **ratio** of each stratum to the product is modelled with a
#'   functional demographic model whose scores are forecast with stationary
#'   ARFIMA models, so the differences between the strata stay bounded and
#'   revert towards their recent level.
#'
#' Forecasts of the strata are products of the two. This is the original
#' method of `demography::coherentfdm()`; [product_ratio()] is a Lee-Carter
#' version of it, and [li_lee()] a related model with a common Lee-Carter
#' factor.
#'
#' Draws simulate the product and each ratio (uncertainty of the scores and
#' resampled model errors, as `simulate()` for functional forecasts in
#' demography) and multiply them; the ratios of each path are rescaled so
#' that they multiply to one, as in `demography::forecast.fdmpr()`.
#' Predictive draws add the extra-Poisson and Poisson variation of the raw
#' deaths.
#'
#' With `common_by` (for example `"sex"`), each level of that stratum
#' variable forms its own group with its own product (for example the
#' quintiles of each sex). If the data carry a hierarchy (see
#' [mt_aggregate_strata()]), the product is computed from the bottom strata;
#' aggregates get their own ratio to it. Every group needs at least two
#' strata ([hyndman_ullah()] models a single stratum).
#'
#' @param order1,order2 Number of principal components of the product and
#'   of the ratios.
#' @param common_by Stratum variable(s) defining groups with their own
#'   product.
#' @param smooth Smooth the rates of each stratum and year over age before
#'   the decomposition (see [hyndman_ullah()]).
#' @param smooth_args List of arguments of `demography::smooth.demogdata()`.
#' @param K Number of most recent years used to forecast the ratios.
#' @param drange Range of the fractional differencing parameter of the
#'   ARFIMA models of the ratios (`c(0, 0.5)`: stationary).
#' @return A model specification.
#' @references
#' Hyndman RJ, Booth H, Yasmeen F (2013). Coherent mortality forecasting: the
#' product-ratio method with functional time series models. *Demography*
#' 50(1): 261-283. \doi{10.1007/s13524-012-0145-5}
#' @examples
#' if (requireNamespace("demography", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 1990:2019)
#'   fit <- mt_fit(coherent_fdm(order1 = 2, order2 = 1), x)
#'   fc <- mt_forecast(fit, h = 10, n_draws = 50, seed = 1)
#' }
#' @export
coherent_fdm <- function(order1 = 6, order2 = 6, common_by = NULL, smooth = TRUE, smooth_args = list(), K = 100,
                         drange = c(0, 0.5)) {
  need_pkg("demography", "for coherent_fdm()")
  check_flag(smooth, "smooth")
  if (!is.list(smooth_args)) mt_abort("`smooth_args` must be a list.", class = "argument")
  label <- if (is.null(common_by)) "Coherent fdm" else
    sprintf("Coherent fdm (product by %s)", paste(common_by, collapse = ", "))
  new_mt_spec("coherent_fdm", list(order1 = order1, order2 = order2, common_by = common_by, smooth = smooth,
                                   smooth_args = smooth_args, K = K, drange = drange),
              joint = TRUE, min_periods = 8, label = label)
}

#' @export
fit_model.mt_spec_coherent_fdm <- function(spec, D, E, ages, periods, strata = NULL, hierarchy = NULL, ...) {
  p <- spec$params
  ns <- dim(D)[3]
  group <- rep(1L, ns)
  if (!is.null(p$common_by)) {
    if (is.null(strata) || !all(p$common_by %in% names(strata))) {
      mt_abort("`common_by` must name stratum variables.", class = "argument")
    }
    key <- do.call(paste, c(lapply(p$common_by, function(v) strata[[v]]), sep = "."))
    group <- match(key, unique(key))
  }
  pool <- rep(TRUE, ns)
  if (!is.null(hierarchy) && !is.null(strata)) pool <- hierarchy_structure(strata, hierarchy)$is_bottom
  pool_size <- vapply(seq_len(max(group)), function(g) {
    if (any(group == g & pool)) sum(group == g & pool) else sum(group == g)
  }, 1L)
  if (any(pool_size < 2)) {
    mt_abort(paste("coherent_fdm() needs at least two strata in every group: the product is their geometric mean,",
                   "and the ratio of a single stratum to it is always one. Use hyndman_ullah() for a single",
                   "stratum."), class = "model_data")
  }
  rates <- lapply(seq_len(ns), function(s) {
    dd <- fdm_demogdata(pmax(D[, , s], 0.5) / E[, , s], E[, , s], ages, periods, "s")
    if (isTRUE(p$smooth)) dd <- fdm_smooth(dd, ages$lower, p$smooth_args)
    as.matrix(dd$rate[[1]])
  })
  order1 <- min(p$order1, dim(D)[1] - 1)
  order2 <- min(p$order2, dim(D)[1] - 1)
  groups <- lapply(seq_len(max(group)), function(g) {
    use <- which(group == g & pool)
    if (!length(use)) use <- which(group == g)
    prod <- exp(Reduce(`+`, lapply(rates[use], log)) / length(use))
    Eg <- Reduce(`+`, lapply(use, function(s) E[, , s]))
    fit_p <- suppressWarnings(demography::fdm(fdm_demogdata(prod, Eg, ages, periods, "product"),
                                              series = "product", order = order1))
    list(fit = fit_p, pool = use)
  })
  ratios <- lapply(seq_len(ns), function(s) {
    g <- groups[[group[s]]]
    prod_obs <- exp(Reduce(`+`, lapply(rates[g$pool], log)) / length(g$pool))
    r <- rates[[s]] / prod_obs
    r[r > 1e9] <- NA
    suppressWarnings(demography::fdm(fdm_demogdata(r, E[, , s], ages, periods, "ratio"), series = "ratio",
                                     order = order2))
  })
  fitted <- array(NA_real_, dim(D))
  for (s in seq_len(ns)) {
    fitted[, , s] <- exp(as.matrix(groups[[group[s]]]$fit$fitted$y) + as.matrix(ratios[[s]]$fitted$y))
  }
  list(groups = groups, group = group, ratios = ratios, K = p$K, drange = p$drange, fitted = fitted)
}

# Evaluate quietly: the model searches of the forecast package print the
# errors of rejected candidate models through try()
quiet_try <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(expr), type = "message")
  out
}

# A demogdata object of one series
fdm_demogdata <- function(rates, pop, ages, periods, name) {
  dimnames(rates) <- dimnames(pop) <- list(ages$lower, periods)
  demography::demogdata(data = rates, pop = pop, ages = ages$lower, years = periods, type = "mortality",
                        label = "mortensor", name = name)
}

#' @export
simulate_model.mt_spec_coherent_fdm <- function(spec, fit, h, n_draws, ...) {
  ns <- length(fit$ratios)
  sim_p <- lapply(fit$groups, function(g) {
    fc <- quiet_try(forecast::forecast(g$fit, h = h, method = "arima"))
    stats::simulate(fc, nsim = n_draws)
  })
  sim_r <- lapply(fit$ratios, function(r) {
    ny <- length(r$year)
    r$weights <- 0 * r$weights + 1
    if (fit$K < ny) r$weights[seq_len(ny - fit$K)] <- 0
    fc <- quiet_try(forecast::forecast(r, h = h, method = "arfima", estim = "mle", drange = fit$drange))
    stats::simulate(fc, nsim = n_draws)
  })
  na <- dim(sim_p[[1]])[1]
  out <- array(NA_real_, c(na, h, ns, n_draws))
  for (g in seq_along(fit$groups)) {
    pool <- fit$groups[[g]]$pool
    # rescale the ratios of each path so that those of the group multiply to one
    adj <- exp(Reduce(`+`, lapply(sim_r[pool], function(a) log(array(a, c(na, h, n_draws))))) / length(pool))
    for (s in which(fit$group == g)) {
      out[, , s, ] <- array(sim_p[[g]], c(na, h, n_draws)) * array(sim_r[[s]], c(na, h, n_draws)) / adj
    }
  }
  out
}
