#' Bootstrap uncertainty of preprocessing and parameters
#'
#' Wraps a model so that its forecast draws include the uncertainty of the
#' estimated parameters and, optionally, of the preprocessing (smoothing,
#' ungrouping). The model is refitted on `n_boot` resampled datasets and the
#' forecast draws are pooled across refits:
#'
#' * `type = "grid"`: deaths on the grid (after preprocessing) are drawn
#'   from Poisson distributions with the fitted rates of the model (or the
#'   preprocessed rates when the model has no fitted rates); the model is
#'   refitted without preprocessing. This is the semi-parametric bootstrap of
#'   Brouhns et al. (2005), generalised to any model.
#' * `type = "data"`: the observed deaths, at their original resolution, are
#'   drawn from Poisson distributions with the observed counts as means, and
#'   preprocessing and model are refitted. This propagates the uncertainty of
#'   ungrouping and smoothing, for any smoother.
#'
#' @param spec A model specification.
#' @param n_boot Number of bootstrap refits.
#' @param type `"grid"` or `"data"`.
#' @param seed Random seed for the resampling.
#' @return A model specification.
#' @references
#' Brouhns N, Denuit M, Van Keilegom I (2005). Bootstrapping the Poisson
#' log-bilinear model for mortality forecasting. *Scandinavian Actuarial
#' Journal* 2005(3): 212-224.
#' @examples
#' x <- mt_example(deaths = "single", years = 2000:2019)
#' fit <- mt_fit(bootstrapped(rwd(), n_boot = 5), x)
#' fc <- mt_forecast(fit, h = 3, n_draws = 50, seed = 1)
#' @export
bootstrapped <- function(spec, n_boot = 50, type = c("grid", "data"), seed = 1) {
  check_class(spec, "mt_spec", "spec")
  type <- match.arg(type)
  check_number(n_boot, "n_boot", lower = 1, integerish = TRUE)
  new_mt_spec("bootstrap", list(base = spec, n_boot = n_boot, type = type, seed = seed),
              joint = spec$joint, min_periods = spec$min_periods, label = sprintf("%s, bootstrap", spec$label),
              smooth = spec$smooth)
}

#' @export
mt_fit.mt_spec_bootstrap <- function(spec, data, prep = NULL, cores = 1, ...) {
  p <- spec$params
  base <- mt_fit(p$base, data, prep = prep, cores = cores, ...)
  fits <- with_seed(p$seed, lapply(seq_len(p$n_boot), function(i) {
    bd <- if (p$type == "grid") resample_grid(base) else resample_data(data)
    tryCatch(mt_fit(p$base, bd, prep = if (p$type == "data") prep else NULL, cores = cores, ...),
             error = function(e) NULL)
  }))
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) mt_abort("Every bootstrap refit failed.", class = "model_fit")
  base$diagnostics <- rbind(base$diagnostics,
                            data.table(stratum = "all", status = "ok",
                                       message = sprintf("%d of %d bootstrap refits succeeded", length(fits),
                                                         p$n_boot), seconds = NA_real_), fill = TRUE)
  # the fitted base spec goes on the stack of inner specs (mt_parameters(),
  # observation noise), as for the other wrappers
  wrap_fit(base, spec, "mt_fit_bootstrap", list(boot = fits))
}

resample_grid <- function(fit) {
  m <- fit$observed
  if (isTRUE(fit$spec$smooth)) {
    # fitted rates of any fit (ensembles: weighted geometric mean of the components)
    f <- tryCatch(stats::fitted(fit), error = function(e) NULL)
    for (s in seq_len(dim(m)[3])) {
      fs <- if (identical(dim(f), dim(m))) f[, , s]
      if (!is.null(fs) && !anyNA(fs)) m[, , s] <- fs
    }
  }
  D <- array(stats::rpois(length(m), m * fit$train_E), dim(m))
  arrays_to_mt_data(D, fit$train_E, fit$age, fit$periods, fit$strata, fit$omega)
}

resample_data <- function(data) {
  data$obs$deaths <- copy(data$obs$deaths)
  data$obs$deaths[, value := stats::rpois(.N, value)]
  data
}

arrays_to_mt_data <- function(D, E, age, periods, strata, omega) {
  d <- dim(D)
  long <- data.table(stratum = rep(seq_len(d[3]), each = d[1] * d[2]), age_lo = rep(age$lower, d[2] * d[3]),
                     age_hi = rep(age$upper, d[2] * d[3]), period_lo = rep(rep(periods, each = d[1]), d[3]),
                     period_hi = rep(rep(periods + 1, each = d[1]), d[3]), deaths = as.vector(D),
                     exposure = as.vector(E))
  by <- setdiff(names(strata), c("stratum", "label"))
  long <- merge(strata[, c("stratum", by), with = FALSE], long, by = "stratum")
  if (!length(by)) long[, stratum := NULL]
  x <- suppressWarnings(mt_data(deaths = long, exposure = long, by = if (length(by)) by else NULL,
                                age = c("age_lo", "age_hi"), period = c("period_lo", "period_hi"), omega = omega))
  mt_set_grid(x, age = age)
}

#' @export
mt_forecast.mt_fit_bootstrap <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                                         exposure = NULL, ...) {
  boots <- fit$boot
  nb <- length(boots)
  alloc <- tabulate(rep_len(seq_len(nb), n_draws), nb)
  out <- with_seed(seed, {
    parts <- lapply(seq_len(nb), function(i) {
      if (alloc[i] == 0) return(NULL)
      bf <- boots[[i]]
      mt_forecast(bf, h = if (is.null(to)) h else to - max(bf$periods), n_draws = alloc[i], predictive = predictive,
                  exposure = exposure)$rates
    })
    abind_draws(Filter(Negate(is.null), parts))
  })
  last <- max(fit$periods)
  res <- new_mt_draws(out, fit$age, last + seq_len(dim(out)[2]), fit$strata, model = fit$spec$label, origin = last,
                      fit = fit, seed = seed)
  res$predictive <- predictive
  res
}

abind_draws <- function(parts) {
  d <- dim(parts[[1]])
  n <- sum(vapply(parts, function(p) dim(p)[4], 0))
  out <- array(NA_real_, c(d[1:3], n))
  k <- 0
  for (p in parts) {
    j <- dim(p)[4]
    out[, , , k + seq_len(j)] <- p
    k <- k + j
  }
  out
}
