#' Bayesian stochastic mortality models in Stan (StanMoMo)
#'
#' Wrapper around \pkg{StanMoMo} (Barigou, Goffard, Loisel and Salhi 2023):
#' Bayesian versions of the Lee-Carter (`"lc"`), Renshaw-Haberman (`"rh"`),
#' age-period-cohort (`"apc"`), Cairns-Blake-Dowd (`"cbd"`) and M6 (`"m6"`)
#' models, fitted with Stan to the deaths and exposures of each stratum
#' with a Poisson or negative binomial likelihood (log link). Period
#' indexes follow random walks with drift (CBD and M6: correlated), cohort
#' effects AR(2) processes.
#'
#' StanMoMo forecasts inside the fit, so the model is fitted with forecasts
#' up to `h_max` years. Draws are posterior draws of the expected rates,
#' rebuilt from the posterior parameters and forecast indexes; they carry
#' the uncertainty of the parameters and of the future indexes. Fitted rates
#' are posterior means.
#'
#' Predictive draws (`mt_forecast(predictive = TRUE)`) are draws of the
#' model's posterior predictive distribution of the deaths, divided by the
#' exposures. For `family = "poisson"` the deaths are Poisson given the rates
#' of each posterior draw and the exposures. For `family = "nb"` they are
#' negative binomial (Delwarde, Denuit and Partrat 2007) with mean rate x
#' exposure and size `phi`, the dispersion parameter of the same posterior
#' draw: variance `mu + mu^2 / phi`, Stan's `neg_binomial_2` parameterisation
#' that StanMoMo uses. StanMoMo samples `aux`, the inverse dispersion, with a
#' half-normal prior, and reports `phi = 1 / aux`. No other observation noise
#' is added. If intervals are too narrow on real data, that is a property of
#' the model; compare the negative binomial family, [recalibrated()] or an
#' [ensemble()]. Wrappers that add the observation noise themselves
#' ([recalibrated()], [anchored()], [reconciled()]) pair each draw with the
#' dispersion of its posterior draw where the draws keep their order, and
#' otherwise draw the dispersion from its posterior.
#'
#' The fit warns when the chains have not converged (largest R-hat of the
#' model parameters above 1.05) or had divergent transitions; then use more
#' `iter` or `chains`. Cohort models need single years of age. StanMoMo
#' needs whole numbers: deaths and exposures are rounded (a warning reports
#' rounding of exposures by more than 1%). Fits take from seconds to
#' minutes per stratum; `cores` runs the chains in parallel.
#'
#' @param model `"lc"`, `"rh"`, `"apc"`, `"cbd"` or `"m6"`.
#' @param family `"poisson"` or `"nb"` (negative binomial).
#' @param h_max Longest forecast horizon the fit supports (at least 2).
#' @param chains,iter Number of Markov chains and iterations per chain (half
#'   of them warm-up).
#' @param cores Cores for the chains.
#' @param seed Seed of the sampler; `NULL` draws one from R's random number
#'   generator, so that `set.seed()` makes fits reproducible.
#' @return A model specification.
#' @references
#' Barigou K, Goffard P-O, Loisel S, Salhi Y (2023). Bayesian model
#' averaging for mortality forecasting using leave-future-out validation.
#' *International Journal of Forecasting* 39(2): 674-690.
#' \doi{10.1016/j.ijforecast.2022.01.011}
#'
#' Delwarde A, Denuit M, Partrat C (2007). Negative binomial version of the
#' Lee-Carter model for mortality forecasting. *Applied Stochastic Models in
#' Business and Industry* 23(5): 385-401.
#' @seealso [stmomo()] for the same models fitted by maximum likelihood.
#' @examples
#' \donttest{
#' if (requireNamespace("StanMoMo", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single", years = 2000:2019)
#'   fit <- mt_fit(stanmomo("lc", chains = 1, iter = 400, h_max = 10), x)
#'   fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
#' }
#' }
#' @export
stanmomo <- function(model = c("lc", "rh", "apc", "cbd", "m6"), family = c("poisson", "nb"), h_max = 50,
                     chains = 4, iter = 2000, cores = 1, seed = NULL) {
  need_pkg("StanMoMo", "for stanmomo()")
  need_pkg("rstan", "for stanmomo()")
  model <- match.arg(model)
  family <- match.arg(family)
  check_number(h_max, "h_max", lower = 2, integerish = TRUE)
  check_number(chains, "chains", lower = 1, integerish = TRUE)
  check_number(iter, "iter", lower = 20, integerish = TRUE)
  names <- c(lc = "Lee-Carter", rh = "Renshaw-Haberman", apc = "APC", cbd = "CBD", m6 = "M6")
  new_mt_spec("stanmomo", list(model = model, family = family, h_max = h_max, chains = chains, iter = iter,
                               cores = cores, seed = seed),
              min_periods = 8,
              label = sprintf("StanMoMo %s%s", names[[model]], if (family == "nb") " (negative binomial)" else ""))
}

#' @export
fit_model.mt_spec_stanmomo <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  cohort <- p$model %in% c("rh", "apc", "m6")
  if (cohort && any(utils::head(ages$upper - ages$lower, -1) != 1)) {
    mt_abort("StanMoMo cohort models need single years of age.", class = "model_data")
  }
  Ei <- pmax(round(E), 1)
  change <- max(abs(Ei / E - 1))
  if (change > 0.01) {
    warning(sprintf("exposures rounded to whole numbers for StanMoMo (largest relative change %.1f%%).",
                    100 * change), call. = FALSE)
  }
  args <- list(death = round(D), exposure = Ei, forecast = p$h_max, family = p$family, chains = p$chains,
               iter = p$iter, cores = p$cores, refresh = 0,
               seed = p$seed %||% sample.int(.Machine$integer.max, 1))
  if (p$model %in% c("cbd", "m6")) args$age <- ages$lower
  fit <- do.call(getExportedValue("StanMoMo", paste0(p$model, "_stan")), args)
  if (!identical(fit@mode, 0L)) mt_abort("StanMoMo sampling failed; see the messages of rstan.", class = "model_fit")
  pars <- switch(p$model, lc = c("a", "b", "k", "k_p"), rh = c("a", "b", "k", "g", "k_p", "gf"),
                 apc = c("a", "k", "g", "k_p", "gf"), cbd = c("k", "k2", "k_p", "k2_p"),
                 m6 = c("k", "k2", "g", "k_p", "k2_p", "gf"))
  # negative binomial: StanMoMo's phi = 1 / aux, the size of neg_binomial_2
  if (p$family == "nb") pars <- c(pars, "phi")
  post <- rstan::extract(fit, pars = pars)
  core <- setdiff(pars, c("k_p", "k2_p", "gf"))
  rhat <- rstan::summary(fit, pars = core)$summary[, "Rhat"]
  rhat_max <- max(rhat, na.rm = TRUE)
  divergent <- rstan::get_num_divergent(fit)
  if (rhat_max > 1.05) {
    warning(sprintf("StanMoMo chains have not converged (largest R-hat %.2f); use more `iter` or `chains`.",
                    rhat_max), call. = FALSE)
  }
  if (divergent > 0) warning(sprintf("StanMoMo: %d divergent transitions.", divergent), call. = FALSE)
  na <- nrow(D)
  nt <- ncol(D)
  # fitted rates: posterior mean of the expected rates over up to 500 draws
  n_post <- nrow(post$k)
  use <- unique(round(seq(1, n_post, length.out = min(500, n_post))))
  eta <- lapply(use, function(i) stanmomo_eta(p$model, post, i, ages$lower, na, nt, fitted = TRUE))
  fitted <- Reduce(`+`, lapply(eta, exp)) / length(use)
  list(post = post, model = p$model, age = ages$lower, nt = nt, h_max = p$h_max, rhat = rhat_max,
       divergent = divergent, fitted = fitted, nb_size = if (p$family == "nb") as.numeric(post$phi))
}

# Linear predictor of one posterior draw i: for the training years
# (`fitted = TRUE`) or for the forecast years 1..h_max
stanmomo_eta <- function(model, post, i, age, na, nt, fitted = FALSE, h = NULL) {
  nh <- if (fitted) nt else h
  kk <- if (fitted) post$k[i, ] else post$k_p[i, seq_len(nh)]
  tt <- if (fitted) seq_len(nt) else nt + seq_len(nh)
  # cohort index of age x and year t: t - x + J (years counted from the first training year)
  coh <- function() {
    gg <- if (fitted) post$g[i, ] else post$gf[i, ]
    matrix(gg[outer(-seq_len(na), tt, `+`) + na], na, nh)
  }
  switch(model,
         lc = post$a[i, ] + outer(post$b[i, ], kk),
         rh = post$a[i, ] + outer(post$b[i, ], kk) + coh(),
         apc = post$a[i, ] + matrix(kk, na, nh, byrow = TRUE) + coh(),
         cbd = , m6 = {
           k2 <- if (fitted) post$k2[i, ] else post$k2_p[i, seq_len(nh)]
           eta <- matrix(kk, na, nh, byrow = TRUE) + outer(age - mean(age), k2)
           if (model == "m6") eta + coh() else eta
         })
}

#' @export
simulate_model.mt_spec_stanmomo <- function(spec, fit, h, n_draws, ...) {
  if (h > fit$h_max) {
    mt_abort(sprintf("The StanMoMo fit forecasts up to %d years; refit with a larger `h_max`.", fit$h_max),
             class = "argument")
  }
  na <- length(fit$age)
  available <- nrow(fit$post$k_p)
  idx <- if (n_draws <= available) sample.int(available, n_draws) else sample.int(available, n_draws, replace = TRUE)
  out <- array(NA_real_, c(na, h, n_draws))
  for (j in seq_len(n_draws)) {
    out[, , j] <- exp(stanmomo_eta(fit$model, fit$post, idx[j], fit$age, na, fit$nt, h = h))
  }
  # negative binomial family: the dispersion of each posterior draw, for
  # predictive draws of the deaths (see predictive_rates())
  if (!is.null(fit$nb_size)) attr(out, "nb_size") <- fit$nb_size[idx]
  out
}
