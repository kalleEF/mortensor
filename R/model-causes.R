#' Coherent multi-cause forecasts with cause shares
#'
#' Forecasts all-cause mortality with any model (`total`) and the shares of
#' the causes in all-cause mortality with a compositional tensor model, so
#' that cause-specific forecasts add up to the all-cause forecast in every
#' draw. This is the Tucker4-ILR method of mortensor 0.1:
#'
#' 1. shares of each cause (and of a residual cause, the part of all-cause
#'    mortality not covered by the listed causes) per age, year and group;
#' 2. isometric log-ratio (ILR) coordinates of the shares (orthonormal
#'    Helmert basis);
#' 3. Tucker decomposition of the array `[age, year, group, coordinate]`
#'    (HOOI) with ranks `ranks`; the year factors are forecast as a
#'    multivariate random walk with drift (bootstrap of joint changes) and
#'    the in-sample residual variation is added to the draws;
#' 4. inverse ILR gives share draws, which multiply the all-cause draws.
#'
#' Shares of rare causes are noisy; smoothing the rates first
#' (`prep = mt_prep(step_smooth())`) kept the intervals calibrated in the
#' multi-cause vignette, where they were too wide without it.
#'
#' The data need a stratum variable `var` with the all-cause level
#' `total_level` and the cause levels; other stratum variables (for example
#' sex) define the groups.
#'
#' @param total Model specification for all-cause mortality.
#' @param var Stratum variable with the causes.
#' @param total_level Level of `var` for all causes.
#' @param ranks Tucker ranks for age, year, group and ILR coordinate.
#' @param residual Add a residual cause so that listed causes may cover less
#'   than all deaths: `TRUE`, `FALSE`, or `"auto"` (added when the listed
#'   causes leave more than 0.1% of all-cause mortality unexplained in the
#'   median cell; a residual that is zero almost everywhere would only add a
#'   degenerate coordinate).
#' @param residual_noise Add in-sample residual variation of the ILR
#'   coordinates to the draws.
#' @return A model specification.
#' @references
#' Egozcue JJ, Pawlowsky-Glahn V, Mateu-Figueras G, Barcelo-Vidal C (2003).
#' Isometric logratio transformations for compositional data analysis.
#' *Mathematical Geology* 35(3): 279-300.
#' @examples
#' x <- mt_example_causes(years = 2010:2019)
#' fit <- mt_fit(cause_shares(total = lee_carter(), ranks = c(2, 2, 2, 2)), x)
#' fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
#' @export
cause_shares <- function(total = lee_carter(), var = "cause", total_level = "all", ranks = c(3, 3, 2, 2),
                         residual = "auto", residual_noise = TRUE) {
  check_class(total, "mt_spec", "total")
  if (length(ranks) != 4) mt_abort("`ranks` needs four values.", class = "argument")
  new_mt_spec("cause_shares", list(total = total, var = var, total_level = total_level, ranks = ranks,
                                   residual = residual, residual_noise = residual_noise),
              min_periods = max(total$min_periods, 6), label = sprintf("Cause shares (%s)", total$label))
}

#' @export
mt_fit.mt_spec_cause_shares <- function(spec, data, prep = NULL, cores = 1, ...) {
  check_class(data, "mt_data", "data")
  p <- spec$params
  d <- mt_prep_apply(prep, data)
  md <- model_inputs(d, spec)
  st <- d$strata
  if (!p$var %in% names(st)) mt_abort(sprintf("`%s` is not a stratum variable.", p$var), class = "argument")
  others <- setdiff(names(st), c("stratum", "label", p$var))
  gkey <- rep("all", nrow(st))
  if (length(others)) gkey <- do.call(paste, c(lapply(others, function(o) st[[o]]), sep = "."))
  groups <- unique(gkey)
  causes <- setdiff(unique(st[[p$var]]), p$total_level)
  if (!length(causes)) mt_abort("No cause levels besides the total.", class = "argument")
  tot_id <- vapply(groups, function(g) {
    i <- st$stratum[gkey == g & st[[p$var]] == p$total_level]
    if (length(i) != 1) mt_abort(sprintf("Group '%s' has no all-cause stratum.", g), class = "model_data")
    i
  }, 1L)
  find_part <- function(cc, g) {
    i <- st$stratum[gkey == g & st[[p$var]] == cc]
    if (length(i) != 1) mt_abort(sprintf("Group '%s' lacks cause '%s'.", g, cc), class = "model_data")
    i
  }
  part_id <- sapply(causes, function(cc) vapply(groups, function(g) find_part(cc, g), 1L))
  part_id <- matrix(part_id, nrow = length(groups))
  # all-cause model on the total strata only
  dt <- d
  dt$obs <- lapply(d$obs, function(o) {
    o <- o[o$stratum %in% tot_id]
    o[, stratum := match(stratum, tot_id)]
    o[]
  })
  dt$strata <- st[match(tot_id, st$stratum)][, stratum := seq_len(.N)][]
  fit_total <- mt_fit(p$total, dt, cores = cores)
  # shares
  m <- md$D / md$E
  na <- dim(m)[1]
  nt <- dim(m)[2]
  ng <- length(groups)
  K <- length(causes)
  S <- array(NA_real_, c(na, nt, ng, K + 1))
  for (g in seq_len(ng)) {
    tot <- pmax(m[, , tot_id[g]], 1e-12)
    for (k in seq_len(K)) S[, , g, k] <- m[, , part_id[g, k]] / tot
  }
  sum_parts <- apply(S[, , , seq_len(K), drop = FALSE], 1:3, sum)
  excess <- mean(sum_parts > 1 + 1e-9)
  use_residual <- if (identical(p$residual, "auto")) stats::median(1 - sum_parts) > 1e-3 else isTRUE(p$residual)
  if (use_residual) S[, , , K + 1] <- pmax(1 - sum_parts, 0) else S <- S[, , , seq_len(K), drop = FALSE]
  S <- pmax(S, 1e-6)
  S <- sweep(S, 1:3, apply(S, 1:3, sum), "/")
  V <- helmert_basis(dim(S)[4])
  Z <- ilr_array(S, V)
  ranks <- pmax(pmin(p$ranks, dim(Z) - c(1, 1, 0, 0), dim(Z)), 1)
  dec <- hooi(Z, ranks)
  Zhat <- tucker_reconstruct(dec$G, dec$U)
  resid <- Z - Zhat
  B <- dec$U[[2]]
  dB <- diff(B)
  share_fit <- list(dec = dec, V = V, resid_sd = apply(resid, c(1, 3, 4), stats::sd), B_last = B[nrow(B), ],
                    B_mu = colMeans(dB), B_resid = t(dB) - colMeans(dB), residual_noise = p$residual_noise)
  # in-sample fitted rates of all strata
  S_hat <- ilr_inverse_array(Zhat, V)
  tot_fitted <- total_fitted(fit_total)
  fitted <- array(NA_real_, dim(m))
  for (g in seq_len(ng)) {
    fitted[, , tot_id[g]] <- tot_fitted[, , g]
    for (k in seq_len(K)) fitted[, , part_id[g, k]] <- tot_fitted[, , g] * S_hat[, , g, k]
  }
  diagnostics <- rbind(cbind(data.table(component = "total"), fit_total$diagnostics),
                       data.table(component = "shares", stratum = "all", status = "ok",
                                  message = if (excess > 0) sprintf("causes exceed all-cause rate in %.1f%% of cells",
                                                                    100 * excess) else "",
                                  seconds = NA_real_), fill = TRUE)
  fits <- list(total = fit_total, shares = share_fit, tot_id = tot_id, part_id = part_id, K = K, fitted = fitted)
  structure(list(spec = spec, prep = prep, fits = fits, diagnostics = diagnostics, age = md$age,
                 periods = md$periods, strata = st, omega = d$omega, observed = m, train_D = md$D, train_E = md$E,
                 exposure_last = md$E[, nt, , drop = FALSE],
                 noise_sd = extra_noise_sd(list(smooth = TRUE, joint = TRUE), list(fitted = fitted), md$D, md$E)),
            class = c("mt_fit_cause_shares", "mt_fit"))
}

total_fitted <- function(fit) {
  if (fit$spec$joint) return(fit$fits$fitted)
  arr <- array(NA_real_, dim(fit$observed))
  for (s in seq_along(fit$fits)) {
    f <- fit$fits[[s]]$fitted
    arr[, , s] <- if (is.null(f)) fit$observed[, , s] else f
  }
  arr
}

#' @export
mt_forecast.mt_fit_cause_shares <- function(fit, h = 10, to = NULL, n_draws = 500, seed = NULL, predictive = FALSE,
                                            exposure = NULL, ...) {
  last <- max(fit$periods)
  if (!is.null(to)) h <- to - last
  f <- fit$fits
  sf <- f$shares
  arr <- with_seed(seed, {
    tot <- mt_forecast(f$total, h = h, n_draws = n_draws)$rates
    U <- sf$dec$U
    na <- nrow(U[[1]])
    ng <- nrow(U[[3]])
    nc <- nrow(U[[4]])
    Bp <- rw_paths(sf$B_last, sf$B_mu, sf$B_resid, h, n_draws, "bootstrap", exp_out = FALSE)
    out <- array(NA_real_, c(na, h, length(fit$strata$stratum), n_draws))
    for (j in seq_len(n_draws)) {
      Uj <- U
      Uj[[2]] <- t(matrix(Bp[, , j], ncol = h))
      Z <- tucker_reconstruct(sf$dec$G, Uj)
      if (sf$residual_noise) {
        sdarr <- aperm(array(sf$resid_sd, c(na, ng, nc, h)), c(1, 4, 2, 3))
        Z <- Z + stats::rnorm(length(Z), 0, as.vector(sdarr))
      }
      Sj <- ilr_inverse_array(Z, sf$V)
      for (g in seq_len(ng)) {
        out[, , f$tot_id[g], j] <- tot[, , g, j]
        for (k in seq_len(f$K)) out[, , f$part_id[g, k], j] <- tot[, , g, j] * Sj[, , g, k]
      }
    }
    if (predictive) out <- predictive_rates(out, fit, exposure, list(smooth = TRUE))
    out
  })
  res <- new_mt_draws(arr, fit$age, last + seq_len(h), fit$strata, model = fit$spec$label, origin = last,
                      fit = fit, seed = seed)
  res$predictive <- predictive
  res
}

# Orthonormal basis of the clr plane (normalised Helmert contrasts), D x (D - 1)
helmert_basis <- function(D) {
  H <- stats::contr.helmert(D)
  sweep(H, 2, sqrt(colSums(H^2)), "/")
}

ilr_array <- function(S, V) {
  d <- dim(S)
  L <- log(matrix(S, ncol = d[4]))
  clr <- L - rowMeans(L)
  array(clr %*% V, c(d[-4], ncol(V)))
}

ilr_inverse_array <- function(Z, V) {
  d <- dim(Z)
  clr <- matrix(Z, ncol = d[4]) %*% t(V)
  E <- exp(clr - apply(clr, 1, max))
  array(E / rowSums(E), c(d[-4], nrow(V)))
}

#' Simulated deaths by cause
#'
#' Extends [mt_example()] with deaths by cause (`"cvd"`, `"cancer"`,
#' `"other"`) and the all-cause total (`"all"`), single years of age. Cause
#' shares change smoothly with age and time; deaths are split multinomially,
#' so causes add up to all causes exactly.
#'
#' @inheritParams mt_example
#' @return An [mt_data] object with strata `sex` x `cause`.
#' @examples
#' mt_example_causes(years = 2010:2019)
#' @export
mt_example_causes <- function(years = 1990:2019, size = 1e5, seed = 1) {
  x <- mt_example(years = years, deaths = "single", size = size, seed = seed)
  D <- mt_obs(x, "deaths")
  E <- mt_obs(x, "exposure")
  with_seed(seed + 1, {
    t0 <- (D$period_lo - min(years)) / 10
    a <- (pmin(D$age_lo, 100) - 30) / 10
    eta_cvd <- -1.5 + 0.35 * a - 0.25 * t0
    eta_can <- 0.3 - 0.15 * a + 0.05 * t0
    p <- cbind(exp(eta_cvd), exp(eta_can), 1)
    p <- p / rowSums(p)
    split <- t(vapply(seq_len(nrow(D)), function(i) as.numeric(stats::rmultinom(1, D$value[i], p[i, ])), numeric(3)))
  })
  causes <- c("cvd", "cancer", "other")
  base <- D[, list(sex, period_lo, age_lo, age_hi)]
  long <- rbind(cbind(base, cause = "all", deaths = D$value),
                rbindlist(lapply(1:3, function(k) cbind(base, cause = causes[k], deaths = split[, k]))))
  ex <- E[, list(sex, period_lo, age_lo, age_hi, exposure = value)]
  ex <- rbindlist(lapply(c("all", causes), function(cc) cbind(ex, cause = cc)))
  mt_data(deaths = long, exposure = ex, by = c("sex", "cause"), age = c("age_lo", "age_hi"), period = "period_lo")
}
