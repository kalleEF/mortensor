#' Tucker3 tensor model across strata
#'
#' Joint model of the array of log rates `[age, period, stratum]` (or of
#' centred log-ratio transformed life-table deaths, `transform = "coda"`)
#' through a Tucker3 decomposition `X = G x1 A x2 B x3 C` with ranks
#' `ranks = c(age, period, stratum)`, estimated by higher-order orthogonal
#' iteration (HOOI; De Lathauwer, De Moor and Vandewalle 2000). The period
#' factors (columns of `B`) are forecast as a multivariate random walk with
#' drift (bootstrap of joint year-to-year changes, which keeps their
#' correlation) or with ARIMA/ETS models per column (`ts = "arima"`,
#' `"ets"`, package \pkg{forecast}).
#'
#' With `transform = "coda"`, life-table deaths (from age `x0`, the first
#' grid age) form a composition over the age cells including the open
#' interval; the open interval's death rate is forecast as an auxiliary
#' random walk with drift on its log, so the open age group is kept.
#'
#' Ranks are reduced to the array dimensions where needed, so the model also
#' runs with one stratum.
#'
#' @param ranks Tucker ranks for age, period and stratum.
#' @param transform `"log"` or `"coda"`.
#' @param ts Forecast method for the period factors: `"rwd"`, `"arima"`,
#'   `"ets"`.
#' @param residual_noise Add the in-sample residual variation (per age and
#'   stratum, normal on the transformed scale) to the draws. It represents
#'   structure that the low-rank model does not capture. On the calibration
#'   benchmark of the validation pipeline (independent period trends in 10
#'   strata) 95% predictive intervals covered 0.95 on average with it and
#'   0.87 without it.
#' @param maxit,tol Settings of the HOOI iterations.
#' @return A model specification.
#' @references
#' De Lathauwer L, De Moor B, Vandewalle J (2000). On the best rank-1 and
#' rank-(R1, R2, ..., RN) approximation of higher-order tensors. *SIAM Journal
#' on Matrix Analysis and Applications* 21(4): 1324-1342.
#' \doi{10.1137/S0895479898346995}
#'
#' Kolda TG, Bader BW (2009). Tensor decompositions and applications. *SIAM
#' Review* 51(3): 455-500. \doi{10.1137/07070111X}
#' @examples
#' x <- mt_example(deaths = "single")
#' fit <- mt_fit(tucker3(ranks = c(3, 3, 2)), x)
#' fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
#' @export
tucker3 <- function(ranks = c(3, 3, 2), transform = c("log", "coda"), ts = c("rwd", "arima", "ets"),
                    residual_noise = TRUE, maxit = 200, tol = 1e-8) {
  transform <- match.arg(transform)
  ts <- match.arg(ts)
  if (ts != "rwd") need_pkg("forecast", sprintf("for ts = \"%s\"", ts))
  if (length(ranks) != 3) mt_abort("`ranks` needs three values.", class = "argument")
  new_mt_spec("tucker3", list(ranks = ranks, transform = transform, ts = ts, residual_noise = residual_noise,
                              maxit = maxit, tol = tol),
              joint = TRUE, min_periods = 6,
              label = sprintf("Tucker3-%s", if (transform == "log") "Direct" else "CoDa"))
}

#' @export
fit_model.mt_spec_tucker3 <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  m <- pmax(D, 0.5) / E
  X <- if (p$transform == "log") log(m) else coda_forward(m, ages)
  dims <- dim(X)
  ranks <- pmin(p$ranks, dims - c(1, 1, 0), dims)
  ranks <- pmax(ranks, 1)
  dec <- hooi(X, ranks, p$maxit, p$tol)
  dims <- dim(X)
  Xhat <- tucker_reconstruct(dec$G, dec$U)
  resid <- X - Xhat
  out <- list(dec = dec, ranks = ranks, transform = p$transform, ts = p$ts, ages = ages,
              resid_sd = apply(resid, c(1, 3), stats::sd), residual_noise = p$residual_noise,
              fit_pct = 100 * (1 - sum(resid^2) / sum((X - mean(X))^2)))
  B <- dec$U[[2]]
  out$B_last <- B[nrow(B), ]
  dB <- diff(B)
  out$B_mu <- colMeans(dB)
  out$B_resid <- t(dB) - out$B_mu
  if (p$ts != "rwd") out$B_models <- fit_series_models(B, p$ts)
  if (p$transform == "coda") {
    open <- is.infinite(ages$upper[length(ages)])
    if (open) {
      lo <- log(m[length(ages), , , drop = FALSE])
      lo <- matrix(lo, dims[2], dims[3])
      out$open_last <- lo[dims[2], ]
      dlo <- diff(lo)
      out$open_mu <- colMeans(dlo)
      out$open_resid <- t(dlo) - out$open_mu
    }
    out$fitted <- coda_backward(Xhat, ages, open_rates = if (open) m[length(ages), , ] else NULL)
  } else {
    out$fitted <- exp(Xhat)
  }
  out
}

#' @export
simulate_model.mt_spec_tucker3 <- function(spec, fit, h, n_draws, ...) {
  dec <- fit$dec
  na <- nrow(dec$U[[1]])
  ns <- nrow(dec$U[[3]])
  out <- array(NA_real_, c(length(fit$ages), h, ns, n_draws))
  Bp <- if (fit$ts == "rwd") {
    rw_paths(fit$B_last, fit$B_mu, fit$B_resid, h, n_draws, "bootstrap", exp_out = FALSE)
  } else {
    simulate_series_models(fit$B_models, h, n_draws)
  }
  op <- NULL
  if (!is.null(fit$open_last)) op <- rw_paths(fit$open_last, fit$open_mu, fit$open_resid, h, n_draws, "bootstrap")
  for (j in seq_len(n_draws)) {
    U <- dec$U
    U[[2]] <- t(matrix(Bp[, , j], ncol = h))
    Xf <- tucker_reconstruct(dec$G, U)
    if (fit$residual_noise) {
      sdarr <- aperm(array(fit$resid_sd, c(na, ns, h)), c(1, 3, 2))
      Xf <- Xf + stats::rnorm(length(Xf), 0, as.vector(sdarr))
    }
    out[, , , j] <- if (fit$transform == "log") exp(Xf) else
      coda_backward(Xf, fit$ages, open_rates = if (is.null(op)) NULL else t(matrix(op[, , j], ncol = h)))
  }
  out
}

# Higher-order orthogonal iteration for an N-way array
hooi <- function(X, ranks, maxit = 200, tol = 1e-8) {
  N <- length(dim(X))
  U <- lapply(seq_len(N), function(n) {
    s <- svd(unfold(X, n), nu = ranks[n], nv = 0)
    s$u[, seq_len(ranks[n]), drop = FALSE]
  })
  fit_old <- 0
  normX <- sum(X^2)
  for (it in seq_len(maxit)) {
    for (n in seq_len(N)) {
      Y <- X
      for (m in setdiff(seq_len(N), n)) Y <- mode_product(Y, t(U[[m]]), m)
      s <- svd(unfold(Y, n), nu = ranks[n], nv = 0)
      U[[n]] <- s$u[, seq_len(ranks[n]), drop = FALSE]
    }
    G <- X
    for (m in seq_len(N)) G <- mode_product(G, t(U[[m]]), m)
    fit <- sum(G^2) / normX
    if (abs(fit - fit_old) < tol) break
    fit_old <- fit
  }
  list(G = G, U = U, iterations = it)
}

unfold <- function(X, n) {
  d <- dim(X)
  perm <- c(n, setdiff(seq_along(d), n))
  matrix(aperm(X, perm), d[n])
}

mode_product <- function(X, M, n) {
  d <- dim(X)
  Y <- M %*% unfold(X, n)
  d2 <- d
  d2[n] <- nrow(M)
  perm <- c(n, setdiff(seq_along(d), n))
  aperm(array(Y, d2[perm]), order(perm))
}

tucker_reconstruct <- function(G, U) {
  X <- G
  for (m in seq_along(U)) X <- mode_product(X, U[[m]], m)
  X
}

# CoDa: rates -> life-table deaths (composition over ages) -> clr. For a
# closed age range the survivors at the upper bound are an extra component.
coda_forward <- function(m, ages) {
  d <- dim(m)
  M <- matrix(m, d[1])
  lt <- lt_columns(M, ages)
  dx <- lt$dx
  if (!is.infinite(ages$upper[d[1]])) dx <- rbind(dx, lt$lx[d[1], ] - lt$dx[d[1], ])
  z <- log(pmax(dx, 1e-12))
  z <- sweep(z, 2, colMeans(z))
  array(z, c(nrow(z), d[-1]))
}

# inverse clr -> life-table deaths -> rates; open interval rate from the auxiliary series
coda_backward <- function(Z, ages, open_rates = NULL) {
  d <- dim(Z)
  na <- length(ages)
  Zm <- matrix(Z, d[1])
  dx <- exp(Zm)
  dx <- sweep(dx, 2, colSums(dx), "/")
  dx <- dx[seq_len(na), , drop = FALSE]
  lx <- rbind(1, 1 - apply(dx, 2, cumsum))[seq_len(na), , drop = FALSE]
  n <- ages$upper - ages$lower
  open <- is.infinite(n[na])
  nn <- matrix(n, na, ncol(dx))
  Lx <- nn * (lx - dx / 2)
  m <- dx / Lx
  if (open) {
    if (is.null(open_rates)) mt_abort("The open interval needs an auxiliary rate series.", class = "model")
    m[na, ] <- as.vector(open_rates)
  }
  array(pmax(m, 1e-12), c(na, d[-1]))
}

# Per-column ARIMA or ETS models (forecast package) for the period factors
fit_series_models <- function(B, ts) {
  lapply(seq_len(ncol(B)), function(r) {
    y <- stats::ts(B[, r])
    if (ts == "arima") forecast::auto.arima(y) else forecast::ets(y, model = "ZZN")
  })
}

simulate_series_models <- function(models, h, n_draws) {
  R <- length(models)
  out <- array(NA_real_, c(R, h, n_draws))
  for (r in seq_len(R)) {
    for (j in seq_len(n_draws)) {
      out[r, , j] <- as.numeric(stats::simulate(models[[r]], nsim = h, future = TRUE, bootstrap = TRUE))
    }
  }
  out
}
