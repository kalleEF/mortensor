#' Tucker3 tensor model across strata (three-way Lee-Carter model)
#'
#' Joint model of the array of log rates `[age, period, stratum]`: the
#' three-way extension of the Lee-Carter model of Russolillo, Giordano and
#' Haberman (2011). Status: published, with the deviations listed below;
#' section and equation numbers refer to that paper.
#'
#' **Model (Section 2, eq. 7-8).** As in Lee-Carter, the mean log rate over
#' the years of each age and stratum, `a[x, i]`, is removed (`centre =
#' TRUE`), and the centred array is decomposed by a Tucker3 decomposition
#' `X - a = G x1 A x2 B x3 C` with ranks `ranks = c(age, period, stratum)`:
#' column-wise orthonormal age factors `A`, period factors `B` and stratum
#' factors `C`, linked by the core array `G`, estimated by least squares.
#' Ranks `c(1, 1, 1)` give eq. (7), `log m[x, t, i] - a[x, i] = b[x] k[t]
#' c[i]` (one period index for all strata); the default `c(2, 2, 2)` is the
#' model of the application (Section 3.2, ten European countries). The
#' squared core elements over the sum of squares of the centred array are
#' the shares of variation of each combination of components (Table 4);
#' their sum is the fit (`fit$fits$fit_pct`, `core_pct`).
#'
#' **Forecasts.** The paper proposes a time-series model for the period
#' index, as in Lee-Carter, and leaves its form open (Section 3.2.2). Here
#' the period factors (columns of `B`) follow a multivariate random walk
#' with drift, as `k[t]` of Lee and Carter (1992); draws resample the joint
#' year-to-year changes, which keeps their correlation, and include the
#' uncertainty of the drift (`ts = "rwd"`). ARIMA or ETS models per column
#' (`ts = "arima"`, `"ets"`, package \pkg{forecast}) are alternatives. The
#' rates are rebuilt with `a`, `A`, `C` and `G` fixed. Draws carry the
#' uncertainty of the forecast period factors only. With fewer period
#' factors than independent trends in the data, the strata share their
#' forecast uncertainty and intervals are too narrow: in an example on
#' simulated Lee-Carter data with Poisson deaths and an own random walk in
#' each of 10 strata, 95% intervals covered about 74% with the default
#' ranks and 95% with `ranks = c(1, 10, 10)`.
#' Predictive draws (`mt_forecast(predictive = TRUE)`) add Poisson variation
#' of the deaths given the rates and the exposures, and no other observation
#' noise. If intervals are too narrow on real data, that is a property of
#' the model (see [recalibrated()], [ensemble()], or overdispersed models
#' such as `stanmomo(family = "nb")`).
#'
#' **Deviations from the paper.**
#'
#' * Estimation: the paper fits the model by alternating least squares
#'   (footnote 3; Kiers' Tucker3 program for Matlab, Appendix); here by
#'   higher-order orthogonal iteration (HOOI; De Lathauwer, De Moor and
#'   Vandewalle 2000), an alternating least-squares algorithm for the same
#'   criterion. The paper also rotates the solution to simple structure,
#'   which leaves the fitted values unchanged; signs are fixed so that the
#'   largest entry of each factor is positive.
#' * Centring: the text and eq. (7) centre over the years; the Matlab session
#'   in the Appendix also centres across the age groups. The package follows
#'   eq. (7), because the term removed by centring across ages varies over
#'   the years and the paper gives no forecast for it.
#' * Package-specific (reference unknown): the forecast model of the period
#'   factors (see above), the choice of the ranks by an information
#'   criterion (`ranks = "bic"`, `"aic"`), and the replacement of zero deaths
#'   by 0.5 before taking logarithms. The paper's data end with the closed
#'   age group 95-99; an open age interval is treated like the other age
#'   groups.
#'
#' **Deprecated:** `transform = "coda"` (centred log-ratio transformed
#' life-table deaths in place of log rates, with an auxiliary random walk
#' for the rate of the open age interval) is a combination of this package
#' without a published source; it still works and warns (class
#' `mt_warning_deprecated`). Use [coda()] (`coherent = TRUE` for several
#' strata), the published compositional model.
#'
#' Ranks are reduced to the array dimensions where needed, so the model also
#' runs with one stratum or one age group. The model is experimental: compare
#' it with established models such as [li_lee()] in a backtest before using
#' it.
#'
#' @param ranks Tucker ranks for age, period and stratum (default `c(2, 2,
#'   2)`, Russolillo et al. 2011, Section 3.2), or `"bic"` / `"aic"` to
#'   choose them by an information criterion on the training data
#'   (package-specific: `n log(RSS / n)` plus the penalty times the number
#'   of parameters of the factor matrices and the core; ranks up to 6 for
#'   age and period and up to 4 for strata). For selection by forecast
#'   accuracy use [tuned()].
#' @param transform `"log"`, or the deprecated `"coda"`.
#' @param ts Forecast method for the period factors: `"rwd"`, `"arima"`,
#'   `"ets"`.
#' @param residual_noise Add the in-sample residual variation (per age and
#'   stratum, normal on the transformed scale) to the draws. The residuals
#'   contain the observation noise of the training years, which predictive
#'   draws add again as Poisson variation, so the noise is then counted
#'   twice. Default `FALSE` (earlier development versions: `TRUE`).
#' @param centre Remove the mean over the years of each age and stratum
#'   before the decomposition and add it back afterwards (Russolillo,
#'   Giordano and Haberman 2011). With `FALSE` the decomposition describes
#'   the log rates themselves.
#' @param maxit,tol Settings of the HOOI iterations.
#' @return A model specification.
#' @references
#' Russolillo M, Giordano G, Haberman S (2011). Extending the Lee-Carter
#' model: a three-way decomposition. *Scandinavian Actuarial Journal*
#' 2011(2): 96-117. \doi{10.1080/03461231003611933}
#'
#' Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
#' *Journal of the American Statistical Association* 87(419): 659-671.
#'
#' De Lathauwer L, De Moor B, Vandewalle J (2000). On the best rank-1 and
#' rank-(R1, R2, ..., RN) approximation of higher-order tensors. *SIAM Journal
#' on Matrix Analysis and Applications* 21(4): 1324-1342.
#' \doi{10.1137/S0895479898346995}
#'
#' Kolda TG, Bader BW (2009). Tensor decompositions and applications. *SIAM
#' Review* 51(3): 455-500. \doi{10.1137/07070111X}
#' @examples
#' x <- mt_example(deaths = "single")
#' fit <- mt_fit(tucker3(), x)
#' fit$fits$fit_pct
#' fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
#' @export
tucker3 <- function(ranks = c(2, 2, 2), transform = c("log", "coda"), ts = c("rwd", "arima", "ets"),
                    residual_noise = FALSE, centre = TRUE, maxit = 200, tol = 1e-8) {
  transform <- match.arg(transform)
  ts <- match.arg(ts)
  check_flag(residual_noise, "residual_noise")
  check_flag(centre, "centre")
  if (transform == "coda") {
    mt_warn(paste("tucker3(transform = \"coda\") is deprecated: this combination of a Tucker3 decomposition with",
                  "compositional life-table deaths has no published source. Use coda() (coherent = TRUE for",
                  "several strata)."), class = "deprecated")
  }
  if (ts != "rwd") need_pkg("forecast", sprintf("for ts = \"%s\"", ts))
  if (is.character(ranks)) {
    check_choice(ranks, c("bic", "aic"), "ranks")
  } else if (length(ranks) != 3) {
    mt_abort("`ranks` needs three values or \"bic\"/\"aic\".", class = "argument")
  }
  new_mt_spec("tucker3", list(ranks = ranks, transform = transform, ts = ts, residual_noise = residual_noise,
                              centre = centre, maxit = maxit, tol = tol),
              joint = TRUE, min_periods = 6,
              label = sprintf("Tucker3-%s%s", if (transform == "log") "Direct" else "CoDa",
                              if (is.character(ranks)) sprintf(" (%s ranks)", toupper(ranks)) else ""))
}

#' @export
fit_model.mt_spec_tucker3 <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  m <- pmax(D, 0.5) / E
  X <- if (p$transform == "log") log(m) else coda_forward(m, ages)
  dims <- dim(X)
  # centring (Russolillo, Giordano and Haberman 2011): mean over the years of
  # each age and stratum, as a[x] of Lee-Carter
  alpha <- if (isTRUE(p$centre)) apply(X, c(1, 3), mean) else NULL
  Xc <- if (is.null(alpha)) X else sweep(X, c(1, 3), alpha)
  rk <- if (is.character(p$ranks)) tucker_select_ranks(Xc, p$ranks)$ranks else p$ranks
  ranks <- pmin(rk, dims - c(1, 1, 0), dims)
  ranks <- pmax(ranks, 1)
  dec <- hooi(Xc, ranks, p$maxit, p$tol)
  Xhat <- tucker_add_centre(tucker_reconstruct(dec$G, dec$U), alpha)
  resid <- X - Xhat
  # shares of the sum of squares of the decomposed array (Russolillo et al.
  # 2011, Table 4): squared core elements; their sum is the fit
  ss <- sum(Xc^2)
  out <- list(dec = dec, ranks = ranks, transform = p$transform, ts = p$ts, ages = ages, alpha = alpha,
              resid_sd = apply(resid, c(1, 3), stats::sd), residual_noise = isTRUE(p$residual_noise),
              core_pct = 100 * dec$G^2 / ss, fit_pct = 100 * sum(dec$G^2) / ss)
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
    Xf <- tucker_add_centre(tucker_reconstruct(dec$G, U), fit$alpha)
    if (fit$residual_noise) {
      sdarr <- aperm(array(fit$resid_sd, c(na, ns, h)), c(1, 3, 2))
      Xf <- Xf + stats::rnorm(length(Xf), 0, as.vector(sdarr))
    }
    out[, , , j] <- if (fit$transform == "log") exp(Xf) else
      coda_backward(Xf, fit$ages, open_rates = if (is.null(op)) NULL else t(matrix(op[, , j], ncol = h)))
  }
  out
}

# Add the means a[x, i] of the centred model to an array [age, period, stratum]
tucker_add_centre <- function(X, alpha) {
  if (is.null(alpha)) X else sweep(X, c(1, 3), alpha, "+")
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
  # SVD signs depend on the LAPACK build; fix them so that the entry with the
  # largest absolute value of every factor is positive (the core absorbs the
  # signs, so the reconstruction is unchanged)
  for (n in seq_len(N)) {
    sg <- apply(U[[n]], 2, function(u) if (u[which.max(abs(u))] < 0) -1 else 1)
    if (any(sg < 0)) {
      U[[n]] <- sweep(U[[n]], 2, sg, "*")
      G <- mode_product(G, diag(sg, length(sg)), n)
    }
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

# Rank selection by information criterion (mortensor 0.1, tucker3_info_criterion)
tucker_select_ranks <- function(X, criterion = "bic", max_ranks = c(6, 6, 4)) {
  dims <- dim(X)
  grid <- expand.grid(a = seq_len(min(max_ranks[1], dims[1] - 1)), b = seq_len(min(max_ranks[2], dims[2] - 1)),
                      c = seq_len(min(max_ranks[3], dims[3])))
  n <- prod(dims)
  ic <- apply(grid, 1, function(r) {
    if (r[3] > r[1] * r[2] || r[1] > r[2] * r[3] || r[2] > r[1] * r[3]) return(Inf)
    dec <- hooi(X, r, maxit = 50, tol = 1e-6)
    rss <- sum((X - tucker_reconstruct(dec$G, dec$U))^2)
    k <- sum(r * dims) + prod(r)
    n * log(rss / n) + (if (criterion == "bic") log(n) else 2) * k
  })
  i <- which.min(ic)
  list(ranks = as.numeric(grid[i, ]), table = cbind(grid, ic = ic))
}
