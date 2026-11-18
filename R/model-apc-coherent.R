# Coherent joint age-period-cohort model: the native sampler behind
# apc_bayes(coherent = TRUE).
#
# Rewrite of bamp_coherent() and predict_coherent() of the development
# branch "coherent_forecasting" of bamp (version 2.2.0 in its DESCRIPTION,
# commit 4363a78 of 2026-07-13; file R/coherent_joint.R), GPL-3. bamp:
# Volker Schmid (author, maintainer), Florian Geressen, Leonhard Held, Evi
# Rainer and Chris Kypridemos (contributors); the coherent models of the
# branch are by Chris Kypridemos. When a bamp
# release on CRAN contains bamp_coherent(), apc_bayes(coherent = TRUE) is to
# call it and this sampler is to be deprecated (dev/coherence-plan.md).
#
# The model is that of the branch for any number of strata S >= 2:
#   logit p[s, i, j] = mu + a[s] + theta[s, i] + phi[j] + psi[k] + d[s, j]
#                      (+ c[s, k] with a cohort deviation),
# binomial deaths given the exposures, shared period phi and cohort psi
# effects with random-walk priors scaled to unit generalised variance
# (Sørbye and Rue 2014), stratum age effects theta[s, ] with one random-walk
# precision, and stratum deviations d[s, ] from the shared period effect,
# stationary AR(1) series with the same rho and innovation precision.
#
# Changes from the branch:
# * Sampler: Metropolis-Hastings with a Gaussian (GMRF) approximation of the
#   conditional posterior of all effects jointly (one Newton step from the
#   current state; Rue and Held 2005, Section 4.4), with sparse Cholesky
#   factors, instead of Polya-Gamma Gibbs steps whose Polya-Gamma draws are
#   normal approximations. The chain targets the exact posterior.
# * The same parameterisation for every S: the deviations sum to zero over
#   the strata in every period (the "equal populations" variant of Cairns et
#   al. 2011, Section 4.5, extended to S strata), or the deviation of a
#   reference stratum is zero (their "population 1 dominant" variant,
#   Section 4.4). The branch uses a time constraint for S = 2 and the sum
#   over strata for S >= 3.
# * rho is sampled for every S (the branch: S = 2 only, fixed for S >= 3),
#   as is the AR(1) coefficient of the optional cohort deviation (the
#   branch: fixed, S = 2 only).
# * The zero-slope constraint on the period effect is imposed only when all
#   three effects are second-order random walks, where the linear trends are
#   not identified; the branch imposes it whenever the period effect is RW2.
# * Exposures are used as they are (no rounding); cells without exposure are
#   left out of the likelihood.
# * The precisions of the deviations have Gamma(1, 5e-4) priors, as the
#   period and cohort precisions; the branch's Gamma(1, 0.05) dominates the
#   data and sets the width of the forecast gaps.
# * A vague normal prior (precision 1e-6) on every effect, as the ridge of
#   the branch; the precision matrix is made positive definite in the
#   directions removed by the constraints by adding A'A, which leaves the
#   densities on the constrained subspace unchanged.

# Structure matrix of an intrinsic random walk of order `order` on `n`
# points, scaled to unit geometric mean of the marginal variances of its
# generalised inverse (Sørbye and Rue 2014)
coh_rw <- function(n, order) {
  Dm <- diff(diag(n), differences = order)
  K <- crossprod(Dm)
  e <- eigen(K, symmetric = TRUE)
  keep <- e$values > max(e$values) * 1e-9
  V <- e$vectors[, keep, drop = FALSE]
  scale <- exp(mean(log(rowSums(V^2 %*% diag(1 / e$values[keep], sum(keep))))))
  list(K = Matrix::Matrix(K * scale, sparse = TRUE), scale = scale, rank = n - order, order = order)
}

# Quadratic form of the AR(1) precision for the rows of Y (one series per
# row)
coh_ar1_quad <- function(Y, rho) {
  n <- ncol(Y)
  sum((Y[, -1, drop = FALSE] - rho * Y[, -n, drop = FALSE])^2) + (1 - rho^2) * sum(Y[, 1]^2)
}

coh_softplus <- function(z) pmax(z, 0) + log1p(exp(-abs(z)))

# Default hyperparameters: Gamma(shape, rate) priors of the precisions and a
# uniform prior of rho on (0, 1). Age, period and cohort as in
# bamp_coherent() of the branch; the deviations take the rate of the period
# and cohort precisions (5e-4). The branch's rate 0.05 dominates the data:
# for England 2001-2019 (women and men, ages 30-89) the data add about 6e-4
# to the rate of the deviation precision, so its posterior and the width of
# the forecast gaps are set by the prior (dev/coherent-apc-check.md)
coh_hyper <- function() {
  list(age = c(1, 0.5), period = c(1, 5e-4), cohort = c(1, 5e-4), dev = c(1, 5e-4), dev_cohort = c(1, 5e-4),
       rho = c(1, 1))
}

# The model of one group of strata: design, priors and constraints.
# D, E: arrays [age group, period, stratum]; ord: random-walk orders of
# age, period and cohort; ref: index of the reference stratum or NULL
coh_setup <- function(D, E, ppa, ord, cohort_dev = FALSE, ref = NULL, hyper = coh_hyper()) {
  I <- dim(D)[1]
  J <- dim(D)[2]
  S <- dim(D)[3]
  M <- as.integer(ppa)
  K <- M * (I - 1L) + J
  len <- c(mu = 1L, a = S, th = S * I, ph = J, ps = K, d = S * J, c = if (cohort_dev) S * K else 0L)
  end <- cumsum(len)
  pos <- lapply(names(len), function(nm) end[[nm]] - len[[nm]] + seq_len(len[[nm]]))
  names(pos) <- names(len)
  i_mu <- pos$mu
  i_a <- pos$a
  i_th <- pos$th
  i_ph <- pos$ph
  i_ps <- pos$ps
  i_d <- pos$d
  i_c <- pos$c
  P <- end[["c"]]
  ok <- which(is.finite(E) & E > 0 & is.finite(D))
  if (any(D[ok] > E[ok] + 1e-8)) {
    mt_abort("Deaths exceed exposures in some cells; the binomial model needs rates below 1.", class = "model_data")
  }
  ai <- arrayInd(ok, dim(D))
  ci <- ai[, 1]
  cj <- ai[, 2]
  cs <- ai[, 3]
  ck <- M * (I - ci) + cj
  cols <- cbind(i_mu, i_a[cs], i_th[(cs - 1L) * I + ci], i_ph[cj], i_ps[ck], i_d[(cs - 1L) * J + cj],
                if (cohort_dev) i_c[(cs - 1L) * K + ck])
  n <- length(ok)
  X <- Matrix::sparseMatrix(i = rep(seq_len(n), ncol(cols)), j = as.vector(cols), x = 1, dims = c(n, P))
  row <- function(idx, val = 1) {
    r <- numeric(P)
    r[idx] <- val
    r
  }
  dev_rows <- function(idx, n) {
    lapply(seq_len(n), function(t) row(if (is.null(ref)) idx[(seq_len(S) - 1L) * n + t] else idx[(ref - 1L) * n + t]))
  }
  rows <- c(list(row(if (is.null(ref)) i_a else i_a[ref])),
            lapply(seq_len(S), function(s) row(i_th[(s - 1L) * I + seq_len(I)])),
            list(row(i_ph), row(i_ps)),
            # linear trends of age, period and cohort are not identified when
            # all three are second-order random walks: pin the period slope
            if (all(ord == 2L)) list(row(i_ph, seq_len(J) - (J + 1) / 2)),
            dev_rows(i_d, J),
            if (cohort_dev) dev_rows(i_c, K))
  A <- do.call(rbind, rows)
  st <- list(I = I, J = J, S = S, M = M, K = K, P = P, X = X, tX = Matrix::t(X), y = D[ok], N = E[ok], cells = ok,
             ord = ord, idx = list(mu = i_mu, a = i_a, th = i_th, ph = i_ph, ps = i_ps, d = i_d, c = i_c),
             Ka = coh_rw(I, ord[1]), Kp = coh_rw(J, ord[2]), Kc = coh_rw(K, ord[3]), A = A, tA = t(A),
             cohort_dev = cohort_dev, ref = ref, hyper = hyper, tau_a = 0.01, ridge = 1e-6, dims = dim(D))
  coh_pattern(st, cols)
}

# Fixed sparsity pattern of the precision matrices H = X'WX + Q (upper
# triangle) and the pieces from which their values are computed in every
# iteration without building new sparse matrices: `tmap` maps the cell
# weights w to the entries of X'WX, and `basis` holds the entries of the
# prior components, whose coefficients depend on the hyperparameters (the
# AR(1) precision is I + rho^2 T2 - rho T1)
coh_pattern <- function(st, cols) {
  S <- st$S
  P <- st$P
  bd <- function(n_before, B, n_after) Matrix::bdiag(Matrix::Diagonal(n_before, 0), B, Matrix::Diagonal(n_after, 0))
  ix <- st$idx
  before <- function(r) min(r) - 1L
  after <- function(r) P - max(r)
  kron <- function(B) kronecker(Matrix::Diagonal(S), B)
  ar_parts <- function(n) {
    list(I = Matrix::Diagonal(n),
         T2 = Matrix::Diagonal(n, x = c(0, rep(1, n - 2), 0)),
         T1 = Matrix::bandSparse(n, k = 1, symmetric = TRUE, diagonals = list(rep(1, n - 1))))
  }
  ad <- ar_parts(st$J)
  basis <- list(tau_a = bd(before(ix$a), Matrix::Diagonal(S), after(ix$a)),
                kappa = bd(before(ix$th), kron(st$Ka$K), after(ix$th)),
                lambda = bd(before(ix$ph), st$Kp$K, after(ix$ph)),
                nu = bd(before(ix$ps), st$Kc$K, after(ix$ps)),
                d_I = bd(before(ix$d), kron(ad$I), after(ix$d)),
                d_T2 = bd(before(ix$d), kron(ad$T2), after(ix$d)),
                d_T1 = bd(before(ix$d), kron(ad$T1), after(ix$d)),
                fixed = Matrix::Diagonal(P, st$ridge) +
                  Matrix::crossprod(Matrix::Matrix(st$A, sparse = TRUE)))
  if (st$cohort_dev) {
    ac <- ar_parts(st$K)
    basis <- c(basis, list(c_I = bd(before(ix$c), kron(ac$I), after(ix$c)),
                           c_T2 = bd(before(ix$c), kron(ac$T2), after(ix$c)),
                           c_T1 = bd(before(ix$c), kron(ac$T1), after(ix$c))))
  }
  gmrf_pattern(st, cols, basis)
}

# The sparsity pattern (upper triangle) of X'WX plus the prior components
# `basis` (sparse P x P matrices), the entries of each component on it
# (st$basis), and the map from the cell weights w to the entries of X'WX
# (st$tmap: every pair of the columns of a cell, weighted by its w). Shared
# by the coherent and the multi-cause models
gmrf_pattern <- function(st, cols, basis) {
  P <- st$P
  pat <- Matrix::crossprod(st$X) + Reduce(`+`, lapply(basis, function(B) abs(B))) + Matrix::Diagonal(P)
  pat <- methods::as(Matrix::forceSymmetric(pat, uplo = "U"), "CsparseMatrix")
  pj <- rep(seq_len(P), diff(pat@p))
  key <- (pj - 1) * P + pat@i + 1
  entries <- function(B) {
    B <- methods::as(methods::as(B, "CsparseMatrix"), "TsparseMatrix")
    up <- B@i <= B@j
    out <- numeric(length(key))
    k <- match(B@j[up] * P + B@i[up] + 1, key)
    out[k] <- out[k] + B@x[up]
    out
  }
  st$basis <- lapply(basis, entries)
  nc <- ncol(cols)
  pr <- which(upper.tri(diag(nc), diag = TRUE), arr.ind = TRUE)
  a <- pmin(cols[, pr[, 1], drop = FALSE], cols[, pr[, 2], drop = FALSE])
  b <- pmax(cols[, pr[, 1], drop = FALSE], cols[, pr[, 2], drop = FALSE])
  e <- match((b - 1) * P + a, key)
  st$tmap <- Matrix::sparseMatrix(i = e, j = rep(seq_len(nrow(cols)), ncol(a)), x = 1,
                                  dims = c(length(key), nrow(cols)))
  pat@x[] <- 0
  st$pattern <- pat
  st
}

# Prior precision matrix of all effects for hyperparameters h, plus A'A:
# the effects live on the subspace A x = 0, where x' A'A x = 0, so the
# density there is unchanged, and A'A makes the matrix positive definite in
# the directions that the constraints remove (the intrinsic random walks
# and the likelihood leave them flat)
coh_prior <- function(st, h) {
  b <- st$basis
  v <- b$fixed + st$tau_a * b$tau_a + h$kappa * b$kappa + h$lambda * b$lambda + h$nu * b$nu +
    h$lambda_d * (b$d_I + h$rho^2 * b$d_T2 - h$rho * b$d_T1)
  if (st$cohort_dev) v <- v + h$lambda_c * (b$c_I + h$rho_c^2 * b$c_T2 - h$rho_c * b$c_T1)
  Q <- st$pattern
  Q@x <- v
  Q
}

# Gaussian approximation of the conditional posterior of the effects at x
# (one Newton step under the constraints A x = 0): the constrained mean,
# the precision H and what the proposal density needs
coh_approx <- function(st, x, Qp, sym) {
  eta <- as.vector(st$X %*% x)
  p <- stats::plogis(eta)
  w <- st$N * p * (1 - p)
  Qx <- as.vector(Qp %*% x)
  g <- as.vector(st$tX %*% (st$y - st$N * p)) - Qx
  # symmetric storage (update() factorises H H' for a general matrix)
  H <- Qp
  H@x <- Qp@x + as.vector(st$tmap %*% w)
  L <- Matrix::update(sym, H)
  V <- as.matrix(Matrix::solve(L, st$tA, system = "A"))
  R <- chol(st$A %*% V)
  correct <- function(v) v - as.vector(V %*% backsolve(R, forwardsolve(t(R), as.vector(st$A %*% v))))
  m <- correct(x + as.vector(Matrix::solve(L, g, system = "A")))
  # diag() of a Cholesky factor gives the squared diagonal of L: log det(H)
  logdet <- sum(log(Matrix::diag(L))) + 2 * sum(log(diag(R)))
  lp <- sum(st$y * eta - st$N * coh_softplus(eta)) - 0.5 * sum(x * Qx)
  list(m = m, H = H, L = L, correct = correct, logdet = logdet, lp = lp)
}

# Log conditional posterior of the effects (up to a constant)
coh_lp <- function(st, x, Qp) {
  eta <- as.vector(st$X %*% x)
  sum(st$y * eta - st$N * coh_softplus(eta)) - 0.5 * sum(x * as.vector(Qp %*% x))
}

# Starting values: the overall logit, the empirical logits by age of each
# stratum (centred), then Newton steps with step halving to the mode
coh_start <- function(st, Qp, sym) {
  x <- numeric(st$P)
  mu <- stats::qlogis(sum(st$y) / sum(st$N))
  x[st$idx$mu] <- mu
  Dm <- array(0, st$dims)
  Em <- array(0, st$dims)
  Dm[st$cells] <- st$y
  Em[st$cells] <- st$N
  for (s in seq_len(st$S)) {
    r <- (rowSums(Dm[, , s, drop = FALSE]) + 0.5) / (rowSums(Em[, , s, drop = FALSE]) + 1)
    th <- stats::qlogis(pmin(pmax(r, 1e-10), 1 - 1e-10)) - mu
    x[st$idx$th[(s - 1L) * st$I + seq_len(st$I)]] <- th - mean(th)
  }
  lp <- coh_lp(st, x, Qp)
  for (it in 1:50) {
    m <- coh_approx(st, x, Qp, sym)$m
    step <- 1
    repeat {
      xn <- x + step * (m - x)
      ln <- coh_lp(st, xn, Qp)
      if (is.finite(ln) && ln >= lp - 1e-10) break
      step <- step / 2
      if (step < 1e-6) break
    }
    if (step < 1e-6) break
    done <- abs(ln - lp) < 1e-8 * (1 + abs(lp))
    x <- xn
    lp <- ln
    if (done) break
  }
  x
}

# Log density (up to a constant) of the constrained Gaussian approximation
# `ap` at z
coh_logq <- function(ap, z) {
  v <- z - ap$m
  0.5 * ap$logdet - 0.5 * sum(v * as.vector(ap$H %*% v))
}

# One chain. `init`: list of effects x and hyperparameters h; returns the
# kept draws. `fix_hyper = TRUE` keeps the hyperparameters at their starting
# values (tests of the sampler of the effects)
coh_chain <- function(st, n_iter, burn_in, thin, init = NULL, fix_hyper = FALSE) {
  hp <- st$hyper
  h <- init$h %||% list(kappa = 1, lambda = 1, nu = 1, lambda_d = 10, rho = 0.5, lambda_c = 10, rho_c = 0.5)
  Qp <- coh_prior(st, h)
  H0 <- Qp
  H0@x <- Qp@x + as.vector(st$tmap %*% (st$N * 0.01))
  sym <- Matrix::Cholesky(H0, perm = TRUE, LDL = FALSE, super = FALSE)
  # start at the mode for the starting hyperparameters
  x <- init$x %||% coh_start(st, Qp, sym)
  keep <- seq.int(burn_in + thin, n_iter, by = thin)
  nk <- length(keep)
  S <- st$S
  out <- list(x = matrix(NA_real_, nk, st$P), h = matrix(NA_real_, nk, 7,
                                                        dimnames = list(NULL, c("kappa", "lambda", "nu", "lambda_d",
                                                                                "rho", "lambda_c", "rho_c"))))
  acc <- 0L
  acc_rho <- c(0L, 0L)
  sd_rho <- c(0.5, 0.5)
  thS <- function(v) matrix(v[st$idx$th], S, st$I, byrow = TRUE)
  devm <- function(v, idx, len) matrix(v[idx], S, len, byrow = TRUE)
  n_dev <- (S - 1L)
  store <- 0L
  for (it in seq_len(n_iter)) {
    # effects: Metropolis-Hastings with the Gaussian approximation
    cur <- coh_approx(st, x, Qp, sym)
    z <- as.vector(Matrix::solve(cur$L, Matrix::solve(cur$L, stats::rnorm(st$P), system = "Lt"), system = "Pt"))
    xs <- cur$m + cur$correct(z)
    pro <- coh_approx(st, xs, Qp, sym)
    la <- pro$lp - cur$lp + coh_logq(pro, x) - coh_logq(cur, xs)
    if (is.finite(la) && log(stats::runif(1)) < la) {
      x <- xs
      acc <- acc + 1L
    }
    if (fix_hyper) {
      if (it %in% keep) {
        store <- store + 1L
        out$x[store, ] <- x
        out$h[store, ] <- unlist(h[colnames(out$h)])
      }
      next
    }
    # precisions: Gibbs
    Th <- thS(x)
    qa <- sum(vapply(seq_len(S), function(s) sum(Th[s, ] * as.vector(st$Ka$K %*% Th[s, ])), 0))
    h$kappa <- stats::rgamma(1, hp$age[1] + S * st$Ka$rank / 2, hp$age[2] + qa / 2)
    ph <- x[st$idx$ph]
    h$lambda <- stats::rgamma(1, hp$period[1] + st$Kp$rank / 2, hp$period[2] + sum(ph * as.vector(st$Kp$K %*% ph)) / 2)
    ps <- x[st$idx$ps]
    h$nu <- stats::rgamma(1, hp$cohort[1] + st$Kc$rank / 2, hp$cohort[2] + sum(ps * as.vector(st$Kc$K %*% ps)) / 2)
    # deviations: precision (Gibbs) and rho (Metropolis on the logit scale)
    # (the step size of rho adapts towards an acceptance rate of 0.44
    # during the burn-in only)
    upd <- function(Y, rho, hyp, sd) {
      lam <- stats::rgamma(1, hyp[1] + n_dev * ncol(Y) / 2, hyp[2] + coh_ar1_quad(Y, rho) / 2)
      lt <- function(r) {
        n_dev / 2 * log(1 - r^2) - lam / 2 * coh_ar1_quad(Y, r) + stats::dbeta(r, hp$rho[1], hp$rho[2], log = TRUE) +
          log(r * (1 - r))
      }
      rs <- stats::plogis(stats::qlogis(rho) + stats::rnorm(1, 0, sd))
      ok <- log(stats::runif(1)) < lt(rs) - lt(rho)
      list(lam = lam, rho = if (ok) rs else rho, ok = ok,
           sd = if (it <= burn_in) exp(log(sd) + min(0.5, 5 / it) * (ok - 0.44)) else sd)
    }
    r <- upd(devm(x, st$idx$d, st$J), h$rho, hp$dev, sd_rho[1])
    h$lambda_d <- r$lam
    h$rho <- r$rho
    sd_rho[1] <- r$sd
    if (it > burn_in) acc_rho[1] <- acc_rho[1] + r$ok
    if (st$cohort_dev) {
      r <- upd(devm(x, st$idx$c, st$K), h$rho_c, hp$dev_cohort, sd_rho[2])
      h$lambda_c <- r$lam
      h$rho_c <- r$rho
      sd_rho[2] <- r$sd
      if (it > burn_in) acc_rho[2] <- acc_rho[2] + r$ok
    }
    Qp <- coh_prior(st, h)
    if (it %in% keep) {
      store <- store + 1L
      out$x[store, ] <- x
      out$h[store, ] <- unlist(h[colnames(out$h)])
    }
  }
  if (!st$cohort_dev) out$h[, c("lambda_c", "rho_c")] <- NA_real_
  out$accept <- acc / n_iter
  out$accept_rho <- acc_rho / max(1, n_iter - burn_in)
  out$last <- list(x = x, h = h)
  out
}

# Gelman-Rubin potential scale reduction factor (point estimate; Gelman and
# Rubin 1992) of each column, from a list of chains [draw, quantity]
coh_psrf <- function(chains) {
  n <- nrow(chains[[1]])
  if (length(chains) < 2 || n < 4) return(rep(NA_real_, ncol(chains[[1]])))
  means <- vapply(chains, colMeans, numeric(ncol(chains[[1]])))
  vars <- vapply(chains, function(m) apply(m, 2, stats::var), numeric(ncol(chains[[1]])))
  means <- matrix(means, ncol = length(chains))
  vars <- matrix(vars, ncol = length(chains))
  W <- rowMeans(vars)
  B <- n * apply(means, 1, stats::var)
  vhat <- (n - 1) / n * W + B / n
  out <- sqrt(vhat / W)
  out[W <= 0] <- NA_real_
  out
}

# Quantities monitored for convergence: log precisions, rho, the period
# effect and the deviations of the last year, and the mean linear predictor
# of each stratum in the last year
coh_monitor <- function(st, ch) {
  x <- ch$x
  h <- ch$h
  J <- st$J
  S <- st$S
  lp_last <- vapply(seq_len(S), function(s) {
    rowMeans(x[, st$idx$th[(s - 1L) * st$I + seq_len(st$I)], drop = FALSE]) + x[, st$idx$a[s]] +
      x[, st$idx$d[(s - 1L) * J + J]]
  }, numeric(nrow(x)))
  lp_last <- matrix(lp_last, nrow(x)) + x[, st$idx$mu] + x[, st$idx$ph[J]]
  hh <- cbind(log(h[, c("kappa", "lambda", "nu", "lambda_d")]), rho = h[, "rho"])
  if (st$cohort_dev) hh <- cbind(hh, log(h[, "lambda_c", drop = FALSE]), rho_c = h[, "rho_c"])
  colnames(lp_last) <- paste0("level_", seq_len(S))
  cbind(hh, period_last = x[, st$idx$ph[J]], lp_last)
}

# Fit of one group: several chains (parallel streams), refits with longer
# chains while not converged
coh_fit <- function(D, E, ppa, ord, cohort_dev, ref, mcmc, chains = 4L, cores = 1L, max_refits = 1L,
                    hyper = coh_hyper()) {
  st <- coh_setup(D, E, ppa, ord, cohort_dev, ref, hyper)
  o <- mcmc$options
  n_iter <- if (is.numeric(o$number_of_iterations)) o$number_of_iterations else 10000
  burn <- if (is.numeric(o$burn_in)) o$burn_in else 2000
  thin <- if (is.numeric(o$step)) o$step else 10
  msg <- character()
  for (attempt in 0:max_refits) {
    res <- par_lapply(seq_len(chains), function(i) coh_chain(st, n_iter, burn, thin), cores = cores,
                      streams = TRUE)
    bad <- vapply(res, function(r) !is.list(r) || is.null(r$x), logical(1))
    if (any(bad)) mt_abort("A chain of the coherent APC sampler failed.", class = "model_fit")
    rhat <- coh_psrf(lapply(res, function(r) coh_monitor(st, r)))
    names(rhat) <- colnames(coh_monitor(st, res[[1]]))
    conv <- all(is.finite(rhat)) && max(rhat) < 1.1
    if (conv) break
    n_iter <- 2 * n_iter
    burn <- 2 * burn
    msg <- c(msg, sprintf("not converged (attempt %d, largest R-hat %.2f)", attempt + 1, max(rhat, na.rm = TRUE)))
  }
  if (!conv) {
    mt_abort(sprintf(paste("The coherent APC chains did not converge after %d refit(s) (largest R-hat %.2f); use",
                           "longer chains (mcmc_control())."), max_refits, max(rhat, na.rm = TRUE)),
             class = "convergence")
  }
  x <- do.call(rbind, lapply(res, `[[`, "x"))
  h <- do.call(rbind, lapply(res, `[[`, "h"))
  # fitted rates: posterior means over up to 500 draws spread over the chains
  use <- unique(round(seq(1, nrow(x), length.out = min(nrow(x), 500))))
  eta <- as.matrix(Matrix::tcrossprod(x[use, , drop = FALSE], coh_cell_design(st)))
  fitted <- array(colMeans(stats::plogis(eta)), c(st$I, st$J, st$S))
  st$X <- NULL
  list(setup = st, x = x, h = h, psrf = rhat, accept = vapply(res, `[[`, 0, "accept"),
       accept_rho = do.call(rbind, lapply(res, `[[`, "accept_rho")), fitted = fitted,
       message = paste(msg, collapse = "; "))
}

# Design of every cell [age group, period, stratum] (cells without exposure
# included), for fitted rates
coh_cell_design <- function(st) {
  I <- st$I
  J <- st$J
  S <- st$S
  K <- st$K
  g <- expand.grid(i = seq_len(I), j = seq_len(J), s = seq_len(S))
  k <- st$M * (I - g$i) + g$j
  ix <- st$idx
  cols <- cbind(ix$mu, ix$a[g$s], ix$th[(g$s - 1L) * I + g$i], ix$ph[g$j], ix$ps[k], ix$d[(g$s - 1L) * J + g$j],
                if (st$cohort_dev) ix$c[(g$s - 1L) * K + k])
  Matrix::sparseMatrix(i = rep(seq_len(nrow(g)), ncol(cols)), j = as.vector(cols), x = 1,
                       dims = c(nrow(g), st$P))
}

# Random-walk extension of an effect by `h` steps with the scaled
# innovation precision
coh_extend_rw <- function(v, h, prec, order) {
  n <- length(v)
  out <- c(v, numeric(h))
  e <- stats::rnorm(h, 0, 1 / sqrt(prec))
  for (t in seq_len(h)) {
    i <- n + t
    out[i] <- if (order == 1L) out[i - 1] + e[t] else 2 * out[i - 1] - out[i - 2] + e[t]
  }
  out
}

# AR(1) extension of the deviations Y [stratum, time] by `h` steps; the
# deviations of every step are centred over the strata (sum to zero), or the
# reference stays at zero
coh_extend_ar1 <- function(Y, h, rho, prec, ref) {
  S <- nrow(Y)
  n <- ncol(Y)
  out <- cbind(Y, matrix(0, S, h))
  for (t in seq_len(h)) {
    v <- rho * out[, n + t - 1] + stats::rnorm(S, 0, 1 / sqrt(prec))
    if (is.null(ref)) v <- v - mean(v) else v[ref] <- 0
    out[, n + t] <- v
  }
  out
}

# Forecast draws of the rates of a group: array [age group, h, stratum,
# n_draws], from posterior draws chosen at random
coh_simulate <- function(fit, h, n_draws) {
  st <- fit$setup
  I <- st$I
  J <- st$J
  S <- st$S
  M <- st$M
  K <- st$K
  ix <- st$idx
  nd <- nrow(fit$x)
  pick <- if (n_draws <= nd) sample.int(nd, n_draws) else sample.int(nd, n_draws, replace = TRUE)
  out <- array(NA_real_, c(I, h, S, n_draws))
  fut <- J + seq_len(h)
  for (r in seq_len(n_draws)) {
    x <- fit$x[pick[r], ]
    hh <- fit$h[pick[r], ]
    ph <- coh_extend_rw(x[ix$ph], h, hh[["lambda"]] * st$Kp$scale, st$ord[2])
    ps <- coh_extend_rw(x[ix$ps], h, hh[["nu"]] * st$Kc$scale, st$ord[3])
    dv <- coh_extend_ar1(matrix(x[ix$d], S, J, byrow = TRUE), h, hh[["rho"]], hh[["lambda_d"]], st$ref)
    cv <- if (st$cohort_dev) {
      coh_extend_ar1(matrix(x[ix$c], S, K, byrow = TRUE), h, hh[["rho_c"]], hh[["lambda_c"]], st$ref)
    }
    Th <- matrix(x[ix$th], S, I, byrow = TRUE)
    for (s in seq_len(S)) {
      for (i in seq_len(I)) {
        k <- M * (I - i) + fut
        eta <- x[ix$mu] + x[ix$a[s]] + Th[s, i] + ph[fut] + ps[k] + dv[s, fut] + (if (st$cohort_dev) cv[s, k] else 0)
        out[i, , s, r] <- stats::plogis(eta)
      }
    }
  }
  out
}

# Index of the reference stratum of a group (strata: rows of the group), or
# NULL without reference
coh_reference <- function(strata, reference) {
  if (is.null(reference)) return(NULL)
  miss <- setdiff(names(reference), names(strata))
  if (length(miss)) {
    mt_abort(sprintf("`reference` names unknown stratum variables: %s.", paste(miss, collapse = ", ")),
             class = "argument")
  }
  hit <- Reduce(`&`, lapply(names(reference), function(v) strata[[v]] %in% reference[[v]]))
  if (sum(hit) != 1L) {
    mt_abort(sprintf("`reference` must select exactly one stratum in every group (%d selected in group %s).",
                     sum(hit), paste(strata$label, collapse = ", ")), class = "argument")
  }
  which(hit)
}

#' @export
mt_fit.mt_spec_apc_coherent <- function(spec, data, prep = NULL, cores = 1, ...) {
  check_class(data, "mt_data", "data")
  pr <- mt_prep_apply_raw(prep, data)
  d <- pr$data
  md <- model_inputs(d, spec)
  D_obs <- raw_deaths(pr$raw, md) %||% md$D
  p <- spec$params
  src <- if (p$resolution == "observed") data else d
  cells <- apc_cells(src, md$periods)
  shape <- apc_shape(md, cells$age, p$shape_years)
  st <- d$strata
  groups <- if (is.null(p$common_by)) list(all = seq_len(nrow(st))) else strata_groups(st, p$common_by)
  small <- names(groups)[lengths(groups) < 2]
  if (length(small)) {
    mt_abort(sprintf("The coherent model needs at least two strata per group; group(s) %s have one.",
                     paste(small, collapse = ", ")), class = "argument")
  }
  ord <- match(c(p$age, p$period, p$cohort), c("rw1", "rw2"))
  chain_cores <- if (isTRUE(p$mcmc$parallel)) max(1L, min(4L, cores)) else 1L
  seed <- p$seed %||% sample.int(.Machine$integer.max, 1L)
  fitted <- array(NA_real_, dim(md$D))
  groups_fit <- list()
  diag <- list()
  for (g in names(groups)) {
    idx <- groups[[g]]
    ref <- coh_reference(st[idx], p$reference)
    t0 <- proc.time()[["elapsed"]]
    res <- with_seed(seed + match(g, names(groups)), coh_fit(
      cells$D[, , idx, drop = FALSE], cells$E[, , idx, drop = FALSE], cells$ppa, ord, p$cohort_deviation, ref,
      p$mcmc, cores = chain_cores, max_refits = p$max_refits))
    secs <- proc.time()[["elapsed"]] - t0
    res$strata <- idx
    groups_fit[[g]] <- res
    for (q in seq_along(idx)) {
      s <- idx[q]
      fitted[, , s] <- res$fitted[shape$group, , q, drop = FALSE] * shape$r[, s]
    }
    diag[[g]] <- data.table(group = g, strata = paste(st$label[idx], collapse = ", "), status = "ok",
                            message = res$message, seconds = secs, max_rhat = max(res$psrf),
                            accept = mean(res$accept))
  }
  structure(list(spec = spec, prep = prep, fits = list(groups = groups_fit, fitted = fitted, shape = shape),
                 diagnostics = rbindlist(diag), age = md$age, periods = md$periods, strata = st, omega = d$omega,
                 observed = D_obs / md$E, train_D = md$D, raw_D = raw_deaths(pr$raw, md), train_E = md$E,
                 exposure_last = md$E[, length(md$periods), , drop = FALSE], level_totals = data$level_totals,
                 apc_age = cells$age),
            class = "mt_fit")
}

#' @export
simulate_model.mt_spec_apc_coherent <- function(spec, fit, h, n_draws, ...) {
  shape <- fit$shape
  ns <- ncol(shape$r)
  out <- array(NA_real_, c(length(shape$group), h, ns, n_draws))
  for (g in fit$groups) {
    sim <- coh_simulate(g, h, n_draws)
    for (q in seq_along(g$strata)) {
      s <- g$strata[q]
      out[, , s, ] <- sim[shape$group, , q, , drop = FALSE] * shape$r[, s]
    }
  }
  out
}

# mt_parameters() of the coherent APC model: posterior medians and 90%
# intervals of the effects (logit scale), the hyperparameters, and the
# Gelman-Rubin statistics of the monitored quantities, per group
coh_parameters <- function(fit) {
  lab <- fit$strata$label
  rbindlist(lapply(names(fit$fits$groups), function(g) {
    f <- fit$fits$groups[[g]]
    st <- f$setup
    q <- function(M) t(apply(as.matrix(M), 2, stats::quantile, c(0.05, 0.5, 0.95), na.rm = TRUE))
    one <- function(stratum, parameter, M) {
      Q <- q(M)
      data.table(stratum = stratum, parameter = parameter, index = as.character(seq_len(nrow(Q))),
                 value = Q[, 2], lower = Q[, 1], upper = Q[, 3])
    }
    common <- paste0("common (", g, ")")
    sl <- lab[f$strata]
    own <- rbindlist(lapply(seq_along(sl), function(s) {
      rbind(one(sl[s], "age", f$x[, st$idx$th[(s - 1L) * st$I + seq_len(st$I)], drop = FALSE]),
            one(sl[s], "level", f$x[, st$idx$a[s], drop = FALSE]),
            one(sl[s], "deviation", f$x[, st$idx$d[(s - 1L) * st$J + seq_len(st$J)], drop = FALSE]),
            if (st$cohort_dev) one(sl[s], "cohort_deviation", f$x[, st$idx$c[(s - 1L) * st$K + seq_len(st$K)],
                                                                  drop = FALSE]))
    }))
    hy <- f$h[, colSums(is.na(f$h)) == 0, drop = FALSE]
    hq <- q(hy)
    rbind(one(common, "intercept", f$x[, st$idx$mu, drop = FALSE]),
          one(common, "period", f$x[, st$idx$ph, drop = FALSE]),
          one(common, "cohort", f$x[, st$idx$ps, drop = FALSE]),
          own,
          data.table(stratum = common, parameter = colnames(hy), index = "", value = hq[, 2], lower = hq[, 1],
                     upper = hq[, 3]),
          data.table(stratum = common, parameter = "rhat", index = names(f$psrf), value = unname(f$psrf),
                     lower = NA_real_, upper = NA_real_))
  }))
}
