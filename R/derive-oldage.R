# Hazard laws for old ages: the models compared by Feehan (2018, Table 1 and
# Appendix E), fitted by maximum likelihood with binomial deaths of those
# reaching each age (Thatcher, Kannisto and Vaupel 1998, chapter 6; Feehan
# 2018, equation 4). Gompertz and Kannisto keep their own Fisher scoring in
# derive-lifetable.R (tail_fit_nd()); the other laws use the generic fit
# below. Used by life_table(open = ) and mt_open_check().

# Number of parameters of each law
oldage_npar <- c(gompertz = 2L, kannisto = 2L, makeham = 3L, weibull = 2L, beard = 3L, logquad = 3L,
                 logistic = 4L, perks = 4L, lynch_brown = 4L)
oldage_laws <- names(oldage_npar)

# Gauss-Legendre nodes and weights on [0, 1] (Golub-Welsch), computed once
gauss_legendre01 <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      n <- 20L
      i <- seq_len(n - 1L)
      b <- i / sqrt(4 * i^2 - 1)
      J <- matrix(0, n, n)
      J[cbind(i, i + 1L)] <- b
      J[cbind(i + 1L, i)] <- b
      e <- eigen(J, symmetric = TRUE)
      cache <<- list(t = (e$values + 1) / 2, w = e$vectors[1, ]^2)
    }
    cache
  }
})

softplus <- function(z) pmax(z, 0) + log1p(exp(-abs(z)))

# Integrated hazard H over [u, u + 1] and its derivatives with respect to the
# internal parameters, for centred ages u = x - centre (vector) and
# parameters th (matrix [parameter, schedule]): matrices [age, schedule].
# Internal parameters (natural ones of Feehan 2018 in oldage_natural()):
#   makeham     (log alpha_c, beta, log gamma): gamma + alpha_c exp(beta u)
#   weibull     (eta, log beta): exp(eta) (x / centre)^(beta - 1)
#   beard       (log alpha_c, log beta, log delta_c): alpha_c exp(beta u) / (1 + delta_c exp(beta u))
#   logquad     (alpha_c, beta_c, gamma): exp(alpha_c + beta_c u + gamma u^2)
#   logistic    (log alpha_c, log beta, log gamma, log delta_c): gamma + beard
#   perks       (log alpha_c, log beta, log gamma, log delta_c):
#               (gamma + alpha_c exp(beta u)) / (1 + delta_c exp(beta u))
#   lynch_brown (alpha, log beta, log gamma, delta_c): alpha + beta atan(gamma (u - delta_c))
# with alpha_c = alpha exp(beta centre), delta_c = delta exp(beta centre)
oldage_H <- function(law, th, u, centre) {
  nu <- length(u)
  n <- ncol(th)
  P <- function(k) matrix(th[k, ], nu, n, byrow = TRUE)
  U <- matrix(u, nu, n)
  switch(law,
         makeham = {
           a <- P(1)
           b <- P(2)
           small <- abs(b) < 1e-6
           bs <- ifelse(small, 1, b)
           e1 <- ifelse(small, 1 + b / 2, expm1(bs) / bs)
           g1 <- ifelse(small, 0.5 + b / 12, exp(bs) / expm1(bs) - 1 / bs)
           G <- exp(a + b * U) * e1
           cg <- exp(P(3))
           list(H = cg + G, dH = list(G, G * (U + g1), cg))
         },
         weibull = {
           eta <- P(1)
           b <- exp(P(2))
           r0 <- (U + centre) / centre
           r1 <- (U + centre + 1) / centre
           s0 <- r0^b
           s1 <- r1^b
           H <- exp(eta) * centre * (s1 - s0) / b
           list(H = H, dH = list(H, exp(eta) * centre * (s1 * log(r1) - s0 * log(r0) - (s1 - s0) / b)))
         },
         beard = {
           bd <- oldage_beard(P(1), P(2), P(3), U)
           list(H = bd$H, dH = list(bd$da, bd$dlb, bd$dld))
         },
         logistic = {
           bd <- oldage_beard(P(1), P(2), P(4), U)
           cg <- exp(P(3))
           list(H = cg + bd$H, dH = list(bd$da, bd$dlb, cg, bd$dld))
         },
         perks = {
           ld <- P(4)
           b <- exp(P(2))
           bd <- oldage_beard(P(1), P(2), ld, U)
           cg <- exp(P(3))
           # gamma / (1 + delta exp(beta t)) integrates to gamma (1 - S / beta) with
           # S = softplus(z1) - softplus(z0), z = log delta + beta t
           z0 <- ld + b * U
           z1 <- z0 + b
           small <- b < 1e-6
           bs <- ifelse(small, 1, b)
           p0 <- stats::plogis(z0)
           p1 <- stats::plogis(z1)
           pm <- stats::plogis(ld + b * (U + 0.5))
           Sb <- ifelse(small, pm, (softplus(z1) - softplus(z0)) / bs)
           Tg <- cg * (1 - Sb)
           dld <- ifelse(small, -cg * pm * (1 - pm), -cg * (p1 - p0) / bs)
           dlb <- ifelse(small, 0, cg * (Sb - ((U + 1) * p1 - U * p0)))
           list(H = Tg + bd$H, dH = list(bd$da, bd$dlb + dlb, Tg, bd$dld + dld))
         },
         logquad = {
           gl <- gauss_legendre01()
           a <- P(1)
           b <- P(2)
           g <- P(3)
           H <- d1 <- d2 <- 0
           for (j in seq_along(gl$t)) {
             V <- U + gl$t[j]
             mu <- gl$w[j] * exp(a + b * V + g * V^2)
             H <- H + mu
             d1 <- d1 + V * mu
             d2 <- d2 + V^2 * mu
           }
           list(H = H, dH = list(H, d1, d2))
         },
         lynch_brown = {
           al <- P(1)
           b <- exp(P(2))
           g <- exp(P(3))
           k0 <- g * (U - P(4))
           k1 <- k0 + g
           Fk <- function(k) k * atan(k) - 0.5 * log1p(k^2)
           A <- (Fk(k1) - Fk(k0)) / g
           list(H = al + b * A,
                dH = list(matrix(1, nu, n), b * A, b / (2 * g) * (log1p(k1^2) - log1p(k0^2)),
                          -b * (atan(k1) - atan(k0))))
         },
         mt_abort(sprintf("Unknown hazard law \"%s\".", law), class = "lifetable"))
}

# The Beard hazard alpha_c exp(beta u) / (1 + delta_c exp(beta u)) with a =
# log alpha_c, lb = log beta, ld = log delta_c: integrated hazard over [u, u +
# 1] (Feehan 2018, eq. 46) and its derivatives. With z = ld + beta t,
# spr(t) = softplus(z) / delta_c (computed without overflow for small delta)
# and R(t) = exp(beta t) / (1 + exp(z)), H = exp(a) (spr(u + 1) - spr(u)) /
# beta, dH/d ld = exp(a) (R(u + 1) - R(u)) / beta - H and dH/d lb = -H +
# exp(a) ((u + 1) R(u + 1) - u R(u))
oldage_beard <- function(a, lb, ld, U) {
  b <- exp(lb)
  spr <- function(t) {
    z <- ld + b * t
    ifelse(z < -30, exp(b * t) * (1 - exp(z) / 2 + exp(2 * z) / 3), softplus(z) * exp(-ld))
  }
  R <- function(t) exp(b * t - softplus(ld + b * t))
  small <- b < 1e-6
  bs <- ifelse(small, 1, b)
  R0 <- R(U)
  R1 <- R(U + 1)
  ea <- exp(a)
  H <- ea * (spr(U + 1) - spr(U)) / bs
  dld <- ea * (R1 - R0) / bs - H
  dlb <- -H + ea * ((U + 1) * R1 - U * R0)
  if (any(small)) {
    Rm <- R(U + 0.5)
    pm <- stats::plogis(ld + b * (U + 0.5))
    H[small] <- (ea * Rm)[small]
    dld[small] <- (-ea * Rm * pm)[small]
    dlb[small] <- 0
  }
  list(H = H, da = H, dlb = dlb, dld = dld)
}

# Binomial log-likelihood of d deaths among N persons reaching each age, with
# q = 1 - exp(-H), without the binomial coefficients; -Inf for a schedule
# whose hazard is not positive at every age
oldage_loglik <- function(H, N, d) {
  q <- -expm1(-H)
  ll <- colSums(d * log(pmax(q, 1e-300)) - (N - d) * H)
  ll[colSums(!is.finite(H) | H <= 0) > 0] <- -Inf
  ll
}

# Solve the symmetric positive definite systems A[, , j] x = b[j, ] for all j
# (Cholesky, vectorised over j; A is [p, p, n], b is [n, p]). Rows of the
# result are NA where A[, , j] is not positive definite
batch_solve <- function(A, b) {
  p <- dim(A)[1]
  n <- dim(A)[3]
  b <- matrix(b, n, p)
  L <- array(0, c(p, p, n))
  ok <- rep(TRUE, n)
  for (j in seq_len(p)) {
    prev <- seq_len(j - 1L)
    s <- A[j, j, ] - (if (j > 1) colSums(L[j, prev, , drop = FALSE]^2, dims = 2) else 0)
    ok <- ok & is.finite(s) & s > 0
    L[j, j, ] <- sqrt(pmax(s, 1e-300))
    for (i in seq_len(p)[-seq_len(j)]) {
      t <- A[i, j, ] - (if (j > 1) colSums(L[i, prev, , drop = FALSE] * L[j, prev, , drop = FALSE], dims = 2) else 0)
      L[i, j, ] <- t / L[j, j, ]
    }
  }
  y <- matrix(0, n, p)
  for (j in seq_len(p)) {
    s <- b[, j]
    for (k in seq_len(j - 1L)) s <- s - L[j, k, ] * y[, k]
    y[, j] <- s / L[j, j, ]
  }
  x <- matrix(0, n, p)
  for (j in rev(seq_len(p))) {
    s <- y[, j]
    for (k in seq_len(p)[-seq_len(j)]) s <- s - L[k, j, ] * x[, k]
    x[, j] <- s / L[j, j, ]
  }
  x[!ok | rowSums(!is.finite(x)) > 0, ] <- NA_real_
  x
}

# Weighted least squares of each column of Y on the columns of X (common
# design) with weights W [row, column]; returns coefficients [p, column]
oldage_wls <- function(X, Y, W) {
  p <- ncol(X)
  n <- ncol(Y)
  A <- array(0, c(p, p, n))
  B <- matrix(0, n, p)
  for (i in seq_len(p)) {
    B[, i] <- colSums(W * X[, i] * Y)
    for (j in seq_len(p)) A[i, j, ] <- colSums(W * X[, i] * X[, j])
  }
  # a tiny ridge keeps nearly collinear designs solvable
  for (i in seq_len(p)) A[i, i, ] <- A[i, i, ] * (1 + 1e-10) + 1e-300
  t(batch_solve(A, B))
}

# Starting values of the internal parameters (a list of matrices [parameter,
# schedule]) from the crude hazards mt = -log(1 - d / N) at the middle of each
# year, by weighted least squares (weights: deaths). Several starts per law;
# the fit keeps the one with the highest likelihood (package-specific; Feehan
# 2018 uses one heuristic start and ten random starts per model)
oldage_starts <- function(law, u, N, d, centre) {
  n <- ncol(N)
  qt <- pmin(pmax(d / N, 1e-8), 1 - 1e-8)
  mt <- -log1p(-qt)
  W <- pmax(d, 1e-8)
  t <- u + 0.5
  X1 <- cbind(1, t)
  mmin <- apply(mt, 2, min)
  mmax <- apply(mt, 2, max)
  rep_row <- function(v) matrix(v, length(u), n, byrow = TRUE)
  lin <- function(Y) oldage_wls(X1, Y, W)
  # beard-type starts: log(m / (1 - m / P)) linear in age for the plateau P
  beard_start <- function(m, P) {
    Pm <- rep_row(P)
    cf <- lin(log(pmax(m, 1e-12) / pmax(1 - m / Pm, 1e-3)))
    rbind(cf[1, ], log(pmax(cf[2, ], 1e-3)), cf[1, ] - log(P))
  }
  plateaus <- list(2 * mmax, 5 * mmax, rep(30, n))
  switch(law,
         makeham = {
           out <- list()
           for (f in c(1e-3, 0.5)) {
             gam <- f * mmin
             cf <- lin(log(pmax(mt - rep_row(gam), 1e-3 * mt)))
             out[[length(out) + 1L]] <- rbind(cf[1, ], cf[2, ], log(gam))
           }
           out
         },
         weibull = {
           cf <- oldage_wls(cbind(1, log((t + centre) / centre)), log(mt), W)
           list(rbind(cf[1, ], log(pmax(cf[2, ] + 1, 0.05))))
         },
         logquad = {
           cf <- oldage_wls(cbind(1, t, t^2), log(mt), W)
           g <- lin(log(mt))
           list(cf, rbind(g, 0))
         },
         beard = lapply(c(list(1.2 * mmax, rep(1, n)), plateaus), function(P) beard_start(mt, P)),
         logistic = {
           out <- list()
           for (f in c(1e-3, 0.5)) {
             gam <- f * mmin
             mg <- pmax(mt - rep_row(gam), 1e-3 * mt)
             for (P in plateaus) {
               s <- beard_start(mg, P)
               out[[length(out) + 1L]] <- rbind(s[1, ], s[2, ], log(gam), s[3, ])
             }
           }
           out
         },
         perks = {
           out <- list()
           for (f in c(1e-3, 0.5)) {
             gam <- f * mmin
             for (P in plateaus) {
               cf <- lin(log(pmax(mt - rep_row(gam), 1e-3 * mt) / pmax(1 - mt / rep_row(P), 1e-3)))
               out[[length(out) + 1L]] <- rbind(cf[1, ], log(pmax(cf[2, ], 1e-3)), log(gam), cf[1, ] - log(P))
             }
           }
           out
         },
         lynch_brown = {
           # alpha and beta by least squares of the crude hazards for a grid of
           # gamma and delta; the three best combinations with beta > 0 and a
           # positive hazard at every age
           Wl <- N / pmax(mt, 1e-12)
           grid <- expand.grid(g = c(0.02, 0.05, 0.1, 0.2, 0.5),
                               dl = c(min(t), mean(t), max(t), max(t) + 10, max(t) + 25))
           cand <- lapply(seq_len(nrow(grid)), function(i) {
             at <- atan(grid$g[i] * (t - grid$dl[i]))
             cf <- oldage_wls(cbind(1, at), mt, Wl)
             fit <- outer(rep(1, length(t)), cf[1, ]) + outer(at, cf[2, ])
             at0 <- atan(grid$g[i] * (u - grid$dl[i]))
             okc <- cf[2, ] > 0 & colSums(outer(at0, cf[2, ]) + rep_row(cf[1, ]) <= 0) == 0
             sse <- colSums(Wl * (mt - fit)^2)
             sse[!okc | !is.finite(sse)] <- Inf
             list(th = rbind(cf[1, ], log(pmax(cf[2, ], 1e-12)), log(grid$g[i]), grid$dl[i]), sse = sse)
           })
           S <- vapply(cand, `[[`, numeric(n), "sse")
           S <- matrix(S, n)
           lapply(1:3, function(r) {
             pick <- apply(S, 1, function(s) order(s)[r])
             th <- matrix(NA_real_, 4, n)
             for (i in unique(pick)) th[, pick == i] <- cand[[i]]$th[, pick == i]
             th
           })
         },
         mt_abort(sprintf("Unknown hazard law \"%s\".", law), class = "lifetable"))
}

# Maximum-likelihood fit of a law (not Gompertz or Kannisto) to the deaths d
# among N persons reaching ages x (matrices [age, schedule]):
# Levenberg-Marquardt steps on the Fisher information of the binomial
# likelihood, from several starts, keeping the best. Returns the internal
# parameters, the centre, the log-likelihood and whether each fit converged
oldage_fit_generic <- function(law, x, N, d, maxit = 2000L) {
  centre <- mean(x) + 0.5
  u <- x - centre
  p <- oldage_npar[[law]]
  n <- ncol(N)
  best_th <- matrix(NA_real_, p, n)
  best_ll <- rep(-Inf, n)
  best_conv <- rep(FALSE, n)
  for (th0 in oldage_starts(law, u, N, d, centre)) {
    r <- oldage_lm(law, th0, u, centre, N, d, maxit)
    better <- is.finite(r$ll) & r$ll > best_ll
    best_th[, better] <- r$th[, better]
    best_ll[better] <- r$ll[better]
    best_conv[better] <- r$converged[better]
  }
  best_ll[!is.finite(best_ll)] <- NA_real_
  list(law = law, theta = best_th, centre = centre, loglik = best_ll, ages = x, converged = best_conv)
}

# Levenberg-Marquardt iterations for all schedules at once: the step solves
# (I + lambda diag(I)) s = score with the Fisher information I; lambda falls
# after a step that raises the log-likelihood and rises otherwise. Converged
# when the predicted gain of a step is negligible, when three accepted steps
# in a row gain almost nothing (a parameter drifting to the boundary of its
# range, such as a Makeham constant of 0), or when no step raises the
# log-likelihood (its maximum to rounding). Not converged: after maxit
# iterations, or when 50 iterations raised the log-likelihood by less than
# 1e-7 per death while the parameters still move (a supremum that is not
# attained, as for Lynch-Brown laws whose inflection lies beyond the data)
oldage_lm <- function(law, th, u, centre, N, d, maxit = 2000L) {
  p <- nrow(th)
  n <- ncol(th)
  bad_start <- colSums(!is.finite(th)) > 0
  th[, bad_start] <- 0
  ll <- oldage_loglik(oldage_H(law, th, u, centre)$H, N, d)
  ll[bad_start] <- -Inf
  tol <- 1e-14 * pmax(colSums(d), 1e-12)
  lam <- rep(1e-3, n)
  done <- !is.finite(ll)
  conv <- rep(FALSE, n)
  flat <- integer(n)
  ll_before <- ll
  for (it in seq_len(maxit)) {
    if (it %% 50L == 0L) {
      stall <- !done & ll - ll_before < 1e-7 * pmax(colSums(d), 1e-12)
      done[stall] <- TRUE
      ll_before <- ll
    }
    act <- which(!done)
    if (!length(act)) break
    hz <- oldage_H(law, th[, act, drop = FALSE], u, centre)
    Na <- N[, act, drop = FALSE]
    da <- d[, act, drop = FALSE]
    q <- pmax(-expm1(-hz$H), 1e-300)
    r <- (da - Na * q) / q
    w <- Na * exp(-hz$H) / q
    sc <- matrix(vapply(hz$dH, function(D) colSums(r * D), numeric(length(act))), length(act), p)
    Im <- array(0, c(p, p, length(act)))
    for (i in seq_len(p)) for (j in seq_len(i)) {
      Im[i, j, ] <- Im[j, i, ] <- colSums(w * hz$dH[[i]] * hz$dH[[j]])
    }
    pending <- rep(TRUE, length(act))
    for (tries in seq_len(40)) {
      k <- which(pending)
      if (!length(k)) break
      A <- Im[, , k, drop = FALSE]
      for (i in seq_len(p)) A[i, i, ] <- A[i, i, ] * (1 + lam[act[k]]) + 1e-300
      st <- batch_solve(A, sc[k, , drop = FALSE])
      gp <- rowSums(st * sc[k, , drop = FALSE]) / 2
      th_new <- th[, act[k], drop = FALSE] + t(st)
      ok_step <- rowSums(!is.finite(st)) == 0
      ll_new <- rep(-Inf, length(k))
      if (any(ok_step)) {
        ll_new[ok_step] <- oldage_loglik(oldage_H(law, th_new[, ok_step, drop = FALSE], u, centre)$H,
                                         N[, act[k][ok_step], drop = FALSE], d[, act[k][ok_step], drop = FALSE])
      }
      acc <- ok_step & ll_new >= ll[act[k]]
      # converged: a step would gain (almost) nothing
      tiny <- ok_step & gp < tol[act[k]] & lam[act[k]] <= 1e-3
      if (any(acc)) {
        j <- act[k[acc]]
        gain <- ll_new[acc] - ll[j]
        th[, j] <- th_new[, acc, drop = FALSE]
        ll[j] <- ll_new[acc]
        lam[j] <- pmax(lam[j] / 10, 1e-12)
        flat[j] <- ifelse(gain < tol[j], flat[j] + 1L, 0L)
      }
      conv[act[k[tiny]]] <- TRUE
      pending[k[acc | tiny]] <- FALSE
      rej <- k[!(acc | tiny)]
      lam[act[rej]] <- lam[act[rej]] * 10
      # no step raises the log-likelihood: its maximum to rounding
      stuck <- rej[lam[act[rej]] > 1e10]
      conv[act[stuck]] <- TRUE
      pending[stuck] <- FALSE
    }
    conv[act[flat[act] >= 3]] <- TRUE
    done[act] <- conv[act] | done[act]
  }
  ll[bad_start] <- -Inf
  list(th = th, ll = ll, converged = conv)
}

# Natural parameters (alpha, beta, gamma, delta of Feehan 2018, Table 1, for
# age x in years; NA where a law has fewer) from the internal ones
oldage_natural <- function(law, th, centre) {
  n <- ncol(th)
  out <- matrix(NA_real_, 4, n, dimnames = list(c("alpha", "beta", "gamma", "delta"), NULL))
  switch(law,
         makeham = {
           out["beta", ] <- th[2, ]
           out["alpha", ] <- exp(th[1, ] - th[2, ] * centre)
           out["gamma", ] <- exp(th[3, ])
         },
         weibull = {
           out["beta", ] <- exp(th[2, ])
           out["alpha", ] <- exp(th[1, ] - (out["beta", ] - 1) * log(centre))
         },
         beard = {
           out["beta", ] <- exp(th[2, ])
           out["alpha", ] <- exp(th[1, ] - out["beta", ] * centre)
           out["delta", ] <- exp(th[3, ] - out["beta", ] * centre)
         },
         logistic = ,
         perks = {
           out["beta", ] <- exp(th[2, ])
           out["alpha", ] <- exp(th[1, ] - out["beta", ] * centre)
           out["gamma", ] <- exp(th[3, ])
           out["delta", ] <- exp(th[4, ] - out["beta", ] * centre)
         },
         logquad = {
           out["gamma", ] <- th[3, ]
           out["beta", ] <- th[2, ] - 2 * th[3, ] * centre
           out["alpha", ] <- th[1, ] - th[2, ] * centre + th[3, ] * centre^2
         },
         lynch_brown = {
           out["alpha", ] <- th[1, ]
           out["beta", ] <- exp(th[2, ])
           out["gamma", ] <- exp(th[3, ])
           out["delta", ] <- th[4, ] + centre
         })
  out
}

# Internal parameters from natural ones (par: alpha, beta, gamma, delta;
# vector or matrix [4, schedule]), the inverse of oldage_natural()
oldage_internal <- function(law, par, centre) {
  par <- matrix(par, 4)
  a <- par[1, ]
  b <- par[2, ]
  g <- par[3, ]
  dl <- par[4, ]
  switch(law,
         makeham = rbind(log(a) + b * centre, b, log(g)),
         weibull = rbind(log(a) + (b - 1) * log(centre), log(b)),
         beard = rbind(log(a) + b * centre, log(b), log(dl) + b * centre),
         logistic = ,
         perks = rbind(log(a) + b * centre, log(b), log(g), log(dl) + b * centre),
         logquad = rbind(a + b * centre + g * centre^2, b + 2 * g * centre, g),
         lynch_brown = rbind(a, log(b), log(g), dl - centre))
}

# Fit of any law to deaths d among N persons reaching ages x (matrices [age,
# schedule]); schedules with missing, negative or no numbers are not fitted.
# Returns law, natural parameters par [4, schedule], log-likelihood, ages,
# convergence, and what oldage_fit_H() needs
oldage_fit <- function(law, x, N, d) {
  n <- ncol(N)
  ok <- colSums(!is.finite(N) | !is.finite(d) | N < 0 | d < 0) == 0 & colSums(N) > 0
  par <- matrix(NA_real_, 4, n, dimnames = list(c("alpha", "beta", "gamma", "delta"), NULL))
  ll <- rep(NA_real_, n)
  conv <- rep(FALSE, n)
  fit <- list(law = law, ages = x)
  if (law %in% c("gompertz", "kannisto")) {
    fit$alpha <- fit$b <- rep(NA_real_, n)
    if (any(ok)) {
      f <- tail_fit_nd(x, N[, ok, drop = FALSE], d[, ok, drop = FALSE], law)
      fit$alpha[ok] <- f$alpha
      fit$b[ok] <- f$b
      conv[ok] <- f$converged
      ll[ok] <- oldage_loglik(tail_hazard(f$alpha, f$b, x, law)$H, N[, ok, drop = FALSE], d[, ok, drop = FALSE])
    }
    par["alpha", ] <- exp(fit$alpha)
    par["beta", ] <- fit$b
  } else {
    fit$theta <- matrix(NA_real_, oldage_npar[[law]], n)
    fit$centre <- mean(x) + 0.5
    if (any(ok)) {
      f <- oldage_fit_generic(law, x, N[, ok, drop = FALSE], d[, ok, drop = FALSE])
      fit$theta[, ok] <- f$theta
      ll[ok] <- f$loglik
      conv[ok] <- f$converged
    }
    par[] <- oldage_natural(law, fit$theta, fit$centre)
  }
  ll[!is.finite(ll)] <- NA_real_
  fit$par <- par
  fit$loglik <- ll
  fit$converged <- conv
  fit
}

# Integrated hazard of a fitted law over [x, x + 1] for the ages x
# (matrix [age, schedule])
oldage_fit_H <- function(fit, x) {
  if (fit$law %in% c("gompertz", "kannisto")) return(tail_hazard(fit$alpha, fit$b, x, fit$law)$H)
  oldage_H(fit$law, fit$theta, x - fit$centre, fit$centre)$H
}

# Fits of several laws and the choice of each schedule by the lowest AIC
# (Feehan 2018, eq. 6): AIC = -2 loglik + 2 k
oldage_fit_aic <- function(x, N, d, laws = oldage_laws) {
  fits <- lapply(stats::setNames(laws, laws), function(l) oldage_fit(l, x, N, d))
  aic <- vapply(fits, function(f) -2 * f$loglik + 2 * oldage_npar[[f$law]], numeric(ncol(N)))
  aic <- matrix(aic, ncol(N))
  choice <- apply(aic, 1, function(a) if (all(is.na(a))) NA_character_ else laws[which.min(a)])
  list(law = "aic", fits = fits, choice = choice, aic = aic, ages = x,
       converged = vapply(seq_len(ncol(N)), function(j) {
         if (is.na(choice[j])) FALSE else fits[[choice[j]]]$converged[j]
       }, TRUE))
}

# Integrated hazards of the law chosen for each schedule
oldage_aic_H <- function(fit, x) {
  H <- matrix(NA_real_, length(x), length(fit$choice))
  for (l in unique(stats::na.omit(fit$choice))) {
    j <- which(fit$choice == l)
    H[, j] <- oldage_fit_H(fit$fits[[l]], x)[, j, drop = FALSE]
  }
  H
}

#' Check old-age mortality models out of sample
#'
#' Fits the hazard laws of [life_table()] to the deaths at ages `from` to
#' `cut - 1` and scores their extrapolation against the observed deaths at
#' ages `cut` and above, the design of Thatcher, Kannisto and Vaupel (1998,
#' chapter 6): they fit six models at ages 80-98 and compare the predicted
#' probabilities of dying with the observed ones at ages 99 and over. The
#' predictions are based only on the data below the cut, so they are
#' independent of the observations they are compared with.
#'
#' Every law is fitted by maximum likelihood with binomial deaths: of the
#' `N_x` persons who reach age `x`, `d_x` die before `x + 1` with
#' probability `q_x = 1 - exp(-H_x)`, where `H_x` is the integral of the
#' hazard over the year (Thatcher et al. 1998, chapter 6; Feehan 2018,
#' equation 4). The deaths are those of the data and `N = E + d / 2`, the
#' person-years `E` plus half the deaths, as for [life_table()] with
#' `exposure` (Thatcher et al. use the numbers reaching each age directly).
#' At the test ages, with the predicted `q`, Thatcher et al. compare the
#' models by
#'
#' * `nll`: minus the binomial log-likelihood of the observed deaths,
#'   `-(d log q + (N - d) log(1 - q))`, without the binomial coefficients,
#'   which are the same for all models (their Table 6.3); the model with the
#'   lowest sum fits best on the likelihood criterion (`rank`);
#' * `chisq`: `(d - N q)^2 / (N q (1 - q))`, distributed like chi-squared
#'   with one degree of freedom per age for a correct model and accurate
#'   data (their Table 6.4), so that the sum over `n_ages` ages has expected
#'   value `n_ages`.
#'
#' The columns `aic` and `delta_aic` give Akaike's information criterion of
#' the fit below the cut, `-2 loglik + 2 k` with `k` parameters, and its
#' difference to the lowest value of the schedule, as Feehan (2018, section
#' 4.2) compares the models; without the binomial coefficients, which cancel
#' in the differences. Thatcher et al. pool the data of decades and 13
#' countries: pool periods with [mt_regroup()] (`period = `) before the
#' check. Cells whose deaths exceed the numbers at risk (`d > N`, possible
#' for period data with tiny exposures) are left out (package-specific
#' safeguard).
#'
#' @param data An [mt_data] object with deaths and exposures in single years
#'   of age from `from` to `to`.
#' @param cut First age of the test ages; the laws are fitted to the ages
#'   `from` to `cut - 1` (Thatcher et al. 1998: 99, fitting ages 80-98).
#' @param from First age of the fit (default 80, the lowest age of the data
#'   of Thatcher et al.).
#' @param to Last test age (default: the last single year of age of the
#'   data).
#' @param laws Hazard laws to check, see [life_table()] (default: all nine
#'   models of Feehan 2018).
#' @param by_age Return one row per test age instead of the sums.
#' @return A `data.table` with the stratum columns, `period` and `law`,
#'   `n_par` (number of parameters), `loglik`, `aic` and `delta_aic` of the
#'   fit, the out-of-sample sums `nll` and `chisq` over `n_ages` test ages,
#'   `rank` (by `nll` within each stratum and period, 1 = best) and
#'   `converged`. With `by_age = TRUE`, one row per test age with `age`,
#'   `N`, `d`, the observed `q_obs = d / N`, the predicted `q_fit`, `nll`
#'   and `chisq`. No thresholds: the check shows which law extrapolated best
#'   in these data.
#' @references
#' Feehan DM (2018). Separating the signal from the noise: evidence for
#' deceleration in old-age death rates. *Demography* 55(6): 2025-2044.
#' \doi{10.1007/s13524-018-0728-x}
#'
#' Thatcher AR, Kannisto V, Vaupel JW (1998). *The Force of Mortality at Ages
#' 80 to 120*. Odense University Press.
#' @seealso [life_table()] for the laws and their use as a closure of the
#'   open age interval.
#' @examples
#' # expected deaths of a cohort under a Kannisto hazard: the laws fitted at
#' # 80-94 are scored at ages 95-104
#' x <- 80:104
#' H <- log1p(5e-5 * exp(0.11 * (x + 1))) / 0.11 - log1p(5e-5 * exp(0.11 * x)) / 0.11
#' N <- 1e5 * exp(-cumsum(c(0, H[-length(H)])))
#' d <- round(N * (1 - exp(-H)))
#' sim <- mt_data(deaths = data.frame(age = x, year = 2000, deaths = d),
#'                exposure = data.frame(age = x, year = 2000, exposure = N - d / 2))
#' chk <- mt_open_check(sim, cut = 95, laws = c("gompertz", "kannisto", "logquad", "beard"))
#' chk[, list(law, nll, chisq, n_ages, rank)]
#' @export
mt_open_check <- function(data, cut = 99, from = 80, to = NULL,
                          laws = c("gompertz", "kannisto", "makeham", "weibull", "beard", "logquad", "logistic",
                                   "perks", "lynch_brown"),
                          by_age = FALSE) {
  check_class(data, "mt_data", "data")
  check_number(cut, "cut", integerish = TRUE)
  check_flag(by_age, "by_age")
  if (!is.character(laws) || !length(laws) || !all(laws %in% oldage_laws) || anyDuplicated(laws)) {
    mt_abort(sprintf("`laws` must name laws among %s.", paste0("\"", oldage_laws, "\"", collapse = ", ")),
             class = "argument")
  }
  # more fitted ages than parameters
  check_number(from, "from", lower = 1, upper = cut - max(oldage_npar[laws]) - 1, integerish = TRUE)
  g <- mt_grid(data)
  single <- g$age$lower[is.finite(g$age$upper) & abs(g$age$upper - g$age$lower - 1) < 1e-9]
  to <- to %||% max(single)
  check_number(to, "to", lower = cut, integerish = TRUE)
  ages <- from:to
  ia <- match(ages, g$age$lower)
  if (anyNA(ia) || !all(ages %in% single)) {
    mt_abort(sprintf("mt_open_check() needs single years of age from %d to %d in the grid of `data`.", from, to),
             class = "argument")
  }
  D <- mt_array(data, "deaths", coarse = "na")[ia, , , drop = FALSE]
  E <- mt_array(data, "exposure", coarse = "na")[ia, , , drop = FALSE]
  np <- dim(D)[2]
  ns <- dim(D)[3]
  D <- matrix(D, length(ages))
  E <- matrix(E, length(ages))
  N <- E + D / 2
  # cells with more deaths than persons at risk are left out
  drop <- !is.finite(N) | !is.finite(D) | D > N
  N[drop] <- NA_real_
  D[drop] <- NA_real_
  fit_rows <- ages < cut
  test_rows <- which(!fit_rows)
  xf <- ages[fit_rows]
  cols <- data.table(stratum = rep(seq_len(ns), each = np), period = rep(g$period$lower, ns))
  Nt <- N[test_rows, , drop = FALSE]
  dt <- D[test_rows, , drop = FALSE]
  has <- is.finite(Nt) & is.finite(dt) & Nt > 0
  res <- lapply(laws, function(l) {
    f <- oldage_fit(l, xf, N[fit_rows, , drop = FALSE], D[fit_rows, , drop = FALSE])
    qf <- -expm1(-oldage_fit_H(f, ages[test_rows]))
    # minus the binomial log-likelihood and the chi-squared of each test age
    nl <- -(ifelse(dt > 0, dt * log(qf), 0) + ifelse(Nt > dt, (Nt - dt) * log1p(-qf), 0))
    cs <- (dt - Nt * qf)^2 / (Nt * qf * (1 - qf))
    if (by_age) {
      keep <- as.vector(has)
      return(data.table(cols[rep(seq_len(nrow(cols)), each = length(test_rows))], law = l,
                        age = rep(ages[test_rows], nrow(cols)), N = as.vector(Nt), d = as.vector(dt),
                        q_obs = as.vector(dt / Nt), q_fit = as.vector(qf), nll = as.vector(nl),
                        chisq = as.vector(cs))[keep])
    }
    # a predicted probability of 0 or 1 that the deaths contradict gives an
    # infinite sum; a prediction of 1 with all persons dying contributes 0
    cs[is.nan(cs) & has & !is.na(qf)] <- 0
    nl[!has] <- 0
    cs[!has] <- 0
    fitted <- !is.na(f$loglik)
    data.table(cols, law = l, n_par = oldage_npar[[l]], loglik = f$loglik,
               aic = -2 * f$loglik + 2 * oldage_npar[[l]],
               nll = ifelse(fitted, colSums(nl), NA_real_), chisq = ifelse(fitted, colSums(cs), NA_real_),
               n_ages = colSums(has), converged = f$converged)
  })
  out <- rbindlist(res)
  if (!by_age) {
    out[, delta_aic := if (all(is.na(aic))) NA_real_ else aic - min(aic, na.rm = TRUE), by = c("stratum", "period")]
    out[, rank := as.integer(rank(ifelse(is.na(nll), Inf, nll), ties.method = "min")), by = c("stratum", "period")]
    out[is.na(nll), rank := NA_integer_]
    setcolorder(out, c("stratum", "period", "law", "n_par", "loglik", "aic", "delta_aic", "nll", "chisq", "n_ages",
                       "rank", "converged"))
  }
  out <- merge(data$strata, out, by = "stratum", sort = FALSE)
  out[, law := factor(law, laws)]
  setorderv(out, c("stratum", "period", "law", if (by_age) "age"))
  out[, law := as.character(law)]
  out[]
}
