# Two-dimensional P-spline models of log death rates over age and year:
# p_splines() (Currie, Durban and Eilers 2004) and cp_splines() (Camarda
# 2019).
#
# The P-spline engine is adapted from MortalitySmooth 2.3.4 (Mort2Dsmooth(),
# Mort2Dsmooth_estimate(), Mort2Dsmooth_update(), Mort2Dsmooth_se(),
# predict.Mort2Dsmooth(), MortSmooth_bbase(), MortSmooth_BWB() and
# MortSmooth_BcoefB() in R/), GPL (>= 2), by Carlo G. Camarda: B-spline bases
# with equally spaced knots on the range of the data widened by 1% at each
# end, difference penalties, the scoring algorithm in array form (row tensors
# of the marginal bases), the effective dimension and the information
# criteria, the overdispersion parameter, and the forecast by a basis over
# years extended on the same knot spacing with weights 0 in the forecast
# years and standard errors from the inverse of B'WB + P. Changes: B-splines
# from splines::splineDesign() on the same knots; row tensors restricted to
# the pairs of B-splines whose supports overlap; ages at the cell ages of
# cell_ages(); knots about five years of age apart for age groups wider than
# one year; the smoothing parameters by a grid search with steps of 1 on the
# log10 scale followed by a continuous search (MortalitySmooth: two greedy
# grid searches with steps of 2 and 0.5); the Poisson deviance with 0 log 0
# = 0 (MortalitySmooth replaces zero counts by 1e-4); step halving when the
# penalised deviance increases; with overdispersion, the
# forecast keeps the dispersion in the penalty and the covariance is psi2
# (B'WB + psi2 P)^-1 (MortalitySmooth's predict() refits with P and
# multiplies the standard errors by sqrt(psi2)); draws from the normal
# approximation of the coefficients.
#
# cp_splines() is written from the paper (Camarda 2019, sections 2 and 3);
# the author's R code (supplementary material of the paper, no licence
# stated) was used only to check details and results.
#
# The dispersion by age of p_splines(overdispersion = "age") is written
# from Djeundje and Currie (2011, sections 2 to 4 and 5.2); no code of the
# authors was used.

# Bases and array arithmetic ---------------------------------------------------

# B-spline basis of degree `deg` at the points `x`, with `ndx` intervals
# between equally spaced knots on the range of `x` widened by `margin` times
# the range at each end (MortalitySmooth: 0.01; Camarda 2019's code: 0).
# `C` holds the first derivatives of the B-splines. A single point gives a
# constant (one coefficient, no penalty).
ps_basis <- function(x, ndx, deg = 3, margin = 0.01) {
  xl <- min(x)
  xr <- max(x)
  if (xr - xl < 1e-9) {
    return(list(B = matrix(1, length(x), 1), C = matrix(0, length(x), 1), knots = NULL, dx = NA_real_, deg = deg))
  }
  lo <- xl - margin * (xr - xl)
  hi <- xr + margin * (xr - xl)
  dx <- (hi - lo) / ndx
  knots <- lo + seq(-deg, ndx + deg) * dx
  ps_basis_knots(x, knots, dx, deg)
}

ps_basis_knots <- function(x, knots, dx, deg) {
  B <- splines::splineDesign(knots, x, ord = deg + 1, outer.ok = TRUE)
  C <- splines::splineDesign(knots, x, ord = deg + 1, derivs = rep(1L, length(x)), outer.ok = TRUE)
  list(B = B, C = C, knots = knots, dx = dx, deg = deg)
}

# The basis of `b` at points `x` beyond its range, with knots on the same
# spacing as predict.Mort2Dsmooth() places them: from the deg-th of the
# knots lo - deg dx, lo - (deg + 1) dx, ... at or below min(x) to the (deg +
# 1)-th of the knots hi + deg dx, hi + (deg + 1) dx, ... at or above max(x),
# where [lo, hi] is the range of the knots of `b`. The B-splines of `b`
# are columns offset + 1, offset + 2, ... of the result.
ps_extend <- function(b, x) {
  if (is.null(b$knots)) mt_abort("Forecasting needs at least two periods.", class = "model_fit")
  deg <- b$deg
  dx <- b$dx
  lo <- b$knots[deg + 1]
  hi <- b$knots[length(b$knots) - deg]
  steps <- deg:(deg + ceiling((max(x) - min(x)) / dx) + 2 * deg + 2)
  below <- lo - steps * dx
  above <- hi + steps * dx
  first <- below[which(below <= min(x) + 1e-9)[deg]]
  last <- above[which(above >= max(x) - 1e-9)[deg + 1]]
  knots <- first + seq(0, round((last - first) / dx)) * dx
  out <- ps_basis_knots(x, knots, dx, deg)
  out$offset <- as.integer(round((b$knots[1] - first) / dx))
  out
}

# Row tensor of a basis, restricted to the pairs of columns i <= j whose
# supports overlap (the other products are zero): columns B[, i] * B[, j]
rt_pairs <- function(B) {
  nz <- crossprod(B != 0) > 0
  ij <- which(nz & upper.tri(nz, diag = TRUE), arr.ind = TRUE)
  list(i = ij[, 1], j = ij[, 2], RT = B[, ij[, 1], drop = FALSE] * B[, ij[, 2], drop = FALSE], k = ncol(B))
}

# Positions in the matrix B'WB [ca cy, ca cy] (coefficients ordered with age
# fastest) of every product of an age pair and a year pair from rt_pairs():
# four entries each, for both orders of both pairs
glam_map <- function(pa, py) {
  ca <- pa$k
  p <- ca * py$k
  na <- length(pa$i)
  a1 <- rep(pa$i, length(py$i))
  a2 <- rep(pa$j, length(py$i))
  y1 <- rep(py$i, each = na)
  y2 <- rep(py$j, each = na)
  pos <- function(a, y) a + (y - 1) * ca
  idx <- c(pos(a1, y1) + (pos(a2, y2) - 1) * p, pos(a2, y1) + (pos(a1, y2) - 1) * p,
           pos(a1, y2) + (pos(a2, y1) - 1) * p, pos(a2, y2) + (pos(a1, y1) - 1) * p)
  # a pair of equal columns gives the same entry twice
  keep <- !duplicated(idx)
  list(p = p, idx = idx[keep], src = rep(seq_along(a1), 4)[keep])
}

# Inner product (Ba (x) By)' diag(vec(W)) (Ba (x) By) of the Kronecker basis
# from the row tensors of the marginal bases (MortSmooth_BWB()), computed
# for the overlapping pairs only: its nonzero entries, positions `idx` and
# values `val` (added to dense matrices in place by the callers, which saves
# copies of the matrix). `pa`, `py`: rt_pairs(); `map`: glam_map()
glam_bwb <- function(pa, py, W, map = glam_map(pa, py)) {
  M <- crossprod(pa$RT, W %*% py$RT)
  list(idx = map$idx, val = M[map$src])
}

# Ba A By' for the coefficient matrix A (MortSmooth_BcoefB())
glam_eta <- function(Ba, By, alpha) Ba %*% matrix(alpha, ncol(Ba)) %*% t(By)

# Difference penalty of order `pord` over `k` coefficients (zero for one
# coefficient)
diff_penalty <- function(k, pord) {
  if (k <= pord) return(matrix(0, k, k))
  crossprod(diff(diag(k), differences = pord))
}

# Default numbers of intervals between knots: one per five ages (years), as
# MortalitySmooth's floor(n / 5) for single years; for age groups, n is the
# number of single years the cell ages span
ps_ndx <- function(x) max(1L, floor((round(max(x) - min(x)) + 1) / 5))

# Solution of A x = b for a symmetric positive definite A, with its Cholesky
# factor; a generalised inverse where the factorisation fails
chol_solve <- function(A, b) {
  R <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(R)) return(list(x = as.numeric(MASS_ginv(A) %*% b), R = NULL))
  list(x = backsolve(R, forwardsolve(R, b, upper.tri = TRUE, transpose = TRUE)), R = R)
}

# Fitting --------------------------------------------------------------------

# Penalised Poisson scoring algorithm of the two-dimensional P-spline model in
# array form. Y, O and wt are [age, year] matrices of deaths, log exposures
# and weights (0: cell without data, such as a missing cell or a forecast
# year); `ma` and `my` the marginal bases with their row tensors (`pr`,
# rt_pairs()); P the penalty matrix (smoothing parameters and dispersion
# included). `extra(alpha)` gives asymmetric penalties (cp_splines()):
# entries of the system matrix (lists `idx`, `val`) and `p`, added to the
# right-hand side; recomputed at every step, without step halving (as in
# Camarda 2019). `phi` [age, year]: dispersions of the cells of the
# quasi-likelihood (Djeundje and Currie 2011, eq. 2.1 and 3.13), which
# divide the weights and the deviance. Convergence: largest change of a
# coefficient below `tol` (relative to its size where it exceeds 1).
# Returns the coefficients, the fitted log rates, the weighted inner
# product B'WB and the Cholesky factor of B'WB + P (+ extra).
ps2_irls <- function(Y, O, wt, ma, my, P, alpha, maxit = 50, tol = 1e-6, extra = NULL, phi = NULL) {
  map <- glam_map(ma$pr, my$pr)
  obs <- wt > 0
  pw <- if (is.null(phi)) wt else wt / phi
  pdev <- function(eta, a) {
    mu <- exp(O + eta)[obs]
    if (any(!is.finite(mu))) return(Inf)
    dv <- if (is.null(phi)) poisson_dev(Y[obs], mu) else sum(poisson_dev_cells(Y[obs], mu) / phi[obs])
    dv + sum(a * (P %*% a))
  }
  # B'WB + P (+ asymmetric penalties at the coefficients `a`), filled in place
  system <- function(e, a) {
    A <- P
    A[e$idx] <- A[e$idx] + e$val
    ex <- NULL
    if (!is.null(extra)) {
      ex <- extra(a)
      for (k in seq_along(ex$idx)) A[ex$idx[[k]]] <- A[ex$idx[[k]]] + ex$val[[k]]
    }
    list(A = A, p = ex$p)
  }
  eta <- glam_eta(ma$B, my$B, alpha)
  obj <- pdev(eta, alpha)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    mu <- exp(O + eta)
    W <- pw * mu
    z <- eta + (Y - mu) / mu
    z[!obs] <- 0
    rhs <- as.vector(crossprod(ma$B, W * z) %*% my$B)
    s <- system(glam_bwb(ma$pr, my$pr, W, map), alpha)
    if (!is.null(s$p)) rhs <- rhs + s$p
    new <- chol_solve(s$A, rhs)$x
    if (anyNA(new)) mt_abort("The P-spline fit failed (singular system).", class = "model_fit")
    step <- new - alpha
    eta_new <- glam_eta(ma$B, my$B, new)
    if (is.null(extra)) {
      obj_new <- pdev(eta_new, new)
      halvings <- 0L
      while (!(obj_new <= obj + 1e-10 * abs(obj)) && halvings < 30L) {
        step <- step / 2
        new <- alpha + step
        eta_new <- glam_eta(ma$B, my$B, new)
        obj_new <- pdev(eta_new, new)
        halvings <- halvings + 1L
      }
      obj <- obj_new
    }
    if (any(!is.finite(eta_new))) mt_abort("The P-spline fit diverged.", class = "model_fit")
    alpha <- new
    eta <- eta_new
    if (max(abs(step) / pmax(abs(alpha), 1)) < tol) {
      converged <- TRUE
      break
    }
  }
  mu <- exp(O + eta)
  e <- glam_bwb(ma$pr, my$pr, pw * mu, map)
  BWB <- matrix(0, map$p, map$p)
  BWB[e$idx] <- e$val
  R <- tryCatch(chol(system(e, alpha)$A), error = function(e) NULL)
  list(alpha = alpha, eta = eta, mu = mu, BWB = BWB, R = R, converged = converged, iterations = it)
}

# Effective dimension (trace of the hat matrix) of a fit from ps2_irls()
ps2_ed <- function(f, P) {
  Ainv <- if (!is.null(f$R)) chol2inv(f$R) else MASS_ginv(f$BWB + P)
  sum(Ainv * f$BWB)
}

# Poisson deviance of each cell (with 0 log 0 = 0)
poisson_dev_cells <- function(y, mu) 2 * (ifelse(y > 0, y * log(y / mu), 0) - (y - mu))

# Diagonal of the hat matrix [age, year] of a fit from ps2_irls() with cell
# weights W (Djeundje and Currie 2011, eq. 3.17 and 4.21): W b' (B'WB +
# P)^-1 b for the rows b = By[t, ] (x) Ba[x, ] of the Kronecker basis, as
# Ba2 M t(By2) with the row tensors Ba2, By2 of the marginal bases and M the
# inverse rearranged over pairs of age and pairs of year coefficients
ps2_leverage <- function(f, P, Ba, By, W) {
  Ainv <- if (!is.null(f$R)) chol2inv(f$R) else MASS_ginv(f$BWB + P)
  ca <- ncol(Ba)
  cy <- ncol(By)
  rt <- function(B) {
    k <- seq_len(ncol(B))
    B[, rep(k, each = length(k)), drop = FALSE] * B[, rep(k, length(k)), drop = FALSE]
  }
  M <- matrix(aperm(array(Ainv, c(ca, cy, ca, cy)), c(3, 1, 4, 2)), ca * ca, cy * cy)
  W * (rt(Ba) %*% M %*% t(rt(By)))
}

# Marginal bases with row tensors for ages (with the specialised coefficient
# for age 0 of Camarda 2019, section 3.1, when `infant`) and years
ps_age_basis <- function(xa, ndx, deg, margin, infant) {
  if (!infant) {
    b <- ps_basis(xa, ndx, deg, margin)
    b$infant <- FALSE
  } else {
    b0 <- ps_basis(xa[-1], ndx, deg, margin)
    k <- ncol(b0$B)
    B <- rbind(c(1, numeric(k)), cbind(0, b0$B))
    # rate of ageing at age 0: forward difference to the first age of the
    # surface (Camarda 2019's code; per year of age)
    C <- rbind((B[2, ] - c(1, numeric(k))) / (xa[2] - xa[1]), cbind(0, b0$C))
    b <- list(B = B, C = C, knots = b0$knots, dx = b0$dx, deg = deg, infant = TRUE)
  }
  b$pr <- rt_pairs(b$B)
  b
}

# Penalties over ages and years of the coefficient matrix [ca, cy]
# (Kronecker form of Currie, Durban and Eilers 2004); with the specialised
# coefficient for age 0, the first element of the difference matrix over age
# is 0 (Camarda 2019, section 3.1)
ps_penalties <- function(ca, cy, pord, infant = FALSE) {
  Da <- if (ca > pord[1]) diff(diag(ca), differences = pord[1]) else matrix(0, 1, ca)
  if (infant) Da[1, 1] <- 0
  list(Pa = kronecker(diag(cy), crossprod(Da)), Py = kronecker(diff_penalty(cy, pord[2]), diag(ca)))
}

# Starting coefficients: penalised least squares of log((D + 1) / E) on the
# bases (MortalitySmooth), with a plane over age and year in cells without
# data
ps_start <- function(Y, E, wt, Ba, By, xa, xy) {
  eta0 <- log((Y + 1) / E)
  obs <- wt > 0 & is.finite(eta0)
  if (!all(obs)) {
    g <- data.frame(eta = eta0[obs], a = rep(xa, length(xy))[obs], y = rep(xy, each = length(xa))[obs])
    cf <- stats::coef(stats::lm(eta ~ a + y, data = g))
    cf[is.na(cf)] <- 0
    plane <- outer(xa, xy, function(a, y) cf[1] + cf[2] * a + cf[3] * y)
    eta0[!obs] <- plane[!obs]
  }
  ia <- solve(crossprod(Ba) + diag(1e-6, ncol(Ba)), t(Ba))
  iy <- solve(crossprod(By) + diag(1e-6, ncol(By)), t(By))
  as.vector(ia %*% eta0 %*% t(iy))
}

# Data of one stratum for the P-spline models: deaths and log exposures with
# weight 0 for missing cells and cells without exposure
ps_data <- function(D, E) {
  wt <- matrix(as.numeric(!is.na(D) & !is.na(E) & E > 0), nrow(D))
  if (!any(wt > 0)) mt_abort("No cells with deaths and positive exposures.", class = "model_fit")
  Y <- ifelse(wt > 0, D, 0)
  O <- ifelse(wt > 0, log(pmax(E, 1e-300)), 0)
  list(Y = Y, O = O, wt = wt, E = ifelse(wt > 0, E, 1))
}

# Two-dimensional P-spline fit with smoothing parameters given or selected
# by an information criterion. With `overdispersion = TRUE`, the dispersion
# psi2 = deviance / (cells - ed) enters the penalty (B'WB + psi2 P) and the
# criterion (deviance / psi2), and selection and estimation of psi2 are
# iterated (at most 5 times; Mort2Dsmooth(overdispersion = TRUE)). Without,
# psi2 is reported after the fit. With `overdispersion = "age"`, a
# dispersion per age that varies smoothly over age divides the weights and
# the deviance (ps2_dispersion_age(); Djeundje and Currie 2011).
ps2_model <- function(D, E, xa, xy, ndx, deg, pord, margin, infant, lambda, criterion, lambda_range,
                      overdispersion) {
  dat <- ps_data(D, E)
  ma <- ps_age_basis(xa, ndx[1], deg[1], margin, infant)
  my <- ps_basis(xy, ndx[2], deg[2], margin)
  my$pr <- rt_pairs(my$B)
  ca <- ncol(ma$B)
  cy <- ncol(my$B)
  pen <- ps_penalties(ca, cy, pord, infant)
  n_obs <- sum(dat$wt > 0)
  k_ic <- if (criterion == "bic") log(n_obs) else 2
  start <- ps_start(dat$Y, dat$E, dat$wt, ma$B, my$B, xa, xy)
  state <- new.env(parent = emptyenv())
  state$alpha <- start
  obs <- dat$wt > 0
  # fit with smoothing parameters `lam` and the dispersion psi2 or the
  # dispersions `phi` by age; the criterion of the deviance scaled by the
  # dispersion (Djeundje and Currie 2011, eq. 3.16)
  fit_at <- function(lam, psi2, phi = NULL) {
    P <- psi2 * (lam[1] * pen$Pa + lam[2] * pen$Py)
    Phi <- if (!is.null(phi)) matrix(phi, nrow(dat$Y), ncol(dat$Y))
    f <- ps2_irls(dat$Y, dat$O, dat$wt, ma, my, P, state$alpha, phi = Phi)
    if (f$converged) state$alpha <- f$alpha
    f$ed <- ps2_ed(f, P)
    f$dev <- poisson_dev(dat$Y[obs], f$mu[obs])
    sdev <- if (is.null(phi)) f$dev / psi2 else sum(poisson_dev_cells(dat$Y[obs], f$mu[obs]) / Phi[obs])
    f$ic <- sdev + k_ic * f$ed
    f$lambda <- lam
    f$psi2 <- psi2
    f$phi <- phi
    f
  }
  dispersion <- function(f) f$dev / max(n_obs - f$ed, 1)
  at_bound <- FALSE
  if (identical(overdispersion, "age")) {
    f <- ps2_dispersion_age(fit_at, dat, ma, my, pen, pord, infant, lambda, lambda_range, ca > 1)
    at_bound <- f$at_bound
  } else if (!is.null(lambda)) {
    f <- fit_at(lambda, 1)
    if (overdispersion) {
      for (i in 1:5) {
        psi2 <- dispersion(f)
        f2 <- fit_at(lambda, psi2)
        done <- abs(dispersion(f2) - psi2) / psi2 < 1e-3
        f <- f2
        if (done) break
      }
    }
  } else {
    psi2 <- 1
    for (round in seq_len(if (overdispersion) 5L else 1L)) {
      sel <- ps_select_lambda(function(ll) fit_at(10^ll, psi2)$ic, log10(lambda_range), ca > 1)
      f <- fit_at(10^sel$ll, psi2)
      if (!overdispersion) break
      psi2_new <- dispersion(f)
      done <- abs(psi2_new - psi2) / psi2_new < 1e-3
      psi2 <- psi2_new
      if (done) break
    }
    if (overdispersion) f <- fit_at(f$lambda, psi2)
    # a smoothing parameter at an end of the range matters only where the
    # criterion improves beyond it (two decades further) by more than 1;
    # otherwise the fit has reached its limit there (no penalty, or the
    # polynomial of the null space)
    if (any(sel$bound != 0)) {
      at_bound <- fit_at(f$lambda * 100^sel$bound, f$psi2)$ic < f$ic - 1
    }
  }
  if (!f$converged) mt_warn("The P-spline fit did not converge.", class = "model_fit", call = NULL)
  f$psi2_estimate <- if (is.null(f$phi)) dispersion(f) else f$phi
  f$at_bound <- NULL
  c(f, list(data = dat, ma = ma, my = my, pen = pen, xa = xa, xy = xy, at_bound = at_bound, n_obs = n_obs,
            overdispersion = overdispersion))
}

# Dispersion that varies smoothly over age (Djeundje and Currie 2011, Model
# 3 of section 5.2, bias corrected scheme of section 4). From phi = 1, two
# steps alternate: the fit of the rates with the dispersions in the weights
# (eq. 3.13), its smoothing parameters selected by the criterion of the
# scaled deviance (eq. 3.16) unless given; and the dispersions, the
# bias-corrected mean Pearson statistics of the ages (eq. 4.20, with the
# contributions of the ages to the effective dimension, eq. 4.21) smoothed
# over age (ps_disp_smooth()). The iteration stops when no dispersion
# changes by more than 1e-5 relative to its value (section 5.2.2), after at
# most 50 rounds. `fit_at(lambda, psi2, phi)`: the fit of ps2_model().
ps2_dispersion_age <- function(fit_at, dat, ma, my, pen, pord, infant, lambda, lambda_range, two) {
  obs <- dat$wt > 0
  n_x <- rowSums(obs)
  Dp <- ps_penalties(ncol(ma$B), 1, pord, infant)$Pa
  st <- new.env(parent = emptyenv())
  phi_step <- function(f) {
    P <- f$lambda[1] * pen$Pa + f$lambda[2] * pen$Py
    v <- rowSums(ps2_leverage(f, P, ma$B, my$B, dat$wt * f$mu / f$phi))
    pearson <- rowSums(ifelse(obs, (dat$Y - f$mu)^2 / f$mu, 0))
    ok <- n_x > 0 & n_x - v > 0
    st$raw <- ifelse(ok, pearson / (n_x - v), NA_real_)
    st$disp <- ps_disp_smooth(ifelse(ok, st$raw, 1), ifelse(ok, n_x, 0), ma$B, Dp, lambda_range, st$disp$beta,
                              st$disp$tau %||% 1)
    st$disp$phi
  }
  phi <- rep(1, nrow(dat$Y))
  converged <- FALSE
  for (round in seq_len(50)) {
    if (is.null(lambda)) {
      sel <- ps_select_lambda(function(ll) fit_at(10^ll, 1, phi)$ic, log10(lambda_range), two)
      f <- fit_at(10^sel$ll, 1, phi)
    } else {
      f <- fit_at(lambda, 1, phi)
    }
    phi_new <- phi_step(f)
    converged <- max(abs(phi_new - phi) / phi_new) < 1e-5
    phi <- phi_new
    if (converged) break
  }
  if (!converged) mt_warn("The dispersion of the P-spline fit did not converge.", class = "model_fit", call = NULL)
  f <- fit_at(f$lambda, 1, phi)
  # as in ps2_model(): a smoothing parameter at an end of the range
  f$at_bound <- is.null(lambda) && any(sel$bound != 0) &&
    fit_at(f$lambda * 100^sel$bound, 1, phi)$ic < f$ic - 1
  c(f, list(dispersion_raw = st$raw, lambda_dispersion = st$disp$lambda, dispersion_tau = st$disp$tau,
            dispersion_rounds = round))
}

# Penalised quasi-likelihood fit of the dispersions over age (Djeundje and
# Currie 2011, eq. 3.14 and 4.22): responses `d` (bias-corrected mean
# Pearson statistics) with E[d] = phi, var(d) = tau phi^2 / n, link -1 /
# phi = B beta, and the penalty Pd (smoothing parameter included). `n`: the
# numbers of cells of the ages (0 for ages without data, whose dispersion
# the smooth fills in). Scoring with step halving on the penalised deviance;
# the dispersions stay positive at every age.
ps_disp_irls <- function(d, n, tau, B, Pd, beta, maxit = 100, tol = 1e-8) {
  ok <- n > 0
  w <- n / tau
  gamma_dev <- function(phi) 2 * (-log(d[ok] / phi[ok]) + (d[ok] - phi[ok]) / phi[ok])
  obj <- function(b) {
    eta <- as.vector(B %*% b)
    if (any(eta >= 0)) return(Inf)
    sum(w[ok] * gamma_dev(-1 / eta)) + sum(b * (Pd %*% b))
  }
  cur <- obj(beta)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    eta <- as.vector(B %*% beta)
    phi <- -1 / eta
    Wd <- w * phi^2
    z <- eta + ifelse(ok, (d - phi) / phi^2, 0)
    new <- chol_solve(crossprod(B, Wd * B) + Pd, as.vector(crossprod(B, Wd * z)))$x
    if (anyNA(new)) mt_abort("The smoothing of the dispersion failed (singular system).", class = "model_fit")
    step <- new - beta
    val <- obj(new)
    halvings <- 0L
    while (!(val <= cur + 1e-10 * abs(cur)) && halvings < 30L) {
      step <- step / 2
      new <- beta + step
      val <- obj(new)
      halvings <- halvings + 1L
    }
    beta <- new
    cur <- val
    if (max(abs(step) / pmax(abs(beta), 1e-3)) < tol) {
      converged <- TRUE
      break
    }
  }
  phi <- -1 / as.vector(B %*% beta)
  BWB <- crossprod(B, (w * phi^2) * B)
  A <- BWB + Pd
  Ainv <- tryCatch(chol2inv(chol(A)), error = function(e) MASS_ginv(A))
  list(beta = beta, phi = phi, ed = sum(Ainv * BWB), sdev = sum(w[ok] * gamma_dev(phi)), converged = converged)
}

# Smooth dispersions over age from the statistics `d` of the ages with `n`
# cells (Djeundje and Currie 2011, section 4 and eq. 5.28): the basis `B`
# over age of the rates and the difference penalty `Dp`; the smoothing
# parameter minimises the scaled BIC of this fit (deviance divided by tau
# plus log(ages) times its effective dimension) over `lambda_range`; tau is
# the mean Pearson statistic of this fit, sum n (d - phi)^2 / phi^2 divided
# by the number of ages minus the effective dimension (package-specific:
# the paper estimates tau without stating how), and selection and tau
# alternate as the constant dispersion in ps2_model() (at most 5 times,
# relative change 1e-3), from the starting values `beta` and `tau` (those
# of the previous round of ps2_dispersion_age()). Without a penalty (one
# age), phi = d.
ps_disp_smooth <- function(d, n, B, Dp, lambda_range, beta = NULL, tau = 1) {
  ok <- n > 0
  K <- sum(ok)
  st <- new.env(parent = emptyenv())
  st$beta <- beta %||% rep(-1 / stats::weighted.mean(d[ok], n[ok]), ncol(B))
  fit_at <- function(l2, tau) {
    f <- ps_disp_irls(d, n, tau, B, l2 * Dp, st$beta)
    if (f$converged) st$beta <- f$beta
    f$ic <- f$sdev + log(K) * f$ed
    f$lambda <- l2
    f$tau <- tau
    f
  }
  if (!any(Dp != 0)) return(fit_at(0, 1))
  for (round in 1:5) {
    l2 <- 10^ps_select_lambda(function(ll) fit_at(10^ll[2], tau)$ic, log10(lambda_range), FALSE)$ll[2]
    f <- fit_at(l2, tau)
    if (K - f$ed < 1) break
    tau_new <- sum(n[ok] * (d[ok] - f$phi[ok])^2 / f$phi[ok]^2) / (K - f$ed)
    done <- abs(tau_new - tau) / tau_new < 1e-3
    tau <- tau_new
    if (done) break
  }
  fit_at(l2, tau)
}

# Smoothing parameters on the log10 scale minimising `ic`: a greedy search on
# a grid with steps of 1 from the centre of the range (as svcm's
# cleversearch(), used by MortalitySmooth), then a continuous search within
# one step of the grid optimum, one coordinate at a time, twice.
# `two`: FALSE when the age direction has a single coefficient
ps_select_lambda <- function(ic, lr, two = TRUE) {
  cache <- new.env(parent = emptyenv())
  f <- function(ll) {
    key <- paste(sprintf("%.4f", ll), collapse = ",")
    if (!exists(key, envir = cache, inherits = FALSE)) {
      v <- ic(ll)
      assign(key, if (is.finite(v)) v else 1e300, envir = cache)
    }
    get(key, envir = cache, inherits = FALSE)
  }
  grid <- seq(lr[1], lr[2], by = 1)
  if (grid[length(grid)] < lr[2] - 1e-9) grid <- c(grid, lr[2])
  dims <- if (two) 1:2 else 2L
  idx <- rep(ceiling(length(grid) / 2), 2)
  if (!two) idx[1] <- 1L
  best <- f(grid[idx])
  repeat {
    moved <- FALSE
    for (d in dims) {
      for (j in c(idx[d] - 1L, idx[d] + 1L)) {
        if (j < 1L || j > length(grid)) next
        cand <- idx
        cand[d] <- j
        v <- f(grid[cand])
        if (v < best) {
          best <- v
          idx <- cand
          moved <- TRUE
        }
      }
    }
    if (!moved) break
  }
  ll <- grid[idx]
  for (round in 1:2) {
    for (d in dims) {
      rng <- c(max(lr[1], ll[d] - 1), min(lr[2], ll[d] + 1))
      opt <- stats::optimize(function(v) {
        l2 <- ll
        l2[d] <- v
        f(l2)
      }, rng, tol = 0.02)
      if (opt$objective < best) {
        best <- opt$objective
        ll[d] <- opt$minimum
      }
    }
  }
  bound <- (ll > lr[2] - 0.05) - (ll < lr[1] + 0.05)
  if (!two) bound[1] <- 0
  list(ll = ll, bound = bound)
}

# Coefficients over the extended year basis: those of the fit and, for the
# `left` B-splines added below and the ones added above, a continuation with
# differences of order `pord` equal to zero (starting values of the
# forecast fits)
ps_continue <- function(alpha, ca, cy_new, pord, left = 0L) {
  A <- matrix(alpha, ca)
  next_col <- function(A, cols) {
    d <- min(pord, length(cols))
    w <- (-1)^(seq_len(d) + 1) * choose(d, seq_len(d))
    A[, cols[seq_len(d)], drop = FALSE] %*% w
  }
  for (i in seq_len(left)) A <- cbind(next_col(A, seq_len(ncol(A))), A)
  while (ncol(A) < cy_new) A <- cbind(A, next_col(A, rev(seq_len(ncol(A)))))
  as.vector(A)
}

# Forecast fit of a P-spline model: the basis over years extended to `h`
# more years (ps_extend(), as predict.Mort2Dsmooth()), the forecast years
# with weights 0, the coefficients refitted with the smoothing parameters
# (and dispersion, or dispersions by age) of the fit. Returns the extended
# fit with its year basis `my` and penalty `P`.
ps2_extend_fit <- function(fit, h, pord, extra_fun = NULL, tol = 1e-6, maxit = 50, Y = NULL, start = NULL) {
  dat <- fit$data
  ny <- ncol(dat$Y)
  # single calendar years
  my <- ps_extend(fit$my, c(fit$xy, fit$xy[ny] + seq_len(h)))
  my$pr <- rt_pairs(my$B)
  ca <- ncol(fit$ma$B)
  cy <- ncol(my$B)
  pen <- ps_penalties(ca, cy, pord, isTRUE(fit$ma$infant))
  P <- fit$psi2 * (fit$lambda[1] * pen$Pa + fit$lambda[2] * pen$Py)
  pad <- function(M, v) cbind(M, matrix(v, nrow(M), h))
  Yx <- pad(Y %||% dat$Y, 0)
  Ox <- pad(dat$O, 0)
  wx <- pad(dat$wt, 0)
  alpha <- start %||% ps_continue(fit$alpha, ca, cy, pord[2], my$offset)
  extra <- if (!is.null(extra_fun)) extra_fun(my, ny) else NULL
  phi <- if (!is.null(fit$phi)) matrix(fit$phi, nrow(Yx), ncol(Yx))
  f <- ps2_irls(Yx, Ox, wx, fit$ma, my, P, alpha, maxit = maxit, tol = tol, extra = extra, phi = phi)
  f$my <- my
  f$P <- P
  f$forecast_cols <- ny + seq_len(h)
  f
}

# Does the age axis start with the first year of life as its own group
# (Camarda 2019, section 3.1; at least four age groups, as in mt_smooth())?
ps_has_infant <- function(age) age$lower[1] == 0 && age$upper[1] == 1 && length(age) > 3

# Settings of a P-spline specification: one value per direction
ps_settings <- function(lambda, criterion, pord, deg, ndx, lambda_range) {
  two <- function(v, name, lower, integerish) {
    if (!is.numeric(v) || !length(v) %in% 1:2 || anyNA(v) || any(v < lower) ||
        (integerish && any(abs(v - round(v)) > 1e-8))) {
      mt_abort(sprintf("`%s` must be one or two %s of at least %s (age, year).", name,
                       if (integerish) "whole numbers" else "numbers", lower), class = "argument", call = NULL)
    }
    rep_len(v, 2)
  }
  if (!is.null(lambda)) lambda <- two(lambda, "lambda", 0, FALSE)
  if (!is.null(ndx)) ndx <- two(ndx, "ndx", 1, TRUE)
  if (!is.numeric(lambda_range) || length(lambda_range) != 2 || anyNA(lambda_range) || any(lambda_range <= 0) ||
      lambda_range[1] >= lambda_range[2]) {
    mt_abort("`lambda_range` must be two positive increasing numbers.", class = "argument", call = NULL)
  }
  list(lambda = lambda, criterion = criterion, pord = two(pord, "pord", 1, TRUE), deg = two(deg, "deg", 1, TRUE),
       ndx = ndx, lambda_range = lambda_range)
}

# P-spline fit of one stratum (common to p_splines() and cp_splines())
ps_fit_stratum <- function(p, D, E, ages, periods, margin, overdispersion) {
  xa <- cell_ages(ages$lower, ages$upper)
  infant <- isTRUE(p$infant) && ps_has_infant(ages)
  ndx <- p$ndx %||% c(ps_ndx(if (infant) xa[-1] else xa), ps_ndx(periods))
  f <- ps2_model(D, E, xa, periods, ndx, p$deg, p$pord, margin, infant, p$lambda, p$criterion, p$lambda_range,
                 overdispersion)
  if (isTRUE(f$at_bound)) {
    mt_warn(sprintf(paste("A smoothing parameter of the P-spline fit is at the boundary of `lambda_range`",
                          "(lambda = %s); widen the range or fix `lambda`."),
                    paste(signif(f$lambda, 3), collapse = ", ")), class = "lambda_bound", call = NULL)
  }
  list(alpha = f$alpha, lambda = f$lambda, psi2 = f$psi2, dispersion = f$psi2_estimate, ed = f$ed, dev = f$dev,
       ic = f$ic, criterion = p$criterion, data = f$data, ma = f$ma, my = f$my, xa = xa, xy = periods,
       pord = p$pord, ndx = ndx, n_obs = f$n_obs, mu = f$mu, fitted = exp(f$eta), at_bound = f$at_bound,
       converged = f$converged, phi = f$phi, dispersion_raw = f$dispersion_raw,
       lambda_dispersion = f$lambda_dispersion, dispersion_tau = f$dispersion_tau)
}

# Draws of the forecast log rates [age, h, draw] from the normal
# approximation of the coefficients of the extended fit `fx`: mean
# fx$alpha, covariance psi2 (B'WB + P)^-1 (P with psi2)
ps_draws <- function(fx, Ba, psi2, n_draws) {
  ca <- ncol(Ba)
  cy <- ncol(fx$my$B)
  h <- length(fx$forecast_cols)
  Byf <- fx$my$B[fx$forecast_cols, , drop = FALSE]
  p <- ca * cy
  Z <- matrix(stats::rnorm(p * n_draws), p)
  U <- if (!is.null(fx$R)) {
    backsolve(fx$R, Z)
  } else {
    e <- eigen(fx$BWB + fx$P, symmetric = TRUE)
    e$vectors %*% (Z / sqrt(pmax(e$values, max(e$values) * 1e-12)))
  }
  U <- U * sqrt(psi2) + fx$alpha
  # Ba A_k Byf' for every draw k, as two matrix products
  T1 <- Byf %*% matrix(aperm(array(U, c(ca, cy, n_draws)), c(2, 1, 3)), cy)
  T2 <- matrix(aperm(array(T1, c(h, ca, n_draws)), c(2, 1, 3)), ca)
  array(Ba %*% T2, c(nrow(Ba), h, n_draws))
}

#' Two-dimensional P-spline model
#'
#' `p_splines()` smooths the log death rates over age and year with a
#' Poisson P-spline surface and forecasts them by extending the surface over
#' the forecast years (Currie, Durban and Eilers 2004). The code follows the
#' reference implementation, `Mort2Dsmooth()` and `predict.Mort2Dsmooth()`
#' of the R package MortalitySmooth 2.3.4 (Camarda), and reproduces its
#' fitted values, forecasts and standard errors for given smoothing
#' parameters.
#'
#' **Model.** The deaths of age `x` and year `t` are Poisson with mean
#' `E[x, t] m[x, t]`. The log rates are a tensor-product B-spline surface,
#' `log m = Ba A t(By)`, with B-splines of degree `deg` on equally spaced
#' knots over age (`Ba`) and over years (`By`) and a matrix of coefficients
#' `A`. Differences of order `pord` between neighbouring coefficients are
#' penalised, over age with weight `lambda[1]` and over years with weight
#' `lambda[2]`. The coefficients maximise the penalised Poisson
#' log-likelihood (penalised scoring algorithm, computed in array form).
#' The smoothing parameters minimise the BIC (default) or the AIC, with the
#' effective dimension (the trace of the hat matrix) as the number of
#' parameters.
#'
#' **Forecast.** The forecast is part of the fit. The basis over years is
#' extended by B-splines on the same knot spacing beyond the last forecast
#' year (and below the first year, as `predict.Mort2Dsmooth()` places the
#' knots), and the forecast years enter with weight 0: they carry no data,
#' and their coefficients follow from the penalty. With `pord = 2` the
#' coefficients continue about linearly, so each age continues the trend of
#' its last years; the order of the penalty determines the form of the
#' forecast. The model is fitted again on the extended basis with the
#' smoothing parameters of the fit, as `predict.Mort2Dsmooth()` does. The
#' penalty over age also acts on the added coefficients, so the fitted rates
#' of the data years change slightly, and the forecast depends a little on
#' the horizon through the last knot (up to 0.02 on the log scale for
#' England and Wales men, ages 30-100, 10 years ahead).
#'
#' **Uncertainty.** Draws come from the normal approximation of the
#' coefficients, with the covariance `(B'WB + P)^-1` of the penalised
#' likelihood (`W` the fitted deaths, `P` the penalty). MortalitySmooth
#' computes the standard errors of forecast log rates from this covariance.
#' Over the forecast years the penalty acts as a prior: the variance grows
#' with the horizon, the faster the smaller `lambda[2]`. Draws are paths of
#' the whole surface, correlated over ages and years as in the model. They
#' describe the expected rates; predictive draws
#' (`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths.
#' With `overdispersion = TRUE`, the dispersion `psi2` (deviance divided by
#' the number of cells minus the effective dimension) enters the choice of
#' the smoothing parameters and the penalty, as in
#' `Mort2Dsmooth(overdispersion = TRUE)`, and the covariance is
#' `psi2 (B'WB + psi2 P)^-1`. MortalitySmooth's `predict()` instead refits
#' the forecast without `psi2` in the penalty and multiplies the standard
#' errors by `sqrt(psi2)`. For England and Wales men (ages 30-100,
#' 1961-2011), the two forecasts differ by up to 0.08 on the log scale ten
#' years ahead.
#'
#' **Dispersion by age.** `overdispersion = "age"` replaces the Poisson
#' assumption by the first two moments of the deaths, `E[D] = E m` and
#' `var(D) = phi[x] E[D]`, with a dispersion `phi[x]` for each age that
#' varies smoothly over age (Djeundje and Currie 2011, eq. 2.1; Model 3,
#' eq. 5.26). Rates and dispersions are estimated by the bias corrected
#' scheme of the paper (section 4):
#'
#' 1. Start with `phi = 1` at every age.
#' 2. Fit the rates by penalised quasi-likelihood: the scoring algorithm
#'    with the weights of the cells divided by their dispersion (eq. 3.13).
#'    The smoothing parameters minimise the criterion of the scaled
#'    deviance: the deviance of each cell divided by its dispersion, plus
#'    `log(n)` (BIC, eq. 3.16) or 2 (AIC) times the effective dimension.
#' 3. For each age, compute the bias-corrected mean Pearson statistic
#'    `d[x]`: the sum of `(D - mu)^2 / mu` over its years, divided by its
#'    number of cells `n[x]` minus its contribution to the effective
#'    dimension, the sum of its diagonal elements of the hat matrix (eq.
#'    4.20 and 4.21).
#' 4. Smooth the statistics over age by penalised quasi-likelihood with
#'    `E[d] = phi`, `var(d) = tau phi^2 / n[x]` and the link `-1 / phi =
#'    Ba beta` (eq. 4.22): the B-splines over age of the rates and a
#'    difference penalty of order `pord[1]` (eq. 5.28). The smoothing
#'    parameter minimises the scaled BIC of this fit.
#' 5. Repeat steps 2 to 4 until no dispersion changes by more than 1e-5
#'    relative to its value (section 5.2.2), at most 50 times.
#'
#' The forecast fit keeps the dispersions of the ages in the weights, and
#' the draws come from the covariance `(B'WB + P)^-1` with `W` the fitted
#' deaths divided by their dispersion. In the paper's example (CMI assured
#' lives, ages 40-90, section 5.2.3) the dispersion down-weights the fit to
#' the data and the selected surface is smoother: the effective dimension
#' is 46 under the Poisson assumption, 36 with one dispersion and 33 with
#' the dispersion by age, and the confidence intervals of the forecast are
#' narrower. With a single age, `"age"` gives one dispersion, the
#' bias-corrected mean Pearson statistic (eq. 4.20). The fit selects the
#' smoothing parameters in every round, so it takes several times as long
#' as the Poisson fit.
#'
#' **Data.** Ages are placed at the midpoints of their cells, an open
#' interval at its lower bound plus half the width of the interval below it
#' (as in [mt_smooth()]). Cells with missing deaths or exposures, and cells
#' without exposure, get weight 0, so the surface interpolates them.
#' `infant = TRUE` gives the first year of life a row of coefficients of
#' its own, penalised over years but not connected over age with the
#' following ages (Camarda 2019, section 3.1); it applies when the grid
#' starts with age group 0-1 and has at least four age groups.
#'
#' **Package-specific choices (reference unknown).** By default the knots
#' are about five years of age and five calendar years apart, `floor(n /
#' 5)` intervals for `n` ages or years as in MortalitySmooth; for age groups
#' wider than one year, `n` is the number of single years that the cell
#' ages span. The smoothing parameters are searched on a grid with steps of
#' 1 on the log10 scale, followed by a continuous search near the best grid
#' point (MortalitySmooth: two grid searches with steps of 2 and 0.5).
#' With `overdispersion = "age"`: `tau` is the mean Pearson statistic of the
#' smooth of step 4, `sum(n[x] (d[x] - phi[x])^2 / phi[x]^2)` divided by the
#' number of ages minus its effective dimension, alternated with the choice
#' of its smoothing parameter as `psi2` with `overdispersion = TRUE` (at
#' most 5 times; the paper estimates `tau` without stating how); that
#' smoothing parameter is searched over `lambda_range` as those of the
#' rates; ages without data get weight 0 in step 4, so the smooth gives
#' their dispersion; predictive draws add Poisson variation of the deaths,
#' as with `overdispersion = TRUE`, because the quasi-likelihood defines
#' the mean and variance of the deaths but no distribution, and the paper
#' gives no predictive distribution. The paper's full extended
#' quasi-likelihood scheme (mean deviances by age, eq. 2.6) is not offered:
#' it underestimates the dispersion (section 4, Figure 4), and the paper
#' reports its results with the bias corrected scheme.
#'
#' @param lambda Smoothing parameters over age and years, `c(age, year)`;
#'   `NULL` (default) selects them by `criterion`.
#' @param criterion `"bic"` (default) or `"aic"`.
#' @param pord Order of the difference penalties, over age and years (one
#'   value for both).
#' @param deg Degree of the B-splines, over age and years.
#' @param ndx Numbers of intervals between knots over age and years; `NULL`
#'   gives about one knot per five years (see Details).
#' @param overdispersion `FALSE` (default): Poisson deaths. `TRUE`: one
#'   dispersion for all cells in the choice of the smoothing parameters and
#'   in the draws (MortalitySmooth). `"age"`: a dispersion for each age that
#'   varies smoothly over age (Djeundje and Currie 2011; see Details).
#' @param lambda_range Search range of the smoothing parameters (the default
#'   of MortalitySmooth). A selected value at an end of the range gives a
#'   warning (class `mt_warning_lambda_bound`) when the criterion improves
#'   by more than 1 two decades beyond it.
#' @param infant Give the first year of life its own coefficients (Camarda
#'   2019).
#' @return A model specification.
#' @references
#' Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
#' mortality rates. *Statistical Modelling* 4(4): 279-298.
#' \doi{10.1191/1471082X04st080oa}
#'
#' Camarda CG (2015). MortalitySmooth: Smoothing and forecasting Poisson
#' counts with P-splines. R package version 2.3.4 (CRAN archive).
#'
#' Camarda CG (2019). Smooth constrained mortality forecasting.
#' *Demographic Research* 41(38): 1091-1130.
#' \doi{10.4054/DemRes.2019.41.38}
#'
#' Djeundje VAB, Currie ID (2011). Smoothing dispersed counts with
#' applications to mortality data. *Annals of Actuarial Science* 5(1):
#' 33-52. \doi{10.1017/S1748499510000047}
#' @seealso [cp_splines()] for forecasts constrained to past age patterns and
#'   time trends; [mt_smooth()] for smoothing without forecasting.
#' @examples
#' x <- mt_example(deaths = "single", years = 1995:2012)
#' fit <- mt_fit(p_splines(), x)
#' fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
#' mt_parameters(fit)
#' # dispersion that varies smoothly over age (Djeundje and Currie 2011),
#' # smoothing parameters of the rates given
#' fit_age <- mt_fit(p_splines(lambda = c(1, 1000), overdispersion = "age"), x)
#' p <- mt_parameters(fit_age)
#' p[parameter == "dispersion" & stratum == "female"]
#' @export
p_splines <- function(lambda = NULL, criterion = c("bic", "aic"), pord = 2, deg = 3, ndx = NULL,
                      overdispersion = FALSE, lambda_range = c(1e-4, 1e6), infant = FALSE) {
  criterion <- match.arg(criterion)
  if (!identical(overdispersion, "age") && !(is.logical(overdispersion) && length(overdispersion) == 1L &&
                                              !is.na(overdispersion))) {
    mt_abort("`overdispersion` must be TRUE, FALSE or \"age\".", class = "argument")
  }
  check_flag(infant, "infant")
  p <- ps_settings(lambda, criterion, pord, deg, ndx, lambda_range)
  p$overdispersion <- overdispersion
  p$infant <- infant
  od <- if (identical(overdispersion, "age")) ", dispersion by age" else if (overdispersion) ", overdispersed"
  spec <- new_mt_spec("p_splines", p, min_periods = 5,
                      label = paste0("P-splines", od, if (infant) ", infant coefficient" else ""))
  spec$missing_ok <- TRUE
  spec
}

#' @export
fit_model.mt_spec_p_splines <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  ps_fit_stratum(p, D, E, ages, periods, margin = 0.01, overdispersion = p$overdispersion)
}

#' @export
simulate_model.mt_spec_p_splines <- function(spec, fit, h, n_draws, ...) {
  fx <- ps2_extend_fit(fit, h, fit$pord)
  if (!fx$converged) mt_warn("The P-spline forecast fit did not converge.", class = "model_fit", call = NULL)
  exp(ps_draws(fx, fit$ma$B, fit$psi2, n_draws))
}

# CP-splines (Camarda 2019) -------------------------------------------------------

# Rates of ageing and of change [age, year] of a fitted surface (Camarda
# 2019, eq. 14-17): derivatives of the log rates over age and over years,
# from the derivatives of the B-splines
cps_derivatives <- function(alpha, ma, my) {
  A <- matrix(alpha, ncol(ma$B))
  list(ageing = ma$C %*% A %*% t(my$B), change = ma$B %*% A %*% t(my$C))
}

# Bounds per age (Camarda 2019, section 3.2): the central ranges with
# `level_ageing` and `level_change` of the rates of ageing and of change
# over the fitting years (empirical quantiles, R's default type 7)
cps_bounds <- function(alpha, ma, my, level_ageing, level_change) {
  d <- cps_derivatives(alpha, ma, my)
  q <- function(M, level) {
    pr <- c((1 - level) / 2, (1 + level) / 2)
    t(apply(M, 1, stats::quantile, probs = pr, names = FALSE))
  }
  a <- q(d$ageing, level_ageing)
  ch <- q(d$change, level_change)
  list(ageing_lower = a[, 1], ageing_upper = a[, 2], change_lower = ch[, 1], change_upper = ch[, 2])
}

# Asymmetric penalties of the forecast years (Camarda 2019, eq. 18-20):
# where the rate of ageing or of change of an age in a forecast year lies
# outside its bounds, kappa times the squared distance to the bound is
# added to the penalised deviance. Returns a function of the extended year
# basis that returns the function of the coefficients used by ps2_irls().
cps_extra <- function(ma, bounds, kappa) {
  kappa <- rep_len(kappa, 2)
  pCa <- rt_pairs(ma$C)
  na <- nrow(ma$B)
  ca <- ncol(ma$B)
  function(my, ny) {
    fy <- seq.int(ny + 1L, nrow(my$B))
    h <- length(fy)
    Byf <- my$B[fy, , drop = FALSE]
    Cyf <- my$C[fy, , drop = FALSE]
    pBy <- rt_pairs(Byf)
    pCy <- rt_pairs(Cyf)
    map_a <- glam_map(pCa, pBy)
    map_t <- glam_map(ma$pr, pCy)
    lo_a <- matrix(bounds$ageing_lower, na, h)
    up_a <- matrix(bounds$ageing_upper, na, h)
    lo_t <- matrix(bounds$change_lower, na, h)
    up_t <- matrix(bounds$change_upper, na, h)
    function(alpha) {
      Am <- matrix(alpha, ca)
      da <- ma$C %*% Am %*% t(Byf)
      dt <- ma$B %*% Am %*% t(Cyf)
      va <- (da > up_a) + (da < lo_a)
      vt <- (dt > up_t) + (dt < lo_t)
      ga <- (da > up_a) * up_a + (da < lo_a) * lo_a
      gt <- (dt > up_t) * up_t + (dt < lo_t) * lo_t
      ea <- glam_bwb(pCa, pBy, va, map_a)
      et <- glam_bwb(ma$pr, pCy, vt, map_t)
      list(idx = list(ea$idx, et$idx), val = list(kappa[1] * ea$val, kappa[2] * et$val),
           p = kappa[1] * as.vector(crossprod(ma$C, ga) %*% Byf) + kappa[2] * as.vector(crossprod(ma$B, gt) %*% Cyf))
    }
  }
}

# Deaths of the residual bootstrap (Camarda 2019, eq. 21-22): for observed
# deaths y and a resampled deviance residual r, the value mu at which y has
# the deviance residual r, the root of g(u) = exp(u) - y u - c with u =
# log(mu) and c = r^2 / 2 + y - y log(y), below log(y) for r > 0 and above
# for r < 0. g is convex, so Newton steps from a point of the root's side
# where g > 0 (y (e^-k + k - 1) > r^2 / 2 at u = log(y) - k with k = 2 + r^2
# / (2 y); y (e^k - k - 1) > r^2 / 2 at u = log(y) + k with k = |r| /
# sqrt(y) + 1) converge to it monotonically. Without deaths (y = 0) the
# equation gives r^2 / 2 for r <= 0, and its limit 0 for r > 0.
boot_deaths <- function(y, r) {
  out <- y
  z <- y == 0
  out[z] <- (r[z] <= 0) * r[z]^2 / 2
  k <- which(!z & r != 0)
  if (!length(k)) return(out)
  yy <- y[k]
  rr <- r[k]
  cc <- rr^2 / 2 + yy - yy * log(yy)
  pos <- rr > 0
  u <- log(yy) + ifelse(pos, -(2 + rr^2 / (2 * yy)), abs(rr) / sqrt(yy) + 1)
  for (it in 1:100) {
    step <- (exp(u) - yy * u - cc) / (exp(u) - yy)
    u <- u - step
    if (max(abs(step)) < 1e-12) break
  }
  out[k] <- exp(u)
  out
}

# Residual bootstrap of the fit (Camarda 2019, section 3.2.2): deviance
# residuals of the observed cells resampled with replacement, bootstrap
# deaths by boot_deaths(), the surface refitted with the smoothing
# parameters of the fit, and new bounds. Returns the bootstrap deaths
# [age, year, replicate] and the bounds of each replicate.
cps_bootstrap <- function(fit, n_boot, level_ageing, level_change) {
  dat <- fit$data
  obs <- dat$wt > 0
  y <- dat$Y[obs]
  mu <- fit$mu[obs]
  res <- sign(y - mu) * sqrt(2 * pmax(ifelse(y > 0, y * log(y / mu), 0) - (y - mu), 0))
  ca <- ncol(fit$ma$B)
  cy <- ncol(fit$my$B)
  pen <- ps_penalties(ca, cy, fit$pord, isTRUE(fit$ma$infant))
  P <- fit$psi2 * (fit$lambda[1] * pen$Pa + fit$lambda[2] * pen$Py)
  Ys <- array(0, c(dim(dat$Y), n_boot))
  bounds <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    Yb <- dat$Y
    Yb[obs] <- boot_deaths(y, res[sample.int(length(res), length(res), replace = TRUE)])
    fb <- ps2_irls(Yb, dat$O, dat$wt, fit$ma, fit$my, P, fit$alpha)
    bounds[[b]] <- cps_bounds(fb$alpha, fit$ma, fit$my, level_ageing, level_change)
    Ys[, , b] <- Yb
  }
  list(Y = Ys, bounds = bounds)
}

# Constrained forecast fit (Camarda 2019, eq. 19) for deaths Y (default:
# the data) with the given bounds
cps_forecast_fit <- function(fit, h, bounds, Y = NULL, start = NULL) {
  ps2_extend_fit(fit, h, fit$pord, extra_fun = cps_extra(fit$ma, bounds, fit$kappa), tol = 1e-4, maxit = 100,
                 Y = Y, start = start)
}

#' Smooth constrained mortality forecasting (CP-splines)
#'
#' `cp_splines()` is the CP-spline model of Camarda (2019): the
#' two-dimensional P-spline model of [p_splines()], with coefficients of
#' its own for the first year of life, whose forecast is constrained to age
#' patterns and time trends like those of the fitting years. The code is
#' written from the paper; the author's R code (supplementary material) was
#' used to check details and results.
#'
#' The model has five steps (Camarda 2019, sections 3.1 to 3.2.2):
#'
#' 1. A two-dimensional P-spline model is fitted to the data years
#'    ([p_splines()], smoothing parameters by BIC). The first year of life
#'    gets a row of coefficients of its own, penalised over years but not
#'    connected over age with the following ages (section 3.1).
#' 2. At every age and data year, two derivatives of the fitted log rates
#'    are computed from the derivatives of the B-splines: over age, the rate
#'    of ageing, and over years, the rate of change (eq. 14-17).
#' 3. Bounds per age: the central range with `level_ageing` (default 95%)
#'    of the rates of ageing over the data years, and the central range with
#'    `level_change` (default 50%, the interquartile range) of the rates of
#'    change (empirical quantiles).
#' 4. The forecast fit (the basis over years extended over the forecast
#'    years with weight 0, as in [p_splines()]) adds asymmetric penalties:
#'    where a derivative in a forecast year lies outside its bounds, its
#'    squared distance to the bound is penalised with weight `kappa` (eq.
#'    18-20). The cells outside their bounds are found again at every step
#'    of the scoring algorithm. Forecast age patterns thus keep the shape of
#'    the fitted ones, and the change of each age over time stays within the
#'    range of its past change.
#' 5. Uncertainty comes from a residual bootstrap (section 3.2.2). The
#'    deviance residuals of the fit of step 1 are resampled with
#'    replacement. With the observed deaths, each resampled residual gives
#'    bootstrap deaths: the value at which the observed deaths have that
#'    deviance residual (eq. 22). Steps 2 to 4 are repeated on every
#'    bootstrap sample, with the smoothing parameters, levels and `kappa`
#'    of the fit.
#'
#' Each draw is the forecast of one bootstrap sample; the `n_boot` samples
#' are used in turn when `n_draws` is larger. The draws thus carry the
#' uncertainty of the fit only: the forecast has no innovations, and
#' Camarda (2019) notes that the intervals are often narrower than those of
#' other methods. Draws describe the expected rates; predictive draws
#' (`mt_forecast(predictive = TRUE)`) add Poisson variation of the deaths,
#' the observation model of the fit (the paper gives no method for
#' prediction intervals of observed rates). With `n_boot = 0` every draw is
#' the forecast of the data.
#'
#' The grid needs the first year of life as its own age group for step 1
#' to give it coefficients of its own (`infant = TRUE` and at least four age
#' groups); otherwise all ages share the surface. Ages and missing cells
#' are treated as in [p_splines()]. The bounds are ranges over the data
#' years, so the model needs at least 10 of them.
#'
#' **Package-specific choices (reference unknown).** As in [p_splines()]:
#' cell ages, the knot spacing for age groups and the search of the
#' smoothing parameters (the paper minimises the BIC; the author's code
#' searches a grid over ranges set for its data). The knots lie on the range
#' of the data without margins, as in the author's code. The basis over
#' years of the forecast is that of the fit extended on the same knot
#' spacing, as the paper describes, with the knots placed as in
#' [p_splines()]; the author's code builds a new basis over the data and
#' forecast years with one knot per five years (for French women, ages
#' 0-100, 1950-2001, this changes the log rates forecast to 2021 by at most
#' 0.001). The rate
#' of ageing at age 0 is the difference of the log rates of the first two
#' age groups divided by the distance between their cell ages (the author's
#' code: single years of age). Cells without deaths have the deviance
#' residual `-sqrt(2 mu)` (with 0 log 0 = 0), and their bootstrap deaths
#' are `r^2 / 2` for a negative resampled residual `r` and 0 for a positive
#' one, the limit of eq. 22; the author's code gives these cells a residual
#' drawn from the other cells. The scoring algorithm with asymmetric
#' penalties stops when no coefficient changes by more than 1e-4 (relative
#' to its size where that exceeds 1), after at most 100 steps, as in the
#' author's code (relative changes).
#'
#' @inheritParams p_splines
#' @param level_ageing,level_change Levels of the central ranges of the
#'   rates of ageing and of change over the data years that bound the
#'   forecast (95% and 50% in Camarda 2019).
#' @param kappa Weight of the asymmetric penalties, one value or one per
#'   constraint (ageing, change); `1e4` in Camarda (2019).
#' @param n_boot Number of bootstrap samples (1000 in Camarda 2019). Each
#'   needs a forecast fit, so forecasts with many draws take time.
#' @param infant Give the first year of life its own coefficients (Camarda
#'   2019, section 3.1).
#' @return A model specification.
#' @references
#' Camarda CG (2019). Smooth constrained mortality forecasting.
#' *Demographic Research* 41(38): 1091-1130.
#' \doi{10.4054/DemRes.2019.41.38}
#'
#' Currie ID, Durban M, Eilers PHC (2004). Smoothing and forecasting
#' mortality rates. *Statistical Modelling* 4(4): 279-298.
#' \doi{10.1191/1471082X04st080oa}
#' @seealso [p_splines()] for the unconstrained forecast.
#' @examples
#' \donttest{
#' x <- mt_example(deaths = "single", years = 1995:2012)
#' fit <- mt_fit(cp_splines(n_boot = 20), x)
#' fc <- mt_forecast(fit, h = 10, n_draws = 20, seed = 1)
#' mt_parameters(fit)
#' }
#' @export
cp_splines <- function(level_ageing = 0.95, level_change = 0.5, kappa = 1e4, n_boot = 100, lambda = NULL,
                       criterion = c("bic", "aic"), pord = 2, deg = 3, ndx = NULL, lambda_range = c(1e-4, 1e6),
                       infant = TRUE) {
  criterion <- match.arg(criterion)
  check_number(level_ageing, "level_ageing", lower = 0, upper = 1)
  check_number(level_change, "level_change", lower = 0, upper = 1)
  if (!is.numeric(kappa) || !length(kappa) %in% 1:2 || anyNA(kappa) || any(kappa < 0)) {
    mt_abort("`kappa` must be one or two non-negative numbers (ageing, change).", class = "argument")
  }
  check_number(n_boot, "n_boot", lower = 0, integerish = TRUE)
  check_flag(infant, "infant")
  p <- ps_settings(lambda, criterion, pord, deg, ndx, lambda_range)
  p <- c(p, list(level_ageing = level_ageing, level_change = level_change, kappa = rep_len(kappa, 2),
                 n_boot = as.integer(n_boot), infant = infant))
  default <- level_ageing == 0.95 && level_change == 0.5
  spec <- new_mt_spec("cp_splines", p, min_periods = 10,
                      label = paste0("CP-splines",
                                     if (!default) sprintf(" (levels %g%%, %g%%)", 100 * level_ageing,
                                                           100 * level_change) else ""))
  spec$missing_ok <- TRUE
  spec
}

#' @export
fit_model.mt_spec_cp_splines <- function(spec, D, E, ages, periods, ...) {
  p <- spec$params
  fit <- ps_fit_stratum(p, D, E, ages, periods, margin = 0, overdispersion = FALSE)
  fit$kappa <- p$kappa
  fit$bounds <- cps_bounds(fit$alpha, fit$ma, fit$my, p$level_ageing, p$level_change)
  if (p$n_boot > 0) fit$boot <- cps_bootstrap(fit, p$n_boot, p$level_ageing, p$level_change)
  fit
}

#' @export
simulate_model.mt_spec_cp_splines <- function(spec, fit, h, n_draws, ...) {
  main <- cps_forecast_fit(fit, h, fit$bounds)
  fc <- function(f) exp(f$eta[, f$forecast_cols, drop = FALSE])
  out <- array(NA_real_, c(nrow(fit$data$Y), h, n_draws))
  if (is.null(fit$boot)) {
    if (!main$converged) mt_warn("The CP-spline forecast fit did not converge.", class = "model_fit", call = NULL)
    out[] <- fc(main)
    return(out)
  }
  nb <- min(dim(fit$boot$Y)[3], n_draws)
  conv <- logical(nb)
  sims <- lapply(seq_len(nb), function(b) {
    f <- cps_forecast_fit(fit, h, fit$boot$bounds[[b]], Y = fit$boot$Y[, , b], start = main$alpha)
    conv[b] <<- f$converged
    fc(f)
  })
  if (!all(conv)) {
    mt_warn(sprintf("The CP-spline forecast fit did not converge for %d of %d bootstrap samples.", sum(!conv), nb),
            class = "model_fit", call = NULL)
  }
  for (j in seq_len(n_draws)) out[, , j] <- sims[[(j - 1L) %% nb + 1L]]
  out
}
