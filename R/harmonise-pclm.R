# Penalized composite link model (PCLM) ----------------------------------------
#
# Method: Eilers (2007), Rizzi, Gampe and Eilers (2015); for rates over age and
# period, Rizzi et al. (2019). The structure (B-spline basis on the latent
# scale, difference penalty, IRLS, lambda chosen by an information criterion)
# follows the 'ungroup' package 1.6.3 (MIT licence; M. D. Pascariu, S. Rizzi,
# J. Schoeley, M. J. Danko). This implementation is written for mortensor and
# differs in the treatment of rates: the exposure offset enters the composite
# link directly, and an open age interval of the grid is expanded into single
# years up to omega, its exposure split by the PCLM of the exposures (Rizzi,
# Gampe and Eilers 2015; open_weights_pclm(); default) or, as an option, with
# stationary-population weights (package-specific, reference unknown; see
# open_weights()), so no artificial-bin reallocation is applied to rates
# (the source of the level shift found in revamp plan 5.1).
#
# Over age and period, the basis is the tensor product of the bases over
# fine ages and periods, and the products with it are computed in array
# form, as in generalised linear array models (Currie, Durban and Eilers
# 2006), without forming the tensor basis (pclm_tensor(), pclm2_irls()).

#' Fit a penalized composite link model
#'
#' Latent values `gamma = exp(B beta)` on a fine grid are linked to observed
#' counts by `mu = C (o * gamma)`, where `o` is an offset (exposures for
#' rates, interval widths for counts). The coefficients are estimated by
#' penalized Poisson IRLS; the smoothing parameter is chosen by BIC or AIC on a
#' log scale unless fixed.
#'
#' @param y Observed counts.
#' @param C Composite link matrix (observations x fine cells), dense or sparse.
#' @param o Offset per fine cell.
#' @param B Basis matrix (fine cells x coefficients), or a tensor basis
#'   over fine ages and periods from `pclm_tensor()`, whose products are
#'   computed in array form (Currie, Durban and Eilers 2006).
#' @param pen List of penalty matrices (one per smoothing parameter).
#' @param lambda Fixed smoothing parameters, or `NULL` to select them.
#' @param criterion `"bic"` or `"aic"`.
#' @param lambda_range Search range for each smoothing parameter.
#' @param beta Optional starting values.
#' @return A list with `gamma`, `beta`, `mu` (fitted observed counts),
#'   `lambda`, `ed` (effective dimension), `dev`, `ic`, `se_eta`,
#'   `converged` and `at_bound`.
#' @references
#' Currie ID, Durban M, Eilers PHC (2006). Generalized linear array models
#' with applications to multidimensional smoothing. *Journal of the Royal
#' Statistical Society Series B* 68(2): 259-280.
#' \doi{10.1111/j.1467-9868.2006.00543.x}
#' @keywords internal
pclm_fit <- function(y, C, o, B, pen, lambda = NULL, criterion = c("bic", "aic"),
                     lambda_range = c(1e-2, 1e8), beta = NULL) {
  criterion <- match.arg(criterion)
  np <- length(pen)
  if (is.null(beta)) {
    tot <- sum(as.numeric(C %*% o))
    beta <- solve_lsq(B, rep(log(max(sum(y), 0.5) / tot), basis_nrow(B)))
  }
  ic_of <- function(fit) fit$dev + (if (criterion == "bic") log(length(y)) else 2) * fit$ed
  if (!is.null(lambda)) {
    fit <- pclm_irls(y, C, o, B, pen, rep_len(lambda, np), beta, se = TRUE)
    fit$ic <- ic_of(fit)
    fit$at_bound <- FALSE
    return(fit)
  }
  lr <- log10(lambda_range)
  state <- new.env()
  state$beta <- beta
  obj <- function(ll) {
    f <- pclm_irls(y, C, o, B, pen, 10^ll, state$beta)
    if (f$converged) state$beta <- f$beta
    v <- ic_of(f)
    if (is.finite(v)) v else 1e300
  }
  ll <- rep(mean(lr), np)
  if (np == 1L) {
    ll <- stats::optimize(obj, lr, tol = 0.05)$minimum
  } else {
    for (round in 1:2) {
      for (k in seq_len(np)) {
        obj_k <- function(z) {
          l2 <- ll
          l2[k] <- z
          obj(l2)
        }
        ll[k] <- stats::optimize(obj_k, lr, tol = 0.1)$minimum
      }
    }
  }
  fit <- pclm_irls(y, C, o, B, pen, 10^ll, state$beta, se = TRUE)
  fit$ic <- ic_of(fit)
  fit$at_bound <- any(ll < lr[1] + 0.1 | ll > lr[2] - 0.1)
  fit
}

pclm_irls <- function(y, C, o, B, pen, lambda, beta, maxit = 100L, tol = 1e-6, se = FALSE) {
  if (inherits(B, "pclm_tensor")) return(pclm2_irls(y, C, o, B, pen, lambda, beta, maxit, tol, se))
  P <- Reduce(`+`, Map(`*`, lambda, pen))
  pdev <- function(b) {
    eta <- as.numeric(B %*% b)
    if (any(!is.finite(eta)) || max(eta) > 50) return(Inf)
    mu <- pmax(as.numeric(C %*% (o * exp(eta))), 1e-12)
    poisson_dev(y, mu) + sum(b * as.numeric(P %*% b))
  }
  obj_old <- pdev(beta)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    eta <- as.numeric(B %*% beta)
    gam <- exp(eta)
    mu <- pmax(as.numeric(C %*% (o * gam)), 1e-12)
    X <- as.matrix(C %*% (B * (o * gam)))
    XtWX <- crossprod(X, X / mu)
    A <- XtWX + P
    rhs <- crossprod(X, (y - mu) / mu) + XtWX %*% beta
    beta_new <- tryCatch({
      R <- chol(A)
      as.numeric(backsolve(R, forwardsolve(t(R), rhs)))
    }, error = function(e) tryCatch(as.numeric(solve(A, rhs)), error = function(e2) NULL))
    if (is.null(beta_new) || anyNA(beta_new)) break
    step <- beta_new - beta
    obj_new <- pdev(beta + step)
    halvings <- 0L
    while (!(obj_new <= obj_old + 1e-10 * abs(obj_old)) && halvings < 25L) {
      step <- step / 2
      obj_new <- pdev(beta + step)
      halvings <- halvings + 1L
    }
    if (!is.finite(obj_new)) break
    beta <- beta + step
    if (max(abs(step)) < tol * (1 + max(abs(beta))) || abs(obj_old - obj_new) < 1e-8 * (1 + abs(obj_new))) {
      converged <- TRUE
      obj_old <- obj_new
      break
    }
    obj_old <- obj_new
  }
  eta <- as.numeric(B %*% beta)
  gam <- exp(eta)
  mu <- pmax(as.numeric(C %*% (o * gam)), 1e-12)
  X <- as.matrix(C %*% (B * (o * gam)))
  XtWX <- crossprod(X, X / mu)
  A <- XtWX + P
  Ainv <- tryCatch(chol2inv(chol(A)), error = function(e) MASS_ginv(A))
  ed <- sum(Ainv * XtWX)
  se_eta <- NULL
  if (se) {
    BV <- as.matrix(B %*% Ainv)
    se_eta <- sqrt(pmax(rowSums(BV * as.matrix(B)), 0))
  }
  list(gamma = gam, beta = beta, mu = mu, lambda = lambda, ed = ed, dev = poisson_dev(y, mu),
       se_eta = se_eta, converged = converged, iterations = it)
}

poisson_dev <- function(y, mu) {
  r <- ifelse(y > 0, y * log(y / mu), 0)
  2 * sum(r - (y - mu))
}

# Minimum-norm least squares (starting values); the basis can have more
# columns than there are cells, for example when smoothing 5-year age groups
solve_lsq <- function(B, z) {
  if (inherits(B, "pclm_tensor")) return(pclm2_lsq(B, z))
  B <- as.matrix(B)
  as.numeric(MASS_ginv(crossprod(B)) %*% crossprod(B, z))
}

MASS_ginv <- function(A) {
  s <- svd(A)
  pos <- s$d > max(s$d) * 1e-10
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

# PCLM over age and period in array form ---------------------------------------
#
# The latent log rates of the fine cells [fine age, period] are eta = Ba
# Theta By' for the bases Ba over fine ages and By over periods, the rows of
# the tensor basis kronecker(By, Ba) (coefficients with age fastest). As in
# generalised linear array models (Currie, Durban and Eilers 2006), the
# tensor basis is never formed. The composite link sums cells within
# observations; an observation within one period t has the row by_t (x) a
# of the working matrix X, with a = Ba' diag(o * gamma) c for its cells c
# in period t, so that X'WX is the sum over periods of (by_t by_t') (x)
# (A_t' W_t A_t), computed from row tensors (glam_map(), as glam_bwb()).
# Observations over several periods (deaths in periods coarser than the
# grid) are added as rows of X. X'WX and the penalties are banded in blocks
# of the period coefficients; the system is kept, solved and inverted within
# the band (band_chol(), band_inverse()). The model, the iterations and
# their stopping rules are those of pclm_irls(); only the arithmetic
# differs.

# Tensor basis kronecker(By, Ba) restricted to the fine cells `cells`
# (positions in the [fine age, period] matrix); `cache` keeps the array
# structure of the composite link between the fits of a search
pclm_tensor <- function(Ba, By, cells = seq_len(nrow(Ba) * nrow(By))) {
  structure(list(Ba = Ba, By = By, cells = cells, cache = new.env(parent = emptyenv())), class = "pclm_tensor")
}

# Rows of a basis (a matrix or a tensor basis) for the cells `keep`
# (logical), the number of rows, and the linear predictor B beta
basis_rows <- function(B, keep) {
  if (!inherits(B, "pclm_tensor")) return(B[keep, , drop = FALSE])
  pclm_tensor(B$Ba, B$By, B$cells[keep])
}

basis_nrow <- function(B) if (inherits(B, "pclm_tensor")) length(B$cells) else nrow(B)

basis_eta <- function(B, beta) {
  if (inherits(B, "pclm_tensor")) glam_eta(B$Ba, B$By, beta)[B$cells] else as.numeric(B %*% beta)
}

# solve_lsq() for a tensor basis: B'B and B'z in array form. With all cells,
# B'B = (By'By) (x) (Ba'Ba), whose singular values are the products of the
# eigenvalues of the two factors; when MASS_ginv() keeps all of them, the
# solution is the inverse of the factors applied to B'z
pclm2_lsq <- function(B, z) {
  Ba <- B$Ba
  By <- B$By
  Z <- matrix(0, nrow(Ba), nrow(By))
  Z[B$cells] <- z
  btz <- crossprod(Ba, Z) %*% By
  if (length(B$cells) == length(Z)) {
    Ga <- crossprod(Ba)
    Gy <- crossprod(By)
    ea <- eigen(Ga, symmetric = TRUE, only.values = TRUE)$values
    ey <- eigen(Gy, symmetric = TRUE, only.values = TRUE)$values
    if (min(ea) > 0 && min(ey) > 0 && min(ea) * min(ey) > max(ea) * max(ey) * 1e-10) {
      return(as.vector(solve(Ga, t(solve(Gy, t(btz))))))
    }
  }
  pa <- rt_pairs(Ba)
  py <- rt_pairs(By)
  map <- glam_map(pa, py)
  W <- matrix(0, nrow(Ba), nrow(By))
  W[B$cells] <- 1
  G <- matrix(0, map$p, map$p)
  e <- glam_bwb(pa, py, W, map)
  G[e$idx] <- e$val
  as.numeric(MASS_ginv(G) %*% as.vector(btz))
}

# Array structure of the composite link C [observation, cell] of a tensor
# basis: the pairs of an observation and a period it covers (`po`, `pt`),
# the matrix `Cp` [pair, cell] of the cells of each pair, the pairs of age
# coefficients that meet in an observation within one period (`pa`), the
# row tensors of By (`py`), the positions of their products in the band of
# X'WX (`bidx`, `bsrc`), the pairs of observations over several periods,
# the band (`band`, in blocks of the period coefficients) of X'WX and of
# the penalties, and the penalties within the band (`pen`)
pclm2_setup <- function(C, B, pen) {
  env <- B$cache
  if (!is.null(env$st) && identical(env$C, C) && identical(env$pen, pen)) return(env$st)
  Ba <- B$Ba
  By <- B$By
  nf <- nrow(Ba)
  nt <- nrow(By)
  ka <- ncol(Ba)
  ky <- ncol(By)
  ac <- (B$cells - 1L) %% nf + 1L
  tc <- (B$cells - 1L) %/% nf + 1L
  tr <- Matrix::mat2triplet(methods::as(methods::as(C, "CsparseMatrix"), "generalMatrix"))
  nz <- tr$x != 0
  key <- (tr$i[nz] - 1L) * nt + tc[tr$j[nz]]
  ukey <- sort(unique(key))
  po <- (ukey - 1L) %/% nt + 1L
  pt <- (ukey - 1L) %% nt + 1L
  Cp <- Matrix::sparseMatrix(i = match(key, ukey), j = tr$j[nz], x = tr$x[nz],
                             dims = c(length(ukey), length(B$cells)))
  Bc <- Ba[ac, , drop = FALSE]
  one <- tabulate(po, nrow(C))[po] == 1L
  s <- which(one)
  m <- which(!one)
  py <- rt_pairs(By)
  st <- list(Ba = Ba, By = By, cells = B$cells, nf = nf, ka = ka, ky = ky, Cp = Cp, Bc = Bc, po = po, pt = pt,
             s = s, m = m, py = py)
  # bandwidth: B-splines over periods that overlap, observations over
  # several periods, and the penalties
  h <- max(py$j - py$i)
  if (length(m)) {
    st$mobs <- unique(po[m])
    st$mi <- match(po[m], st$mobs)
    st$jj <- rep(seq_len(ky), each = ka)
    st$kk <- rep(seq_len(ka), ky)
    nzy <- By != 0
    lo <- apply(nzy, 1, function(v) min(which(v)))
    hi <- apply(nzy, 1, function(v) max(which(v)))
    h <- max(h, tapply(hi[pt[m]], st$mi, max) - tapply(lo[pt[m]], st$mi, min))
  }
  pz <- which(Reduce(`|`, lapply(pen, function(P) as.matrix(P) != 0)), arr.ind = TRUE)
  if (length(pz)) h <- max(h, abs((pz[, 1] - 1L) %/% ka - (pz[, 2] - 1L) %/% ka))
  st$band <- band_layout(ky, h, ka)
  st$pen <- lapply(pen, band_from_dense, band = st$band)
  if (length(s)) {
    pat <- as.matrix((Cp[s, , drop = FALSE] != 0) %*% (Bc != 0)) > 0
    ij <- which(crossprod(pat) > 0 & upper.tri(diag(ka), diag = TRUE), arr.ind = TRUE)
    st$pa <- list(i = ij[, 1], j = ij[, 2], k = ka)
    map <- glam_map(st$pa, py)
    K <- ka * ky
    r <- (map$idx - 1L) %% K
    cl <- (map$idx - 1L) %/% K
    jr <- r %/% ka + 1L
    jc <- cl %/% ka + 1L
    up <- jc >= jr
    # position in block row jr, block jc - jr
    st$bidx <- r[up] %% ka + 1L + ((jc[up] - jr[up]) * ka + cl[up] %% ka) * ka +
      (jr[up] - 1L) * ka^2 * (st$band$h + 1L)
    st$bsrc <- map$src[up]
  }
  env$C <- C
  env$pen <- pen
  env$st <- st
  st
}

# Fitted values and X'WX (within the band) at the coefficients `beta`, with
# the working matrix `A` [pair, age coefficient]
pclm2_system <- function(st, C, o, beta) {
  gam <- exp(glam_eta(st$Ba, st$By, beta)[st$cells])
  og <- o * gam
  mu <- pmax(as.numeric(C %*% og), 1e-12)
  A <- as.matrix(st$Cp %*% (og * st$Bc))
  w <- 1 / mu
  bd <- st$band
  XtWX <- array(0, c(st$ka, (bd$h + 1L) * st$ka, bd$nb))
  if (length(st$s)) {
    s <- st$s
    RT <- A[s, st$pa$i, drop = FALSE] * A[s, st$pa$j, drop = FALSE]
    Mt <- rowsum(w[st$po[s]] * RT, st$pt[s])
    M <- crossprod(Mt, st$py$RT[as.integer(rownames(Mt)), , drop = FALSE])
    XtWX[st$bidx] <- M[st$bsrc]
  }
  if (length(st$m)) {
    m <- st$m
    Xm <- rowsum(st$By[st$pt[m], st$jj, drop = FALSE] * A[m, st$kk, drop = FALSE], st$mi)
    wX <- w[st$mobs] * Xm
    for (j in seq_len(bd$nb)) {
      cols <- seq_len((min(bd$h, bd$nb - j) + 1L) * st$ka)
      XtWX[, cols, j] <- XtWX[, cols, j] + crossprod(Xm[, band_ix(j, st$ka), drop = FALSE],
                                                     wX[, (j - 1L) * st$ka + cols, drop = FALSE])
    }
  }
  list(gamma = gam, mu = mu, A = A, XtWX = XtWX)
}

# X'r for a vector r over observations
pclm2_xtr <- function(st, A, r) as.vector(crossprod(A, r[st$po] * st$By[st$pt, , drop = FALSE]))

# pclm_irls() for a tensor basis
pclm2_irls <- function(y, C, o, B, pen, lambda, beta, maxit = 100L, tol = 1e-6, se = FALSE) {
  st <- pclm2_setup(C, B, pen)
  bd <- st$band
  P <- Reduce(`+`, Map(`*`, lambda, st$pen))
  pdev <- function(b) {
    eta <- glam_eta(st$Ba, st$By, b)[st$cells]
    if (any(!is.finite(eta)) || max(eta) > 50) return(Inf)
    mu <- pmax(as.numeric(C %*% (o * exp(eta))), 1e-12)
    poisson_dev(y, mu) + band_quad(P, bd, b)
  }
  obj_old <- pdev(beta)
  converged <- FALSE
  for (it in seq_len(maxit)) {
    s <- pclm2_system(st, C, o, beta)
    A <- s$XtWX + P
    rhs <- pclm2_xtr(st, s$A, (y - s$mu) / s$mu) + band_mv(s$XtWX, bd, beta)
    f <- band_chol(A, bd)
    beta_new <- if (!is.null(f)) band_solve(f, bd, rhs) else
      tryCatch(as.numeric(solve(band_dense(A, bd), rhs)), error = function(e) NULL)
    if (is.null(beta_new) || anyNA(beta_new)) break
    step <- beta_new - beta
    obj_new <- pdev(beta + step)
    halvings <- 0L
    while (!(obj_new <= obj_old + 1e-10 * abs(obj_old)) && halvings < 25L) {
      step <- step / 2
      obj_new <- pdev(beta + step)
      halvings <- halvings + 1L
    }
    if (!is.finite(obj_new)) break
    beta <- beta + step
    if (max(abs(step)) < tol * (1 + max(abs(beta))) || abs(obj_old - obj_new) < 1e-8 * (1 + abs(obj_new))) {
      converged <- TRUE
      obj_old <- obj_new
      break
    }
    obj_old <- obj_new
  }
  s <- pclm2_system(st, C, o, beta)
  A <- s$XtWX + P
  f <- band_chol(A, bd)
  S <- if (!is.null(f)) {
    band_inverse(f, bd)
  } else {
    Ad <- band_dense(A, bd)
    band_from_dense(tryCatch(chol2inv(chol(Ad)), error = function(e) MASS_ginv(Ad)), bd)
  }
  ed <- sum(S * s$XtWX * bd$w)
  se_eta <- NULL
  if (se) se_eta <- sqrt(pmax(pclm2_var_eta(st, S), 0))
  list(gamma = s$gamma, beta = beta, mu = s$mu, lambda = lambda, ed = ed, dev = poisson_dev(y, s$mu),
       se_eta = se_eta, converged = converged, iterations = it)
}

# Variances of eta for the cells, the diagonal of B S B' for the rows b =
# by_t (x) ba_x of the tensor basis: the sum over the pairs (j, k) of
# B-splines over periods that overlap of by_t[j] by_t[k] ba_x' S[j, k] ba_x
pclm2_var_eta <- function(st, S) {
  py <- st$py
  Q <- vapply(seq_along(py$i), function(k) {
    rowSums((st$Ba %*% S[, band_cols(py$j[k] - py$i[k], st$ka), py$i[k]]) * st$Ba)
  }, numeric(st$nf))
  Q <- matrix(Q, st$nf)
  V <- Q %*% t(py$RT * rep(ifelse(py$i == py$j, 1, 2), each = nrow(py$RT)))
  V[st$cells]
}

# Block-banded symmetric matrices ----------------------------------------------
#
# A symmetric matrix of nb x nb blocks of size m whose blocks more than h
# off the diagonal are zero is kept as the array [m, (h + 1) m, nb] of its
# block rows: block row j holds the blocks [j, j], [j, j + 1], ...,
# [j, j + h] (zero beyond the last block). band_chol() gives the upper
# triangular Cholesky factor A = U'U in the same form (it has the same
# band), or NULL when A is not positive definite; band_solve() solves A x =
# b with it. band_inverse() gives the blocks of the inverse S = A^-1 within
# the band from U S = U^-T, block row j from the last up: U[j, j] S[j, i] =
# -sum_k U[j, k] S[k, i] for i > j, and S[j, j] = (U[j, j]' U[j, j])^-1 -
# U[j, j]^-1 sum_k U[j, k] S[k, j], k = j + 1, ..., j + h; these blocks give
# the trace of the hat matrix and the variances of eta.

# Layout: `w` weights the blocks of a block row in sums over the whole
# matrix (1 on the diagonal, 2 off it)
band_layout <- function(nb, h, m) {
  h <- as.integer(max(0L, min(h, nb - 1L)))
  list(nb = nb, h = h, m = m, w = rep(rep(c(1, rep(2, h)), each = m * m), nb))
}

band_seq <- function(from, to) if (to >= from) seq.int(from, to) else integer(0)

# Positions of block j in a vector, and of block [j, j + d] in block row j
band_ix <- function(j, m) (j - 1L) * m + seq_len(m)

band_cols <- function(d, m) d * m + seq_len(m)

band_chol <- function(A, band) {
  m <- band$m
  h <- band$h
  U <- array(0, dim(A))
  for (j in seq_len(band$nb)) {
    Ar <- A[, , j]
    for (k in band_seq(max(1L, j - h), j - 1L)) {
      # blocks [k, j], ..., [k, k + h] of row k update blocks [j, j], ..., [j, k + h]
      cols <- ((j - k) * m + 1L):((h + 1L) * m)
      upd <- seq_along(cols)
      Ar[, upd] <- Ar[, upd] - crossprod(U[, band_cols(j - k, m), k], U[, cols, k])
    }
    R <- tryCatch(chol(Ar[, seq_len(m)]), error = function(e) NULL)
    if (is.null(R)) return(NULL)
    U[, seq_len(m), j] <- R
    if (h > 0L) U[, -seq_len(m), j] <- backsolve(R, Ar[, -seq_len(m), drop = FALSE], transpose = TRUE)
  }
  U
}

band_solve <- function(U, band, b) {
  m <- band$m
  h <- band$h
  z <- c(b, numeric(h * m))
  y <- numeric(length(z))
  for (j in seq_len(band$nb)) {
    ix <- band_ix(j, m)
    y[ix] <- backsolve(U[, seq_len(m), j], z[ix], transpose = TRUE)
    if (h > 0L) {
      nx <- j * m + seq_len(h * m)
      z[nx] <- z[nx] - crossprod(U[, -seq_len(m), j], y[ix])
    }
  }
  x <- numeric(length(z))
  for (j in rev(seq_len(band$nb))) {
    ix <- band_ix(j, m)
    r <- y[ix]
    if (h > 0L) r <- r - U[, -seq_len(m), j] %*% x[j * m + seq_len(h * m)]
    x[ix] <- backsolve(U[, seq_len(m), j], r)
  }
  x[seq_len(band$nb * m)]
}

band_inverse <- function(U, band) {
  m <- band$m
  h <- band$h
  nb <- band$nb
  S <- array(0, dim(U))
  for (j in rev(seq_len(nb))) {
    R <- U[, seq_len(m), j]
    Sjj <- chol2inv(R)
    if (h > 0L && j < nb) {
      # the blocks S[k, i], k, i = j + 1, ..., j + h (zero beyond the last block)
      W <- matrix(0, h * m, h * m)
      for (a in seq_len(min(h, nb - j))) {
        cols <- ((a - 1L) * m + 1L):(h * m)
        W[band_ix(a, m), cols] <- S[, seq_along(cols), j + a]
      }
      W[lower.tri(W)] <- t(W)[lower.tri(W)]
      Uo <- U[, -seq_len(m), j]
      So <- -backsolve(R, Uo %*% W)
      S[, -seq_len(m), j] <- So
      Sjj <- Sjj - backsolve(R, Uo %*% t(So))
    }
    S[, seq_len(m), j] <- (Sjj + t(Sjj)) / 2
  }
  S
}

# A x and x'A x for a banded matrix A: block row j adds A[j, j + d] x[j +
# d] to block j and A[j, j + d]' x[j] to block j + d (d > 0)
band_mv <- function(A, band, x) {
  m <- band$m
  span <- (band$h + 1L) * m
  xp <- c(x, numeric(band$h * m))
  y <- numeric(length(xp))
  for (j in seq_len(band$nb)) {
    ix <- band_ix(j, m)
    win <- (j - 1L) * m + seq_len(span)
    Aj <- A[, , j]
    y[ix] <- y[ix] + Aj %*% xp[win]
    if (band$h > 0L) {
      off <- win[-seq_len(m)]
      y[off] <- y[off] + crossprod(Aj[, -seq_len(m), drop = FALSE], xp[ix])
    }
  }
  y[seq_len(band$nb * m)]
}

band_quad <- function(A, band, x) sum(x * band_mv(A, band, x))

# Dense matrix of a banded matrix, and the blocks within the band of a
# dense matrix
band_dense <- function(A, band) {
  m <- band$m
  D <- matrix(0, m * band$nb, m * band$nb)
  for (j in seq_len(band$nb)) {
    cols <- seq_len((min(band$h, band$nb - j) + 1L) * m)
    D[band_ix(j, m), (j - 1L) * m + cols] <- A[, cols, j]
  }
  D[lower.tri(D)] <- t(D)[lower.tri(D)]
  D
}

band_from_dense <- function(D, band) {
  m <- band$m
  D <- as.matrix(D)
  A <- array(0, c(m, (band$h + 1L) * m, band$nb))
  for (j in seq_len(band$nb)) {
    cols <- seq_len((min(band$h, band$nb - j) + 1L) * m)
    A[, cols, j] <- D[band_ix(j, m), (j - 1L) * m + cols]
  }
  A
}

#' B-spline basis and difference penalty on cell positions
#' @keywords internal
pclm_basis <- function(x, lo, hi, kr = 2, deg = 3, pord = 2) {
  span <- hi - lo
  nseg <- max(2L, ceiling(span / kr))
  dx <- span / nseg
  knots <- seq(lo - deg * dx, hi + deg * dx, by = dx)
  B <- splines::splineDesign(knots, x, ord = deg + 1, outer.ok = TRUE)
  K <- ncol(B)
  D <- diff(diag(K), differences = min(pord, K - 1))
  list(B = B, P = crossprod(D))
}

#' Fine latent axis: grid cells, with an open last cell expanded to omega
#' @keywords internal
fine_axis <- function(age, omega) {
  n <- length(age)
  if (!is_open(age)) {
    return(list(lower = age$lower, upper = age$upper, grid = seq_len(n), open = FALSE))
  }
  lo_open <- age$lower[n]
  if (omega <= lo_open + 1) {
    return(list(lower = age$lower, upper = c(age$upper[-n], omega), grid = seq_len(n), open = FALSE))
  }
  ext <- seq(lo_open, omega - 1)
  list(lower = c(age$lower[-n], ext), upper = c(age$upper[-n], ext + 1),
       grid = c(seq_len(n - 1L), rep(n, length(ext))), open = TRUE)
}

# 0/1 matrix linking fine cells to observed intervals
link_matrix <- function(lo, hi, fine_lo, fine_hi) {
  C <- matrix(0, length(lo), length(fine_lo))
  for (i in seq_along(lo)) {
    C[i, fine_lo >= lo[i] - 1e-9 & fine_hi <= hi[i] + 1e-9] <- 1
  }
  C
}

# Stationary-population weights within the open interval of the grid
# (mt_ungroup(open_split = "stationary")): the exposure of the open interval
# is split over its single years in proportion to l(x) = exp(-cumsum(m))
# implied by the fitted rates (package-specific, reference unknown); see
# ungroup_rate_block() for the iteration
open_weights <- function(rate) {
  l <- exp(-c(0, cumsum(rate[-length(rate)])))
  l / sum(l)
}
