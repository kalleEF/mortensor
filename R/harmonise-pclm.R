# Penalized composite link model (PCLM) ----------------------------------------
#
# Method: Eilers (2007), Rizzi, Gampe and Eilers (2015). The structure (B-spline
# basis on the latent scale, difference penalty, IRLS, lambda chosen by an
# information criterion) follows the 'ungroup' package 1.6.3 (MIT licence;
# M. D. Pascariu, S. Rizzi, J. Schoeley, M. J. Danko). This implementation is
# written for mortensor and differs in the treatment of rates: the exposure
# offset enters the composite link directly and the open age interval is
# expanded into single years up to omega with stationary-population weights,
# so no artificial-bin reallocation is applied to rates (the source of the
# level shift found in revamp plan 5.1).

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
#' @param B Basis matrix (fine cells x coefficients).
#' @param pen List of penalty matrices (one per smoothing parameter).
#' @param lambda Fixed smoothing parameters, or `NULL` to select them.
#' @param criterion `"bic"` or `"aic"`.
#' @param lambda_range Search range for each smoothing parameter.
#' @param beta Optional starting values.
#' @return A list with `gamma`, `beta`, `mu` (fitted observed counts),
#'   `lambda`, `ed` (effective dimension), `dev`, `ic`, `se_eta`,
#'   `converged` and `at_bound`.
#' @keywords internal
pclm_fit <- function(y, C, o, B, pen, lambda = NULL, criterion = c("bic", "aic"),
                     lambda_range = c(1e-2, 1e8), beta = NULL) {
  criterion <- match.arg(criterion)
  np <- length(pen)
  if (is.null(beta)) {
    tot <- sum(as.numeric(C %*% o))
    beta <- solve_lsq(B, rep(log(max(sum(y), 0.5) / tot), nrow(B)))
  }
  ic_of <- function(fit) fit$dev + (if (criterion == "bic") log(length(y)) else 2) * fit$ed
  if (!is.null(lambda)) {
    fit <- pclm_irls(y, C, o, B, pen, rep_len(lambda, np), beta)
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
        ll[k] <- stats::optimize(function(z) { l2 <- ll; l2[k] <- z; obj(l2) }, lr, tol = 0.05)$minimum
      }
    }
  }
  fit <- pclm_irls(y, C, o, B, pen, 10^ll, state$beta)
  fit$ic <- ic_of(fit)
  fit$at_bound <- any(ll < lr[1] + 0.1 | ll > lr[2] - 0.1)
  fit
}

pclm_irls <- function(y, C, o, B, pen, lambda, beta, maxit = 100L, tol = 1e-6) {
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
    beta_new <- tryCatch(as.numeric(solve(A, rhs)), error = function(e) NULL)
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
  BV <- as.matrix(B %*% Ainv)
  se_eta <- sqrt(pmax(rowSums(BV * as.matrix(B)), 0))
  list(gamma = gam, beta = beta, mu = mu, lambda = lambda, ed = ed, dev = poisson_dev(y, mu),
       se_eta = se_eta, converged = converged, iterations = it)
}

Matrix_diag <- function(x) Matrix::Diagonal(x = x)

poisson_dev <- function(y, mu) {
  r <- ifelse(y > 0, y * log(y / mu), 0)
  2 * sum(r - (y - mu))
}

solve_lsq <- function(B, z) {
  B <- as.matrix(B)
  as.numeric(qr.solve(crossprod(B) + diag(1e-8, ncol(B)), crossprod(B, z)))
}

MASS_ginv <- function(A) {
  s <- svd(A)
  pos <- s$d > max(s$d) * 1e-10
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
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

# Stationary-population weights within the open interval
open_weights <- function(rate) {
  l <- exp(-c(0, cumsum(rate[-length(rate)])))
  l / sum(l)
}
