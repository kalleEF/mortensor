# Filtered jump-off of the Lee-Carter model (De Jong and Tickle 2006).
#
# The state-space form of the Lee-Carter model (their eq. 6), with the
# age effects a and b of a Lee-Carter fit:
#   y[t] = a + b alpha[t] + eps[t],          eps[t] ~ N(0, sigma^2 I),
#   alpha[t + 1] = alpha[t] + delta + w[t],  w[t] ~ N(0, lambda^2 sigma^2 / b'b),
# where y[t] are the log rates of year t and lambda is the paper's
# signal-to-noise ratio, defined with b'b = 1 (alpha here has the scale of
# the package's k, sum(b) = 1). Given a and b, y[t] splits into the per-year
# index z[t] = b'(y[t] - a) / b'b = alpha[t] + e[t], e[t] ~ N(0, sigma^2 / b'b),
# and the residuals r[t] = y[t] - a - b z[t], orthogonal to b and
# independent of z, which carry sigma^2 only. With the level alpha[1] and
# the drift delta as unknown constants, z ~ N(X beta, sigma^2 / b'b V) with
# X = [1, t - 1], beta = (alpha[1], delta) and V = I + lambda^2 M,
# M[t, s] = min(t, s) - 1 (eq. 8 of the paper). The log-likelihood profiled
# over beta (generalised least squares) and sigma^2,
#   -(n p / 2) log((S_r + b'b q) / (n p)) - log det(V) / 2,
# with S_r the sum of squared residuals r and q the GLS residual sum of
# squares of z, is maximised over lambda (Section 5, steps 1 and 2, with b
# fixed). The filtered index of the last year, E(alpha[n] | y[1], ..., y[n])
# (eq. 15), is the best linear unbiased prediction
#   x[n]' beta + c' V^-1 (z - X beta),  c[t] = lambda^2 (t - 1),
# which equals the Kalman filter with a diffuse initial state for the level
# and the drift (De Jong 1991); its mean squared error includes the
# uncertainty of beta.
#
# y: log rates [age, year] of the ages with deaths; a, b: their Lee-Carter
# age effects. `lambda`: a fixed signal-to-noise ratio instead of its
# estimate (0: no state noise; large values: little observation noise).
lc_filter <- function(y, a, b, lambda = NULL) {
  y <- as.matrix(y)
  n <- ncol(y)
  p <- nrow(y)
  if (n < 3) mt_abort("The filtered jump-off needs at least three years.", class = "model_fit")
  bb <- sum(b^2)
  if (!is.finite(bb) || bb <= 0) mt_abort("The filtered jump-off needs an age pattern b with non-zero values.",
                                         class = "model_fit")
  z <- colSums(b * (y - a)) / bb
  s_r <- sum((y - a - outer(b, z))^2)
  tt <- seq_len(n) - 1
  X <- cbind(1, tt)
  M <- outer(tt, tt, pmin)
  gls <- function(lam2) {
    R <- chol(diag(n) + lam2 * M)
    Xw <- backsolve(R, X, transpose = TRUE)
    zw <- backsolve(R, z, transpose = TRUE)
    beta <- qr.coef(qr(Xw), zw)
    list(R = R, Xw = Xw, beta = beta, q = sum((zw - Xw %*% beta)^2), ldet = 2 * sum(log(diag(R))))
  }
  profile <- function(lam2) {
    g <- gls(lam2)
    -0.5 * (n * p * log((s_r + bb * g$q) / (n * p)) + g$ldet)
  }
  if (is.null(lambda)) {
    # maximum over log(lambda^2): a grid, then a one-dimensional search
    # around the best grid point; lambda = 0 (no state noise) is a candidate
    grid <- seq(-12, 12, by = 1)
    ll <- vapply(exp(grid), profile, 0)
    i <- which.max(ll)
    lo <- grid[max(i - 1, 1)]
    hi <- grid[min(i + 1, length(grid))]
    opt <- stats::optimize(function(u) profile(exp(u)), c(lo, hi), maximum = TRUE, tol = 1e-8)
    lam2 <- if (profile(0) >= opt$objective) 0 else exp(opt$maximum)
  } else {
    check_number(lambda, "lambda", lower = 0)
    lam2 <- lambda^2
  }
  g <- gls(lam2)
  sigma2 <- (s_r + bb * g$q) / (n * p)
  # V^-1 c and V^-1 (z - X beta) from the Cholesky factor
  vinv <- function(v) backsolve(g$R, backsolve(g$R, v, transpose = TRUE))
  cc <- lam2 * tt
  vc <- vinv(cc)
  xn <- c(1, n - 1)
  level <- sum(xn * g$beta) + sum(vc * (z - X %*% g$beta))
  u <- xn - crossprod(X, vc)
  xvx <- crossprod(g$Xw)
  mse <- sigma2 / bb * (lam2 * (n - 1) - sum(cc * vc) + as.numeric(crossprod(u, solve(xvx, u))))
  list(level = level, se = sqrt(max(mse, 0)), drift = unname(g$beta[2]), lambda = sqrt(lam2),
       sigma = sqrt(sigma2), z = z, loglik = profile(lam2))
}
