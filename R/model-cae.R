#' Common age effect model for several populations
#'
#' Joint model of the strata as populations that share their age effects and
#' keep their own period effects: the common age effect (CAE) model of
#' Kleinow (2015). Status: published, with the package-specific forecast and
#' choices listed below; section and equation numbers refer to that paper.
#'
#' **Model (Sections 2 and 3.1, eq. 1-3).** The log death rates `y[x, t, i]
#' = log(D[x, t, i] / E[x, t, i])` of each population `i` are centred over
#' the years: `m[x, t, i] = y[x, t, i] - a[x, i]`, with `a[x, i]` the mean
#' over the years. The centred rates follow
#' `m[, , i] = beta kappa[, , i] + e[, , i]`. The age effects `beta` (ages x
#' `order`) are common to all populations; the period effects `kappa[, , i]`
#' (`order` x years) belong to each population. With one population the
#' model is the Lee-Carter model with `order` age-period terms, estimated by
#' a singular value decomposition; for `order = 1` it gives the fitted rates
#' and forecasts of `lee_carter(method = "svd", adjust = "none")`.
#'
#' **Identification (Section 2).** The columns of `beta` have length 1 and
#' are orthogonal to each other. The period effects of each population sum
#' to zero over the years, because the rates are centred. The sign of each
#' column of `beta` is chosen so that its sum is positive (package-specific).
#'
#' **Estimation (Section 3.2).** With `Q[i] = m[, , i] t(m[, , i])`, `beta`
#' is the orthogonal matrix that minimises the sum over populations of the
#' squared off-diagonal elements of `t(beta) Q[i] beta`: the paper's
#' least-squares modification of the common principal component analysis of
#' Flury (1984). Every population has the same weight whatever its size, as
#' in the paper, which prefers this to maximum likelihood because large
#' populations would dominate the age effects. The first `order` columns
#' are kept, and `kappa[, , i] = t(beta) m[, , i]`. The paper requires more
#' years than age groups (Section 2), and so does `cae()`: with fewer years,
#' use wider age groups ([mt_regroup()]) or fewer ages ([mt_subset()]);
#' Kleinow fitted ages 18-52 and 53-87 separately.
#'
#' **Forecasts (package-specific).** Kleinow (2015) does not forecast. Here
#' the period effects of all populations and orders follow one multivariate
#' random walk with drift, `kappa[t] = kappa[t - 1] + mu + e[t]`, with normal
#' innovations whose covariance matrix is estimated from the year-to-year
#' changes (divided by the number of changes minus one), so the innovations
#' are correlated between populations and between orders. This is the
#' random walk with drift of Lee and Carter (1992) for every index, and the
#' multivariate random walk with drift that StMoMo uses for several period
#' indexes (Villegas, Kaishev and Millossovich 2018), here over all
#' populations together. Draws include the innovations and, as in
#' [lee_carter()], the uncertainty of the drift (normal, with the innovation
#' covariance divided by the number of changes; StMoMo keeps the drift
#' fixed). `a` and `beta` are fixed. Forecasts start from the fitted rates of
#' the last year ([anchored()] starts them from the observed rates). Each
#' population has its own drift, so the forecasts of different populations
#' can diverge; [li_lee()] and [coherent_fdm()] keep them together.
#' Predictive draws (`mt_forecast(predictive = TRUE)`) add Poisson variation
#' of the deaths given the rates and the exposures (the paper makes no
#' assumption about the distribution of the errors).
#'
#' **Other package-specific choices.**
#'
#' * The criterion is minimised by Jacobi rotations (Cardoso and Souloumiac
#'   1996), each of which minimises it in one plane, started from the
#'   eigenvectors of the sum of the `Q[i]`, until no rotation of more than
#'   1e-12 radians lowers the criterion by more than 1e-20 times the sum of
#'   squares of the matrices (at most 500 sweeps over all pairs of
#'   columns; a warning of class `mt_warning_cae` reports when this limit
#'   is reached). The paper uses
#'   the F-G algorithm (Flury and Constantine 1985) with the modification of
#'   Clarkson (1988). The criterion can have several local minima; the
#'   estimate is the one reached from this start.
#' * The columns of `beta` are ordered by the sum over populations of the
#'   diagonal elements of `t(beta) Q[i] beta`. The paper orders them by the
#'   diagonal of one population and notes that the order can differ between
#'   populations, with no general solution (Section 3.2).
#' * Zero deaths are replaced by 0.5 before taking logarithms.
#' * `common_by` (for example `"sex"`) gives each group of strata its own
#'   common age effects, estimated from the strata of the group. The paper
#'   names this extension and leaves it for further research (Section
#'   3.1). The period effects of all strata still form one random walk, with
#'   innovations correlated across groups.
#' * If the data carry a hierarchy (see [mt_aggregate_strata()]), `beta` is
#'   estimated from the bottom strata only, and the aggregates get their
#'   period effects by projection on it, as the other strata (as in
#'   [li_lee()]).
#'
#' **Choice of `order`.** The fit holds the mean squared error of the
#' centred log rates (eq. 4) and the approximate BIC of eq. (5), `N
#' log(MSE) + k log(N)` with `N` the number of rates and `k = order (K T +
#' n)` parameters for `K` populations, `T` years and `n` ages (one `beta` per
#' group of `common_by`); [mt_parameters()] reports them (rows `mse` and
#' `bic`). In the paper's application (men in ten countries), `order = 2`
#' had the lowest BIC in both age ranges (Table 2). Compare orders also in a
#' backtest.
#'
#' @param order Number of common age effects (`p` in the paper).
#' @param common_by Stratum variable(s) defining groups of strata with their
#'   own common age effects (package-specific, see Details).
#' @return A model specification. The fit (`fit$fits`) holds the common age
#'   effects `beta` (one matrix per group), the centring `a` `[age,
#'   stratum]`, the period effects `kappa` `[order, year, stratum]`, the
#'   drift `drift` and innovation covariance `sigma` of the stacked period
#'   effects (order within stratum), `mse` and `bic`.
#' @references
#' Kleinow T (2015). A common age effect model for the mortality of multiple
#' populations. *Insurance: Mathematics and Economics* 63: 147-152.
#' \doi{10.1016/j.insmatheco.2015.03.023}
#'
#' Flury BN (1984). Common principal components in k groups. *Journal of
#' the American Statistical Association* 79(388): 892-898.
#' \doi{10.1080/01621459.1984.10477108}
#'
#' Flury BN, Constantine G (1985). Algorithm AS 211: the F-G
#' diagonalization algorithm. *Applied Statistics* 34(2): 177-183.
#' \doi{10.2307/2347375}
#'
#' Clarkson DB (1988). Remark AS R71: a remark on algorithm AS 211. The F-G
#' diagonalization algorithm. *Applied Statistics* 37(1): 147-151.
#' \doi{10.2307/2347513}
#'
#' Cardoso J-F, Souloumiac A (1996). Jacobi angles for simultaneous
#' diagonalization. *SIAM Journal on Matrix Analysis and Applications*
#' 17(1): 161-164. \doi{10.1137/S0895479893259546}
#'
#' Lee RD, Carter LR (1992). Modeling and forecasting U.S. mortality.
#' *Journal of the American Statistical Association* 87(419): 659-671.
#'
#' Villegas AM, Kaishev VK, Millossovich P (2018). StMoMo: an R package for
#' stochastic mortality modeling. *Journal of Statistical Software* 84(3):
#' 1-38. \doi{10.18637/jss.v084.i03}
#' @seealso [li_lee()] and [coherent_fdm()] for coherent joint models;
#'   [coupled()] to give separate forecasts the dependence of a joint model.
#' @examples
#' # more years (30) than age groups (9): ages 50 and over in 5-year groups
#' x <- mt_subset(mt_regroup(mt_example(), age = "deaths"), ages = c(50, Inf))
#' fit <- mt_fit(cae(order = 2), x)
#' p <- mt_parameters(fit)
#' p[parameter %in% c("mse", "bic")]
#' fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
#' @export
cae <- function(order = 1, common_by = NULL) {
  check_number(order, "order", lower = 1, integerish = TRUE)
  if (!is.null(common_by) && (!is.character(common_by) || !length(common_by))) {
    mt_abort("`common_by` must name stratum variables.", class = "argument")
  }
  new_mt_spec("cae", list(order = as.integer(order), common_by = common_by), joint = TRUE, min_periods = 5,
              label = paste0("CAE (order ", order, if (!is.null(common_by)) {
                paste0(", common by ", paste(common_by, collapse = ", "))
              }, ")"))
}

#' @export
fit_model.mt_spec_cae <- function(spec, D, E, ages, periods, strata = NULL, hierarchy = NULL, ...) {
  p <- spec$params$order
  d <- dim(D)
  na <- d[1]
  nt <- d[2]
  ns <- d[3]
  if (nt <= na) {
    mt_abort(sprintf(paste("cae() needs more years than age groups (Kleinow 2015, Section 2); the data have %d",
                           "years and %d age groups. Use wider age groups (mt_regroup()) or fewer ages",
                           "(mt_subset(ages = ))."), nt, na), class = "model_data")
  }
  if (p > na) mt_abort(sprintf("`order` (%d) exceeds the number of age groups (%d).", p, na), class = "argument")
  y <- log_rates(D, E)
  a <- matrix(apply(y, c(1, 3), mean), na, ns)
  m <- sweep(y, c(1, 3), a)
  group <- rep(1L, ns)
  if (!is.null(spec$params$common_by)) {
    cb <- spec$params$common_by
    if (is.null(strata) || !all(cb %in% names(strata))) {
      mt_abort("`common_by` must name stratum variables.", class = "argument")
    }
    key <- do.call(paste, c(lapply(cb, function(v) strata[[v]]), sep = "."))
    group <- match(key, unique(key))
  }
  pool <- rep(TRUE, ns)
  if (!is.null(hierarchy) && !is.null(strata)) pool <- hierarchy_structure(strata, hierarchy)$is_bottom
  ng <- max(group)
  members <- lapply(seq_len(ng), function(g) {
    use <- which(group == g & pool)
    if (!length(use)) use <- which(group == g)
    use
  })
  cpca <- lapply(members, function(use) {
    Q <- lapply(use, function(s) tcrossprod(matrix(m[, , s], na, nt)))
    cpca_ls(Q)
  })
  if (any(!vapply(cpca, `[[`, logical(1), "converged"))) {
    mt_warn("CAE: the common principal components did not converge within 500 sweeps.", class = "cae")
  }
  beta <- lapply(cpca, function(cp) cp$beta[, seq_len(p), drop = FALSE])
  kappa <- array(NA_real_, c(p, nt, ns), dimnames = list(NULL, as.character(periods), NULL))
  fitted <- array(NA_real_, d)
  for (s in seq_len(ns)) {
    b <- beta[[group[s]]]
    kappa[, , s] <- crossprod(b, matrix(m[, , s], na, nt))
    fitted[, , s] <- exp(a[, s] + b %*% matrix(kappa[, , s], p))
  }
  # Approximate BIC of Kleinow (2015, eq. 4-5) over the strata that define
  # the age effects: N log(MSE) + k log(N), k = sum over groups of p (K T + n)
  used <- unlist(members)
  res <- log(fitted[, , used, drop = FALSE]) - y[, , used, drop = FALSE]
  n_obs <- length(res)
  mse <- sum(res^2) / n_obs
  k <- sum(vapply(members, function(u) p * (length(u) * nt + na), 0))
  bic <- n_obs * log(mse) + k * log(n_obs)
  # multivariate random walk with drift of the stacked period effects (the
  # orders of stratum 1, then of stratum 2, ...): K [order and stratum, year]
  K <- matrix(aperm(kappa, c(1, 3, 2)), p * ns)
  dK <- K[, -1, drop = FALSE] - K[, -nt, drop = FALSE]
  drift <- rowMeans(dK)
  R <- dK - drift
  sigma <- tcrossprod(R) / max(nt - 2, 1)
  list(beta = beta, a = a, kappa = kappa, group = group, members = members, drift = drift, sigma = sigma,
       root = mvn_root(sigma), last = K[, nt], n_changes = nt - 1, fitted = fitted, mse = mse, bic = bic,
       lambda = lapply(cpca, function(cp) cp$lambda[seq_len(p), , drop = FALSE]),
       sweeps = vapply(cpca, `[[`, 0, "sweeps"))
}

#' @export
simulate_model.mt_spec_cae <- function(spec, fit, h, n_draws, ...) {
  na <- nrow(fit$a)
  ns <- ncol(fit$a)
  p <- dim(fit$kappa)[1]
  nd <- length(fit$last)
  # cumulative sums over the horizon as a product with a triangular matrix
  U <- 1 * upper.tri(diag(h), diag = TRUE)
  out <- array(NA_real_, c(na, h, ns, n_draws))
  for (j in seq_len(n_draws)) {
    drift <- fit$drift + as.vector(fit$root %*% stats::rnorm(nd)) / sqrt(fit$n_changes)
    eps <- fit$root %*% matrix(stats::rnorm(nd * h), nd, h)
    path <- fit$last + (eps + drift) %*% U
    for (s in seq_len(ns)) {
      rows <- (s - 1L) * p + seq_len(p)
      out[, , s, j] <- exp(fit$a[, s] + fit$beta[[fit$group[s]]] %*% path[rows, , drop = FALSE])
    }
  }
  out
}

# Least-squares common principal components (Kleinow 2015, Section 3.2):
# the orthogonal matrix B that minimises the sum over the matrices Q[[i]] of
# the squared off-diagonal elements of t(B) Q[[i]] B, by Jacobi rotations
# (Cardoso and Souloumiac 1996). A rotation by `theta` in the plane of
# columns p and q changes the difference of the diagonal elements of each
# matrix to h1 cos(2 theta) + h2 sin(2 theta), with h1 = A[p, p] - A[q, q]
# and h2 = 2 A[p, q]; the sum of their squares, and so the diagonal part of
# the criterion, is largest for (cos(2 theta), sin(2 theta)) along the first
# eigenvector of sum(h h'), theta in (-pi/4, pi/4]. Started from the
# eigenvectors of the sum of the matrices (exact for one matrix). Columns
# are ordered by the sum over the matrices of their diagonal elements
# (`lambda`, [column, matrix]), with signs so that each column sums to a
# non-negative value
cpca_ls <- function(Q, tol = 1e-12, max_sweeps = 500) {
  n <- nrow(Q[[1]])
  k <- length(Q)
  V <- eigen(Reduce(`+`, Q), symmetric = TRUE)$vectors
  A <- array(unlist(lapply(Q, function(q) crossprod(V, q %*% V))), c(n, n, k))
  scale <- sum(A^2)
  sweeps <- 0L
  converged <- TRUE
  if (n > 1 && k > 1) {
    converged <- FALSE
    for (sw in seq_len(max_sweeps)) {
      sweeps <- sw
      rotated <- FALSE
      for (p in seq_len(n - 1L)) {
        for (q in (p + 1L):n) {
          h1 <- A[p, p, ] - A[q, q, ]
          h2 <- A[p, q, ] + A[q, p, ]
          g11 <- sum(h1^2)
          g22 <- sum(h2^2)
          g12 <- sum(h1 * h2)
          # gain of the diagonal part; a plane in which the matrices are
          # already diagonal with equal elements has no gain and an
          # undefined angle
          gain <- (g22 - g11) / 2 + sqrt(((g11 - g22) / 2)^2 + g12^2)
          theta <- 0.25 * atan2(2 * g12, g11 - g22)
          if (abs(theta) <= tol || gain <= 1e-20 * scale) next
          rotated <- TRUE
          cs <- cos(theta)
          sn <- sin(theta)
          Ap <- A[p, , , drop = FALSE]
          A[p, , ] <- cs * Ap + sn * A[q, , , drop = FALSE]
          A[q, , ] <- -sn * Ap + cs * A[q, , , drop = FALSE]
          Ap <- A[, p, , drop = FALSE]
          A[, p, ] <- cs * Ap + sn * A[, q, , drop = FALSE]
          A[, q, ] <- -sn * Ap + cs * A[, q, , drop = FALSE]
          vp <- V[, p]
          V[, p] <- cs * vp + sn * V[, q]
          V[, q] <- -sn * vp + cs * V[, q]
        }
      }
      if (!rotated) {
        converged <- TRUE
        break
      }
    }
  }
  lambda <- matrix(vapply(seq_len(k), function(i) diag(matrix(A[, , i], n)), numeric(n)), n)
  o <- order(rowSums(lambda), decreasing = TRUE)
  V <- V[, o, drop = FALSE]
  lambda <- lambda[o, , drop = FALSE]
  sg <- ifelse(colSums(V) < 0, -1, 1)
  V <- sweep(V, 2, sg, "*")
  off <- sum(vapply(seq_len(k), function(i) {
    a <- crossprod(V, Q[[i]] %*% V)
    sum(a^2) - sum(diag(a)^2)
  }, 0))
  list(beta = V, lambda = lambda, criterion = off, sweeps = sweeps, converged = converged)
}

# A square root L of a covariance matrix S (L t(L) = S): the transposed
# Cholesky factor, or from the eigen decomposition (negative eigenvalues of
# rounding set to zero) where S is singular, as with more indexes than
# year-to-year changes
mvn_root <- function(S) {
  R <- tryCatch(chol(S), error = function(e) NULL)
  if (!is.null(R)) return(t(R))
  e <- eigen(S, symmetric = TRUE)
  e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values))
}
