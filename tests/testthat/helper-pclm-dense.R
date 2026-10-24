# Dense reference of the two-dimensional PCLM ---------------------------------
#
# The implementation of mt_ungroup(dims = 2) and mt_smooth(dims = 2) in
# mortensor 0.2.0: pclm_fit(), pclm_irls() and solve_lsq() of
# R/harmonise-pclm.R, ungroup_rates() and ungroup_rate_block() of
# R/harmonise-ungroup.R, copied without changes except the prefix ref_, the
# namespace of the package functions they call and the wrapping of two long
# lines. It forms the tensor basis kronecker(B_period, B_age) and the
# products with it as dense matrices. The package computes the same model
# with the arithmetic of generalised linear array models; the tests in
# test-harmonise-pclm2d.R compare the two with the functions of the last
# section.

ref_pclm_fit <- function(y, C, o, B, pen, lambda = NULL, criterion = c("bic", "aic"),
                     lambda_range = c(1e-2, 1e8), beta = NULL) {
  criterion <- match.arg(criterion)
  np <- length(pen)
  if (is.null(beta)) {
    tot <- sum(as.numeric(C %*% o))
    beta <- ref_solve_lsq(B, rep(log(max(sum(y), 0.5) / tot), nrow(B)))
  }
  ic_of <- function(fit) fit$dev + (if (criterion == "bic") log(length(y)) else 2) * fit$ed
  if (!is.null(lambda)) {
    fit <- ref_pclm_irls(y, C, o, B, pen, rep_len(lambda, np), beta, se = TRUE)
    fit$ic <- ic_of(fit)
    fit$at_bound <- FALSE
    return(fit)
  }
  lr <- log10(lambda_range)
  state <- new.env()
  state$beta <- beta
  obj <- function(ll) {
    f <- ref_pclm_irls(y, C, o, B, pen, 10^ll, state$beta)
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
  fit <- ref_pclm_irls(y, C, o, B, pen, 10^ll, state$beta, se = TRUE)
  fit$ic <- ic_of(fit)
  fit$at_bound <- any(ll < lr[1] + 0.1 | ll > lr[2] - 0.1)
  fit
}

ref_pclm_irls <- function(y, C, o, B, pen, lambda, beta, maxit = 100L, tol = 1e-6, se = FALSE) {
  P <- Reduce(`+`, Map(`*`, lambda, pen))
  pdev <- function(b) {
    eta <- as.numeric(B %*% b)
    if (any(!is.finite(eta)) || max(eta) > 50) return(Inf)
    mu <- pmax(as.numeric(C %*% (o * exp(eta))), 1e-12)
    mortensor:::poisson_dev(y, mu) + sum(b * as.numeric(P %*% b))
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
  Ainv <- tryCatch(chol2inv(chol(A)), error = function(e) mortensor:::MASS_ginv(A))
  ed <- sum(Ainv * XtWX)
  se_eta <- NULL
  if (se) {
    BV <- as.matrix(B %*% Ainv)
    se_eta <- sqrt(pmax(rowSums(BV * as.matrix(B)), 0))
  }
  list(gamma = gam, beta = beta, mu = mu, lambda = lambda, ed = ed, dev = mortensor:::poisson_dev(y, mu),
       se_eta = se_eta, converged = converged, iterations = it)
}

ref_solve_lsq <- function(B, z) {
  B <- as.matrix(B)
  as.numeric(mortensor:::MASS_ginv(crossprod(B)) %*% crossprod(B, z))
}

ref_ungroup_rates <- function(D, E, g, omega, ctl, dims) {
  idx <- mortensor:::obs_grid_index(D, g, "deaths")
  D <- cbind(D, idx)[!is.na(a1)]
  Earr <- mortensor:::obs_to_matrix(E, g)
  fa <- mortensor:::fine_axis(g$age, omega)
  xa <- (fa$lower + fa$upper) / 2
  n_open <- sum(fa$grid == length(g$age) & fa$open)
  bas <- mortensor:::pclm_basis(xa, min(fa$lower), max(fa$upper), ctl$kr, ctl$deg, ctl$pord)
  sep0 <- isTRUE(ctl$infant) && g$age$lower[1] == 0 && g$age$upper[1] == 1 && length(g$age) > 3
  if (sep0) {
    k1 <- fa$grid > 1L
    fa1 <- list(lower = fa$lower[k1], upper = fa$upper[k1], grid = fa$grid[k1], open = fa$open)
    bas1 <- mortensor:::pclm_basis(xa[k1], min(fa1$lower), max(fa1$upper), ctl$kr, ctl$deg, ctl$pord)
  }
  # rows of age 0 set aside: only when age 0 is never part of a wider cell
  infant <- function(o) {
    inf <- sep0 & o$a1 == 1L & o$a2 == 1L & o$p1 == o$p2
    if (all(inf) || any(o$a1 == 1L & !inf)) inf[] <- FALSE
    inf
  }
  out <- list()
  dg <- list()
  for (s in unique(D$stratum)) {
    Ds <- D[stratum == s]
    if (dims == 1) {
      if (any(Ds$p1 != Ds$p2)) {
        mortensor:::mt_abort("Deaths have periods coarser than the grid; use dims = 2.", class = "ungroup")
      }
      for (p in sort(unique(Ds$p1))) {
        o <- Ds[p1 == p]
        e_grid <- Earr[, p, s]
        inf <- infant(o)
        if (any(inf)) {
          r <- ref_ungroup_rate_block(o[!inf], list(e_grid), fa1, bas1$B, list(bas1$P), ctl, omega, n_open)
        } else {
          r <- ref_ungroup_rate_block(o, list(e_grid), fa, bas$B, list(bas$P), ctl, omega, n_open)
        }
        out[[length(out) + 1L]] <- data.table::data.table(stratum = s, a = c(o$a1[inf], r$a), p = p,
                                              value = c(o$value[inf], r$deaths),
                                              provenance = c(o$provenance[inf], r$provenance))
        dg[[length(dg) + 1L]] <- data.table::data.table(measure = "deaths", stratum = s, period = g$period$labels[p],
                                            method = paste0(ctl$method_label %||% "pclm_rate", "_1d"),
                                            lambda_age = r$fit$lambda[1],
                                            lambda_period = NA_real_, ed = r$fit$ed,
                                            max_rel_diff = r$rel, at_bound = isTRUE(r$fit$at_bound),
                                            open_iterations = r$open_iterations, open_converged = r$open_converged)
      }
    } else {
      ps <- seq(min(Ds$p1), max(Ds$p2))
      xp <- (g$period$lower[ps] + g$period$upper[ps]) / 2
      bp <- mortensor:::pclm_basis(xp, g$period$lower[min(ps)], g$period$upper[max(ps)], ctl$kr, ctl$deg, ctl$pord)
      inf <- infant(Ds)
      fb <- if (any(inf)) list(fa = fa1, bas = bas1) else list(fa = fa, bas = bas)
      B2 <- kronecker(bp$B, fb$bas$B)
      pen <- list(kronecker(diag(ncol(bp$B)), fb$bas$P), kronecker(bp$P, diag(ncol(fb$bas$B))))
      r <- ref_ungroup_rate_block(Ds[!inf], lapply(ps, function(p) Earr[, p, s]), fb$fa, B2, pen, ctl, omega, n_open,
                              ps)
      out[[length(out) + 1L]] <- data.table::data.table(stratum = s, a = c(Ds$a1[inf], r$a), p = c(Ds$p1[inf], r$p),
                                            value = c(Ds$value[inf], r$deaths),
                                            provenance = c(Ds$provenance[inf], r$provenance))
      dg[[length(dg) + 1L]] <- data.table::data.table(measure = "deaths", stratum = s, period = "all",
                                          method = paste0(ctl$method_label %||% "pclm_rate", "_2d"),
                                          lambda_age = r$fit$lambda[1],
                                          lambda_period = r$fit$lambda[2],
                                          ed = r$fit$ed, max_rel_diff = r$rel, at_bound = isTRUE(r$fit$at_bound),
                                          open_iterations = r$open_iterations, open_converged = r$open_converged)
    }
  }
  list(obs = mortensor:::grid_table(data.table::rbindlist(out), g), diag = data.table::rbindlist(dg, fill = TRUE))
}

ref_ungroup_rate_block <- function(o, e_list, fa, B, pen, ctl, omega, n_open, ps = NULL) {
  ps <- ps %||% o$p1[1]
  nf <- length(fa$lower)
  ng <- max(fa$grid)
  # composite link: rows = observed cells; columns = fine cells x periods
  rows <- lapply(seq_len(nrow(o)), function(i) {
    cols_a <- which(fa$lower >= o$age_lo[i] - 1e-9 & fa$upper <= pmin(o$age_hi[i], omega) + 1e-9)
    cols_p <- match(seq(o$p1[i], o$p2[i]), ps)
    as.vector(outer(cols_a, (cols_p - 1L) * nf, `+`))
  })
  C <- Matrix::sparseMatrix(i = rep(seq_along(rows), lengths(rows)), j = unlist(rows), x = 1,
                            dims = c(nrow(o), nf * length(ps)))
  e_grid <- do.call(cbind, e_list)
  if (anyNA(e_grid)) {
    mortensor:::mt_abort("Exposures are missing for grid cells needed by the ungrouping.", class = "ungroup")
  }
  open_idx <- which(fa$grid == ng & fa$open)
  w <- matrix(1, nf, length(ps))
  stationary <- length(open_idx) && identical(ctl$open_split, "stationary")
  if (length(open_idx)) {
    if (stationary) {
      # starting weights: rates rising by 10% per year of age in the open interval
      w[open_idx, ] <- mortensor:::open_weights(sum(o$value) / sum(e_grid) * exp(0.1 * seq_along(open_idx)))
    } else {
      # the exposure of the open interval split by the PCLM of the exposures
      # of the grid cells (Rizzi, Gampe and Eilers 2015)
      for (j in seq_along(ps)) w[open_idx, j] <- mortensor:::open_weights_pclm(e_grid[, j], fa, ctl)
    }
  }
  fit_w <- function(w, lambda, beta) {
    off <- as.vector(e_grid[fa$grid, , drop = FALSE] * w)
    cov <- Matrix::colSums(C) > 0 & off > 0
    fit <- ref_pclm_fit(o$value, C[, cov, drop = FALSE], off[cov], B[cov, , drop = FALSE], pen, lambda,
                    ctl$criterion, ctl$lambda_range, beta = beta)
    list(fit = fit, off = off, cov = cov)
  }
  r <- fit_w(w, ctl$lambda, NULL)
  open_it <- NA_integer_
  open_conv <- NA
  if (stationary) {
    # the exposure of the open interval is split with the stationary weights
    # of the fitted rates: weights and fit are iterated to a fixed point, with
    # the smoothing parameters fixed; they are then selected again with the
    # converged weights, until the selection no longer changes
    open_it <- 0L
    open_conv <- TRUE
    select <- is.null(ctl$lambda)
    at_bound <- isTRUE(r$fit$at_bound)
    for (round in seq_len(if (select) 3L else 1L)) {
      lam <- r$fit$lambda
      repeat {
        rate_m <- matrix(exp(as.numeric(B %*% r$fit$beta)), nf, length(ps))
        w_new <- w
        for (j in seq_along(ps)) w_new[open_idx, j] <- mortensor:::open_weights(rate_m[open_idx, j])
        if (max(abs(w_new - w)) < 1e-6) break
        if (open_it >= 50L) {
          open_conv <- FALSE
          break
        }
        open_it <- open_it + 1L
        w <- w_new
        r <- fit_w(w, lam, r$fit$beta)
      }
      if (!select || !open_conv || round == 3L) break
      r2 <- fit_w(w, NULL, r$fit$beta)
      at_bound <- isTRUE(r2$fit$at_bound)
      if (all(abs(log10(r2$fit$lambda / lam)) < 0.01)) break
      r <- r2
    }
    r$fit$at_bound <- at_bound
  }
  fit <- r$fit
  off <- r$off
  cov <- r$cov
  deaths_f <- rep(0, length(off))
  deaths_f[cov] <- off[cov] * fit$gamma
  rel <- max(abs(fit$mu / pmax(o$value, 1e-12) - 1)[o$value > 0], 0)
  if (ctl$preserve_totals) {
    deaths_f[cov] <- mortensor:::rescale_to_groups(deaths_f[cov], C[, cov, drop = FALSE], o$value)
  }
  dm <- matrix(deaths_f, nf, length(ps))
  dgrid <- rowsum(dm, fa$grid)
  covered <- rowsum(matrix(as.numeric(Matrix::colSums(C) > 0), nf, length(ps)), fa$grid) > 0
  # grid cells of the fine axis (from age 1 when age 0 is set aside)
  ga <- sort(unique(fa$grid))
  a <- rep(ga, length(ps))
  p <- rep(ps, each = length(ga))
  keep <- as.vector(covered)
  prov <- rep(ctl$provenance %||% "pclm", length(a))
  vals <- as.vector(dgrid)
  if (!isFALSE(ctl$keep_observed)) {
    exact <- o[a1 == a2 & p1 == p2]
    m <- match(paste(exact$a1, exact$p1), paste(a, p))
    vals[m] <- exact$value
    prov[m] <- "observed"
  }
  list(a = a[keep], p = p[keep], deaths = vals[keep], provenance = prov[keep], fit = fit, rel = rel,
       open_iterations = open_it, open_converged = open_conv)
}

# Comparison of the package with the reference --------------------------------

# largest relative difference: per element (missing in both: no
# difference), and relative to the largest element (for coefficients, which
# can be close to 0)
rel_diff <- function(a, b) {
  if (!identical(is.na(a), is.na(b))) return(Inf)
  ok <- !is.na(b)
  max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-300), 0)
}
norm_diff <- function(a, b) max(abs(a - b)) / max(abs(b))

# the data and diagnostics of a 2D fit of the package and of the reference
expect_same_2d <- function(new, ref, what, tol = 1e-8) {
  dn <- mt_obs(new, "deaths")
  dr <- mt_obs(ref, "deaths")
  testthat::expect_identical(dn[, c("stratum", "age_lo", "age_hi", "period_lo", "period_hi", "provenance")],
                   dr[, c("stratum", "age_lo", "age_hi", "period_lo", "period_hi", "provenance")])
  testthat::expect_lt(rel_diff(dn$value, dr$value), tol)
  gn <- attr(new, what)
  gr <- attr(ref, what)
  testthat::expect_identical(gn[, c("stratum", "method", "at_bound")], gr[, c("stratum", "method", "at_bound")])
  for (v in c("lambda_age", "lambda_period", "ed", "max_rel_diff")) testthat::expect_lt(rel_diff(gn[[v]], gr[[v]]), tol)
}

# the fit objects of pclm_fit(): coefficients, latent rates, fitted counts,
# standard errors of the linear predictor, smoothing parameters and
# effective dimension
expect_same_fit <- function(fn, fr, tol = 1e-8) {
  testthat::expect_lt(norm_diff(fn$beta, fr$beta), tol)
  for (v in c("gamma", "mu", "se_eta", "lambda", "ed", "dev")) testthat::expect_lt(rel_diff(fn[[v]], fr[[v]]), tol)
  testthat::expect_identical(fn$converged, fr$converged)
  testthat::expect_identical(fn$iterations, fr$iterations)
}

# inputs of ungroup_rate_block() for the 2D PCLM of the deaths of stratum s,
# as ungroup_rates() builds them
block_inputs <- function(x, s = 1L) {
  g <- mt_grid(x)
  D <- mt_obs(x, "deaths")
  D <- cbind(D, mortensor:::obs_grid_index(D, g, "deaths"))
  D <- D[D$stratum == s & !is.na(D$a1)]
  fa <- mortensor:::fine_axis(g$age, x$omega)
  ps <- seq(min(D$p1), max(D$p2))
  ba <- mortensor:::pclm_basis((fa$lower + fa$upper) / 2, min(fa$lower), max(fa$upper))
  bp <- mortensor:::pclm_basis((g$period$lower[ps] + g$period$upper[ps]) / 2, g$period$lower[min(ps)],
                               g$period$upper[max(ps)])
  E <- mortensor:::obs_to_matrix(mt_obs(x, "exposure"), g)
  list(o = D, e = lapply(ps, function(p) E[, p, s]), fa = fa, ba = ba, bp = bp, ps = ps, omega = x$omega,
       pen = list(kronecker(diag(ncol(bp$B)), ba$P), kronecker(bp$P, diag(ncol(ba$B)))))
}

pclm_ctl <- function(lambda, ...) {
  utils::modifyList(list(lambda = lambda, criterion = "bic", kr = 2, deg = 3, pord = 2, lambda_range = c(1e-2, 1e8),
                         preserve_totals = FALSE, open_split = "pclm"), list(...))
}

# one block of the package (tensor basis) and of the reference (dense basis)
compare_block <- function(b, ctl) {
  new <- mortensor:::ungroup_rate_block(b$o, b$e, b$fa, mortensor:::pclm_tensor(b$ba$B, b$bp$B), b$pen, ctl,
                                        b$omega, NULL, b$ps)
  ref <- ref_ungroup_rate_block(b$o, b$e, b$fa, kronecker(b$bp$B, b$ba$B), b$pen, ctl, b$omega, NULL, b$ps)
  expect_same_fit(new$fit, ref$fit)
  testthat::expect_lt(rel_diff(new$deaths, ref$deaths), 1e-8)
  invisible(list(new = new, ref = ref))
}
