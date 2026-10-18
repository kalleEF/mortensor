# Numerical core of the compositional (CoDa) models coda(), coda(coherent =
# TRUE) and coda_causes(), written from Bergeron-Boucher, Canudas-Romo,
# Oeppen and Vaupel (2017) and Oeppen (2008). The supplementary R code of
# Bergeron-Boucher et al. (2017) (CoDaForecast.R, Demographic Research)
# states no licence; it is not copied here, and its results serve as test
# references (tests/testthat/test-coda.R).

# Layout of the composition ----------------------------------------------------

# Parts of the composition of life-table deaths for an age axis, with the
# method for the old ages (`open`: "auto", "kannisto" or "regression"):
# * "kannisto": the rates from age 80 are those of a Kannisto model fitted
#   per year to the single ages 80-100, extended to age 120, where the life
#   table closes (Bergeron-Boucher et al. 2017, Appendix C); parts are the age
#   groups below 80 and the single ages 80-120;
# * "regression": an open interval without the Kannisto extension; parts are
#   the age groups (the deaths of the open interval are the survivors at its
#   lower bound), and the log rate of the open interval is regressed on the
#   indexes (package-specific);
# * "closed": a closed age range without the Kannisto extension; parts are
#   the age groups and the survivors at the upper bound.
coda_layout <- function(ages, open = "auto") {
  na <- length(ages)
  lo <- ages$lower
  up <- ages$upper
  is_open <- is.infinite(up[na])
  kan <- coda_kannisto_ages(ages)
  if (open == "kannisto" && is.null(kan)) {
    mt_abort(paste("open = \"kannisto\" needs single-year age groups from 80 to at least 89 and an open interval",
                   "(or upper bound) at or below age 120; use open = \"regression\" for other age groups."),
             class = "model_data", call = NULL)
  }
  method <- if (open != "regression" && !is.null(kan)) "kannisto" else if (is_open) "regression" else "closed"
  out <- list(ages = ages, method = method)
  if (method == "kannisto") {
    below <- which(lo < 80)
    ext <- 80:120
    out$below <- below
    out$fit <- kan$fit
    out$width <- c(up[below] - lo[below], rep(1, length(ext)))
    out$parts <- c(ages$labels[below], as.character(ext))
    # part of each closed grid age from 80, and the parts of the open interval
    out$grid80 <- kan$grid80
    out$part80 <- length(below) + match(lo[kan$grid80], ext)
    out$open_parts <- if (is_open) length(below) + which(ext >= lo[na]) else integer()
  } else {
    out$width <- c(up - lo, if (method == "closed") 1)
    out$parts <- c(ages$labels, if (method == "closed") "survivors")
  }
  out
}

# Grid ages for the Kannisto extension: single-year closed age groups from 80
# (at least ten), the last bound at or below 120; NULL where it does not apply.
# `fit`: the age groups 80-100 to which the model is fitted.
coda_kannisto_ages <- function(ages) {
  na <- length(ages)
  lo <- ages$lower
  up <- ages$upper
  last <- if (is.infinite(up[na])) lo[na] else up[na]
  closed80 <- which(lo >= 80 & is.finite(up))
  if (!80 %in% lo || length(closed80) < 10 || last > 120) return(NULL)
  if (any(up[closed80] - lo[closed80] != 1) || any(lo[closed80] != round(lo[closed80]))) return(NULL)
  list(fit = closed80[lo[closed80] <= 100], grid80 = closed80)
}

# Rates of the composition [model age, year] of one stratum, and the log rate
# of the open interval for method "regression". Zero deaths are replaced
# multiplicatively (coda_replace_zeros()); with the Kannisto extension, the
# ages from 80 come from the Kannisto model.
coda_schedule <- function(D, E, layout, what) {
  D <- matrix(D, length(layout$ages))
  E <- matrix(E, length(layout$ages))
  if (layout$method == "kannisto") {
    below <- layout$below
    mb <- coda_replace_zeros(D[below, , drop = FALSE], what) / E[below, , drop = FALSE]
    x <- layout$ages$lower[layout$fit]
    kan <- coda_kannisto(D[layout$fit, , drop = FALSE], E[layout$fit, , drop = FALSE], x, what)
    return(list(m = rbind(mb, kan$m), open_log = NULL, kannisto = kan$par))
  }
  m <- coda_replace_zeros(D, what) / E
  list(m = m, open_log = if (layout$method == "regression") log(m[nrow(m), ]) else NULL)
}

# clr of the life-table deaths [part, year] of rate schedules M [model age,
# year], radix 1 at the first age ([life_table()] conventions; with the
# Kannisto extension all survivors die in the last single year, 120)
coda_clr <- function(M, layout) {
  M <- matrix(M, if (layout$method == "kannisto") length(layout$width) else length(layout$ages))
  if (layout$method == "kannisto") {
    nn <- matrix(layout$width, nrow(M), ncol(M))
    Q <- lt_qx(M, nn, nn / 2)
    Q[nrow(Q), ] <- 1
    L <- rbind(1, apply(1 - Q, 2, cumprod))[seq_len(nrow(M)), , drop = FALSE]
    dx <- L * Q
  } else {
    lt <- lt_columns(M, layout$ages)
    dx <- lt$dx
    if (layout$method == "closed") dx <- rbind(dx, lt$lx[nrow(M), ] - lt$dx[nrow(M), ])
  }
  z <- log(pmax(dx, 1e-12))
  sweep(z, 2, colMeans(z))
}

# Rates on the age axis [age, column] from clr life-table deaths Z [part,
# ...]: m = d / L for each part; with the Kannisto extension, the open
# interval gets the deaths over the person-years of its parts; with method
# "regression", `open_log` gives the log rate of the open interval
coda_layout_rates <- function(Z, layout, open_log = NULL) {
  ages <- layout$ages
  na <- length(ages)
  P <- exp(matrix(Z, length(layout$parts)))
  P <- sweep(P, 2, colSums(P), "/")
  l <- apply(P, 2, function(p) rev(cumsum(rev(p))))
  if (is.null(dim(l))) l <- matrix(l, nrow(P))
  L <- layout$width * (l - P / 2)
  out <- matrix(NA_real_, na, ncol(P))
  if (layout$method == "kannisto") {
    nb <- length(layout$below)
    out[layout$below, ] <- P[seq_len(nb), , drop = FALSE] / L[seq_len(nb), , drop = FALSE]
    out[layout$grid80, ] <- P[layout$part80, , drop = FALSE] / L[layout$part80, , drop = FALSE]
    if (length(layout$open_parts)) {
      op <- layout$open_parts
      out[na, ] <- colSums(P[op, , drop = FALSE]) / colSums(L[op, , drop = FALSE])
    }
  } else {
    m <- P / L
    out[] <- m[seq_len(na), , drop = FALSE]
    if (layout$method == "regression") {
      if (is.null(open_log)) mt_abort("The open interval needs its rates.", class = "model")
      out[na, ] <- exp(as.vector(open_log))
    }
  }
  dim(out) <- c(na, dim(Z)[-1])
  pmax(out, 1e-12)
}

# Zero deaths ------------------------------------------------------------------

# Multiplicative replacement of zero deaths [age, year] (Bergeron-Boucher et
# al. 2017, Appendix C): in a year with z zero cells, each gets the share
# delta of the year's deaths and the other shares are multiplied by
# 1 - z delta, so that the total is kept; delta is half the smallest positive
# count over all ages and years divided by the deaths of the year. Package
# safeguard: z delta is at most 1/2.
coda_replace_zeros <- function(D, what) {
  if (!length(D)) return(D)
  pos <- is.finite(D) & D > 0
  if (!any(pos)) mt_abort(sprintf("No deaths in %s in the training years.", what), class = "model_fit")
  if (all(pos)) return(D)
  share <- mean(!pos)
  if (share > 0.1) {
    mt_warn(sprintf(paste("%.0f%% of the cells of %s have no deaths; they were replaced multiplicatively by half the",
                          "smallest count. Wider age groups (mt_regroup()) or smoothing (step_smooth()) avoid this."),
                    100 * share, what), class = "sparse")
  }
  tot <- colSums(D)
  if (any(tot <= 0)) {
    mt_abort(sprintf(paste("No deaths in %s in %d training year(s): the composition of deaths needs deaths in every",
                           "year. Wider age groups (mt_regroup()) or fewer strata avoid this."), what, sum(tot <= 0)),
             class = "model_fit")
  }
  z <- colSums(!pos)
  delta <- pmin(min(D[pos]) / 2 / tot, 0.5 / pmax(z, 1))
  nr <- nrow(D)
  P <- D / rep(tot, each = nr)
  P <- ifelse(pos, P * rep(1 - z * delta, each = nr), rep(delta, each = nr))
  P * rep(tot, each = nr)
}

# Kannisto model ---------------------------------------------------------------

# Kannisto model (Thatcher, Kannisto and Vaupel 1998) of each year (column):
# hazard mu(x) = a exp(b (x - 80)) / (1 + a exp(b (x - 80))), fitted by Poisson
# maximum likelihood to the deaths D and exposures E of the single ages x
# (Bergeron-Boucher et al. 2017, Appendix C); rates at ages 80-120 [41, year].
# Years without a finite fit (no deaths) take the fit to all years together
# (package safeguard, with a warning).
coda_kannisto <- function(D, E, x, what = "the stratum") {
  par <- vapply(seq_len(ncol(D)), function(t) kannisto_ml(D[, t], E[, t], x - 80), numeric(2))
  bad <- !is.finite(colSums(par))
  if (any(bad)) {
    pooled <- kannisto_ml(rowSums(D), rowSums(E), x - 80)
    if (!all(is.finite(pooled))) {
      mt_abort(sprintf("The Kannisto model cannot be fitted to the deaths at ages 80 and over of %s; use open = %s.",
                       what, "\"regression\""), class = "model_fit")
    }
    par[, bad] <- pooled
    mt_warn(sprintf(paste("The Kannisto model could not be fitted to the deaths at ages 80 and over of %s in %d",
                          "year(s); the fit to all training years together is used for them."), what, sum(bad)),
            class = "sparse")
  }
  list(m = stats::plogis(outer(0:40, par[2, ]) + rep(par[1, ], each = 41)), par = par)
}

# Poisson maximum likelihood of the Kannisto model, mu = plogis(c + b t), by
# Fisher scoring with step halving; c(log a, b), or NA without a finite fit
kannisto_ml <- function(D, E, t) {
  ok <- is.finite(D) & is.finite(E) & E > 0
  D <- D[ok]
  E <- E[ok]
  t <- t[ok]
  if (sum(D) <= 0 || length(D) < 3) return(c(NA_real_, NA_real_))
  X <- cbind(1, t)
  ll <- function(th) {
    mu <- stats::plogis(drop(X %*% th))
    sum(D * log(mu) - E * mu)
  }
  y <- stats::qlogis(pmin(pmax((D + 0.5) / E, 1e-6), 0.99))
  th <- tryCatch(qr.coef(qr(X), y), error = function(e) c(-3, 0.1))
  if (anyNA(th)) th <- c(-3, 0.1)
  l0 <- ll(th)
  for (it in seq_len(200)) {
    mu <- stats::plogis(drop(X %*% th))
    g <- crossprod(X, (1 - mu) * (D - E * mu))
    H <- crossprod(X, X * (E * mu * (1 - mu)^2))
    step <- tryCatch(drop(solve(H, g)), error = function(e) NULL)
    if (is.null(step) || !all(is.finite(step))) return(c(NA_real_, NA_real_))
    s <- 1
    repeat {
      new <- th + s * step
      l1 <- ll(new)
      if ((is.finite(l1) && l1 >= l0 - 1e-10 * abs(l0)) || s < 1e-8) break
      s <- s / 2
    }
    th <- new
    if (!is.finite(l1)) return(c(NA_real_, NA_real_))
    done <- abs(l1 - l0) <= 1e-12 * (1 + abs(l1)) && max(abs(s * step)) < 1e-8
    l0 <- l1
    if (done) return(unname(th))
  }
  unname(th)
}

# CoDa Lee-Carter fit ----------------------------------------------------------

# CoDa model of clr life-table deaths Z [part, year] (Bergeron-Boucher et al.
# 2017, eq. 5-6; Oeppen 2008, Section 3): the centre `alpha` (clr of the
# closed geometric mean over the years), the first `rank` components of the
# singular value decomposition of the centred matrix (`b` [part, rank], `k`
# [year, rank]) and an ARIMA model of each index. `common`: clr of a common
# component [part, year] removed before the decomposition (coherent model).
# `adjust`: the first forecast step of an ARIMA(0, 1, q) index with drift is
# corrected to the drift (Bergeron-Boucher et al. 2017, eq. 20).
coda_core <- function(Z, rank = 1, order = c(0, 1, 1), drift = TRUE, common = NULL, adjust = TRUE) {
  alpha <- rowMeans(Z)
  H <- Z - alpha
  if (!is.null(common)) H <- H - common
  rank <- min(rank, nrow(H) - 1, ncol(H) - 1)
  sv <- coda_svd(H, rank)
  list(alpha = alpha, b = sv$b, k = sv$k, explained = sv$explained, rank = rank,
       ts = lapply(seq_len(rank), function(j) index_arima_fit(sv$k[, j], order, drift)),
       order = order, drift = drift, adjust = adjust, resid = H - sv$b %*% t(sv$k), Z_last = Z[, ncol(Z)])
}

# First `rank` components of the SVD of a centred matrix C [part, year]:
# b = right singular vectors (parts), k = singular value times the left
# singular vectors (years). Signs: aligned with `ref` (columns of b) where
# given, otherwise each index increases over the years (Bergeron-Boucher et
# al. 2017, p. 533). `explained`: share of the squared singular values of
# each component.
coda_svd <- function(C, rank = 1, ref = NULL) {
  s <- svd(t(C), nu = rank, nv = rank)
  b <- s$v[, seq_len(rank), drop = FALSE]
  k <- s$u[, seq_len(rank), drop = FALSE] * rep(s$d[seq_len(rank)], each = nrow(s$u))
  for (j in seq_len(rank)) {
    flip <- if (!is.null(ref)) sum(b[, j] * ref[, j]) < 0 else k[nrow(k), j] < k[1, j]
    if (isTRUE(flip)) {
      b[, j] <- -b[, j]
      k[, j] <- -k[, j]
    }
  }
  list(b = b, k = k, explained = s$d[seq_len(rank)]^2 / max(sum(s$d^2), .Machine$double.xmin))
}

# Least-squares loadings of the log rate of the open interval on the
# indexes (columns of X), the open-interval analogue of b[x] (package-specific)
coda_open_fit <- function(lo, X) {
  X <- cbind(1, X)
  cf <- qr.coef(qr(X), lo)
  cf[is.na(cf)] <- 0
  res <- lo - X %*% cf
  list(a = cf[[1]], b = unname(cf[-1]), last = lo[length(lo)], y = lo,
       sd = sqrt(sum(res^2) / max(length(lo) - ncol(X), 1)))
}

# Point forecasts of the indexes [h, rank]: means of the ARIMA models, the
# first step of an ARIMA(0, 1, q) index with drift corrected to the drift
# (eq. 20 of Bergeron-Boucher et al. 2017) where `adjust`
coda_index_mean <- function(core, h) {
  matrix(vapply(seq_len(core$rank), function(j) {
    f <- core$ts[[j]]
    kf <- as.vector(index_arima_paths(f, h, 1, innovations = "none"))
    if (isTRUE(core$adjust) && f$order[1] == 0 && f$order[2] == 1 && f$order[3] > 0 && f$drift) {
      kf <- kf + f$const - (kf[1] - f$k[length(f$k)])
    }
    kf
  }, numeric(h)), h, core$rank)
}

# Simulated paths of the components, sum_j b*_j k*_j(T + s) [part, h, draw]
# (Bergeron-Boucher et al. 2017, Appendix B): the residuals of the fit are
# resampled cell by cell and added to the fitted components, the components
# are estimated again, ARIMA models of the same order are fitted to the new
# indexes, and their future paths are simulated with resampled innovations;
# ceiling(sqrt(n_draws)) refits share the paths. `open_y`: log rate of the
# open interval (or its part not explained by the common component),
# regressed on the indexes of each refit [h, draw]. The order of the draws is
# random.
coda_boot <- function(core, h, n_draws, open_y = NULL) {
  r <- core$rank
  np <- nrow(core$b)
  fitted <- core$b %*% t(core$k)
  R <- core$resid
  out <- array(0, c(np, h, n_draws))
  open <- if (!is.null(open_y)) matrix(0, h, n_draws) else NULL
  n_rep <- ceiling(sqrt(n_draws))
  n_k <- ceiling(n_draws / n_rep)
  j0 <- 0L
  while (j0 < n_draws) {
    m <- min(n_k, n_draws - j0)
    idx <- j0 + seq_len(m)
    Hs <- fitted + matrix(R[sample.int(length(R), length(R), replace = TRUE)], nrow(R))
    sv <- coda_svd(Hs, r, ref = core$b)
    paths <- lapply(seq_len(r), function(j) {
      index_arima_paths(index_arima_fit(sv$k[, j], core$order, core$drift), h, m)
    })
    for (j in seq_len(r)) out[, , idx] <- out[, , idx, drop = FALSE] + outer(sv$b[, j], paths[[j]])
    if (!is.null(open_y)) {
      cf <- coda_open_fit(open_y, sv$k)
      open[, idx] <- cf$a
      for (j in seq_len(r)) open[, idx] <- open[, idx] + cf$b[j] * paths[[j]]
    }
    j0 <- j0 + m
  }
  perm <- sample.int(n_draws)
  list(Z = out[, , perm, drop = FALSE], open = if (!is.null(open)) open[, perm, drop = FALSE])
}

# Centring of simulated clr deaths Z [part, h, draw] (and log open rates
# [h, draw]) on a point forecast (clr `P1` and log open rate `o1` of the
# first forecast year): every draw is perturbed by the ratio of the point
# forecast to the per-part median of the closed draws of the first forecast
# year, and the log open rates are shifted by the difference to their median
# (the jump-off of the simulations in the code of Bergeron-Boucher et al.
# 2017; their eq. 19 shifts the bounds of e0)
coda_centre <- function(Z, P1, open = NULL, o1 = NULL) {
  C1 <- exp(matrix(Z[, 1, ], dim(Z)[1]))
  C1 <- sweep(C1, 2, colSums(C1), "/")
  med <- apply(C1, 1, stats::median)
  off <- (P1 - mean(P1)) - (log(med) - mean(log(med)))
  Z <- Z + off
  if (!is.null(open)) open <- open + (o1 - stats::median(open[1, ]))
  list(Z = Z, open = open)
}

# ARIMA models of the indexes ----------------------------------------------

# ARIMA(p, d, q) model of a period index by maximum likelihood: stats::arima()
# as called by forecast::Arima(method = "ML"), with a drift for d = 1
# (`drift`) and a mean for d = 0. Where the fit fails or the series is too
# short, the model without ARMA terms (random walk with or without drift,
# mean, or random walk of the changes for d = 2; `fallback = TRUE`).
index_arima_fit <- function(k, order, drift = TRUE) {
  n <- length(k)
  d <- order[2]
  use_drift <- d == 1 && isTRUE(drift)
  npar <- order[1] + order[3] + use_drift + (d == 0)
  w <- if (d > 0) diff(k, differences = d) else k
  fit <- NULL
  if (sum(order[-2]) > 0 && length(w) > npar + 2 && stats::sd(w) > 1e-10 * (1 + mean(abs(w)))) {
    xreg <- if (use_drift) cbind(drift = seq_len(n)) else NULL
    for (method in c("ML", "CSS-ML")) {
      fit <- tryCatch(suppressWarnings(stats::arima(k, order = order, xreg = xreg, include.mean = d == 0,
                                                    method = method)), error = function(e) NULL)
      if (!is.null(fit)) break
    }
  }
  if (is.null(fit)) {
    const <- if (use_drift || d == 0) mean(w) else 0
    u <- w - const
    return(list(order = c(0, d, 0), drift = use_drift, phi = numeric(), theta = numeric(), const = const, k = k,
                innov = u, pool = u - mean(u), sigma = sqrt(sum(u^2) / max(length(u) - (const != 0), 1)),
                coef = c(if (use_drift) c(drift = const), if (d == 0) c(intercept = const)),
                fallback = sum(order[-2]) > 0))
  }
  cf <- stats::coef(fit)
  p <- order[1]
  q <- order[3]
  res <- as.numeric(stats::residuals(fit))
  list(order = order, drift = use_drift,
       phi = if (p) unname(cf[paste0("ar", seq_len(p))]) else numeric(),
       theta = if (q) unname(cf[paste0("ma", seq_len(q))]) else numeric(),
       const = if (use_drift) unname(cf[["drift"]]) else if (d == 0) unname(cf[["intercept"]]) else 0,
       k = k, innov = res[(d + 1):n], pool = res - mean(res), sigma = sqrt(sum(res^2) / max(n - d - length(cf), 1)),
       coef = cf, fallback = FALSE)
}

# Paths [h, n] of an index from its ARIMA model, conditional on the observed
# series: recursion of the ARMA model of the differenced series around its
# constant, with innovations resampled from the centred residuals (as
# simulate(future = TRUE, bootstrap = TRUE) of the forecast package) or set
# to zero (`"none"`: the mean forecast given the past innovations)
index_arima_paths <- function(f, h, n, innovations = c("bootstrap", "none")) {
  innovations <- match.arg(innovations)
  d <- f$order[2]
  w <- if (d > 0) diff(f$k, differences = d) else f$k
  u <- w - f$const
  p <- length(f$phi)
  q <- length(f$theta)
  U <- matrix(rev(utils::tail(u, p)), p, n)
  Ep <- matrix(rev(utils::tail(f$innov, q)), q, n)
  lastk <- f$k[length(f$k)]
  lastd <- if (d == 2) diff(utils::tail(f$k, 2)) else 0
  kk <- rep(lastk, n)
  dk <- rep(lastd, n)
  out <- matrix(NA_real_, h, n)
  for (s in seq_len(h)) {
    e <- if (innovations == "bootstrap") f$pool[sample.int(length(f$pool), n, replace = TRUE)] else numeric(n)
    un <- e
    if (p) un <- un + colSums(f$phi * U)
    if (q) un <- un + colSums(f$theta * Ep)
    wn <- un + f$const
    if (d == 0) {
      kk <- wn
    } else if (d == 1) {
      kk <- kk + wn
    } else {
      dk <- dk + wn
      kk <- kk + dk
    }
    out[s, ] <- kk
    if (p) U <- rbind(un, U[-p, , drop = FALSE])
    if (q) Ep <- rbind(e, Ep[-q, , drop = FALSE])
  }
  out
}

# Checks of an ARIMA order given to a constructor
check_arima_order <- function(order, name) {
  if (!is.numeric(order) || length(order) != 3 || anyNA(order) || any(order < 0) || any(order != round(order)) ||
      order[2] > 2 || any(order[-2] > 3)) {
    mt_abort(sprintf("`%s` must be an ARIMA order c(p, d, q) of whole numbers with p, q <= 3 and d <= 2.", name),
             class = "argument", call = NULL)
  }
  as.numeric(order)
}

# Label of an ARIMA model of an index
arima_label <- function(order, drift) {
  sprintf("ARIMA(%s)%s", paste(order, collapse = ","), if (order[2] == 1 && isTRUE(drift)) " with drift" else "")
}
