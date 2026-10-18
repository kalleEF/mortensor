# Score-optimal probabilistic forecast reconciliation
#
# Panagiotelis A, Gamakumara P, Athanasopoulos G, Hyndman RJ (2023).
# Probabilistic forecast reconciliation: properties, evaluation and score
# optimisation. European Journal of Operational Research 306(2): 693-706.
#
# Reconciled draws are S (d + G y) of base draws y; d and G minimise the
# total energy score of past forecasts. The optimiser is adapted from
# ProbReco 0.1.2 (scoreopt(), scoreopt.control() and total_score(); GPL-3;
# author Anastasios Panagiotelis): the parameter layout c(d, vec(G)), the
# Adam updates with their defaults, the start at least squares (OLS)
# reconciliation and the stopping rule. ProbReco differentiates the Monte
# Carlo energy score with Stan Math; here its gradient is written out
# analytically, so no compiled code is needed. The training forecasts come
# from inner rolling origins within the training data (as in recalibrated()).
# Differences from ProbReco (documented in reconciled()): the two samples of
# each step are drawn from disjoint sets of a fixed set of inner draws
# (scoreopt_pairs()), and outputs and inputs are rescaled before Adam.

# Settings of the training (`control` of reconciled()): draws per inner
# forecast, draws per optimisation step (Q of ProbReco), and the Adam
# settings of ProbReco::scoreopt.control() (maxIter is `max_iter`)
scoreopt_control <- function(control = list()) {
  def <- list(n_draws = 100, draws_per_step = 50, max_iter = 500, eta = 0.001, beta1 = 0.9, beta2 = 0.999,
              tol = 1e-4, epsilon = 1e-8)
  if (!is.list(control) || (length(control) && (is.null(names(control)) || any(!nzchar(names(control)))))) {
    mt_abort("`control` must be a named list, for example list(max_iter = 200).", class = "argument")
  }
  bad <- setdiff(names(control), names(def))
  if (length(bad)) {
    mt_abort(sprintf("Unknown element(s) of `control`: %s; known are %s.", paste(bad, collapse = ", "),
                     paste(names(def), collapse = ", ")), class = "argument")
  }
  ctl <- utils::modifyList(def, control)
  check_number(ctl$n_draws, "control$n_draws", lower = 10, integerish = TRUE)
  check_number(ctl$draws_per_step, "control$draws_per_step", lower = 2, integerish = TRUE)
  check_number(ctl$max_iter, "control$max_iter", lower = 0, integerish = TRUE)
  check_number(ctl$eta, "control$eta", lower = 0)
  check_number(ctl$beta1, "control$beta1", lower = 0, upper = 0.9999)
  check_number(ctl$beta2, "control$beta2", lower = 0, upper = 0.9999)
  check_number(ctl$tol, "control$tol", lower = 0)
  check_number(ctl$epsilon, "control$epsilon", lower = 1e-300)
  if (ctl$eta > 0 && ctl$max_iter > 0 && ctl$tol >= ctl$eta) {
    mt_warn(paste("`control$tol` is not below `control$eta`: the steps of Adam are at most about `eta`, so the",
                  "optimisation stops after the first step. Give a smaller `tol` (0 runs all `max_iter` steps)."),
            class = "argument")
  }
  ctl
}

# Training of d and G for each age from inner rolling-origin forecasts of
# the training data. `data` holds every node of the hierarchy, `fit` is the
# base fit on all training data and `hs` its hierarchy structure. Returns d
# [bottom, age], G [bottom, node, age] and a table per age
scoreopt_train <- function(p, data, prep, fit, hs, cores = 1) {
  ctl <- scoreopt_control(p$control)
  g <- mt_grid(data)
  last <- max(g$period$lower)
  first_ok <- min(g$period$lower) + p$base$min_periods - 1
  origins <- calibration_origins(last, p$years, p$n_origins)
  if (any(origins < first_ok)) {
    origins <- origins[origins >= first_ok]
    if (!length(origins)) {
      mt_abort(sprintf(paste("Score-optimal reconciliation on the last %d years needs at least %d training years;",
                             "give a shorter `years`."), p$years, p$base$min_periods + p$years), class = "model_data")
    }
    mt_warn(sprintf("Only %d training origin(s) of the reconciliation leave enough training years.", length(origins)),
            class = "reconcile")
  }
  additive <- p$hierarchy$additive
  out <- with_seed(p$seed, {
    pieces <- par_lapply(origins, function(o) {
      tryCatch(scoreopt_origin(p, data, o, last, prep, ctl$n_draws, hs, fit$strata, additive),
               error = function(e) e)
    }, cores, streams = TRUE)
    bad <- vapply(pieces, inherits, logical(1), what = "error")
    if (any(bad)) {
      mt_warn(sprintf("Training origin(s) %s of the reconciliation failed: %s", paste(origins[bad], collapse = ", "),
                      conditionMessage(pieces[[which(bad)[1]]])), class = "reconcile")
    }
    pieces <- Filter(Negate(is.null), pieces[!bad])
    if (!length(pieces)) mt_abort("Every training origin of the reconciliation failed.", class = "reconcile")
    ages <- seq_along(fit$age)
    res <- par_lapply(ages, function(a) scoreopt_age(pieces, a, hs, ctl), cores, streams = TRUE)
    list(res = res, origins = vapply(pieces, `[[`, 0, "origin"))
  })
  res <- out$res
  m <- length(hs$bottom)
  n <- nrow(hs$S)
  na <- length(res)
  G <- array(vapply(res, `[[`, numeric(m * n), "G"), c(m, n, na),
             dimnames = list(bottom = fit$strata$label[hs$bottom], node = fit$strata$label, age = fit$age$labels))
  d <- matrix(vapply(res, `[[`, numeric(m), "d"), m, na, dimnames = dimnames(G)[c(1, 3)])
  tab <- data.table(age = fit$age$labels, forecasts = vapply(res, `[[`, 0, "periods"),
                    scale = vapply(res, `[[`, 0, "scale"), iterations = vapply(res, `[[`, 0, "iterations"),
                    energy_ols = vapply(res, `[[`, 0, "es_start"), energy_opt = vapply(res, `[[`, 0, "es_end"))
  untrained <- tab$forecasts < 2
  if (any(untrained)) {
    mt_warn(sprintf(paste("Age(s) %s have fewer than 2 training forecasts with observed values for every stratum;",
                          "they are reconciled by OLS."), paste(utils::head(tab$age[untrained], 6), collapse = ", ")),
            class = "reconcile")
  }
  list(d = d, G = G, table = tab, origins = out$origins, control = ctl)
}

# One inner origin: base draws (reordered as in the forecast with
# dependence = "empirical") and observed values of every node for the next `horizon` years,
# on the deaths scale (additive = "deaths") or the rate scale. `noise`:
# observation noise of the bottom strata per draw, as the forecast adds it
scoreopt_origin <- function(p, data, origin, last, prep, n_draws, hs, strata, additive) {
  train <- mt_window(data, to = origin)
  fit <- mt_fit(p$base, train, prep = prep)
  if (!identical(fit$strata$label, strata$label)) {
    mt_abort("The inner fit has other strata than the final fit.", class = "reconcile")
  }
  last_fit <- max(fit$periods)
  to <- min(last_fit + p$horizon, last)
  if (to <= last_fit) return(NULL)
  test <- mt_window(data, from = last_fit + 1, to = to)
  test <- mt_set_grid(test, age = fit$age)
  E <- test_exposures(test, fit)
  fc <- suppressWarnings(mt_forecast(fit, to = to, n_draws = n_draws, exposure = E))
  arr <- fc$rates
  dn <- dimnames(arr)
  # exposures and observed deaths in the order of the forecast's periods and strata
  E <- E[, match(dn[[2]], dimnames(E)[[2]]), match(strata$label, dimnames(E)[[3]]), drop = FALSE]
  component <- fc$component
  if (identical(p$dependence, "empirical")) {
    cr <- copula_reorder(arr, fit$raw_D %||% fit$train_D, fit$train_E)
    arr <- cr$arr
    component <- permute_component(component, cr$perm)
  }
  b <- hs$bottom
  fc$component <- component
  noisy <- add_observation_noise(arr, fit, fc, b, E)
  # observed values; deaths filled from totals (mt_fill_suppressed()) are not observations
  t2 <- test
  t2$obs$deaths <- t2$obs$deaths[is.na(provenance) | provenance != "estimated"]
  D <- mt_array(t2, "deaths", coarse = "na")
  D <- D[, match(dn[[2]], dimnames(D)[[2]]), match(strata$label, dimnames(D)[[3]]), drop = FALSE]
  by_deaths <- identical(additive, "deaths")
  noise <- (noisy - arr)[, , b, , drop = FALSE]
  if (by_deaths) {
    X <- arr * as.vector(E)
    noise <- noise * as.vector(E[, , b, drop = FALSE])
    Y <- D
  } else {
    X <- arr
    Y <- D / E
  }
  if (!any(noise != 0, na.rm = TRUE)) noise <- NULL
  list(origin = origin, Y = Y, X = X, noise = noise)
}

# Training forecasts of one age (every origin and horizon with observed
# values and draws of all nodes) and the optimisation of d and G
scoreopt_age <- function(pieces, a, hs, ctl) {
  S <- hs$S
  n <- nrow(S)
  m <- ncol(S)
  Y <- list()
  X <- list()
  Z <- list()
  for (pc in pieces) {
    for (t in seq_len(dim(pc$Y)[2])) {
      y <- pc$Y[a, t, ]
      x <- matrix(pc$X[a, t, , ], n)
      if (anyNA(y) || !all(is.finite(x))) next
      Y[[length(Y) + 1L]] <- y
      X[[length(X) + 1L]] <- x
      Z[[length(Z) + 1L]] <- if (!is.null(pc$noise)) matrix(pc$noise[a, t, , ], m) else matrix(0, m, ncol(x))
    }
  }
  G0 <- solve(crossprod(S), t(S))
  nt <- length(Y)
  if (nt < 2) {
    return(list(d = rep(0, m), G = as.vector(G0), periods = nt, scale = NA_real_, iterations = 0L,
                es_start = NA_real_, es_end = NA_real_))
  }
  nd <- vapply(X, ncol, 1L)
  off <- cumsum(c(0L, nd))[seq_len(nt)]
  Y <- do.call(cbind, Y)
  X <- do.call(cbind, X)
  Z <- do.call(cbind, Z)
  # scales: of the output (the mean observed bottom value) and of each node
  # as an input (its mean observed value). With b = sc (d' + G' (x / sx)),
  # d = sc d' and G = sc G' diag(1 / sx), the energy score is the same up
  # to the factor sc, so the optimum does not change; the steps of Adam then
  # move d and every column of G by comparable amounts at all ages and for
  # nodes of different sizes (totals hold many bottom series)
  sc <- mean(abs(Y[hs$bottom, ]))
  if (!is.finite(sc) || sc <= 0) sc <- 1
  sx <- rowMeans(abs(Y))
  sx[!is.finite(sx) | sx <= 0] <- sc
  Y <- Y / sc
  X <- X / sx
  Z <- if (any(Z != 0)) Z / sc
  start <- c(rep(0, m), as.vector(sweep(G0, 2, sx, "*") / sc))
  # columns of the base draws: two samples of `draws_per_step` draws of each
  # training forecast per step, from disjoint sets of its draws
  per <- rep(seq_len(nt), each = ctl$draws_per_step)
  pick <- function(i1, i2, per) {
    list(x1 = X[, i1, drop = FALSE], x2 = X[, i2, drop = FALSE], e1 = if (!is.null(Z)) Z[, i1, drop = FALSE],
         e2 = if (!is.null(Z)) Z[, i2, drop = FALSE], per = per)
  }
  draw <- function() {
    ij <- lapply(nd, scoreopt_pairs, q = ctl$draws_per_step)
    pick(off[per] + unlist(lapply(ij, `[[`, 1L)), off[per] + unlist(lapply(ij, `[[`, 2L)), per)
  }
  opt <- scoreopt_sgd(Y, draw, S, start, ctl)
  # in-sample energy score per training forecast at the start (OLS) and at
  # the end, with every draw paired with the next one (fixed, not sampled)
  all_per <- rep(seq_len(nt), nd)
  k <- seq_len(sum(nd)) - off[all_per] - 1L
  fixed <- pick(seq_len(sum(nd)), off[all_per] + (k + 1L) %% nd[all_per] + 1L, all_per)
  es <- function(par) sc * scoreopt_score(Y, fixed, S, par, gradient = FALSE)$value / nt
  list(d = sc * opt$d, G = as.vector(sc * sweep(opt$G, 2, sx, "/")), periods = nt, scale = sc,
       iterations = opt$iterations, es_start = es(start), es_end = es(c(opt$d, as.vector(opt$G))))
}

# Two samples of q of the n draws of a training forecast for one step of the
# optimisation, from disjoint sets of draws, so that no draw is paired with
# itself (ProbReco simulates the two samples independently): without
# replacement when n >= 2 q, otherwise with replacement from two disjoint
# random halves of the draws
scoreopt_pairs <- function(n, q) {
  if (n >= 2 * q) {
    s <- sample.int(n, 2 * q)
    return(list(s[seq_len(q)], s[q + seq_len(q)]))
  }
  p <- sample.int(n)
  k <- n %/% 2
  list(p[sample.int(k, q, replace = TRUE)], p[k + sample.int(k, q, replace = TRUE)])
}

# Adam on the Monte Carlo total energy score, as ProbReco::scoreopt().
# `draw()` returns the base draws of one step: list(x1, x2 = two
# sets of draws [node, column] (scoreopt_pairs()), e1, e2 = observation noise of the
# bottom strata [bottom, column] or NULL, per = training forecast of each
# column)
scoreopt_sgd <- function(Y, draw, S, Gvec, control) {
  n <- nrow(S)
  m <- ncol(S)
  mom <- 0
  v <- 0
  i <- 1L
  dif <- Inf
  val <- NA_real_
  pre <- scoreopt_pre(Y, S)
  while (i <= control$max_iter && any(dif > control$tol)) {
    sv <- scoreopt_score(Y, draw(), S, Gvec, pre = pre)
    mom <- control$beta1 * mom + (1 - control$beta1) * sv$grad
    v <- control$beta2 * v + (1 - control$beta2) * sv$grad^2
    step <- control$eta * (mom / (1 - control$beta1^i)) / (sqrt(v / (1 - control$beta2^i)) + control$epsilon)
    Gvec <- Gvec - step
    dif <- abs(step)
    val <- sv$value
    i <- i + 1L
  }
  list(d = Gvec[seq_len(m)], G = matrix(Gvec[-seq_len(m)], m, n), value = val, iterations = i - 1L)
}

# Monte Carlo estimate of the total energy score of S (d + G x + e) and its
# gradient with respect to c(d, vec(G)) (ProbReco's total_score() with
# alpha = 1), from the draws `z` of draw():
#   sum_t 1/Q sum_q ( ||y_t - S (d + G x_tq + e_tq)|| - 1/2 ||S (G (x_tq - x*_tq) + e_tq - e*_tq)|| )
# The norms are computed with S'S and S'y (`pre`), which saves most of the
# products with S. Pairs at distance zero contribute nothing to the gradient
scoreopt_score <- function(Y, z, S, Gvec, gradient = TRUE, pre = scoreopt_pre(Y, S)) {
  n <- nrow(S)
  m <- ncol(S)
  k <- ncol(z$x1)
  w <- 1 / tabulate(z$per)[z$per]
  d <- Gvec[seq_len(m)]
  G <- matrix(Gvec[-seq_len(m)], m, n)
  XD <- cbind(z$x1, z$x1 - z$x2)
  GX <- G %*% XD
  B <- GX[, seq_len(k), drop = FALSE] + d
  V <- GX[, k + seq_len(k), drop = FALSE]
  if (!is.null(z$e1)) {
    B <- B + z$e1
    V <- V + (z$e1 - z$e2)
  }
  SB <- pre$StS %*% cbind(B, V)
  SV <- SB[, k + seq_len(k), drop = FALSE]
  SB <- SB[, seq_len(k), drop = FALSE]
  SY <- pre$StY[, z$per, drop = FALSE]
  nr <- sqrt(pmax(pre$yy[z$per] - 2 * colSums(SY * B) + colSums(B * SB), 0))
  nu <- sqrt(pmax(colSums(V * SV), 0))
  value <- sum(w * nr) - 0.5 * sum(w * nu)
  if (!gradient) return(list(value = value))
  A <- cbind((SY - SB) * rep(ifelse(nr > 0, w / nr, 0), each = m),
             SV * rep(ifelse(nu > 0, 0.5 * w / nu, 0), each = m))
  list(value = value, grad = c(-rowSums(A[, seq_len(k), drop = FALSE]), -as.vector(tcrossprod(A, XD))))
}

scoreopt_pre <- function(Y, S) list(StS = crossprod(S), StY = crossprod(S, Y), yy = colSums(Y^2))

# Forecast: y_tilde = S (d + G y) per age, year and draw on the scale of
# reconcile_array(); bottom values below zero are set to zero and the
# aggregates are their sums. Returns what reconcile_array() returns
scoreopt_apply <- function(arr, hs, additive, so, exposure, d) {
  S <- hs$S
  b <- hs$bottom
  agg_rows <- which(!hs$is_bottom)
  E <- expand_exposure(exposure, d, dimnames(arr), hs)
  if (additive == "deaths") {
    Sx <- hs$Sx %||% S
    for (t in seq_len(d[2])) E[, t, ] <- matrix(E[, t, b], d[1]) %*% t(Sx)
    Y <- arr * as.vector(E)
  } else {
    Y <- arr
  }
  if (dim(so$G)[3] != d[1]) mt_abort("The reconciliation was trained for other ages.", class = "reconcile")
  out <- Y
  gap_before <- 0
  zero_share <- 0
  draw_zero <- logical(d[4])
  for (a in seq_len(d[1])) {
    G <- matrix(so$G[, , a], length(b))
    dv <- so$d[, a]
    for (t in seq_len(d[2])) {
      Yc <- matrix(Y[a, t, , ], nrow = d[3])
      if (length(agg_rows)) {
        tot_b <- S[agg_rows, , drop = FALSE] %*% Yc[b, , drop = FALSE]
        gap_before <- max(gap_before, max(abs(Yc[agg_rows, , drop = FALSE] - tot_b) / pmax(abs(tot_b), 1e-12)))
      }
      Yb <- G %*% Yc + dv
      neg <- !is.na(Yb) & Yb < 0
      zero_share <- zero_share + sum(neg)
      draw_zero <- draw_zero | colSums(neg) > 0
      Yb[neg] <- 0
      out[a, t, , ] <- S %*% Yb
    }
  }
  zero_share <- zero_share / (d[1] * d[2] * d[4] * length(b))
  if (mean(draw_zero) > 0.05) {
    mt_warn(sprintf(paste("In %.0f%% of the draws some bottom forecasts were negative after score-optimal",
                          "reconciliation and were set to zero; their distributions are truncated. Check the",
                          "training of the reconciliation (fit$reconcile$scoreopt$table) or use another method."),
                    100 * mean(draw_zero)), class = "reconcile")
  }
  rates <- if (additive == "deaths") out / as.vector(E) else out
  rates <- pmax(rates, 1e-12)
  list(rates = rates, exposure = E,
       info = list(method = "energy_score", max_relative_gap_before = gap_before, share_set_to_zero = zero_share,
                   mean_shrinkage = NA_real_, share_draws_with_zeros = mean(draw_zero),
                   max_relative_change_fixed = NA_real_))
}
