# Multi-cause age-period-cohort model: native sampler behind apc_causes().
#
# Rewrite of bamp_multicause() and predict_multicause() of the development
# branch "coherent_forecasting" of bamp (version 2.2.0 in its DESCRIPTION,
# commit 4363a78 of 2026-07-13; file R/coherent_cause.R), GPL-3. bamp:
# Volker Schmid (author, maintainer), Florian Geressen, Leonhard Held, Evi
# Rainer and Chris Kypridemos (contributors); the multi-cause model is by
# Chris Kypridemos. When a bamp release on CRAN contains bamp_multicause(),
# apc_causes() is to call it and this sampler is to be deprecated
# (dev/coherence-plan.md).
#
# Model of each group of strata (the causes of one sex, say): the deaths of
# the causes split the all-cause deaths of a cell by stick-breaking, a
# sequence of binomial models of the conditional shares (Linderman, Johnson
# and Adams 2015 for the representation): the deaths of cause c are
# binomial out of the deaths R[c] not yet assigned (R[1]: all deaths; R[c]
# is R[c - 1] minus the deaths of cause c - 1), with
# logit pi[c, i, j] = mu[c] + theta[c, i] + phi[c, j] + psi[c, k].
# The period effects of the shares share a matrix-normal random-walk prior
# with precision Omega (x) K_period, Omega Wishart; the cohort effects the
# same with their own Omega. The all-cause rates follow a single
# binomial-logit APC model with the same priors. Cause rates are shares
# times all-cause rates, so they add up to the all-cause rate in every draw.
#
# Changes from the branch:
# * The sampler of apc_bayes(coherent = TRUE) (model-apc-coherent.R): exact
#   Metropolis-Hastings with a Gaussian approximation of all effects in
#   place of Polya-Gamma Gibbs steps with normal approximations.
# * One zero-slope constraint on the period effect of each share, and only
#   when age, period and cohort are all second-order random walks; the
#   branch pins the period and the cohort slope when both are RW2, which
#   removes an identified direction.
# * The all-cause model is fitted with the same native sampler (the branch
#   calls bamp()), its period and cohort precisions with the Wishart priors
#   of dimension 1 (Gamma priors).
# * The latent-factor coupling (`factor`) of the branch is not taken over.

mc_hyper <- function() list(age = c(1, 0.5), omega = c(2, 1e-4), omega_c = c(2, 1e-4))

# The model of C fields (shares, or the all-cause total) of one group. Y, N:
# arrays [age group, period, field]
mc_setup <- function(Y, N, ppa, ord, hyper = mc_hyper()) {
  I <- dim(Y)[1]
  J <- dim(Y)[2]
  C <- dim(Y)[3]
  M <- as.integer(ppa)
  K <- M * (I - 1L) + J
  len <- c(mu = C, th = C * I, ph = C * J, ps = C * K)
  end <- cumsum(len)
  pos <- lapply(names(len), function(nm) end[[nm]] - len[[nm]] + seq_len(len[[nm]]))
  names(pos) <- names(len)
  P <- end[["ps"]]
  ok <- which(is.finite(N) & N > 0 & is.finite(Y))
  if (any(Y[ok] > N[ok] + 1e-8)) mt_abort("Deaths exceed their denominators in some cells.", class = "model_data")
  ai <- arrayInd(ok, dim(Y))
  ci <- ai[, 1]
  cj <- ai[, 2]
  cf <- ai[, 3]
  ck <- M * (I - ci) + cj
  cols <- cbind(pos$mu[cf], pos$th[(cf - 1L) * I + ci], pos$ph[(cf - 1L) * J + cj], pos$ps[(cf - 1L) * K + ck])
  X <- Matrix::sparseMatrix(i = rep(seq_along(ok), ncol(cols)), j = as.vector(cols), x = 1, dims = c(length(ok), P))
  row <- function(idx, val = 1) {
    r <- numeric(P)
    r[idx] <- val
    r
  }
  rows <- unlist(lapply(seq_len(C), function(f) {
    c(list(row(pos$th[(f - 1L) * I + seq_len(I)]), row(pos$ph[(f - 1L) * J + seq_len(J)]),
           row(pos$ps[(f - 1L) * K + seq_len(K)])),
      # the unidentified linear trend when all three are RW2: pin the period slope
      if (all(ord == 2L)) list(row(pos$ph[(f - 1L) * J + seq_len(J)], seq_len(J) - (J + 1) / 2)))
  }), recursive = FALSE)
  A <- do.call(rbind, rows)
  st <- list(I = I, J = J, C = C, M = M, K = K, P = P, X = X, tX = Matrix::t(X), y = Y[ok], N = N[ok], cells = ok,
             ord = ord, idx = pos, Ka = coh_rw(I, ord[1]), Kp = coh_rw(J, ord[2]), Kc = coh_rw(K, ord[3]), A = A,
             tA = t(A), hyper = hyper, ridge = 1e-6, dims = dim(Y))
  zero <- function() Matrix::sparseMatrix(i = integer(), j = integer(), x = numeric(), dims = c(P, P))
  place <- function(idx, B) {
    Z <- zero()
    Z[idx, idx] <- B
    Z
  }
  basis <- list(fixed = Matrix::Diagonal(P, st$ridge) + Matrix::crossprod(Matrix::Matrix(A, sparse = TRUE)))
  for (f in seq_len(C)) basis[[paste0("kappa", f)]] <- place(pos$th[(f - 1L) * I + seq_len(I)], st$Ka$K)
  pairs <- which(upper.tri(diag(C), diag = TRUE), arr.ind = TRUE)
  for (r in seq_len(nrow(pairs))) {
    E <- matrix(0, C, C)
    E[pairs[r, 1], pairs[r, 2]] <- E[pairs[r, 2], pairs[r, 1]] <- 1
    basis[[paste0("om", r)]] <- place(pos$ph, kronecker(Matrix::Matrix(E, sparse = TRUE), st$Kp$K))
    basis[[paste0("oc", r)]] <- place(pos$ps, kronecker(Matrix::Matrix(E, sparse = TRUE), st$Kc$K))
  }
  st$pairs <- pairs
  gmrf_pattern(st, cols, basis)
}

# Prior precision for hyperparameters h (kappa: C precisions of the age
# effects; Om, Oc: C x C precisions of the period and cohort effects), plus
# A'A (see coh_prior())
mc_prior <- function(st, h) {
  b <- st$basis
  v <- b$fixed
  for (f in seq_len(st$C)) v <- v + h$kappa[f] * b[[paste0("kappa", f)]]
  for (r in seq_len(nrow(st$pairs))) {
    v <- v + h$Om[st$pairs[r, 1], st$pairs[r, 2]] * b[[paste0("om", r)]] +
      h$Oc[st$pairs[r, 1], st$pairs[r, 2]] * b[[paste0("oc", r)]]
  }
  Q <- st$pattern
  Q@x <- v
  Q
}

# Starting values: the overall logit and the empirical logits by age of
# each field (centred), then the mode
mc_start <- function(st, Qp, sym) {
  x <- numeric(st$P)
  Y <- array(0, st$dims)
  N <- array(0, st$dims)
  Y[st$cells] <- st$y
  N[st$cells] <- st$N
  for (f in seq_len(st$C)) {
    mu <- stats::qlogis(min(max((sum(Y[, , f]) + 0.5) / (sum(N[, , f]) + 1), 1e-10), 1 - 1e-10))
    x[st$idx$mu[f]] <- mu
    r <- (rowSums(Y[, , f, drop = FALSE]) + 0.5) / (rowSums(N[, , f, drop = FALSE]) + 1)
    th <- stats::qlogis(pmin(pmax(r, 1e-10), 1 - 1e-10)) - mu
    x[st$idx$th[(f - 1L) * st$I + seq_len(st$I)]] <- th - mean(th)
  }
  gmrf_mode(st, x, Qp, sym)
}

# Mode of the conditional posterior of the effects for fixed
# hyperparameters: Newton steps under the constraints, with step halving
gmrf_mode <- function(st, x, Qp, sym) {
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

# Wishart draw of the precision Omega of matrix-normal effects Z [time,
# field] with structure matrix K: prior Wishart(C + df_extra, (v0 I)^-1)
mc_omega <- function(Z, K, rank, hyp) {
  C <- ncol(Z)
  S <- solve(diag(hyp[2], C) + as.matrix(Matrix::crossprod(Z, K %*% Z)))
  matrix(stats::rWishart(1, C + hyp[1] + rank, (S + t(S)) / 2)[, , 1], C, C)
}

# One chain of the multi-cause sampler
mc_chain <- function(st, n_iter, burn_in, thin) {
  hp <- st$hyper
  C <- st$C
  h <- list(kappa = rep(1, C), Om = diag(C), Oc = diag(C))
  Qp <- mc_prior(st, h)
  H0 <- Qp
  H0@x <- Qp@x + as.vector(st$tmap %*% (st$N * 0.01 + 1e-3))
  sym <- Matrix::Cholesky(H0, perm = TRUE, LDL = FALSE, super = FALSE)
  x <- mc_start(st, Qp, sym)
  keep <- seq.int(burn_in + thin, n_iter, by = thin)
  nk <- length(keep)
  out <- list(x = matrix(NA_real_, nk, st$P), kappa = matrix(NA_real_, nk, C), Om = array(NA_real_, c(nk, C, C)),
              Oc = array(NA_real_, c(nk, C, C)))
  acc <- 0L
  store <- 0L
  for (it in seq_len(n_iter)) {
    cur <- coh_approx(st, x, Qp, sym)
    z <- as.vector(Matrix::solve(cur$L, Matrix::solve(cur$L, stats::rnorm(st$P), system = "Lt"), system = "Pt"))
    xs <- cur$m + cur$correct(z)
    pro <- coh_approx(st, xs, Qp, sym)
    la <- pro$lp - cur$lp + coh_logq(pro, x) - coh_logq(cur, xs)
    if (is.finite(la) && log(stats::runif(1)) < la) {
      x <- xs
      acc <- acc + 1L
    }
    for (f in seq_len(C)) {
      th <- x[st$idx$th[(f - 1L) * st$I + seq_len(st$I)]]
      h$kappa[f] <- stats::rgamma(1, hp$age[1] + st$Ka$rank / 2, hp$age[2] + sum(th * as.vector(st$Ka$K %*% th)) / 2)
    }
    h$Om <- mc_omega(matrix(x[st$idx$ph], st$J, C), st$Kp$K, st$Kp$rank, hp$omega)
    h$Oc <- mc_omega(matrix(x[st$idx$ps], st$K, C), st$Kc$K, st$Kc$rank, hp$omega_c)
    Qp <- mc_prior(st, h)
    if (it %in% keep) {
      store <- store + 1L
      out$x[store, ] <- x
      out$kappa[store, ] <- h$kappa
      out$Om[store, , ] <- h$Om
      out$Oc[store, , ] <- h$Oc
    }
  }
  out$accept <- acc / n_iter
  out
}

# Monitored quantities: log precisions of the age effects, log diagonals of
# Omega, the last period effect and the mean linear predictor of the last
# year of each field
mc_monitor <- function(st, ch) {
  C <- st$C
  J <- st$J
  x <- ch$x
  lp <- vapply(seq_len(C), function(f) {
    x[, st$idx$mu[f]] + x[, st$idx$ph[(f - 1L) * J + J]] +
      rowMeans(x[, st$idx$th[(f - 1L) * st$I + seq_len(st$I)], drop = FALSE])
  }, numeric(nrow(x)))
  om <- vapply(seq_len(C), function(f) log(ch$Om[, f, f]), numeric(nrow(x)))
  oc <- vapply(seq_len(C), function(f) log(ch$Oc[, f, f]), numeric(nrow(x)))
  out <- cbind(log(ch$kappa), matrix(om, nrow(x)), matrix(oc, nrow(x)), matrix(lp, nrow(x)))
  colnames(out) <- c(paste0("kappa_", seq_len(C)), paste0("omega_period_", seq_len(C)),
                     paste0("omega_cohort_", seq_len(C)), paste0("level_", seq_len(C)))
  out
}

# Fit of the fields of one group: chains on parallel streams, refits with
# longer chains while not converged; fitted probabilities from up to 500
# draws
mc_fit <- function(Y, N, ppa, ord, mcmc, chains = 4L, cores = 1L, max_refits = 1L, hyper = mc_hyper()) {
  st <- mc_setup(Y, N, ppa, ord, hyper)
  o <- mcmc$options
  n_iter <- if (is.numeric(o$number_of_iterations)) o$number_of_iterations else 10000
  burn <- if (is.numeric(o$burn_in)) o$burn_in else 2000
  thin <- if (is.numeric(o$step)) o$step else 10
  msg <- character()
  for (attempt in 0:max_refits) {
    res <- par_lapply(seq_len(chains), function(i) mc_chain(st, n_iter, burn, thin), cores = cores, streams = TRUE)
    if (any(vapply(res, function(r) !is.list(r) || is.null(r$x), logical(1)))) {
      mt_abort("A chain of the multi-cause sampler failed.", class = "model_fit")
    }
    mon <- lapply(res, function(r) mc_monitor(st, r))
    rhat <- stats::setNames(coh_psrf(mon), colnames(mon[[1]]))
    conv <- all(is.finite(rhat)) && max(rhat) < 1.1
    if (conv) break
    n_iter <- 2 * n_iter
    burn <- 2 * burn
    msg <- c(msg, sprintf("not converged (attempt %d, largest R-hat %.2f)", attempt + 1, max(rhat, na.rm = TRUE)))
  }
  if (!conv) {
    mt_abort(sprintf(paste("The multi-cause chains did not converge after %d refit(s) (largest R-hat %.2f); use",
                           "longer chains (mcmc_control())."), max_refits, max(rhat, na.rm = TRUE)),
             class = "convergence")
  }
  bind3 <- function(nm) {
    a <- lapply(res, `[[`, nm)
    out <- array(NA_real_, c(sum(vapply(a, function(z) dim(z)[1], 0)), st$C, st$C))
    off <- 0
    for (z in a) {
      out[off + seq_len(dim(z)[1]), , ] <- z
      off <- off + dim(z)[1]
    }
    out
  }
  st$X <- NULL
  st$tX <- NULL
  list(setup = st, x = do.call(rbind, lapply(res, `[[`, "x")), kappa = do.call(rbind, lapply(res, `[[`, "kappa")),
       Om = bind3("Om"), Oc = bind3("Oc"), psrf = rhat, accept = vapply(res, `[[`, 0, "accept"),
       message = paste(msg, collapse = "; "))
}

# Draws of the probabilities of every field in the years `years` (indexes
# beyond J are forecast): array [age group, years, field, draws], from the
# posterior draws `pick`
mc_draws <- function(fit, h, pick) {
  st <- fit$setup
  I <- st$I
  J <- st$J
  C <- st$C
  M <- st$M
  K <- st$K
  ix <- st$idx
  nd <- length(pick)
  yrs <- seq_len(J + h)
  out <- array(NA_real_, c(I, J + h, C, nd))
  ext <- function(Z, prec, order) {
    if (h == 0) return(Z)
    R <- chol(solve(prec))
    n <- nrow(Z)
    Z <- rbind(Z, matrix(0, h, ncol(Z)))
    for (t in n + seq_len(h)) {
      e <- as.vector(stats::rnorm(ncol(Z)) %*% R)
      Z[t, ] <- if (order == 1L) Z[t - 1, ] + e else 2 * Z[t - 1, ] - Z[t - 2, ] + e
    }
    Z
  }
  for (r in seq_len(nd)) {
    x <- fit$x[pick[r], ]
    Ph <- ext(matrix(x[ix$ph], J, C), st$Kp$scale * fit$Om[pick[r], , ], st$ord[2])
    Ps <- ext(matrix(x[ix$ps], K, C), st$Kc$scale * fit$Oc[pick[r], , ], st$ord[3])
    for (f in seq_len(C)) {
      th <- x[ix$th[(f - 1L) * I + seq_len(I)]]
      for (i in seq_len(I)) {
        out[i, , f, r] <- stats::plogis(x[ix$mu[f]] + th[i] + Ph[yrs, f] + Ps[M * (I - i) + yrs, f])
      }
    }
  }
  out
}

# Shares of the parts from the conditional probabilities of the
# stick-breaking fields: array [age, year, part, draw] with one more part
stick_shares <- function(p) {
  d <- dim(p)
  out <- array(NA_real_, c(d[1], d[2], d[3] + 1L, d[4]))
  rest <- array(1, d[-3])
  for (f in seq_len(d[3])) {
    out[, , f, ] <- rest * p[, , f, ]
    rest <- rest * (1 - p[, , f, ])
  }
  out[, , d[3] + 1L, ] <- rest
  out
}

#' Multi-cause age-period-cohort model with coupled cause trends
#'
#' Forecasts the deaths of several causes together with all causes, so that
#' the cause-specific rates add up to the all-cause rate in every draw and
#' the trends of the causes move together as estimated. Package-specific
#' (reference unknown): no published source of this combination is known to
#' us.
#'
#' **Model.** In each group of strata that differ only in `var` (for example
#' the causes of one sex), the all-cause deaths follow a binomial-logit
#' age-period-cohort model with random-walk priors (as [apc_bayes()]). The
#' cause deaths split the all-cause deaths of each cell by stick-breaking: a
#' sequence of binomial-logit models of the share of each cause among the
#' deaths not yet assigned, `logit pi[c, i, j] = mu[c] + theta[c, i] +
#' phi[c, j] + psi[c, k]`, the causes ordered by their number of deaths and
#' the last cause taking the rest (Linderman, Johnson and Adams 2015 for the
#' stick-breaking representation of a multinomial). The period effects of
#' the shares have a joint random-walk prior whose precision across causes
#' is a matrix with a Wishart prior, so that period trends of different
#' causes can be correlated, also negatively (one cause replacing another);
#' the cohort effects likewise. Cause rates are shares times the all-cause
#' rate: they add up, and hazards built from them are additive. Within an
#' age group, every cause takes the age pattern of the all-cause rate on
#' the modelling grid.
#'
#' **Sources.** The model is the multi-cause model of the development
#' branch "coherent_forecasting" of bamp (`bamp_multicause()`, by Chris
#' Kypridemos), rewritten here because the branch is not released; when a
#' bamp release contains it, `apc_causes()` will call bamp and the native
#' sampler will be deprecated. Changes: the exact sampler of
#' `apc_bayes(coherent = TRUE)`; one zero-slope constraint on the period
#' effect of each share, only when all effects are second-order random
#' walks (the branch also pins the cohort slope, which removes an
#' identified trend); the all-cause model fitted with the same sampler, its
#' period and cohort precisions with Gamma priors (Wishart priors of
#' dimension 1); the latent-factor coupling of the branch is not taken over.
#' Priors: Gamma(1, 0.5) for the precisions of the age effects, Wishart
#' with `C + 2` degrees of freedom and scale `(1e-4 I)^-1` for the period
#' and cohort precisions of `C` fields, as in the branch.
#'
#' **Sampler and convergence** as [apc_bayes()] with `coherent = TRUE`: four
#' chains, the Gelman-Rubin statistic of the precisions of the age effects,
#' the diagonals of the Wishart precisions, the last period effect and the
#' level of every field in the last year.
#'
#' @param var Stratum variable with the causes.
#' @param total_level Level of `var` for all causes.
#' @param residual Add a residual cause, all deaths not covered by the
#'   listed causes: `TRUE`, `FALSE`, or `"auto"` (added when the listed
#'   causes leave more than 0.1% of all-cause deaths unexplained in the
#'   median cell). The residual cause is modelled but not returned.
#' @inheritParams apc_bayes
#' @return A model specification.
#' @references
#' Linderman SW, Johnson MJ, Adams RP (2015). Dependent multinomial models
#' made easy: stick-breaking with the Pólya-gamma augmentation.
#' arXiv:1506.05843.
#'
#' Schmid VJ, Held L (2007). Bayesian age-period-cohort modeling and
#' prediction: BAMP. *Journal of Statistical Software* 21(8): 1-15.
#' \doi{10.18637/jss.v021.i08}
#' @seealso [coda_causes()] for a published model of causes of death, and
#'   [reconciled()] with a hierarchy of causes.
#' @examples
#' \donttest{
#' x <- mt_example_causes(years = 2005:2019)
#' fit <- mt_fit(apc_causes(mcmc = mcmc_control("test"), seed = 1), x)
#' fc <- mt_forecast(fit, h = 5, n_draws = 100, seed = 1)
#' }
#' @export
apc_causes <- function(var = "cause", total_level = "all", residual = "auto", age = "rw2", period = "rw1",
                       cohort = "rw1", mcmc = mcmc_control(), resolution = c("observed", "grid"), shape_years = 5,
                       max_refits = 1, seed = NULL) {
  if (!is.character(var) || length(var) != 1L) mt_abort("`var` must be one stratum variable.", class = "argument")
  if (length(total_level) != 1L) mt_abort("`total_level` must be one level.", class = "argument")
  if (!(isTRUE(residual) || isFALSE(residual) || identical(residual, "auto"))) {
    mt_abort("`residual` must be TRUE, FALSE or \"auto\".", class = "argument")
  }
  for (nm in c("age", "period", "cohort")) {
    if (!get(nm) %in% c("rw1", "rw2")) mt_abort(sprintf("`%s` must be \"rw1\" or \"rw2\".", nm), class = "argument")
  }
  new_mt_spec("apc_causes", list(var = var, total_level = total_level, residual = residual, age = age, period = period,
                                 cohort = cohort, mcmc = mcmc, resolution = match.arg(resolution),
                                 shape_years = shape_years, max_refits = max_refits, seed = seed),
              joint = TRUE, min_periods = 5, label = "Bayesian APC by cause")
}

#' @export
mt_fit.mt_spec_apc_causes <- function(spec, data, prep = NULL, cores = 1, ...) {
  check_class(data, "mt_data", "data")
  p <- spec$params
  pr <- mt_prep_apply_raw(prep, data)
  d <- pr$data
  md <- model_inputs(d, spec)
  st <- d$strata
  if (!p$var %in% names(st)) mt_abort(sprintf("`var`: \"%s\" is not a stratum variable.", p$var), class = "argument")
  src <- if (p$resolution == "observed") data else d
  cells <- apc_cells(src, md$periods)
  shape <- apc_shape(md, cells$age, p$shape_years)
  others <- setdiff(names(st), c("stratum", "label", p$var))
  key <- if (length(others)) do.call(paste, c(unname(as.list(st[, others, with = FALSE])), sep = ".")) else
    rep("all", nrow(st))
  groups <- split(seq_len(nrow(st)), factor(key, levels = unique(key)))
  ord <- match(c(p$age, p$period, p$cohort), c("rw1", "rw2"))
  chain_cores <- if (isTRUE(p$mcmc$parallel)) max(1L, min(4L, cores)) else 1L
  seed <- p$seed %||% sample.int(.Machine$integer.max, 1L)
  fitted <- array(NA_real_, dim(md$D))
  out <- list()
  diag <- list()
  for (g in names(groups)) {
    idx <- groups[[g]]
    tot <- idx[st[[p$var]][idx] == p$total_level]
    cs <- setdiff(idx, tot)
    if (length(tot) != 1L) {
      mt_abort(sprintf("Group %s needs one stratum with %s = \"%s\".", g, p$var, p$total_level), class = "model_data")
    }
    Dt <- cells$D[, , tot]
    Et <- cells$E[, , tot]
    Yc <- cells$D[, , cs, drop = FALSE]
    rest <- Dt - apply(Yc, c(1, 2), sum)
    if (any(rest < -0.5)) mt_abort(sprintf("In group %s the causes exceed all causes.", g), class = "model_data")
    rest <- pmax(rest, 0)
    add <- isTRUE(p$residual) || (identical(p$residual, "auto") && stats::median(rest / pmax(Dt, 1)) > 0.001)
    parts <- if (add) array(c(Yc, rest), c(dim(Yc)[1:2], length(cs) + 1L)) else Yc
    part_strata <- c(cs, if (add) NA_integer_)
    if (dim(parts)[3] < 2L) {
      mt_abort(sprintf("Group %s needs at least two causes (with the residual).", g), class = "argument")
    }
    # stick-breaking order: most deaths first; the last part takes the rest
    o <- order(apply(parts, 3, sum), decreasing = TRUE)
    parts <- parts[, , o, drop = FALSE]
    part_strata <- part_strata[o]
    nf <- dim(parts)[3] - 1L
    R <- array(0, c(dim(parts)[1:2], nf))
    R[, , 1] <- Dt
    for (f in seq_len(nf - 1L)) R[, , f + 1L] <- R[, , f] - parts[, , f]
    t0 <- proc.time()[["elapsed"]]
    gi <- match(g, names(groups))
    ftot <- with_seed(seed + 2L * gi, mc_fit(array(Dt, c(dim(Dt), 1L)), array(Et, c(dim(Et), 1L)), cells$ppa, ord,
                                            p$mcmc, cores = chain_cores, max_refits = p$max_refits))
    fsh <- with_seed(seed + 2L * gi + 1L, mc_fit(parts[, , seq_len(nf), drop = FALSE], R, cells$ppa, ord, p$mcmc,
                                                cores = chain_cores, max_refits = p$max_refits))
    secs <- proc.time()[["elapsed"]] - t0
    res <- list(total = ftot, shares = fsh, total_stratum = tot, part_strata = part_strata)
    # fitted rates: posterior means of the all-cause rate and of the shares
    # (independent posteriors, draws paired at random)
    nd <- min(500L, nrow(ftot$x), nrow(fsh$x))
    fr <- mc_cause_rates(res, 0L, sample.int(nrow(ftot$x), nd), sample.int(nrow(fsh$x), nd))
    for (q in seq_len(dim(fr)[3])) {
      s <- c(tot, part_strata)[q]
      if (is.na(s)) next
      fitted[, , s] <- apply(fr[, , q, , drop = FALSE], c(1, 2), mean)[shape$group, , drop = FALSE] *
        shape$r[, tot]
    }
    out[[g]] <- res
    diag[[g]] <- data.table(group = g, strata = paste(st$label[idx], collapse = ", "), status = "ok",
                            message = paste(Filter(nzchar, c(ftot$message, fsh$message)), collapse = "; "),
                            seconds = secs,
                            max_rhat = max(c(ftot$psrf, fsh$psrf)), accept = mean(c(ftot$accept, fsh$accept)))
  }
  structure(list(spec = spec, prep = prep, fits = list(groups = out, fitted = fitted, shape = shape),
                 diagnostics = rbindlist(diag), age = md$age, periods = md$periods, strata = st, omega = d$omega,
                 observed = (raw_deaths(pr$raw, md) %||% md$D) / md$E, train_D = md$D,
                 raw_D = raw_deaths(pr$raw, md), train_E = md$E,
                 exposure_last = md$E[, length(md$periods), , drop = FALSE], level_totals = data$level_totals,
                 apc_age = cells$age),
            class = "mt_fit")
}

# Rates of all causes and of every part of a group: array [age group, year,
# 1 + parts, draws], the first slice all causes (years: the fitted years
# for h = 0, else the h forecast years)
mc_cause_rates <- function(res, h, pick_tot, pick_sh) {
  tot <- mc_draws(res$total, h, pick_tot)
  sh <- stick_shares(mc_draws(res$shares, h, pick_sh))
  yrs <- if (h == 0) seq_len(dim(tot)[2]) else dim(tot)[2] - h + seq_len(h)
  tot <- tot[, yrs, 1, , drop = FALSE]
  sh <- sh[, yrs, , , drop = FALSE]
  out <- array(NA_real_, c(dim(sh)[1:2], dim(sh)[3] + 1L, dim(sh)[4]))
  out[, , 1, ] <- tot[, , 1, ]
  for (q in seq_len(dim(sh)[3])) out[, , q + 1L, ] <- sh[, , q, ] * tot[, , 1, ]
  out
}

#' @export
simulate_model.mt_spec_apc_causes <- function(spec, fit, h, n_draws, ...) {
  shape <- fit$shape
  ns <- ncol(shape$r)
  out <- array(NA_real_, c(length(shape$group), h, ns, n_draws))
  for (res in fit$groups) {
    pick <- function(f) sample.int(nrow(f$x), n_draws, replace = n_draws > nrow(f$x))
    r <- mc_cause_rates(res, h, pick(res$total), pick(res$shares))
    tot <- res$total_stratum
    total_grid <- r[shape$group, , 1, , drop = FALSE] * shape$r[, tot]
    out[, , tot, ] <- total_grid
    # each cause: its share of the all-cause rate on the grid
    for (q in seq_along(res$part_strata)) {
      s <- res$part_strata[q]
      if (is.na(s)) next
      share <- r[, , q + 1L, , drop = FALSE] / r[, , 1, , drop = FALSE]
      out[, , s, ] <- share[shape$group, , 1, , drop = FALSE] * total_grid
    }
  }
  out
}

# mt_parameters() of apc_causes(): posterior medians and 90% intervals of
# the effects of the all-cause model and of each share, the correlations of
# the period and cohort effects of the shares across causes (from the
# posterior mean of the inverse Wishart precisions), and the Gelman-Rubin
# statistics
mc_parameters <- function(fit) {
  lab <- fit$strata$label
  rbindlist(lapply(names(fit$fits$groups), function(g) {
    res <- fit$fits$groups[[g]]
    q <- function(M) t(apply(as.matrix(M), 2, stats::quantile, c(0.05, 0.5, 0.95)))
    one <- function(who, parameter, M) {
      Q <- q(M)
      data.table(stratum = who, parameter = parameter, index = as.character(seq_len(nrow(Q))), value = Q[, 2],
                 lower = Q[, 1], upper = Q[, 3])
    }
    field <- function(f, st, who, fit1) {
      rbind(one(who, "intercept", fit1$x[, st$idx$mu[f], drop = FALSE]),
            one(who, "age", fit1$x[, st$idx$th[(f - 1L) * st$I + seq_len(st$I)], drop = FALSE]),
            one(who, "period", fit1$x[, st$idx$ph[(f - 1L) * st$J + seq_len(st$J)], drop = FALSE]),
            one(who, "cohort", fit1$x[, st$idx$ps[(f - 1L) * st$K + seq_len(st$K)], drop = FALSE]))
    }
    names_parts <- ifelse(is.na(res$part_strata), paste0("residual (", g, ")"), lab[res$part_strata])
    sh <- res$shares
    cors <- function(O) {
      S <- apply(O, 1, function(m) stats::cov2cor(solve(m)))
      matrix(rowMeans(matrix(S, ncol = dim(O)[1])), dim(O)[2])
    }
    nf <- sh$setup$C
    cp <- cors(sh$Om)
    cc <- cors(sh$Oc)
    pr <- which(upper.tri(cp), arr.ind = TRUE)
    rbind(field(1L, res$total$setup, lab[res$total_stratum], res$total),
          rbindlist(lapply(seq_len(nf), function(f) field(f, sh$setup, paste("share", names_parts[f]), sh))),
          if (nrow(pr)) data.table(stratum = paste0("shares (", g, ")"), parameter = "period_correlation",
                                   index = paste(names_parts[pr[, 1]], names_parts[pr[, 2]], sep = ":"),
                                   value = cp[pr], lower = NA_real_, upper = NA_real_),
          if (nrow(pr)) data.table(stratum = paste0("shares (", g, ")"), parameter = "cohort_correlation",
                                   index = paste(names_parts[pr[, 1]], names_parts[pr[, 2]], sep = ":"),
                                   value = cc[pr], lower = NA_real_, upper = NA_real_),
          data.table(stratum = paste0("all causes (", g, ")"), parameter = "rhat", index = names(res$total$psrf),
                     value = unname(res$total$psrf), lower = NA_real_, upper = NA_real_),
          data.table(stratum = paste0("shares (", g, ")"), parameter = "rhat", index = names(sh$psrf),
                     value = unname(sh$psrf), lower = NA_real_, upper = NA_real_))
  }))
}
