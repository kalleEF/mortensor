#' Life tables from death rates on any age axis
#'
#' Period life table for single or grouped ages with an open last interval.
#' Interval widths come from the age axis, so abridged and single-year
#' schedules are handled by the same code. The relations are the standard
#' ones of Preston, Heuveline and Guillot (2001):
#'
#' * `q_x = n m / (1 + (n - a) m)` for a closed interval of width `n`; where
#'   this exceeds 1 (rates above about `2 / n`, far beyond observed values),
#'   the constant-hazard value `1 - exp(-n m)` is used, so that survivors
#'   stay positive;
#' * `a_x` defaults to `n / 2` for every closed interval, the simple
#'   convention of deaths spread evenly over the interval (supply `ax` to
#'   override). At age 0, `n / 2` overstates the average time lived by
#'   infants who die; the Human Mortality Database uses the age-0 values of
#'   Andreev and Kingkade (2015), which depend on sex and on the infant
#'   death rate and are not implemented here. The effect on e0 at low
#'   infant mortality is below 0.01 years;
#' * open interval: `"reciprocal"` (default) uses `L = l / m` (Preston,
#'   Heuveline and Guillot 2001). `"gompertz"` and `"kannisto"` close the
#'   table with the Gompertz hazard `mu(x) = a exp(b x)` or the Kannisto
#'   hazard `mu(x) = a exp(b x) / (1 + a exp(b x))`, whose logit is linear in
#'   age and which levels off at old ages (Thatcher, Kannisto and Vaupel
#'   1998, chapter 2). They are fitted as in Thatcher, Kannisto and Vaupel
#'   (1998, chapter 6), by maximum likelihood with binomial deaths: of the
#'   `N` persons who reach age `x`, `d` die before `x + 1` with probability
#'   `q = 1 - exp(-H)`, where `H` is the integral of the hazard over the
#'   year. The fit uses the last `k` closed single-year ages. `N` and `d`
#'   are the survivors and deaths of the life table itself, the numbers of a
#'   stationary population with these rates, unless `exposure` gives the
#'   observed person-years `E` (then `d = m E` and `N = E + (1 - a) d`).
#'   Thatcher et al. fit observed numbers, at ages 80 to 98; the life-table
#'   numbers are a package-specific choice for schedules without counts
#'   (reference unknown). For the cohorts born in 1871-1880 the life-table
#'   numbers are the observed ones, and the fit reproduces the published
#'   parameters (their Appendix D) from their probabilities of dying (Table
#'   6.1). Earlier versions of mortensor fitted least squares to the log or
#'   logit rates. The fitted hazard is applied year by year from the start
#'   of the open interval to `omega` (survival `exp(-H)` and person-years
#'   `l (1 - exp(-H)) / H` in each year), and the survivors at `omega` are
#'   closed with `L = l / H` at the integrated hazard of the following year;
#' * a rate of 0 in the open interval (no deaths, as in small populations)
#'   would give an infinite life expectancy: the last positive rate of the
#'   closed intervals is used for the open interval instead (column `mx`
#'   keeps the given rate). Only a schedule without any positive rate gives
#'   an infinite value.
#' @references
#' Andreev EM, Kingkade WW (2015). Average age at death in infancy and
#' infant mortality level: reconsidering the Coale-Demeny formulas at
#' current levels of low mortality. *Demographic Research* 33: 363-390.
#'
#' Preston SH, Heuveline P, Guillot M (2001). *Demography: Measuring and
#' Modeling Population Processes*. Blackwell.
#'
#' Thatcher AR, Kannisto V, Vaupel JW (1998). *The Force of Mortality at Ages
#' 80 to 120*. Odense University Press.
#'
#' @param mx Death rates, one per interval of `age`.
#' @param age An [mt_axis] (default: single ages from 0 with an open last
#'   interval).
#' @param ax Optional average years lived in the interval by those who die.
#' @param open `"reciprocal"`, `"gompertz"` or `"kannisto"`.
#' @param k Number of ages used for the Gompertz or Kannisto fit.
#' @param omega Age up to which a Gompertz or Kannisto tail is applied
#'   (default 110, the default `omega` of [mt_data()]); ignored with
#'   `open = "reciprocal"`.
#' @param radix Starting number of survivors.
#' @param exposure Optional person-years of the age intervals, which weight
#'   the maximum-likelihood fit of a Gompertz or Kannisto closure with the
#'   observed numbers at risk; by default the numbers of the life table.
#' @return A `data.table` with `age_lo`, `age_hi`, `n`, `mx`, `ax`, `qx`,
#'   `lx`, `dx`, `Lx`, `Tx`, `ex`. With a Gompertz or Kannisto closure, the
#'   attribute `"tail"` holds the fitted `a` and `b`, the ages of the fit and
#'   whether it converged.
#' @examples
#' m <- 3e-4 + 2e-5 * exp(0.1 * (30:100))
#' lt <- life_table(m, single_ages(30, 100))
#' lt[age_lo %in% c(30, 65, 100)]
#' # Kannisto closure fitted by maximum likelihood to ages 90-99
#' lk <- life_table(m, single_ages(30, 100), open = "kannisto")
#' attr(lk, "tail")
#' @export
life_table <- function(mx, age = NULL, ax = NULL, open = c("reciprocal", "gompertz", "kannisto"), k = 10,
                       radix = 1e5, omega = 110, exposure = NULL) {
  open <- match.arg(open)
  age <- age %||% single_ages(0, length(mx) - 1)
  check_class(age, "mt_axis", "age")
  if (length(mx) != length(age)) mt_abort("`mx` needs one value per age interval.", class = "lifetable")
  check_number(omega, "omega", lower = 1)
  check_number(k, "k", lower = 3, integerish = TRUE)
  if (!is.null(exposure) && (!is.numeric(exposure) || length(exposure) != length(mx))) {
    mt_abort("`exposure` needs one value per age interval.", class = "lifetable")
  }
  r <- lt_columns(matrix(mx, ncol = 1), age, ax, open, k, omega, exposure)
  out <- data.table(age_lo = age$lower, age_hi = age$upper, n = age$upper - age$lower, mx = mx, ax = r$ax[, 1],
                    qx = r$qx[, 1], lx = radix * r$lx[, 1], dx = radix * r$dx[, 1], Lx = radix * r$Lx[, 1],
                    Tx = radix * r$Tx[, 1], ex = r$ex[, 1])
  if (!is.null(r$tail)) {
    setattr(out, "tail", data.table(law = open, a = exp(r$tail$alpha), b = r$tail$b, age_from = min(r$tail$ages),
                                    age_to = max(r$tail$ages) + 1, converged = r$tail$converged))
  }
  out[]
}

# Column-wise life tables: M is a matrix [age, schedule]; `exposure`
# (person-years, a vector per age or a matrix like M) weights the fit of a
# Gompertz or Kannisto closure
lt_columns <- function(M, age, ax = NULL, open = "reciprocal", k = 10, omega = 110, exposure = NULL) {
  na <- nrow(M)
  n <- age$upper - age$lower
  is_open <- is.infinite(n[na])
  nn <- matrix(n, na, ncol(M))
  A <- if (is.null(ax)) nn / 2 else matrix(ax, na, ncol(M))
  Q <- lt_qx(M, nn, A)
  if (is_open) {
    m_open <- open_rates(M)
    Q[na, ] <- 1
    A[na, ] <- 1 / m_open
  }
  P <- 1 - Q
  L <- rbind(1, apply(P, 2, cumprod))[seq_len(na), , drop = FALSE]
  Dx <- L * Q
  Lx <- nn * (L - Dx) + A * Dx
  fit <- NULL
  if (is_open) {
    if (open == "reciprocal") {
      Lx[na, ] <- L[na, ] / m_open
    } else {
      if (!is.null(exposure) && !length(exposure) %in% c(na, length(M))) {
        mt_abort("`exposure` needs one value per age interval (and schedule).", class = "lifetable")
      }
      fit <- tail_fit(M, L, Q, A, age, k, open, exposure)
      if (any(!fit$converged & !is.na(fit$alpha))) {
        mt_warn(sprintf("The %s fit of the open interval did not converge for %d of %d schedules.", open,
                        sum(!fit$converged & !is.na(fit$alpha)), ncol(M)), class = "lifetable")
      }
      Lx[na, ] <- L[na, ] * tail_open_ex(fit, age$lower[na], omega, open)
    }
  }
  Tx <- apply(Lx, 2, function(z) rev(cumsum(rev(z))))
  if (is.null(dim(Tx))) Tx <- matrix(Tx, nrow = na)
  list(ax = A, qx = Q, lx = L, dx = Dx, Lx = Lx, Tx = Tx, ex = Tx / L, tail = fit)
}

# Rates of the open interval used in the life table: a rate of 0 (no
# deaths) would make life expectancy infinite and is replaced by the last
# positive rate of the closed intervals
open_rates <- function(M) {
  na <- nrow(M)
  m <- M[na, ]
  zero <- !is.na(m) & m <= 0
  if (any(zero) && na > 1) {
    m[zero] <- apply(M[-na, zero, drop = FALSE], 2, function(v) {
      v <- v[!is.na(v) & v > 0]
      if (length(v)) v[length(v)] else 0
    })
  }
  m
}

# Probabilities of dying in closed intervals from rates; where the linear
# formula exceeds 1, the constant-hazard value 1 - exp(-n m). Missing rates
# give missing probabilities (and life expectancies), not certain death.
lt_qx <- function(M, nn, A) {
  Q <- nn * M / (1 + (nn - A) * M)
  hi <- is.finite(Q) & Q > 1 | is.finite(M) & M > 0 & !is.finite(Q)
  Q[hi] <- 1 - exp(-nn[hi] * M[hi])
  Q[!is.finite(Q) & !is.na(M)] <- 1
  Q
}

# Gompertz hazard mu(x) = a exp(b x) or Kannisto hazard mu(x) = a exp(b x) /
# (1 + a exp(b x)) (Thatcher, Kannisto and Vaupel 1998, chapter 2), with
# alpha = log(a): the integrated hazard H over [x, x + 1] and its
# derivatives with respect to alpha and b, matrices [age, schedule] for the
# ages x and the parameters of each schedule (vectors alpha, b)
tail_hazard <- function(alpha, b, x, law) {
  nx <- length(x)
  A <- matrix(alpha, nx, length(alpha), byrow = TRUE)
  Bm <- matrix(b, nx, length(b), byrow = TRUE)
  X <- matrix(x, nx, length(alpha))
  small <- abs(Bm) < 1e-6
  Bs <- ifelse(small, 1, Bm)
  if (law == "gompertz") {
    # H = a exp(b x) (exp(b) - 1) / b; d log H / d b = x + exp(b) / (exp(b) - 1) - 1 / b
    e1 <- ifelse(small, 1 + Bm / 2, expm1(Bs) / Bs)
    g1 <- ifelse(small, 0.5 + Bm / 12, exp(Bs) / expm1(Bs) - 1 / Bs)
    H <- exp(A + Bm * X) * e1
    return(list(H = H, da = H, db = H * (X + g1)))
  }
  # Kannisto: H = (s(x + 1) - s(x)) / b with s(t) = log(1 + exp(alpha + b t))
  softplus <- function(z) pmax(z, 0) + log1p(exp(-abs(z)))
  z0 <- A + Bm * X
  z1 <- z0 + Bm
  p0 <- stats::plogis(z0)
  p1 <- stats::plogis(z1)
  H <- (softplus(z1) - softplus(z0)) / Bs
  da <- (p1 - p0) / Bs
  db <- (-H + (X + 1) * p1 - X * p0) / Bs
  if (any(small)) {
    pm <- stats::plogis(A + Bm * (X + 0.5))
    H[small] <- pm[small]
    da[small] <- (pm * (1 - pm))[small]
    db[small] <- ((X + 0.5) * pm * (1 - pm))[small]
  }
  list(H = H, da = da, db = db)
}

# Maximum-likelihood fit of a Gompertz or Kannisto hazard to the last k
# single-year ages below the open interval, as in Thatcher, Kannisto and
# Vaupel (1998, chapter 6): of the N_x persons who reach age x, d_x die
# before x + 1, binomially with q_x = 1 - exp(-H_x). N_x and d_x are the
# survivors and deaths of the life table, or, with exposures E (person-years),
# the deaths D = m E and N = E + (1 - a) D. Fisher scoring from least squares
# on the log (logit) rates, with step halving, for all schedules (columns)
# at once. Returns alpha = log(a), b, the ages and whether each fit converged
tail_fit <- function(M, L, Q, A, age, k, law, exposure = NULL) {
  na <- nrow(M)
  single <- which(age$upper - age$lower == 1)
  use <- tail(single[single < na], k)
  if (length(use) < 3) {
    mt_abort(sprintf("The %s closure needs at least 3 single-year ages below the open interval.",
                     if (law == "kannisto") "Kannisto" else "Gompertz"), class = "lifetable")
  }
  x <- age$lower[use]
  if (is.null(exposure)) {
    N <- L[use, , drop = FALSE]
    d <- N * Q[use, , drop = FALSE]
  } else {
    E <- matrix(exposure, na, ncol(M))[use, , drop = FALSE]
    d <- M[use, , drop = FALSE] * E
    N <- E + (1 - A[use, , drop = FALSE]) * d
  }
  nc <- ncol(M)
  alpha <- b <- rep(NA_real_, nc)
  converged <- rep(FALSE, nc)
  ok <- colSums(!is.finite(N) | !is.finite(d) | N < 0 | d < 0) == 0 & colSums(N) > 0
  if (!any(ok)) return(list(alpha = alpha, b = b, ages = x, converged = converged))
  # start: least squares on the log (logit) rates at the middle of each year
  m <- pmax(M[use, ok, drop = FALSE], 1e-12)
  y <- if (law == "kannisto") stats::qlogis(pmin(m, 0.99)) else log(m)
  Xs <- cbind(1, x + 0.5)
  cf <- solve(crossprod(Xs), crossprod(Xs, y))
  al <- cf[1, ]
  bb <- cf[2, ]
  N <- N[, ok, drop = FALSE]
  d <- d[, ok, drop = FALSE]
  loglik <- function(H, j = seq_len(ncol(d))) {
    colSums(d[, j, drop = FALSE] * log(pmax(-expm1(-H), 1e-300)) - (N[, j, drop = FALSE] - d[, j, drop = FALSE]) * H)
  }
  ll <- loglik(tail_hazard(al, bb, x, law)$H)
  conv <- done <- rep(FALSE, length(al))
  for (it in seq_len(100)) {
    act <- which(!done & is.finite(ll))
    if (!length(act)) break
    th <- tail_hazard(al[act], bb[act], x, law)
    q <- pmax(-expm1(-th$H), 1e-300)
    w <- N[, act, drop = FALSE] * exp(-th$H) / q
    r <- (d[, act, drop = FALSE] - N[, act, drop = FALSE] * q) / q
    ua <- colSums(r * th$da)
    ub <- colSums(r * th$db)
    iaa <- colSums(w * th$da^2)
    iab <- colSums(w * th$da * th$db)
    ibb <- colSums(w * th$db^2)
    det <- iaa * ibb - iab^2
    step_a <- (ibb * ua - iab * ub) / det
    step_b <- (iaa * ub - iab * ua) / det
    bad <- !is.finite(step_a) | !is.finite(step_b)
    step_a[bad] <- 0
    step_b[bad] <- 0
    # step halving until the log-likelihood does not fall
    f <- rep(1, length(act))
    for (half in 0:30) {
      ll_new <- loglik(tail_hazard(al[act] + f * step_a, bb[act] + f * step_b, x, law)$H, act)
      worse <- !(ll_new >= ll[act] - 1e-12 * abs(ll[act]))
      if (!any(worse)) break
      f[worse] <- f[worse] / 2
    }
    acc <- !worse
    al[act[acc]] <- al[act[acc]] + f[acc] * step_a[acc]
    bb[act[acc]] <- bb[act[acc]] + f[acc] * step_b[acc]
    ll[act[acc]] <- ll_new[acc]
    # converged when the step is negligible, or when no step along the
    # scoring direction raises the log-likelihood (its maximum to rounding)
    conv[act] <- !bad & ((acc & pmax(abs(f * step_a), abs(f * step_b)) < 1e-10) | !acc)
    done[act] <- conv[act] | bad
  }
  alpha[ok] <- al
  b[ok] <- bb
  converged[ok] <- conv
  list(alpha = alpha, b = b, ages = x, converged = converged)
}

# Life expectancy at the start of the open interval under the fitted tail:
# the integrated hazard of each year from the start of the open interval to
# omega, with person-years L = l (1 - p) / H (constant hazard within the
# year), and the survivors at omega closed with L = l / H at the integrated
# hazard of the year from omega
tail_open_ex <- function(fit, lo, omega, law) {
  ages <- seq(lo, max(lo, ceiling(omega) - 1))
  H <- tail_hazard(fit$alpha, fit$b, c(ages, max(ages) + 1), law)$H
  n <- length(ages)
  Hy <- pmax(H[seq_len(n), , drop = FALSE], 1e-12)
  p <- exp(-Hy)
  l <- rbind(1, matrix(apply(p, 2, cumprod), nrow = n))
  Lx <- l[seq_len(n), , drop = FALSE] * (1 - p) / Hy
  colSums(Lx) + l[n + 1, ] / pmax(H[n + 1, ], 1e-12)
}

#' Life expectancy from rates or forecast draws
#'
#' Life expectancy of the period life table of [life_table()] (Preston,
#' Heuveline and Guillot 2001), for draws path by path.
#'
#' @param x A matrix or array of rates with age in the first dimension, an
#'   [mt_draws] object, or an [mt_data] object (observed rates with Poisson
#'   intervals, see [mt_observed()]).
#' @param age An [mt_axis] (taken from `x` for draws).
#' @param at Age at which life expectancy is computed (start of an
#'   interval; default the first age).
#' @param probs Quantiles reported for draws.
#' @param paths For draws, return the values of every path, an array
#'   `[period, stratum, draw]`.
#' @param ... Passed to the life-table calculation (`ax`, `open`, `k`, `omega`). A
#'   rate of 0 in the open interval is replaced as in [life_table()].
#' @return For arrays, an array without the age dimension. For draws, a
#'   `data.table` with the mean and quantiles per stratum and period, computed
#'   path by path.
#' @references
#' Preston SH, Heuveline P, Guillot M (2001). *Demography: Measuring and
#' Modeling Population Processes*. Blackwell.
#' @examples
#' x <- mt_example(deaths = "single")
#' m <- mt_array(x, "deaths") / mt_array(x, "exposure")
#' life_expectancy(m, mt_grid(x)$age, at = 65)[1:3, ]
#' @export
life_expectancy <- function(x, age = NULL, at = NULL, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), paths = FALSE,
                            ...) {
  if (inherits(x, "mt_data")) x <- mt_observed(x)
  if (inherits(x, "mt_draws")) {
    if (paths) {
      e <- life_expectancy(x$rates, x$age, at = at, ...)
      dimnames(e) <- list(period = x$periods, stratum = x$strata$label, draw = NULL)
      return(e)
    }
    return(life_expectancy_draws(x, at, probs, ...))
  }
  check_class(age, "mt_axis", "age")
  d <- dim(x) %||% length(x)
  M <- matrix(x, nrow = d[1])
  i <- if (is.null(at)) 1L else match(at, age$lower)
  if (is.na(i)) mt_abort("`at` must be the start of an age interval.", class = "lifetable")
  r <- lt_columns(M, age, ...)
  e <- r$ex[i, ]
  if (length(d) > 1) array(e, d[-1], dimnames = dimnames(x)[-1]) else e
}
