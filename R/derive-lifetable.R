#' Life tables from death rates on any age axis
#'
#' Period life table for single or grouped ages with an open last interval.
#' Interval widths come from the age axis, so abridged and single-year
#' schedules are handled by the same code.
#'
#' * `q_x = n m / (1 + (n - a) m)`; where this exceeds 1 (rates above
#'   about `2 / n`, far beyond observed values), the constant-hazard value
#'   `1 - exp(-n m)` is used, so that survivors stay positive;
#' * `a_x` defaults to `n / 2` for every closed interval (supply `ax` to
#'   override, for example for age 0, where `n / 2` overstates the average
#'   time lived by infants who die; the effect on e0 at low infant mortality
#'   is below 0.01 years);
#' * open interval: `"reciprocal"` uses `L = l / m`; `"gompertz"` fits a
#'   Gompertz line to the log rates of the last `k` closed single-year ages and
#'   integrates it from the start of the open interval to age 130;
#'   `"kannisto"` does the same with the Kannisto model, a line for the logit
#'   of the rates, whose hazard levels off at old ages (Thatcher, Kannisto and
#'   Vaupel 1998).
#' @references
#' Thatcher AR, Kannisto V, Vaupel JW (1998). *The Force of Mortality at Ages
#' 80 to 120*. Odense University Press.
#'
#' @param mx Death rates, one per interval of `age`.
#' @param age An [mt_axis] (default: single ages from 0 with an open last
#'   interval).
#' @param ax Optional average years lived in the interval by those who die.
#' @param open `"reciprocal"`, `"gompertz"` or `"kannisto"`.
#' @param k Number of ages used for the Gompertz fit.
#' @param radix Starting number of survivors.
#' @return A `data.table` with `age_lo`, `age_hi`, `n`, `mx`, `ax`, `qx`,
#'   `lx`, `dx`, `Lx`, `Tx`, `ex`.
#' @examples
#' m <- 3e-4 + 2e-5 * exp(0.1 * (30:100))
#' lt <- life_table(m, single_ages(30, 100))
#' lt[age_lo %in% c(30, 65, 100)]
#' @export
life_table <- function(mx, age = NULL, ax = NULL, open = c("reciprocal", "gompertz", "kannisto"), k = 10,
                       radix = 1e5) {
  open <- match.arg(open)
  age <- age %||% single_ages(0, length(mx) - 1)
  check_class(age, "mt_axis", "age")
  if (length(mx) != length(age)) mt_abort("`mx` needs one value per age interval.", class = "lifetable")
  r <- lt_columns(matrix(mx, ncol = 1), age, ax, open, k)
  out <- data.table(age_lo = age$lower, age_hi = age$upper, n = age$upper - age$lower, mx = mx, ax = r$ax[, 1],
                    qx = r$qx[, 1], lx = radix * r$lx[, 1], dx = radix * r$dx[, 1], Lx = radix * r$Lx[, 1],
                    Tx = radix * r$Tx[, 1], ex = r$ex[, 1])
  out[]
}

# Column-wise life tables: M is a matrix [age, schedule]
lt_columns <- function(M, age, ax = NULL, open = "reciprocal", k = 10) {
  na <- nrow(M)
  n <- age$upper - age$lower
  is_open <- is.infinite(n[na])
  nn <- matrix(n, na, ncol(M))
  A <- if (is.null(ax)) nn / 2 else matrix(ax, na, ncol(M))
  Q <- lt_qx(M, nn, A)
  if (is_open) {
    Q[na, ] <- 1
    A[na, ] <- 1 / M[na, ]
  }
  P <- 1 - Q
  L <- rbind(1, apply(P, 2, cumprod))[seq_len(na), , drop = FALSE]
  Dx <- L * Q
  Lx <- nn * (L - Dx) + A * Dx
  if (is_open) {
    Lx[na, ] <- if (open == "reciprocal") L[na, ] / M[na, ] else L[na, ] * gompertz_open_ex(M, age, k, open)
  }
  Tx <- apply(Lx, 2, function(z) rev(cumsum(rev(z))))
  if (is.null(dim(Tx))) Tx <- matrix(Tx, nrow = na)
  list(ax = A, qx = Q, lx = L, dx = Dx, Lx = Lx, Tx = Tx, ex = Tx / L)
}

# Probabilities of dying in closed intervals from rates; where the linear
# formula exceeds 1, the constant-hazard value 1 - exp(-n m)
lt_qx <- function(M, nn, A) {
  Q <- nn * M / (1 + (nn - A) * M)
  hi <- is.finite(Q) & Q > 1 | is.finite(M) & M > 0 & !is.finite(Q)
  Q[hi] <- 1 - exp(-nn[hi] * M[hi])
  Q[!is.finite(Q)] <- 1
  Q
}

# Life expectancy at the start of the open interval under a Gompertz extension
gompertz_open_ex <- function(M, age, k, law = "gompertz") {
  na <- nrow(M)
  single <- which(age$upper - age$lower == 1)
  use <- tail(single[single < na], k)
  if (length(use) < 3) mt_abort("The Gompertz extension needs at least 3 single-year ages below the open interval.",
                                class = "lifetable")
  x <- age$lower[use] + 0.5
  X <- cbind(1, x)
  Y <- M[use, , drop = FALSE]
  Y <- if (law == "kannisto") log(pmin(Y, 0.99) / (1 - pmin(Y, 0.99))) else log(Y)
  coef <- solve(crossprod(X), crossprod(X, Y))
  ages <- seq(age$lower[na], 129) + 0.5
  eta <- cbind(1, ages) %*% coef
  mm <- if (law == "kannisto") exp(eta) / (1 + exp(eta)) else exp(eta)
  p <- exp(-mm)
  l <- rbind(1, apply(p, 2, cumprod))[seq_along(ages), , drop = FALSE]
  Lx <- l * (1 - p) / pmax(mm, 1e-12)
  colSums(Lx)
}

#' Life expectancy from rates or forecast draws
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
#' @param ... Passed to the life-table calculation (`ax`, `open`, `k`).
#' @return For arrays, an array without the age dimension. For draws, a
#'   `data.table` with the mean and quantiles per stratum and period, computed
#'   path by path.
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
