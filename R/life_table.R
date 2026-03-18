#' Convert mortality rates to a death distribution (dx)
#'
#' Applies a standard life-table calculation to convert age-specific mortality
#' rates (mx) to a death distribution (dx). Uses the classic nax approximation
#' for grouped ages.
#'
#' @param mx Numeric vector of age-specific mortality rates.
#' @param ages Integer vector of age-group start values.
#' @param width Integer age-group width (default 5).
#' @return Numeric vector of dx values (sums to ~1).
#' @export
compute_dx <- function(mx, ages, width = 5L) {
  n <- length(mx)
  lx <- numeric(n + 1L)
  dx <- numeric(n)
  lx[1] <- 1.0

  for (i in seq_len(n - 1L)) {
    ax <- width / 2
    qx <- (width * mx[i]) / (1 + (width - ax) * mx[i])
    qx <- min(max(qx, 0), 1)
    dx[i] <- lx[i] * qx
    lx[i + 1] <- lx[i] - dx[i]
  }
  # Open-ended interval: all remaining die

  dx[n] <- lx[n]

  # Normalize to a proper distribution
  total <- sum(dx)
  if (total > 0) dx <- dx / total
  # Replace any zeros with tiny value (for CLR)
  dx[dx <= 0] <- 1e-15
  dx <- dx / sum(dx)
  dx
}

#' Convert a death distribution back to mortality rates
#'
#' Inverse of [compute_dx()]. Uses the life-table relationship and
#' Gompertz extrapolation for the open age interval.
#'
#' @param dx Numeric vector of death distribution values.
#' @param ages Integer vector of age-group start values.
#' @param width Integer age-group width (default 5).
#' @return Numeric vector of age-specific mortality rates (mx).
#' @export
dx_to_mx <- function(dx, ages, width = 5L) {
  n   <- length(dx)
  dx  <- dx / sum(dx)
  lx  <- numeric(n + 1L)
  lx[1] <- 1.0
  mx  <- numeric(n)

  for (i in seq_len(n - 1L)) {
    lx[i + 1] <- lx[i] - dx[i]
    if (lx[i] < 1e-15) {
      mx[i] <- mx[max(i - 1, 1)]
      next
    }
    qx <- dx[i] / lx[i]
    qx <- min(max(qx, 1e-15), 1 - 1e-15)
    ax <- width / 2
    mx[i] <- qx / (width - (width - ax) * qx)
  }

  # Open interval: Gompertz extrapolation from last two closed groups
  if (n >= 3 && mx[n - 1] > 1e-15 && mx[n - 2] > 1e-15) {
    log_r <- log(mx[n - 1]) - log(mx[n - 2])
    mx[n] <- mx[n - 1] * exp(log_r)
  } else {
    mx[n] <- mx[max(n - 1, 1)]
  }

  mx <- pmax(mx, 1e-12)
  mx
}

#' Compute a standard abridged life table
#'
#' Converts age-specific mortality rates (mx) into a full life table with
#' columns for qx, px, lx, dx, Lx, Tx, and ex. Handles both grouped and
#' single-year age data, with special treatment of the open age interval.
#'
#' @param mx Numeric vector of age-specific mortality rates.
#' @param ages Numeric vector of age-group start values.
#' @param width Integer or numeric vector of age-group widths. If scalar,
#'   applied to all groups except the open interval. Default 5.
#' @param radix Numeric radix (initial cohort size). Default 100000.
#' @param open_age_method Character; method for the open age interval's ax
#'   value: \code{"kannisto"} (default) fits a Kannisto model to extrapolate,
#'   \code{"gompertz"} uses Gompertz extrapolation, or \code{"reciprocal"}
#'   uses \code{1/mx} for the open interval.
#' @return A \code{data.frame} with columns: \code{age}, \code{n} (width),
#'   \code{mx}, \code{qx}, \code{px}, \code{lx}, \code{dx}, \code{Lx},
#'   \code{Tx}, \code{ex}.
#' @export
life_table <- function(mx, ages, width = 5, radix = 100000,
                       open_age_method = c("reciprocal", "kannisto", "gompertz")) {
  open_age_method <- match.arg(open_age_method)
  n_ages <- length(mx)
  stopifnot(length(ages) == n_ages, n_ages >= 2)


  # Expand width to a vector
  if (length(width) == 1L) {
    n_vec <- rep(width, n_ages)
  } else {
    stopifnot(length(width) == n_ages)
    n_vec <- width
  }
  # Last group is open interval — width is Inf conceptually
  n_vec[n_ages] <- NA_real_

  # --- qx from mx ---
  qx <- numeric(n_ages)
  for (i in seq_len(n_ages - 1L)) {
    ax <- n_vec[i] / 2
    qx[i] <- (n_vec[i] * mx[i]) / (1 + (n_vec[i] - ax) * mx[i])
    qx[i] <- min(max(qx[i], 0), 1)
  }
  qx[n_ages] <- 1.0  # everyone dies in the open interval

  # --- lx, dx ---
  lx <- numeric(n_ages)
  dx <- numeric(n_ages)
  lx[1] <- radix
  for (i in seq_len(n_ages - 1L)) {
    dx[i] <- lx[i] * qx[i]
    lx[i + 1] <- lx[i] - dx[i]
  }
  dx[n_ages] <- lx[n_ages]

  px <- 1 - qx

  # --- Lx (person-years lived) ---
  Lx <- numeric(n_ages)
  for (i in seq_len(n_ages - 1L)) {
    ax <- n_vec[i] / 2
    Lx[i] <- n_vec[i] * lx[i + 1] + ax * dx[i]
  }
  # Open interval: Lx = lx / mx (or via parametric extrapolation)
  if (mx[n_ages] > 0) {
    if (open_age_method == "reciprocal") {
      Lx[n_ages] <- lx[n_ages] / mx[n_ages]
    } else if (open_age_method == "kannisto") {
      Lx[n_ages] <- .open_interval_kannisto(mx, ages, lx[n_ages])
    } else {
      Lx[n_ages] <- .open_interval_gompertz(mx, ages, lx[n_ages])
    }
  } else {
    Lx[n_ages] <- 0
  }

  # --- Tx, ex ---
  Tx <- numeric(n_ages)
  Tx[n_ages] <- Lx[n_ages]
  for (i in (n_ages - 1):1) {
    Tx[i] <- Tx[i + 1] + Lx[i]
  }
  ex <- ifelse(lx > 0, Tx / lx, 0)

  data.frame(
    age = ages, n = n_vec, mx = mx, qx = qx, px = px,
    lx = lx, dx = dx, Lx = Lx, Tx = Tx, ex = ex
  )
}

#' Open interval person-years via Kannisto extrapolation
#' @keywords internal
.open_interval_kannisto <- function(mx, ages, lx_open) {
  # Kannisto: mu(x) = a*exp(b*x) / (1 + a*exp(b*x))
  # Fit from last 3 closed age groups
  n <- length(mx)
  if (n < 4) return(if (mx[n] > 0) lx_open / mx[n] else 0)

  idx <- (n - 3):(n - 1)
  x <- ages[idx]
  y <- mx[idx]
  y <- pmax(y, 1e-10)

  # Approximate: log(mu/(1-mu)) ~ log(a) + b*x  (logistic link)
  logit_y <- log(y / pmax(1 - y, 1e-10))
  fit <- tryCatch(stats::lm(logit_y ~ x), error = function(e) NULL)
  if (is.null(fit)) return(if (mx[n] > 0) lx_open / mx[n] else 0)

  a <- exp(stats::coef(fit)[1])
  b <- stats::coef(fit)[2]

  # Integrate survival from open age onward (numerical)
  x_open <- ages[n]
  ages_ext <- seq(x_open, x_open + 50, by = 1)
  mu <- a * exp(b * ages_ext) / (1 + a * exp(b * ages_ext))
  mu <- pmin(mu, 1)
  cum_haz <- cumsum(mu)
  surv <- exp(-cum_haz)
  lx_open * sum(surv)
}

#' Open interval person-years via Gompertz extrapolation
#' @keywords internal
.open_interval_gompertz <- function(mx, ages, lx_open) {
  n <- length(mx)
  if (n < 3) return(if (mx[n] > 0) lx_open / mx[n] else 0)

  # Gompertz: log(mu(x)) = alpha + beta*x
  idx <- (n - 2):(n - 1)
  x <- ages[idx]
  log_mx <- log(pmax(mx[idx], 1e-10))
  beta <- (log_mx[2] - log_mx[1]) / (x[2] - x[1])
  alpha <- log_mx[1] - beta * x[1]

  x_open <- ages[n]
  ages_ext <- seq(x_open, x_open + 50, by = 1)
  mu <- exp(alpha + beta * ages_ext)
  cum_haz <- cumsum(mu)
  surv <- exp(-cum_haz)
  lx_open * sum(surv)
}

#' Compute life expectancy
#'
#' Generic function for computing life expectancy from mortality data or
#' forecasts. Dispatches to methods for \code{mort_data} and
#' \code{mort_forecast} objects.
#'
#' @param object A \code{mort_data} or \code{mort_forecast} object.
#' @param ... Additional arguments passed to methods.
#' @return A \code{data.table} with life expectancy values.
#' @export
life_expectancy <- function(object, ...) {
  UseMethod("life_expectancy")
}

#' @describeIn life_expectancy Compute life expectancy for observed mortality data
#'
#' Computes period life expectancy at the youngest age for all groups and years
#' in a \code{mort_data} object.
#'
#' @param width Age-group width. Default 5.
#' @param open_age_method Method for open interval; see \code{\link{life_table}}.
#' @return A \code{data.table} with columns: \code{year}, \code{group},
#'   \code{ex} (life expectancy at the youngest observed age).
#' @export
life_expectancy.mort_data <- function(object, width = 5,
                                      open_age_method = "reciprocal", ...) {
  mx_arr <- object$mx
  ages   <- object$ages
  years  <- object$years
  groups <- object$groups
  dims   <- dim(mx_arr)

  .group_label <- function(groups, g) {
    if (is.data.frame(groups))
      paste(as.character(groups[g, ]), collapse = "_")
    else if (!is.null(groups))
      as.character(groups)[g]
    else
      paste0("Group_", g)
  }

  results <- vector("list", dims[2] * dims[3])
  k <- 1L
  for (g in seq_len(dims[3])) {
    for (t in seq_len(dims[2])) {
      lt <- life_table(mx_arr[, t, g], ages, width = width,
                       open_age_method = open_age_method)
      row <- data.table::data.table(
        year  = years[t],
        group = .group_label(groups, g),
        ex    = lt$ex[1]
      )
      if (is.data.frame(groups))
        for (col_name in names(groups))
          row[[col_name]] <- groups[[col_name]][g]
      results[[k]] <- row
      k <- k + 1L
    }
  }
  data.table::rbindlist(results)
}

#' @describeIn life_expectancy Compute life expectancy for mortality forecasts
#'
#' Computes period life expectancy from forecast mortality rates, including
#' prediction interval bounds if available.
#'
#' @return A \code{data.table} with columns: \code{year}, \code{group},
#'   \code{ex}, and optionally \code{ex_lower}, \code{ex_upper}.
#' @export
life_expectancy.mort_forecast <- function(object, width = 5,
                                           open_age_method = "reciprocal", ...) {
  mx_point <- object$mx_point
  ages     <- object$ages
  fc_years <- object$forecast_years
  groups   <- object$groups
  dims     <- dim(mx_point)
  has_pi   <- !is.null(object$mx_lower) && !is.null(object$mx_upper)

  .group_label <- function(groups, g) {
    if (is.data.frame(groups))
      paste(as.character(groups[g, ]), collapse = "_")
    else if (!is.null(groups))
      as.character(groups)[g]
    else
      paste0("Group_", g)
  }

  results <- vector("list", dims[2] * dims[3])
  k <- 1L
  for (g in seq_len(dims[3])) {
    for (t in seq_len(dims[2])) {
      lt <- life_table(mx_point[, t, g], ages, width = width,
                       open_age_method = open_age_method)
      row <- data.table::data.table(
        year  = fc_years[t],
        group = .group_label(groups, g),
        ex    = lt$ex[1]
      )
      if (has_pi) {
        lt_lo <- life_table(object$mx_upper[, t, g], ages, width = width,
                            open_age_method = open_age_method)
        lt_hi <- life_table(object$mx_lower[, t, g], ages, width = width,
                            open_age_method = open_age_method)
        row$ex_lower <- lt_lo$ex[1]  # higher mx -> lower ex
        row$ex_upper <- lt_hi$ex[1]  # lower mx -> higher ex
      }
      if (is.data.frame(groups))
        for (col_name in names(groups))
          row[[col_name]] <- groups[[col_name]][g]
      results[[k]] <- row
      k <- k + 1L
    }
  }
  data.table::rbindlist(results)
}

#' Inverse CLR transform to death distribution
#'
#' Converts CLR (centered log-ratio) coordinates back to a dx composition
#' using the softmax (inverse CLR) transform.
#'
#' @param clr_vec Numeric vector of CLR coordinates.
#' @return Numeric vector of dx values (positive, sums to 1).
#' @export
inv_clr_to_dx <- function(clr_vec) {
  e <- exp(clr_vec)
  dx <- e / sum(e)
  dx[dx <= 0] <- 1e-15
  dx / sum(dx)
}
