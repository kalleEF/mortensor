#' Build a 3D mortality array from a data.table
#'
#' Constructs a `[n_ages x n_years x n_groups]` array from long-format data.
#'
#' @param dt A `data.table` with columns for year, sex, quintile/deprivation,
#'   age_group, and the value column to extract.
#' @param value_col Character name of the column containing mortality rates.
#' @param ages Integer vector of age-group start values.
#' @param years Integer vector of years.
#' @param groups A data.frame with one row per group and columns matching
#'   stratification variables in `dt` (e.g. `sex`, `quintile`).
#' @param floor Numeric floor for zero/NA cells (default 1e-12).
#'   Set to NULL to skip floor replacement.
#' @return A 3D array `[n_ages x n_years x n_groups]`.
#' @export
build_array <- function(dt, value_col, ages, years, groups, floor = 1e-12) {
  n_a <- length(ages)
  n_t <- length(years)
  if (!is.data.frame(groups)) groups <- as.data.frame(groups)
  n_g <- nrow(groups)

  arr <- array(NA_real_, dim = c(n_a, n_t, n_g))
  group_cols <- names(groups)

  for (k in seq_len(n_g)) {
    for (j in seq_along(years)) {
      # Build matching condition
      mask <- dt[["year"]] == years[j]
      for (gc in group_cols)
        mask <- mask & dt[[gc]] == groups[[gc]][k]
      sub <- dt[mask]
      data.table::setorder(sub, age_group)
      if (nrow(sub) == n_a) arr[, j, k] <- sub[[value_col]]
    }
  }

  if (!is.null(floor) && grepl("mx", value_col, fixed = TRUE)) {
    bad <- is.na(arr) | arr <= 0
    n_zero <- sum(bad, na.rm = TRUE)
    arr[bad] <- floor
    if (n_zero > 0)
      message("build_array: replaced ", n_zero, " zero/NA cells with ", floor)
  }
  arr
}

#' Reconstruct a 3D array from Tucker3 factors
#'
#' Computes `X = A %*% G %*% kronecker(C, B)^T` where G is the unfolded
#' core tensor, then reshapes to `[n_a x n_t x n_g]`.
#'
#' @param A Mode-1 factor matrix (ages).
#' @param B Mode-2 factor matrix (years).
#' @param C Mode-3 factor matrix (groups).
#' @param G_mat Mode-1 unfolded core tensor.
#' @return 3D array `[nrow(A) x nrow(B) x nrow(C)]`.
#' @export
reconstruct_array <- function(A, B, C, G_mat) {
  X_mat <- A %*% G_mat %*% t(kronecker(C, B))
  array(X_mat, dim = c(nrow(A), nrow(B), nrow(C)))
}
