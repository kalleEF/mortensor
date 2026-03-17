#' Read mortality data from a long-format table
#'
#' The primary entry point for loading mortality data into \pkg{mortensor}.
#' Accepts any long-format `data.frame` or `data.table` containing mortality
#' rates (or death counts + population) with columns for year, age, and
#' optionally one or more stratification variables (e.g. sex, region,
#' deprivation).
#'
#' @section Expected input format:
#' A long-format table with one row per year x age x group combination.
#' Required columns (names configurable):
#' \describe{
#'   \item{Year column}{Integer calendar year (e.g. 2001, 2002).
#'     Specify via \code{year_col}.}
#'   \item{Age column}{Age-group identifier, either a label like \code{"30-34"}
#'     or an integer start age like \code{30}. Specify via \code{age_col}.}
#'   \item{Rate column or deaths + population columns}{Either a single
#'     column of pre-computed mortality rates (specify via \code{mx_col}), or
#'     two columns: one for death counts (\code{deaths_col}) and one for
#'     population exposure (\code{pop_col}), in which case
#'     \code{mx = deaths / population} is computed internally.}
#'   \item{Group columns (optional)}{Any additional columns that define strata
#'     (e.g. \code{"sex"}, \code{"region"}). Specify via \code{group_cols}. If
#'     \code{NULL}, \code{age_col} and \code{year_col} are auto-detected and
#'     all remaining non-value columns are treated as group columns.}
#' }
#'
#' @param dt A `data.frame` or `data.table` in long format.
#' @param mx_col Character name of the mortality-rate column. Either `mx_col`
#'   or both `deaths_col` and `pop_col` must be provided.
#' @param ages Integer vector of age-group start values (e.g.
#'   \code{c(30, 35, 40, 90)}), in the same order used to sort the age column.
#' @param years Integer vector of calendar years to keep (default: all in
#'   `dt`).
#' @param group_cols Character vector of column names defining groups/strata
#'   (e.g. `c("sex", "region")`). If `NULL`, auto-detected as all columns
#'   that are not `year_col`, `age_col`, `mx_col`, `deaths_col`, or
#'   `pop_col`.
#' @param pop_col Optional character name of a population-exposure column.
#'   Required when using `deaths_col`.
#' @param deaths_col Optional character name of a death-count column. If
#'   provided together with `pop_col`, mortality rates are computed as
#'   `deaths / population`.
#' @param cause Optional cause-of-death label (e.g. `"IHD"`).
#' @param age_labels Optional character vector of age-group labels (e.g.
#'   \code{c("30-34", "35-39")}). If \code{NULL}, labels are derived from
#'   \code{ages}.
#' @param year_col Character name of the year column (default `"year"`).
#' @param age_col Character name of the age column (default `"age_group"`).
#' @return A [mort_data] object.
#'
#' @examples
#' \dontrun{
#' # Minimal example - rates already computed, single stratum
#' dt <- data.table::data.table(
#'   year      = rep(2001:2020, each = 13),
#'   age_group = rep(seq(30, 90, 5), 20),
#'   mx        = runif(260, 0.001, 0.05)
#' )
#' md <- read_mortality(dt, mx_col = "mx", ages = seq(30, 90, 5))
#'
#' # With deaths + population, two strata
#' dt2 <- data.table::data.table(
#'   year   = rep(2001:2020, each = 26),
#'   age    = rep(rep(seq(30, 90, 5), each = 2), 20),
#'   sex    = rep(c("Male", "Female"), 260),
#'   deaths = rpois(520, 50),
#'   pop    = rpois(520, 10000)
#' )
#' md2 <- read_mortality(dt2, deaths_col = "deaths", pop_col = "pop",
#'                       ages = seq(30, 90, 5), age_col = "age",
#'                       group_cols = "sex")
#' }
#' @export
read_mortality <- function(dt, mx_col = NULL, ages,
                           years = NULL,
                           group_cols = NULL,
                           pop_col = NULL,
                           deaths_col = NULL,
                           cause = NULL,
                           age_labels = NULL,
                           year_col = "year",
                           age_col = "age_group") {
  if (!data.table::is.data.table(dt))
    dt <- data.table::as.data.table(dt)


  # --- Validate value columns ---
  if (is.null(mx_col) && is.null(deaths_col))
    stop("Provide either 'mx_col' (rate column) or 'deaths_col' + 'pop_col' ",
         "(count + exposure columns).")
  if (!is.null(deaths_col) && is.null(pop_col))
    stop("'pop_col' is required when 'deaths_col' is specified.")
  if (!is.null(mx_col) && !mx_col %in% names(dt))
    stop("Column '", mx_col, "' not found in data. Available: ",
         paste(names(dt), collapse = ", "))
  if (!is.null(deaths_col) && !deaths_col %in% names(dt))
    stop("Column '", deaths_col, "' not found in data.")
  if (!is.null(pop_col) && !pop_col %in% names(dt))
    stop("Column '", pop_col, "' not found in data.")
  if (!year_col %in% names(dt))
    stop("Year column '", year_col, "' not found in data.")
  if (!age_col %in% names(dt))
    stop("Age column '", age_col, "' not found in data.")

  # --- Standardise column names internally ---
  dt <- data.table::copy(dt)
  if (year_col != "year") {
    data.table::setnames(dt, year_col, "year")
    year_col <- "year"
  }
  if (age_col != "age_group") {
    data.table::setnames(dt, age_col, "age_group")
    age_col <- "age_group"
  }

  # --- Compute mx from deaths + pop if needed ---
  if (!is.null(deaths_col) && is.null(mx_col)) {
    mx_col <- ".mx_computed"
    dt[, (mx_col) := as.numeric(get(deaths_col)) /
                      pmax(as.numeric(get(pop_col)), 1)]
  }

  # --- Derive group_cols if not provided ---
  value_cols <- c("year", "age_group", mx_col, deaths_col, pop_col)
  value_cols <- value_cols[!is.null(value_cols)]
  if (is.null(group_cols)) {
    group_cols <- setdiff(names(dt), value_cols)
    if (length(group_cols) == 0L) {
      # No strata — add a dummy group
      dt[, .group := "all"]
      group_cols <- ".group"
    }
  }
  bad_gc <- setdiff(group_cols, names(dt))
  if (length(bad_gc) > 0)
    stop("Group column(s) not found in data: ", paste(bad_gc, collapse = ", "))

  # --- Subset to requested years ---
  if (is.null(years)) years <- sort(unique(dt$year))
  dt <- dt[year %in% years]

  # --- Build groups data.frame ---
  groups <- unique(dt[, group_cols, with = FALSE])
  data.table::setorderv(groups, group_cols)
  groups <- as.data.frame(groups)

  # --- Build 3D arrays ---
  mx_arr <- build_array(dt, mx_col, ages, years, groups)
  pop_arr <- if (!is.null(pop_col)) {
    build_array(dt, pop_col, ages, years, groups, floor = NULL)
  }

  mort_data(
    mx         = mx_arr,
    pop        = pop_arr,
    ages       = ages,
    years      = years,
    groups     = groups,
    cause      = cause,
    age_labels = age_labels
  )
}


#' Create a mort_data object from pre-built arrays
#'
#' Convenience wrapper around [mort_data()] for users who already have a 3D
#' mortality-rate array (e.g. from the Human Mortality Database, the
#' `demography` package, or any other source).
#'
#' @param mx 3D array `[n_ages x n_years x n_groups]` of mortality rates.
#' @param ages Integer vector of age-group start values.
#' @param years Integer vector of calendar years.
#' @param groups A data.frame with one row per group (3rd dimension of `mx`).
#'   Column names define the strata (e.g. `sex`, `region`).
#' @param pop Optional 3D array of population exposure (same dims as `mx`).
#' @param cause Optional cause-of-death label.
#' @param age_labels Optional character vector of age-group labels.
#' @return A [mort_data] object.
#' @export
read_mortality_from_arrays <- function(mx, ages, years, groups,
                                       pop = NULL, cause = NULL,
                                       age_labels = NULL) {
  if (!is.array(mx) || length(dim(mx)) != 3)
    stop("'mx' must be a 3D array [n_ages x n_years x n_groups].")
  if (!is.data.frame(groups)) groups <- as.data.frame(groups)
  if (length(ages) != dim(mx)[1])
    stop("length(ages) = ", length(ages), " but dim(mx)[1] = ", dim(mx)[1])
  if (length(years) != dim(mx)[2])
    stop("length(years) = ", length(years), " but dim(mx)[2] = ", dim(mx)[2])
  if (nrow(groups) != dim(mx)[3])
    stop("nrow(groups) = ", nrow(groups), " but dim(mx)[3] = ", dim(mx)[3])

  mort_data(
    mx     = mx,
    pop    = pop,
    ages   = ages,
    years  = years,
    groups = groups,
    cause  = cause,
    age_labels = age_labels
  )
}
