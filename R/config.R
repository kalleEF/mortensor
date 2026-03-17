#' Default disease columns (ONS UK data)
#'
#' Returns the 19 raw death column names that map to the 9 named disease groups
#' in the UK ONS data. Used for computing the residual ("other causes").
#' These are specific to the ONS data format and are used by the ONS helper
#' functions.
#'
#' @return Character vector of column names.
#' @export
default_disease_cols <- function() {
  c("ihd",
    "breast_cancer", "colorectal", "lung_cancer", "prostate", "other_cancer",
    "dementia",
    "stroke",
    "copd",
    "afib", "heart_failure", "hypertension",
    "t1_diabetes", "t2_diabetes", "diabetes_other",
    "alc_liver", "nonalc_liver",
    "alc_specific", "suicide")
}

#' Default disease ordering for coherent ILR composition
#'
#' Defines the canonical ordering of disease groups for the ILR basis in
#' Tucker4-ILR and WLS reconciliation. The "residual" category is last.
#'
#' @return Character vector of cause short names.
#' @export
default_disease_order <- function() {
  c("ihd", "cancer", "dementia", "stroke", "copd",
    "cvd_other", "diabetes", "liver", "alc_suicide",
    "residual")
}

#' Create default cause-specific configurations
#'
#' Returns a named list of `mort_config` objects for the standard UK disease
#' groups (all-causes, 9 diseases, and residual).
#'
#' @param ... Additional overrides passed to each `mort_config()`.
#' @return Named list of `mort_config` objects.
#' @export
default_configs <- function(...) {
  list(
    all = mort_config(
      cause_cols = "all_causes",
      cause_label = "All Causes",
      cause_short = "all",
      is_cause_specific = FALSE, ...
    ),
    ihd = mort_config(
      cause_cols = "ihd",
      cause_label = "Chronic IHD",
      cause_short = "ihd", ...
    ),
    cancer = mort_config(
      cause_cols = c("breast_cancer", "colorectal", "lung_cancer",
                     "prostate", "other_cancer"),
      cause_label = "Cancer (all)",
      cause_short = "cancer", ...
    ),
    dementia = mort_config(
      cause_cols = "dementia",
      cause_label = "Dementia",
      cause_short = "dementia", ...
    ),
    stroke = mort_config(
      cause_cols = "stroke",
      cause_label = "Stroke",
      cause_short = "stroke", ...
    ),
    copd = mort_config(
      cause_cols = "copd",
      cause_label = "COPD",
      cause_short = "copd", ...
    ),
    cvd_other = mort_config(
      cause_cols = c("afib", "heart_failure", "hypertension"),
      cause_label = "Cardiovascular Other",
      cause_short = "cvd_other", ...
    ),
    diabetes = mort_config(
      cause_cols = c("t1_diabetes", "t2_diabetes", "diabetes_other"),
      cause_label = "Diabetes (all)",
      cause_short = "diabetes", ...
    ),
    liver = mort_config(
      cause_cols = c("alc_liver", "nonalc_liver"),
      cause_label = "Liver Disease",
      cause_short = "liver", ...
    ),
    alc_suicide = mort_config(
      cause_cols = c("alc_specific", "suicide"),
      cause_label = "Alcohol & Suicide",
      cause_short = "alc_suicide", ...
    ),
    residual = mort_config(
      cause_cols = "RESIDUAL",
      cause_label = "Other Causes",
      cause_short = "residual", ...
    )
  )
}
