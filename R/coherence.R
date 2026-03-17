#' Enforce sex-differential coherence on forecasts
#'
#' Adjusts forecasted mortality rates to ensure Female <= Male within
#' each age-group/year/quintile combination. Uses observed sex ratios
#' as the adjustment target.
#'
#' @param fc_dt  data.table with columns: `age_group`, `year`, `sex`,
#'   `quintile`, `mx_point`, `mx_lo95`, `mx_hi95`, `mx_median`.
#' @param obs_dt data.table with column `mx_obs` for observed ratios.
#' @param sex_constraint Logical; enforce Female <= Male?
#' @return Modified `fc_dt` with added columns: `mx_adjusted`,
#'   `mx_lo95_adj`, `mx_hi95_adj`, `mx_median_adj`.
#' @export
enforce_coherence <- function(fc_dt, obs_dt, sex_constraint = TRUE) {
  fc_dt[, mx_adjusted := mx_point]
  data.table::setorder(fc_dt, cause, year, age_group, sex, quintile)
  all_quints <- sort(unique(fc_dt$quintile))

  if (sex_constraint) {
    for (iter in 1:5) {
      n_violations <- 0
      for (ag_val in unique(fc_dt$age_group)) {
        for (yr in unique(fc_dt$year)) {
          for (qq in all_quints) {
            idx_m <- fc_dt[age_group == ag_val & year == yr &
                            sex == "Male" & quintile == qq, which = TRUE]
            idx_f <- fc_dt[age_group == ag_val & year == yr &
                            sex == "Female" & quintile == qq, which = TRUE]
            if (length(idx_m) == 1 && length(idx_f) == 1) {
              mx_m <- fc_dt$mx_adjusted[idx_m]
              mx_f <- fc_dt$mx_adjusted[idx_f]
              if (mx_f >= mx_m) {
                n_violations <- n_violations + 1
                avg <- (mx_m + mx_f) / 2
                obs_r <- obs_dt[age_group == ag_val & quintile == qq,
                                .(r = mean(mx_obs[sex == "Female"], na.rm = TRUE) /
                                      mean(mx_obs[sex == "Male"], na.rm = TRUE))]$r
                obs_r <- min(max(obs_r, 0.02, na.rm = TRUE), 0.98, na.rm = TRUE)
                fc_dt$mx_adjusted[idx_f] <- avg * 2 * obs_r / (1 + obs_r)
                fc_dt$mx_adjusted[idx_m] <- avg * 2 / (1 + obs_r)
              }
            }
          }
        }
      }
      if (n_violations == 0) break
    }
  }

  # Scale PIs proportionally
  fc_dt[, ratio := mx_adjusted / pmax(mx_point, 1e-15)]
  fc_dt[, mx_lo95_adj   := mx_lo95   * ratio]
  fc_dt[, mx_hi95_adj   := mx_hi95   * ratio]
  fc_dt[, mx_median_adj := mx_median * ratio]
  fc_dt[, ratio := NULL]

  fc_dt
}
