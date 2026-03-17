#' Compare all candidate models on a train/test split
#'
#' Runs Tucker3-CoDa, Tucker3-Direct, Lee-Carter, FDM, StMoMo-APC, RWD,
#' ARIMA-perSeries, and ETS-Damped on the training data; evaluates on the
#' holdout test set.  Optionally builds an ensemble via inverse-MASE
#' weighting and NNLS stacking.
#'
#' @param mx_array  3D array \[n_a x n_t x n_g\] of (smoothed) mortality rates.
#' @param pop_array 3D population exposure array.
#' @param years     Integer vector of years.
#' @param ages      Integer vector of age-group starts.
#' @param groups    data.frame with `sex` and `quintile` columns.
#' @param train_idx Integer indices of training years.
#' @param test_idx  Integer indices of test years.
#' @param calib_coda  Output of [calibrate_tucker3()] with `method="coda"`.
#' @param calib_direct Output of [calibrate_tucker3()] with `method="direct"`.
#' @param width     Age-group width (default 5).
#' @param calib_all Output of [calibrate_all_models()] (NULL for defaults).
#' @param mx_eval_ref Optional raw-scale array for MASE evaluation.
#' @return Named list: `results` (data.table), `fc_arrays`,
#'   `ens_best_k`, `ens_best_names`, `ens_best_weights`.
#' @export
compare_models <- function(mx_array, pop_array, years, ages,
                           groups, train_idx, test_idx,
                           calib_coda, calib_direct, width = 5,
                           calib_all = NULL, mx_eval_ref = NULL) {
  n_test    <- length(test_idx)
  arr_train <- mx_array[, train_idx, , drop = FALSE]
  arr_test  <- mx_array[, test_idx,  , drop = FALSE]
  pop_train <- pop_array[, train_idx, , drop = FALSE]
  yrs_train <- years[train_idx]

  arr_test_eval  <- if (!is.null(mx_eval_ref))
    mx_eval_ref[, test_idx,  , drop = FALSE] else arr_test
  arr_train_eval <- if (!is.null(mx_eval_ref))
    mx_eval_ref[, train_idx, , drop = FALSE] else arr_train

  h_total  <- max(test_idx) - max(train_idx)
  fc_steps <- test_idx - max(train_idx)
  extract_fc <- function(fc_arr) fc_arr[, fc_steps, , drop = FALSE]

  HAS_STMOMO <- requireNamespace("StMoMo", quietly = TRUE)

  model_names <- c("Tucker3-CoDa", "Tucker3-Direct", "Lee-Carter", "FDM")
  if (HAS_STMOMO) model_names <- c(model_names, "StMoMo-APC")
  model_names <- c(model_names, "RWD", "ARIMA-perSeries", "ETS-Damped")

  results <- data.table::data.table(
    model = model_names,
    rmse = NA_real_, mae = NA_real_, mase = NA_real_,
    mase_trend = NA_real_, mase_combined = NA_real_,
    complete = TRUE)
  fc_arrays <- list()

  score_fc <- function(model_nm, fc_sub) {
    fc_arrays[[model_nm]] <<- fc_sub
    valid <- !is.na(as.vector(fc_sub))
    if (!any(valid)) return()
    results[model == model_nm,
            `:=`(rmse = compute_rmse_log(as.vector(arr_test_eval)[valid],
                                         as.vector(fc_sub)[valid]),
                 mae  = compute_mae(as.vector(arr_test_eval)[valid],
                                    as.vector(fc_sub)[valid]),
                 mase = compute_mase_array(arr_test_eval, fc_sub, arr_train_eval))]
  }

  # --- Tucker3-CoDa ---
  tryCatch({
    cp <- calib_coda
    fc <- tucker3_coda_forecast(arr_train, yrs_train, ages, h_total,
                                cp$R_A, cp$R_B, cp$R_C, cp$ts_method,
                                fast = TRUE, width = width)
    score_fc("Tucker3-CoDa", extract_fc(fc$mx_forecast))
  }, error = function(e) NULL)

  # --- Tucker3-Direct ---
  tryCatch({
    cp <- calib_direct
    fc <- tucker3_direct_forecast(arr_train, yrs_train, h_total,
                                  cp$R_A, cp$R_B, cp$R_C, cp$ts_method,
                                  fast = TRUE)
    score_fc("Tucker3-Direct", extract_fc(fc$mx_forecast))
  }, error = function(e) NULL)

  # --- Lee-Carter ---
  tryCatch({
    lc_adj <- if (!is.null(calib_all$lc)) calib_all$lc$adjust else "dt"
    lc_cp  <- if (!is.null(calib_all$lc)) calib_all$lc$chooseperiod else FALSE
    lc <- lca_forecast(arr_train, pop_train, yrs_train, ages, groups, h_total,
                       adjust = lc_adj, chooseperiod = lc_cp)
    score_fc("Lee-Carter", extract_fc(lc$mx_forecast))
  }, error = function(e) NULL)

  # --- FDM ---
  tryCatch({
    fdm_ord <- if (!is.null(calib_all$fdm)) calib_all$fdm$order else 6
    fdm_mth <- if (!is.null(calib_all$fdm)) calib_all$fdm$method else "arima"
    fdm_r <- fdm_forecast(arr_train, pop_train, yrs_train, ages, groups, h_total,
                          order = fdm_ord, method = fdm_mth)
    score_fc("FDM", extract_fc(fdm_r$mx_forecast))
  }, error = function(e) NULL)

  # --- StMoMo-APC ---
  if (HAS_STMOMO) {
    tryCatch({
      stm_mt <- if (!is.null(calib_all$stmomo)) calib_all$stmomo$model_type else "APC"
      sm <- stmomo_forecast(arr_train, pop_train, yrs_train, ages, groups,
                            h_total, model_type = stm_mt)
      if (!is.null(sm$mx_forecast))
        score_fc("StMoMo-APC", extract_fc(sm$mx_forecast))
    }, error = function(e) NULL)
  }

  # --- RWD ---
  tryCatch({
    dw <- if (!is.null(calib_all$rwd)) calib_all$rwd$drift_window else 0
    rwd_r <- rwd_forecast(arr_train, yrs_train, h_total, drift_window = dw)
    score_fc("RWD", extract_fc(rwd_r$mx_forecast))
  }, error = function(e) NULL)

  # --- ARIMA-perSeries ---
  tryCatch({
    ar_mp <- if (!is.null(calib_all$arima)) calib_all$arima$max.p else 2
    ar_mq <- if (!is.null(calib_all$arima)) calib_all$arima$max.q else 1
    ar_md <- if (!is.null(calib_all$arima)) calib_all$arima$max.d else 1
    ar_r <- arima_forecast(arr_train, yrs_train, h_total,
                           max.p = ar_mp, max.q = ar_mq, max.d = ar_md)
    score_fc("ARIMA-perSeries", extract_fc(ar_r$mx_forecast))
  }, error = function(e) NULL)

  # --- ETS-Damped ---
  tryCatch({
    ets_m <- if (!is.null(calib_all$ets)) calib_all$ets$ets_model else "AAN"
    ets_d <- if (!is.null(calib_all$ets)) calib_all$ets$damped else TRUE
    ets_r <- ets_forecast(arr_train, yrs_train, h_total,
                          ets_model = ets_m, damped = ets_d)
    score_fc("ETS-Damped", extract_fc(ets_r$mx_forecast))
  }, error = function(e) NULL)

  # Compute combined MASE & completeness flags
  for (mnm in names(fc_arrays)) {
    if (mnm == "Ensemble") next
    fa <- fc_arrays[[mnm]]
    mt <- compute_mase_group_trend(arr_test_eval, fa, arr_train_eval)
    mc <- compute_mase_combined(arr_test_eval, fa, arr_train_eval)
    has_na <- any(is.na(fa))
    results[model == mnm,
            `:=`(mase_trend = mt, mase_combined = mc, complete = !has_na)]
  }

  # ---- Ensemble: inverse-combined-MASE + NNLS stacking ----
  valid_models <- results[!is.na(mase_combined) & is.finite(mase_combined) &
                          complete == TRUE]
  best_ens_combined <- Inf; best_ens_mase <- Inf
  best_ens_arr <- NULL; best_ens_k <- 2L
  best_ens_names <- character(0); best_ens_weights <- NULL

  if (nrow(valid_models) >= 2) {
    for (try_k in 2:min(nrow(valid_models), 5)) {
      top_names_k   <- valid_models[order(mase_combined)][1:try_k, model]
      top_combined_k <- valid_models[order(mase_combined)][1:try_k, mase_combined]
      inv_comb_k    <- 1 / top_combined_k
      weights_k     <- inv_comb_k / sum(inv_comb_k)

      ens_arr_k <- array(0, dim = dim(arr_test_eval))
      w_used_k <- 0
      for (ii in seq_along(top_names_k)) {
        nm <- top_names_k[ii]
        if (!is.null(fc_arrays[[nm]])) {
          fa <- fc_arrays[[nm]]; fa[is.na(fa)] <- 0
          ens_arr_k <- ens_arr_k + weights_k[ii] * fa
          w_used_k <- w_used_k + weights_k[ii]
        }
      }
      if (w_used_k > 0) ens_arr_k <- ens_arr_k / w_used_k
      ens_combined_k <- compute_mase_combined(arr_test_eval, ens_arr_k, arr_train_eval)
      if (is.finite(ens_combined_k) && ens_combined_k < best_ens_combined) {
        best_ens_combined <- ens_combined_k
        best_ens_mase  <- compute_mase_array(arr_test_eval, ens_arr_k, arr_train_eval)
        best_ens_arr   <- ens_arr_k
        best_ens_k     <- try_k
        best_ens_names <- top_names_k
        best_ens_weights <- weights_k
      }
    }

    # NNLS stacking alternative
    if (nrow(valid_models) >= 3) {
      tryCatch({
        target_vec <- as.vector(log(pmax(arr_test_eval, 1e-12)))
        M <- do.call(cbind, lapply(valid_models$model, function(nm) {
          fa <- fc_arrays[[nm]]; fa[is.na(fa)] <- 1e-12
          as.vector(log(pmax(fa, 1e-12)))
        }))
        nnls_fit <- tryCatch(nnls::nnls(M, target_vec), error = function(e) NULL)
        if (is.null(nnls_fit)) {
          nmod <- ncol(M)
          nnls_obj <- function(w) sum((target_vec - M %*% pmax(w, 0))^2)
          opt <- stats::optim(rep(1 / nmod, nmod), nnls_obj, method = "L-BFGS-B",
                              lower = rep(0, nmod))
          stack_w_raw <- pmax(opt$par, 0)
        } else {
          stack_w_raw <- nnls_fit$x
        }
        active <- stack_w_raw > 1e-6
        if (sum(active) >= 2) {
          stack_w_raw[!active] <- 0
          stack_w <- stack_w_raw / sum(stack_w_raw)
          stack_names <- valid_models$model
          stack_arr_t <- array(0, dim = dim(arr_test_eval))
          for (ii in seq_along(stack_names)) {
            if (stack_w[ii] > 0) {
              fa <- fc_arrays[[stack_names[ii]]]; fa[is.na(fa)] <- 0
              stack_arr_t <- stack_arr_t + stack_w[ii] * fa
            }
          }
          stack_comb <- compute_mase_combined(arr_test_eval, stack_arr_t, arr_train_eval)
          if (is.finite(stack_comb) && stack_comb < best_ens_combined) {
            best_ens_combined <- stack_comb
            best_ens_mase  <- compute_mase_array(arr_test_eval, stack_arr_t, arr_train_eval)
            best_ens_arr   <- stack_arr_t
            best_ens_k     <- sum(active)
            best_ens_names <- stack_names[active]
            best_ens_weights <- stack_w[active]
          }
        }
      }, error = function(e) NULL)
    }

    fc_arrays[["Ensemble"]] <- best_ens_arr
    ens_trend <- compute_mase_group_trend(arr_test_eval, best_ens_arr, arr_train_eval)
    ens_row <- data.table::data.table(
      model = "Ensemble", rmse = NA_real_, mae = NA_real_,
      mase = best_ens_mase, mase_trend = ens_trend,
      mase_combined = best_ens_combined, complete = TRUE)
    results <- data.table::rbindlist(list(results, ens_row), fill = TRUE)
  }

  results[, mase_rank := data.table::frank(mase_combined,
                                           ties.method = "min", na.last = TRUE)]
  list(results = results, fc_arrays = fc_arrays,
       ens_best_k = best_ens_k,
       ens_best_names = best_ens_names,
       ens_best_weights = best_ens_weights)
}
