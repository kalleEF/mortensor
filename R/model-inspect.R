#' Parameters of a fitted model
#'
#' Returns the estimated parameters of a fit in one long table, for every
#' model type that has interpretable parameters:
#'
#' * [lee_carter()]: `a` (age pattern), `b` (age pattern of change), `k`
#'   (period index), and the drift and innovation standard deviation of `k`;
#'   with `rotate = TRUE` also the ultimate pattern `b_ultimate`; with
#'   `family = "nb"` also the dispersion `nb_dispersion` (variance `mu + phi
#'   mu^2`);
#' * [li_lee()], [product_ratio()]: the common factor of each group (`B`,
#'   `K`, its drift and standard deviation), and per stratum the deviations
#'   `a`, `b`, `k` and their AR(1) coefficient `phi` and standard deviation;
#' * [coda()]: the centre `a` (geometric mean distribution of life-table
#'   deaths), the loadings `b`, the index `k` with the order (`order_k`),
#'   coefficients (`arima_k`) and innovation standard deviation (`sd_k`) of
#'   its ARIMA model, the share of variation explained by the component, the
#'   parameters of the open age interval (`a_open`, `b_open`) or of the
#'   Kannisto model of each year (`kannisto_a`, `kannisto_b`); the coherent
#'   model gives per group of `common_by` the common `a`, `B`, `K` (with
#'   `order_K`, `arima_K`, `sd_K`) and per stratum `a`, `b`, `k` and the
#'   ARIMA model of `k`. [coda_causes()]: per group, as [coda()], with parts
#'   by age and cause (`b1`, `k1`, `b2`, `k2`, ... with several components);
#' * [coda_multiway()]: the common age patterns `b`, time indexes `k` (with
#'   `order_k`, `arima_k`, `sd_k`), population markers `gamma` and core `g`
#'   (index `p,q,r`), numbered with several components, and the share of
#'   variation explained; per stratum the centre `a`, the parameters of the
#'   open age interval or of the Kannisto model of each year;
#' * [rwd()]: drift and innovation standard deviation per age (and, with
#'   `observation_noise = TRUE`, the variances of level and noise);
#' * [apc_bayes()]: posterior medians and 90% intervals of the age, period
#'   and cohort effects (logit scale, bamp's identification), and the
#'   Gelman-Rubin statistics (`rhat`, values below 1.1 indicate convergence)
#'   of the identified quantities;
#' * [p_splines()], [cp_splines()]: the smoothing parameters over age
#'   (`lambda_age`) and years (`lambda_year`), the effective dimension
#'   (`ed`), the criterion (`bic` or `aic`) and the dispersion (deviance
#'   divided by the number of cells minus the effective dimension); for
#'   [cp_splines()] also the bounds of the rates of ageing and of change per
#'   age (`ageing_lower`, `ageing_upper`, `change_lower`, `change_upper`).
#' * [saint()]: per trend (`"trend <group>"`) the indexes `alpha`, `beta`,
#'   `kappa`, `zeta` by year, the frailty variance `sigma2`, the
#'   pseudo-log-likelihood `loglik`, the `drift` and covariance `cov` of
#'   `alpha` and `beta` (or, with `ecm`, the parameters of the
#'   error-correction model: `a11`, `a21`, `a32`, `a42`, `C_D`, `C_delta`,
#'   `lambda`, `Omega`) and the standard deviations `sd_kappa`, `sd_zeta`;
#'   per stratum with a spread the spread parameters `y1` ... `y5` by year
#'   and their VAR(1) (`ar` or `A`, and `Omega`);
#' * [bohk_rau()]: posterior medians and 90% intervals of the curve
#'   parameters `b1` and `b2` by age and of the variance and hyperparameters,
#'   the least-squares values `theta1`, `theta2` (exponential core model) and
#'   the Gelman-Rubin statistics (`rhat`);
#' * [apci()]: `alpha`, `beta` (by age), `kappa` (by year), `gamma` (by
#'   year of birth; 0 for the cohorts with four or fewer observations), the
#'   mean year `ybar`, the smoothing parameters (`lambda_alpha`,
#'   `lambda_beta`), the effective dimension `ed` and the `bic`, and the
#'   ARIMA models of `kappa` and `gamma` (`order_kappa`, `arima_kappa`,
#'   `sd_kappa`, `order_gamma`, `arima_gamma`, `sd_gamma`);
#' * [nordpred_apc()]: estimates and 90% Wald intervals of the age, period
#'   and cohort effects (cohorts by central year of birth) and of the drift,
#'   the projected drift, the number of periods fitted, the choice of the
#'   recent trend and the p-values of the goodness of fit and of the test for
#'   the recent trend;
#' * [cae()]: the common age effects `beta1`, `beta2`, ... of each group
#'   (rows labelled `common <g> (<strata>)`), and per stratum the centring
#'   `a`, the period effects `kappa1`, `kappa2`, ..., their drift
#'   (`drift_kappa1`, ...), the standard deviation of their innovations
#'   (`sd_kappa1`, ...) and the correlation of their innovations with those
#'   of the same order of every stratum (`cor_kappa1`, ..., indexed by the
#'   other stratum); rows of stratum `"model"` give the mean squared error
#'   of the centred log rates (`mse`) and the approximate BIC of Kleinow
#'   (2015, eq. 5; `bic`).
#'
#' Fits of wrapped packages (StMoMo, demography) keep the original fitted
#' object in `fit$fits[[i]]`, where the methods of those packages apply.
#'
#' Wrappers return the parameters of the model they wrap ([recalibrated()]
#' adds its factors; [coupled()] adds those of the template, with column
#' `component` `"base"` or `"template"`). For an [ensemble()], the rows of
#' each component carry its name in column `component`, and rows with
#' `parameter = "weight"` give the weights; components without
#' interpretable parameters (wrapped packages) contribute their weight
#' only. For [stratified()] and
#' `tuned(by = )`, the rows of each group carry its name in column `group`.
#'
#' @param fit An `mt_fit` object.
#' @return A `data.table` with `stratum`, `parameter`, `index` (age label,
#'   year or position), `value` and, for Bayesian fits, `lower` and `upper`;
#'   for ensembles also `component`, for fits by groups also `group`.
#' @seealso [summary.mt_fit()], [residuals.mt_fit()].
#' @examples
#' x <- mt_example(deaths = "single")
#' p <- mt_parameters(mt_fit(lee_carter(), x))
#' head(p)
#' @export
mt_parameters <- function(fit) {
  check_class(fit, "mt_fit", "fit")
  out <- if (!is.null(fit$components)) {
    ensemble_parameters(fit)
  } else if (!is.null(fit$group_fits)) {
    group_parameters(fit)
  } else {
    sp <- fit$spec
    # wrappers that keep the fit of their base model (the fitted base specs are
    # stacked in inner_specs; tuned() fits the spec it chose)
    inner <- fit$inner_specs
    while (inherits(sp, c("mt_spec_recalibrated", "mt_spec_reconciled", "mt_spec_bootstrap", "mt_spec_anchored",
                          "mt_spec_coupled"))) {
      sp <- if (length(inner)) inner[[1]] else sp$params$base
      inner <- inner[-1]
    }
    base_parameters(fit, sp)
  }
  tpl <- fit$coupling$template
  if (!is.null(tpl)) {
    # coupled(): the parameters of the template, which gives the dependence
    tp <- tryCatch(mt_parameters(tpl), mt_error_argument = function(e) NULL)
    out <- copy(out)
    out[, component := "base"]
    if (!is.null(tp)) out <- rbind(out, tp[, component := "template"], fill = TRUE)
  }
  cal <- fit$calibration
  if (!is.null(cal)) out <- rbind(out, calibration_parameters(cal), fill = TRUE)
  out
}

# Rows of mt_parameters() for the calibration of recalibrated(): factors (one
# row per group and band, or per horizon with index the horizon) or the
# tuning factors of the error method
calibration_parameters <- function(cal) {
  if (!"group" %in% names(cal)) {
    return(data.table(stratum = "calibration: past errors", parameter = "tuning", index = as.character(cal$horizon),
                      value = cal$tuning))
  }
  lab <- if (!is.null(cal$band)) paste0(cal$group, ", ", cal$band) else cal$group
  if ("horizon" %in% names(cal)) {
    return(data.table(stratum = rep(paste("calibration:", lab), 2),
                      parameter = rep(c("scale", "drift"), each = nrow(cal)),
                      index = rep(as.character(cal$horizon), 2), value = c(cal$scale, cal$drift)))
  }
  data.table(stratum = rep(paste("calibration:", lab), 3),
             parameter = rep(c("scale_h1", "slope", "drift"), each = nrow(cal)), index = "",
             value = c(cal$scale_h1, cal$slope, cal$drift))
}

# Parameters of each component of an ensemble (column `component`) and the
# weights; components without a method contribute their weight only
ensemble_parameters <- function(fit) {
  comps <- fit$components
  w <- fit$weights %||% stats::setNames(rep(1 / length(comps), length(comps)), names(comps))
  parts <- lapply(names(comps), function(n) {
    p <- tryCatch(mt_parameters(comps[[n]]), mt_error_argument = function(e) NULL)
    if (!is.null(p)) p[, component := n]
    p
  })
  weights <- data.table(stratum = "ensemble", parameter = "weight", index = names(w), value = unname(w),
                        component = names(w))
  rbindlist(c(parts, list(weights)), fill = TRUE)
}

# Parameters of each group of a fit by groups (stratified(), tuned(by = )),
# with the group's name in column `group`
group_parameters <- function(fit) {
  nm <- names(fit$groups) %||% as.character(seq_along(fit$groups))
  rbindlist(lapply(seq_along(fit$group_fits), function(g) {
    p <- mt_parameters(fit$group_fits[[g]])
    p[, group := nm[g]]
    p
  }), fill = TRUE)
}

base_parameters <- function(fit, sp) {
  lab <- fit$strata$label
  ages <- fit$age$labels
  yrs <- as.character(fit$periods)
  row <- function(stratum, parameter, index, value) {
    data.table(stratum = stratum, parameter = parameter, index = as.character(index), value = as.numeric(value))
  }
  if (inherits(sp, "mt_spec_lee_carter")) {
    return(rbindlist(lapply(seq_along(fit$fits), function(s) {
      f <- fit$fits[[s]]
      if (is.null(f)) return(NULL)
      rbind(row(lab[s], "a", ages, f$a), row(lab[s], "b", ages, f$b), row(lab[s], "k", yrs, f$k),
            row(lab[s], "drift_k", "", f$drift), row(lab[s], "sd_k", "", f$sigma_k),
            if (!is.null(f$dispersion)) row(lab[s], "nb_dispersion", "", f$dispersion),
            if (!is.null(f$b_ultimate)) row(lab[s], "b_ultimate", ages, f$b_ultimate))
    })))
  }
  if (inherits(sp, c("mt_spec_coda", "mt_spec_coda_causes"))) return(coda_parameters(fit, sp))
  if (inherits(sp, "mt_spec_coda_multiway")) return(coda_multiway_parameters(fit))
  if (inherits(sp, c("mt_spec_p_splines", "mt_spec_cp_splines"))) {
    return(rbindlist(lapply(seq_along(fit$fits), function(s) {
      f <- fit$fits[[s]]
      if (is.null(f)) return(NULL)
      b <- f$bounds
      rbind(row(lab[s], "lambda_age", "", f$lambda[1]), row(lab[s], "lambda_year", "", f$lambda[2]),
            row(lab[s], "ed", "", f$ed), row(lab[s], f$criterion, "", f$ic),
            row(lab[s], "dispersion", "", f$dispersion),
            if (!is.null(b)) {
              rbindlist(lapply(c("ageing_lower", "ageing_upper", "change_lower", "change_upper"),
                               function(nm) row(lab[s], nm, ages, b[[nm]])))
            })
    })))
  }
  if (inherits(sp, "mt_spec_rwd")) {
    return(rbindlist(lapply(seq_along(fit$fits), function(s) {
      f <- fit$fits[[s]]
      if (is.null(f)) return(NULL)
      out <- rbind(row(lab[s], "drift", ages, f$mu), row(lab[s], "sd", ages, f$sigma))
      if (!is.null(f$s2eta)) {
        out <- rbind(out, row(lab[s], "var_level", ages, f$s2eta), row(lab[s], "var_noise", ages, f$s2e))
      }
      out
    })))
  }
  if (inherits(sp, "mt_spec_multipop")) {
    f <- fit$fits
    groups <- vapply(seq_along(f$commons), function(g) paste(lab[f$group == g], collapse = ", "), "")
    common <- rbindlist(lapply(seq_along(f$commons), function(g) {
      cm <- f$commons[[g]]
      gl <- sprintf("common %d (%s)", g, groups[g])
      rbind(if (!is.null(cm$a) && f$base_type != "li_lee") row(gl, "A", ages, cm$a),
            row(gl, "B", ages, cm$b), row(gl, "K", yrs, cm$k),
            row(gl, "drift_K", "", cm$drift), row(gl, "sd_K", "", cm$sigma_k))
    }))
    own <- rbindlist(lapply(seq_along(f$strata), function(s) {
      st <- f$strata[[s]]
      r <- ncol(st$b)
      rbind(row(lab[s], "a", ages, st$a),
            rbindlist(lapply(seq_len(r), function(j) {
              rbind(row(lab[s], paste0("b", j), ages, st$b[, j]), row(lab[s], paste0("k", j), yrs, st$k[, j]),
                    row(lab[s], paste0("phi", j), "", st$ar[[j]]$phi),
                    row(lab[s], paste0("sd_k", j), "", st$ar[[j]]$sigma))
            })))
    }))
    return(rbind(common, own))
  }
  if (inherits(sp, "mt_spec_cae")) {
    f <- fit$fits
    p <- dim(f$kappa)[1]
    ns <- length(lab)
    sdv <- sqrt(diag(f$sigma))
    cr <- f$sigma / outer(sdv, sdv)
    common <- rbindlist(lapply(seq_along(f$beta), function(g) {
      gl <- sprintf("common %d (%s)", g, paste(lab[f$group == g], collapse = ", "))
      rbindlist(lapply(seq_len(p), function(j) row(gl, paste0("beta", j), ages, f$beta[[g]][, j])))
    }))
    own <- rbindlist(lapply(seq_len(ns), function(s) {
      rbind(row(lab[s], "a", ages, f$a[, s]),
            rbindlist(lapply(seq_len(p), function(j) {
              i <- (s - 1L) * p + j
              same <- (seq_len(ns) - 1L) * p + j
              rbind(row(lab[s], paste0("kappa", j), yrs, f$kappa[j, , s]),
                    row(lab[s], paste0("drift_kappa", j), "", f$drift[i]),
                    row(lab[s], paste0("sd_kappa", j), "", sdv[i]),
                    row(lab[s], paste0("cor_kappa", j), lab, cr[i, same]))
            })))
    }))
    return(rbind(common, own, row("model", c("mse", "bic"), "", c(f$mse, f$bic))))
  }
  if (inherits(sp, "mt_spec_saint")) return(saint_parameters(fit))
  if (inherits(sp, "mt_spec_bohk_rau")) return(bohk_rau_parameters(fit))
  if (inherits(sp, "mt_spec_apci")) return(apci_parameters(fit))
  if (inherits(sp, "mt_spec_nordpred")) return(nordpred_parameters(fit))
  if (inherits(sp, "mt_spec_apc_bayes")) {
    e <- apc_effects(fit)
    e <- e[, list(stratum, parameter = effect, index = as.character(index), value = q0.5, lower = q0.05,
                  upper = q0.95)]
    r <- rbindlist(lapply(seq_along(fit$fits), function(s) {
      ps <- fit$fits[[s]]$psrf
      if (is.null(ps)) return(NULL)
      data.table(stratum = lab[s], parameter = "rhat", index = names(ps), value = unname(ps), lower = NA_real_,
                 upper = NA_real_)
    }))
    return(rbind(e, r))
  }
  mt_abort(sprintf(paste("mt_parameters() has no method for %s; the fitted objects are in `fit$fits` (for",
                         "wrapped packages, use their own methods)."), sp$label), class = "argument")
}

# bohk_rau(): posterior medians and 90% intervals of the curve parameters by
# age and of the scalar parameters; theta of the exponential core model; the
# Gelman-Rubin statistics (value only)
bohk_rau_parameters <- function(fit) {
  lab <- fit$strata$label
  ages <- fit$age$labels
  q <- function(x) {
    if (is.matrix(x)) t(apply(x, 1, stats::quantile, c(0.5, 0.05, 0.95))) else
      matrix(stats::quantile(x, c(0.5, 0.05, 0.95)), 1)
  }
  row <- function(s, parameter, index, x) {
    v <- q(x)
    data.table(stratum = s, parameter = parameter, index = as.character(index), value = v[, 1], lower = v[, 2],
               upper = v[, 3])
  }
  rbindlist(lapply(seq_along(fit$fits), function(s) {
    f <- fit$fits[[s]]
    if (is.null(f)) return(NULL)
    smp <- f$samples
    scal <- intersect(c("sigma", "sigma1", "sigma2", "mu1", "mu2", "omega1", "omega2", "r"), names(smp))
    rbind(row(lab[s], "b1", ages, smp$b1), row(lab[s], "b2", ages, smp$b2),
          rbindlist(lapply(scal, function(v) row(lab[s], v, "", smp[[v]]))),
          if (!is.null(smp$theta)) {
            data.table(stratum = lab[s], parameter = rep(c("theta1", "theta2"), each = length(ages)),
                       index = rep(ages, 2), value = as.vector(smp$theta), lower = NA_real_, upper = NA_real_)
          },
          data.table(stratum = lab[s], parameter = "rhat", index = names(smp$rhat), value = unname(smp$rhat),
                     lower = NA_real_, upper = NA_real_))
  }))
}

# nordpred_apc(): estimates with 90% Wald intervals of the age, period and
# cohort effects (cohorts labelled by the central year of birth) and the
# drift, the projected drift, and the choices of the fit
nordpred_parameters <- function(fit) {
  lab <- fit$strata$label
  rbindlist(lapply(seq_along(fit$fits), function(s) {
    f <- fit$fits[[s]]
    if (is.null(f)) return(NULL)
    se <- sqrt(diag(f$vcov))
    nm <- names(f$coef)
    na <- f$n_age
    np <- f$n_periods
    eff <- function(prefix, levels, index, parameter) {
      v <- numeric(length(levels))
      e <- numeric(length(levels))
      i <- match(paste0(prefix, levels), nm)
      v[!is.na(i)] <- f$coef[i[!is.na(i)]]
      e[!is.na(i)] <- se[i[!is.na(i)]]
      data.table(stratum = lab[s], parameter = parameter, index = as.character(index), value = v,
                 lower = v - 1.645 * e, upper = v + 1.645 * e)
    }
    # central year of birth of cohort k: the cells of the oldest age group in
    # the first period fitted (k = 1), then 5 years per cohort
    k <- seq_len(na + np - 1L)
    birth <- f$first_years[1] - f$age_lower[f$rows][na] + 5 * (k - 1L)
    i_d <- match("drift", nm)
    drift <- f$coef[[i_d]]
    proj <- if (f$recent) drift - f$coef[[match(paste0("period", np - 1L), nm)]] else drift
    rbind(eff("age", seq_len(na), f$age_labels[f$rows], "age"),
          data.table(stratum = lab[s], parameter = "drift", index = "", value = drift,
                     lower = drift - 1.645 * se[[i_d]], upper = drift + 1.645 * se[[i_d]]),
          eff("period", seq_len(np), f$period_labels, "period"),
          eff("cohort", k, birth, "cohort"),
          data.table(stratum = lab[s], parameter = c("drift_projected", "n_periods", "recent", "gof_p_value",
                                                    "recent_p_value"), index = "",
                     value = c(proj, np, as.numeric(f$recent), f$gof_p, f$p_recent), lower = NA_real_,
                     upper = NA_real_))
  }))
}

#' Fitted rates, residuals and summaries of fitted models
#'
#' `fitted()` returns the in-sample fitted death rates as an array `[age,
#' period, stratum]` (`NA` for models without fitted values, such as the
#' plain random walk). `residuals()` returns deviance residuals (default) or
#' Pearson residuals of the death counts, assuming Poisson deaths given the
#' fitted rates, as a long table; values far beyond 2 in absolute value
#' point to structure that the model misses (for example cohort effects) or
#' to overdispersion. `summary()` describes the model, its sources of
#' uncertainty and the fit; [mt_plot_residuals()] shows the residuals as an
#' age-period heat map.
#'
#' @param object An `mt_fit` object.
#' @param type `"deviance"` or `"pearson"`.
#' @param ... Not used.
#' @return `fitted()`: an array; `residuals()`: a `data.table` with
#'   `stratum`, `age`, `period`, `deaths`, `fitted_deaths`, `residual`;
#'   `summary()`: an object of class `summary.mt_fit`, printed.
#' @examples
#' x <- mt_example(deaths = "single")
#' fit <- mt_fit(lee_carter(), x)
#' summary(fit)
#' head(residuals(fit))
#' @name fit-inspection
NULL

#' @rdname fit-inspection
#' @export
fitted.mt_fit <- function(object, ...) {
  d <- dim(object$train_D)
  out <- array(NA_real_, d, dimnames(object$train_D))
  if (inherits(object, "mt_fit_ensemble")) {
    # weighted geometric mean of the components' fitted rates
    acc <- array(0, d)
    ws <- array(0, d)
    for (n in names(object$components)) {
      f <- tryCatch(stats::fitted(object$components[[n]]), error = function(e) NULL)
      if (is.null(f) || !identical(dim(f), d)) next
      ok <- !is.na(f) & f > 0
      acc[ok] <- acc[ok] + object$weights[[n]] * log(f[ok])
      ws[ok] <- ws[ok] + object$weights[[n]]
    }
    out[ws > 0] <- exp(acc[ws > 0] / ws[ws > 0])
    return(out)
  }
  fits <- object$fits
  if (is.list(fits) && "fitted" %in% names(fits)) {
    if (identical(dim(fits$fitted), d)) out[] <- fits$fitted
  } else if (is.list(fits)) {
    for (s in seq_len(min(length(fits), d[3]))) {
      f <- if (is.list(fits[[s]])) fits[[s]]$fitted
      if (!is.null(f) && length(f) == d[1] * d[2]) out[, , s] <- f
    }
  }
  out
}

#' @rdname fit-inspection
#' @export
residuals.mt_fit <- function(object, type = c("deviance", "pearson"), ...) {
  type <- match.arg(type)
  m <- stats::fitted(object)
  D <- object$raw_D %||% object$train_D
  mu <- m * object$train_E
  r <- if (type == "pearson") {
    (D - mu) / sqrt(mu)
  } else {
    sign(D - mu) * sqrt(2 * pmax(ifelse(D > 0, D * log(D / mu), 0) - (D - mu), 0))
  }
  dn <- dimnames(D)
  out <- data.table(stratum = rep(dn[[3]], each = length(dn[[1]]) * length(dn[[2]])),
                    age = rep(dn[[1]], length(dn[[2]]) * length(dn[[3]])),
                    period = rep(rep(as.numeric(dn[[2]]), each = length(dn[[1]])), length(dn[[3]])),
                    deaths = as.vector(D), fitted_deaths = as.vector(mu), residual = as.vector(r))
  out[!is.na(residual)]
}

#' @rdname fit-inspection
#' @export
summary.mt_fit <- function(object, ...) {
  res <- tryCatch(stats::residuals(object), error = function(e) NULL)
  structure(list(label = object$spec$label, assumptions = model_assumptions(object$spec),
                 strata = nrow(object$strata), ages = object$age$labels, periods = range(object$periods),
                 prep = object$prep, diagnostics = object$diagnostics,
                 deviance = if (!is.null(res) && nrow(res)) sum(res$residual^2) else NA_real_,
                 cells = if (!is.null(res)) nrow(res) else NA_integer_,
                 large_residuals = if (!is.null(res) && nrow(res)) mean(abs(res$residual) > 2) else NA_real_,
                 calibration = object$calibration, tuning = object$tuning_choice,
                 anchor = if (!is.null(object$anchor)) stats::quantile(object$anchor$offset, c(0, 0.5, 1))),
            class = "summary.mt_fit")
}

#' @export
print.summary.mt_fit <- function(x, ...) {
  cat(sprintf("%s\n", x$label))
  cat(sprintf("  data: %d strata; ages %s to %s; periods %d-%d\n", x$strata, x$ages[1], x$ages[length(x$ages)],
              x$periods[1], x$periods[2]))
  if (!is.null(x$prep)) {
    cat(sprintf("  preprocessing: %s\n", paste(vapply(x$prep$steps, function(s) paste0("step_", s$type), ""),
                                               collapse = ", ")))
    over_time <- vapply(x$prep$steps, function(s) {
      identical(s$type, "smooth") && (identical(s$args$along, "period") || identical(s$args$dims, 2))
    }, logical(1))
    if (any(over_time)) {
      cat("    smoothing over time narrows forecast intervals; check coverage in a backtest (recalibrated())\n")
    }
  }
  cat("  model and uncertainty:\n")
  cat(paste0("    ", strwrap(x$assumptions, 76), collapse = "\n"), "\n", sep = "")
  if (is.finite(x$deviance)) {
    cat(sprintf(paste("  fit: Poisson deviance %.1f on %d cells (%.2f per cell); %.1f%% of deviance residuals",
                      "beyond +/-2\n"), x$deviance, x$cells, x$deviance / x$cells, 100 * x$large_residuals))
  }
  if (!is.null(x$anchor)) {
    cat(sprintf("  jump-off offsets (log ratio observed / fitted): min %.3f, median %.3f, max %.3f\n", x$anchor[1],
                x$anchor[2], x$anchor[3]))
  }
  if (!is.null(x$tuning)) {
    tu <- as.data.frame(x$tuning)
    cat("  tuning (chosen settings and inner score):\n")
    cat(paste0("    ", utils::capture.output(print(tu, row.names = FALSE, digits = 4)), collapse = "\n"), "\n",
        sep = "")
  }
  if (!is.null(x$calibration) && !"group" %in% names(x$calibration)) {
    cal <- x$calibration
    cat(sprintf(paste("  calibration: resampled errors of inner forecasts at horizons 1-%d (%s inner origins per",
                      "horizon); tuning factor %s\n"), nrow(cal), paste(unique(range(cal$origins)), collapse = "-"),
                paste(unique(sprintf("%.2f", range(cal$tuning))), collapse = "-")))
  } else if (!is.null(x$calibration) && "horizon" %in% names(x$calibration)) {
    cal <- x$calibration
    for (k in unique(paste(cal$group, cal$band, sep = "\r"))) {
      ck <- cal[paste(group, band, sep = "\r") == k]
      grp <- if (ck$band[1] != "all ages") paste0(ck$group[1], ", ages ", ck$band[1]) else ck$group[1]
      dr <- if (any(ck$drift != 0)) sprintf("; drift %s", paste(sprintf("%+.4f", ck$drift), collapse = " ")) else ""
      cat(sprintf("  calibration (%s): spread factors at horizons 1-%d: %s%s\n", grp, nrow(ck),
                  paste(sprintf("%.2f", ck$scale), collapse = " "), dr))
    }
  } else if (!is.null(x$calibration)) {
    cal <- x$calibration
    grp <- if (!is.null(cal$band) && any(cal$band != "all ages")) paste0(cal$group, ", ages ", cal$band) else cal$group
    ci <- if (!is.null(cal$scale_h1_lo)) sprintf(" (%.2f-%.2f)", cal$scale_h1_lo, cal$scale_h1_hi) else ""
    # slope and drift are exactly 0 when they were not estimated
    hz <- ifelse(cal$slope != 0, sprintf(" at horizon 1, x h^%.2f", cal$slope), " at all horizons")
    dr <- ifelse(cal$drift != 0, sprintf("; drift %+.4f", cal$drift), "")
    cat(sprintf(paste("  calibration (%s): spread factor %.2f%s%s%s; inner CRPS",
                      "%.4f -> %.4f, 95%% coverage %.2f -> %.2f\n"), grp, cal$scale_h1, ci, hz, dr,
                cal$crps_before, cal$crps_after, cal$cov95_before, cal$cov95_after), sep = "")
  }
  bad <- x$diagnostics[status != "ok" | nzchar(message)]
  if (nrow(bad)) {
    cat(sprintf("  problems in %d stratum/strata, for example %s: %s\n", nrow(bad), bad$stratum[1],
                bad$message[1]))
  }
  invisible(x)
}

# Short statement of the model and the sources of uncertainty in its draws
model_assumptions <- function(spec) {
  p <- spec$params
  switch(class(spec)[1],
         mt_spec_lee_carter = paste(
           "log m(x,t) = a(x) + b(x) k(t), estimated by",
           if (identical(p$family, "nb")) {
             paste("negative binomial maximum likelihood (Delwarde, Denuit and Partrat 2007): deaths with mean mu",
                   "and variance mu + phi mu^2, one dispersion phi for all ages and years.")
           } else if (p$method == "poisson") {
             "Poisson maximum likelihood (Brouhns et al. 2002)."
           } else {
             paste0("SVD of the centred log rates (Lee and Carter 1992), ", switch(
               p$adjust %||% "dt",
               dt = "with k(t) refitted so that the fitted total deaths of each year match the observed ones",
               dxt = paste("with k(t) refitted by Poisson maximum likelihood of the deaths by age given a(x) and",
                           "b(x) (Booth, Maindonald and Smith 2002)"),
               e0 = paste("with k(t) refitted so that the fitted life expectancy at the first age matches the",
                          "observed one (Lee and Miller 2001)"),
               none = "with k(t) from the SVD (not refitted)"), ".")
           },
           "k(t) follows a random walk with drift; draws include",
           "its innovations and the uncertainty of the drift. a(x) and b(x) are fixed",
           "(parameter uncertainty with bootstrapped()). Predictive draws add",
           if (identical(p$family, "nb")) {
             "negative binomial variation of the deaths given the rates, with the estimated dispersion."
           } else {
             "Poisson variation of the deaths given the rates."
           },
           "Jump-off:", p$jump_off, "rates.",
           if (isTRUE(p$rotate)) {
             sprintf(paste("b(x) rotates towards an ultimate pattern (flat below 65) as life expectancy at birth",
                           "of each path rises from %g to %g (Li, Lee and Gerland 2013); k(t) is solved again so",
                           "that each path keeps the life expectancy at birth of the unrotated model; rates change",
                           "from the jump-off year."), p$rotation[["e0_start"]], p$rotation[["e0_end"]])
           }),
         mt_spec_p_splines = paste(
           "Two-dimensional P-splines (Currie, Durban and Eilers 2004; MortalitySmooth): log m(x,t) is a B-spline",
           "surface over age and year with difference penalties on its coefficients, fitted by penalised Poisson",
           sprintf("likelihood with smoothing parameters chosen by %s.", toupper(p$criterion)),
           "The forecast extends the basis over the forecast years with weight 0, so the penalty (order",
           sprintf("%d over years) extrapolates the trend of each age.", p$pord[2]),
           "Draws come from the normal approximation of the coefficients (covariance of the penalised",
           if (isTRUE(p$overdispersion)) "likelihood, scaled by the overdispersion);" else "likelihood);",
           "their variance grows with the horizon. Predictive draws add Poisson variation of the deaths."),
         mt_spec_cp_splines = paste(
           "CP-splines (Camarda 2019): two-dimensional P-splines, with coefficients of their own for age 0 where",
           "the grid has that age group; asymmetric penalties keep the forecast rate of ageing of each age within",
           sprintf("the central %g%% range of the fitting years and its rate of change over time within the",
                   100 * p$level_ageing),
           sprintf("central %g%% range. Draws are the forecasts of %d residual-bootstrap samples:",
                   100 * p$level_change, p$n_boot),
           "uncertainty of the fit, without innovations, so intervals are narrow. Predictive draws add Poisson",
           "variation of the deaths."),
         mt_spec_coda = coda_assumptions(spec),
         mt_spec_coda_causes = coda_assumptions(spec),
         mt_spec_coda_multiway = coda_assumptions(spec),
         mt_spec_rwd = if (isTRUE(p$observation_noise)) {
           paste("Each age: log rate = level + noise, level a random walk with drift (local level model;",
                 "deprecated variant). Variances from a moment estimator of this package (Poisson part plus a",
                 "common extra part, clamped, level variance partly pooled over ages). Draws include level and",
                 "drift uncertainty; predictive draws add Poisson variation of the deaths.")
         } else {
           paste0("Each age: log rate a random walk with drift, estimated from the year-to-year changes",
                  if (!is.null(p$drift_window)) sprintf(" of the last %d years", p$drift_window) else "",
                  "; draws resample the observed changes (", p$innovations, ") and the drift. Noise of observed ",
                  "rates accumulates, so intervals are wide where deaths are few.")
         },
         mt_spec_multipop = paste(
           if (p$type == "li_lee") "Li-Lee (Li and Lee 2005): log m(x,t,i) = a(x,i) + B(x) K(t) + b(x,i) k(t,i)." else
             paste("Product-ratio, Lee-Carter version (deprecated; no published source, see coherent_fdm()):",
                   "log m(x,t,i) = common Lee-Carter + a(x,i) + b(x,i) k(t,i)."),
           "The common factor (per group of `common_by`)",
           "follows a random walk with drift; each stratum's k(t,i) an AR(1) without intercept",
           "(|phi| <= 0.99, a bound of this package), independent across strata, so deviations from the common",
           "trend stay bounded. Draws include the innovations and drift uncertainty of K and the AR(1) innovations;",
           "a, b and phi are fixed. Predictive draws add Poisson variation of the deaths given the rates."),
         mt_spec_cae = paste(
           sprintf("Common age effect model (Kleinow 2015) of order %d: log m(x,t,i) = a(x,i) + sum_j", p$order),
           "beta_j(x) kappa_j(t,i); the age effects beta (orthonormal columns) are common to the populations",
           if (!is.null(p$common_by)) "of each group of `common_by`",
           "and estimated by least-squares common principal components of the centred log rates (equal weight",
           "per population); kappa = t(beta) times the centred rates. Forecast (package-specific; the paper",
           "does not forecast): the period effects of all populations follow one multivariate random walk with",
           "drift with correlated normal innovations; draws include the innovations and the drift uncertainty;",
           "a and beta are fixed. Each population keeps its own drift, so forecasts can diverge. Predictive",
           "draws add Poisson variation of the deaths given the rates."),
         mt_spec_coupled = paste(model_assumptions(p$base), sprintf(paste(
           "Coupled (ensemble copula coupling, ECC-R; Schefzik et al. 2013): the draws of each age, year and",
           "stratum keep the values of this model, reordered so that their ranks follow the draws of %s; the",
           "dependence between ages, years and strata is that of the template. Observation noise of predictive",
           "draws is added after the coupling."), p$template$label)),
         mt_spec_saint = paste(
           "SAINT (Jarner and Kryger 2011; Jarner and Jallbjorn 2022): log m(x,t) of each stratum is the log",
           "trend plus a spread y(t)' r(x) with five piecewise linear age regressors (eq. 5.11). Trend of the pooled",
           if (is.null(p$by)) "reference strata:" else
             sprintf("reference strata of each group of %s:", paste(p$by, collapse = ", ")),
           if (isFALSE(p$background)) "Gamma frailty model exp(-sigma2 M) mu0 (no background)," else
             "Gamma frailty model exp(-sigma2 M) mu0 + exp(zeta),",
           "log mu0 = alpha + beta (x - 75) + kappa (x - 75) 1{x < 75}, fitted by pseudo-likelihood (EM estimates,",
           "profile over sigma2);", if (isFALSE(p$background)) "kappa a random walk without drift;" else
             "kappa and zeta random walks without drift;",
           if (is.null(p$ecm)) "alpha and beta a bivariate random walk with drift." else
             sprintf(paste("alpha and beta of the two groups an error-correction model with the jump-off differences",
                           "as equilibrium and the mean change of %s as common drift."), p$ecm),
           "Spread fitted by Poisson maximum likelihood per year; a VAR(1) without mean,",
           if (is.numeric(p$spread_ar)) sprintf("A = %g I,", p$spread_ar) else if (p$spread_ar == "diagonal")
             "A diagonal in [0, 0.99]," else "A unrestricted,",
           "so the spread reverts towards zero. Draws include the innovations of the trend indexes and of the",
           "spread; parameters are fixed. Frailty selection makes rates of improvement change over time.",
           "Predictive draws add Poisson variation of the deaths given the rates. Jump-off: fitted rates."),
         mt_spec_coherent_fdm = paste(
           "Product-ratio method with functional data models (Hyndman, Booth and Yasmeen 2013; demography):",
           "rates smoothed over age; the geometric mean of the strata (per group of `common_by`) and the ratio",
           "of each stratum to it are decomposed into principal components; the product's scores follow ARIMA",
           "models, the ratios' scores stationary ARFIMA models, so deviations from the product stay bounded.",
           "Draws simulate the scores and resample the model errors; predictive draws add Poisson variation",
           "of the deaths given the rates."),
         mt_spec_stanmomo = paste(
           "Bayesian stochastic mortality model fitted with Stan (StanMoMo; Barigou et al. 2023):", p$model,
           "with", p$family, "likelihood and log link; period indexes random walks with drift, cohort effects",
           "AR(2). Draws are posterior draws of the expected rates (parameters and future indexes);",
           if (identical(p$family, "nb")) {
             paste("predictive draws are negative binomial deaths with the dispersion of the same posterior",
                   "draw (posterior predictive distribution).")
           } else {
             "predictive draws add Poisson variation of the deaths given the rates."
           }),
         mt_spec_bohk_rau = paste(
           "Rates of mortality improvement rho(x,t) = 1 - m(x,t+1) / m(x,t) of the",
           if (isTRUE(p$smooth)) "rates smoothed with 2D P-splines" else "observed rates",
           "(Bohk-Ewald and Rau 2017); each age follows",
           if (p$core == "exponential") "a curve exp(b1(x) + b2(x) log t)" else "a linear trend b1(x) + b2(x) t",
           "with normal noise, the curves tied together by a hierarchical prior; posterior and future",
           "improvements sampled by MCMC.",
           if (!is.null(p$limits)) sprintf("Forecast improvements are limited to %g-%g.", p$limits[1], p$limits[2]),
           "Rates start from the last", if (isTRUE(p$smooth)) "smoothed" else "observed",
           "rates and fall by the improvements;",
           if (p$draws == "quantiles") {
             paste("draw j applies the same quantile of the improvements at every age and year (the paper's",
                   "method), so draws are comonotone across ages and years.")
           } else {
             "every draw is one posterior predictive path (package-specific)."
           },
           "Predictive draws add Poisson variation of the deaths given the rates. Convergence to reference",
           "populations is not implemented."),
         mt_spec_apci = paste(
           "APCI model (Richards, Currie, Kleinow and Ritchie 2019): log m(x,y) = alpha(x) + beta(x) (y - ybar) +",
           "kappa(y) + gamma(y - x), Poisson maximum likelihood subject to five constraints on kappa and gamma",
           paste0("(over-constrained), gamma = 0 for cohorts with four or fewer observations",
                  if (isTRUE(p$smooth)) "; alpha and beta P-splines over age with smoothing parameters by BIC." else
                    "."),
           "kappa and gamma follow ARIMA models without a mean",
           if (is.null(p$kappa_order) || is.null(p$gamma_order)) "(orders chosen by AICc and RMSE)" else "",
           "and are projected independently; draws hold their innovations with the ARIMA coefficients, alpha,",
           "beta and the estimated kappa and gamma fixed, so beta(x) (y - ybar) continues linearly. Predictive",
           "draws add Poisson variation of the deaths given the rates."),
         mt_spec_nordpred = paste(
           "Nordpred (Moller et al. 2003): age-period-cohort model of 5-year age groups and 5-year periods (annual",
           "data summed to periods ending with the last year), Poisson deaths with",
           if (p$link == "power5") "rates (age + drift + period + cohort)^5;" else "log link;",
           "the drift (or the slope of the last 10 years) is projected and cut by",
           paste0(round(100 * p$cut_trend), "%", collapse = ", "), "in the future periods. Annual rates are",
           "interpolated between the midpoints of the periods (package-specific). Draws take the coefficients from",
           "their asymptotic normal distribution (package-specific): they hold no variation of future period",
           "effects and no overdispersion, so intervals are narrow. Predictive draws add Poisson variation of the",
           "deaths given the rates."),
         mt_spec_stmomo = if (identical(p$model, "plat")) {
           paste("Reduced Plat model (Plat 2009) as defined in the StMoMo vignette (Villegas et al. 2018): logit",
                 "q(x,t) = a(x) + k1(t) + (xbar - x) k2(t) + g(t - x); period indexes a multivariate random walk",
                 "with drift, cohort index", sprintf("ARIMA(%s)", paste(p$gc_order, collapse = ",")),
                 "with a constant; draws from",
                 "StMoMo's simulate() (index uncertainty, parameters fixed). Predictive draws add Poisson",
                 "variation of the deaths given the rates.")
         } else {
           paste(spec$label, "- see the help page of its constructor for the model and its uncertainty.")
         },
         mt_spec_apc_bayes = paste(
           "Bayesian age-period-cohort model (bamp; Schmid and Held 2007): logit rate = age + period + cohort",
           "effects with random-walk priors, binomial likelihood, each stratum fitted separately. Draws are",
           "posterior samples of the rates (parameters and random walks of period and cohort effects);",
           "predictive draws add Poisson variation of the deaths given the rates (bamp's binomial likelihood",
           "and overdispersion are not carried into them)."),
         mt_spec_ensemble = paste("Mixture of the component models: each draw comes from one component, chosen in",
                                  "proportion to the weights; uncertainty of the components and their spread."),
         mt_spec_reconciled = paste("Base model reconciled draw by draw over the hierarchy on the deaths (or rates)",
                                    switch(spec$params$method,
                                           bottom_up = "scale: aggregates are the sums of the bottom strata.",
                                           proportional = "scale: parts of fixed nodes scaled in proportion.",
                                           energy_score = "scale: map learnt by score-optimal reconciliation.",
                                           paste0("scale: ", spec$params$method, " projection."))),
         mt_spec_recalibrated = if (identical(p$method, "errors")) paste(model_assumptions(p$base), sprintf(paste(
           "Recalibrated by past errors: draws are the median of the base draws plus the errors (all ages and",
           "strata of one origin) of forecasts from %d origins in the last %d training years, resampled per",
           "horizon and scaled to %g%% coverage of these errors; they describe observed rates."),
           p$n_origins, p$years, 100 * p$level)) else paste(model_assumptions(p$base), sprintf(paste(
           "Calibrated: the spread%s of the draws%s is rescaled%s to minimise the %s of %s in forecasts from %d",
           "origins in the last %d training years."), if (p$bias) " and drift" else "",
           if (length(p$ages)) " (per age band)" else "",
           if (identical(p$scale, "per_horizon")) " for each horizon" else "",
           if (identical(p$target, "interval")) sprintf("%g%% interval score", 100 * p$level) else "CRPS",
           switch(p$quantity %||% "rate", rate = "log death rates", life_expectancy = "life expectancy at birth",
                  asr = "the log ASR"), p$n_origins, p$years)),
         mt_spec_anchored = paste(model_assumptions(p$base), sprintf(paste(
           "Anchored: forecasts start at the observed level of the last %d training year(s), with offsets per",
           "age and stratum%s %s."), p$years,
           if (isTRUE(p$shrink)) " (shrunk where deaths are few)" else "",
           if (is.finite(p$decay)) sprintf("fading out over %s years", format(p$decay)) else "kept at all horizons")),
         paste(spec$label, "- see the help page of its constructor for the model and its uncertainty."))
}

#' Residual heat map of a fitted model
#'
#' @param fit An `mt_fit` object.
#' @param strata Labels of the strata to show (default: up to four).
#' @param type `"deviance"` or `"pearson"`.
#' @return A ggplot object: residuals by age and period, one panel per
#'   stratum. Diagonal bands indicate cohort effects, vertical bands period
#'   shocks.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   x <- mt_example(deaths = "single")
#'   mt_plot_residuals(mt_fit(lee_carter(), x))
#' }
#' @export
mt_plot_residuals <- function(fit, strata = NULL, type = c("deviance", "pearson")) {
  need_ggplot()
  type <- match.arg(type)
  r <- stats::residuals(fit, type = type)
  strata <- check_strata_labels(strata %||% utils::head(fit$strata$label, 4), fit$strata$label)
  r <- r[stratum %in% strata]
  r[, age_i := match(age, fit$age$labels)]
  lim <- max(3, stats::quantile(abs(r$residual), 0.99))
  ggplot2::ggplot(r, ggplot2::aes(x = period, y = age_i, fill = pmax(pmin(residual, lim), -lim))) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-lim, lim),
                                  name = paste(type, "residual")) +
    ggplot2::scale_y_continuous(breaks = seq_along(fit$age$labels), labels = fit$age$labels) +
    ggplot2::facet_wrap(~stratum) +
    ggplot2::labs(x = "Period", y = "Age") +
    theme_mortensor(legend = "right")
}

# Check stratum labels for plots and list the available ones
check_strata_labels <- function(strata, available) {
  bad <- setdiff(strata, available)
  if (length(bad)) {
    mt_abort(sprintf("Unknown stratum label(s): %s. Available: %s.", paste(bad, collapse = ", "),
                     paste(utils::head(available, 12), collapse = ", ")), class = "argument")
  }
  strata
}
