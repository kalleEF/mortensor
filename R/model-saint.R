#' SAINT model: trend of a reference population and stationary spreads
#'
#' The SAINT model (Spread Adjusted InterNational Trend) of Jarner and Kryger
#' (2011) in the revised form of Jarner and Jallbjørn (2022). It forecasts the
#' adult mortality of small populations as deviations (spreads) from the
#' trend of a large reference population. The trend carries the long-term
#' development. The spreads are stationary, so the populations do not
#' diverge from the reference.
#'
#' The reference is one or more strata of the data (`reference`), for
#' example the national total or several countries. Their deaths and
#' exposures are pooled into one reference population (Jarner and Kryger
#' 2011, Section 3.1). With `by` (for example `"sex"`), each group gets its
#' own trend; both papers fit women and men separately. All other strata of a
#' group are targets with a spread from the group's trend. A group with a
#' single reference stratum forecasts that stratum by the trend. In a group
#' with several reference strata, each of them is also a target, as Jarner
#' and Kryger (2011, concluding remarks) propose for a group of related
#' populations.
#'
#' @section Version:
#' This is the version of Jarner and Jallbjørn (2022), the model in use at
#' the Danish pension fund ATP. Its authors replaced three components of the
#' original model: the deterministic parametric trend by a stochastic frailty
#' model, so that forecasts describe the uncertainty of the trend (Section
#' 2.2); the trend's age pattern by one that fits high ages better (Section
#' 2.4); and the three spread regressors (level, slope, curvature) by five
#' (Section 5.4). Section and equation numbers below refer to Jarner and
#' Jallbjørn (2022) unless the 2011 paper is named.
#'
#' @section Trend (Sections 4 and 5.1):
#' The reference rates are `mu(t, x) = exp(-sigma2 M(t, x)) mu0(t, x) +
#' exp(zeta_t)` (eqs. 5.1-5.3), with baseline `log mu0(t, x) = alpha_t +
#' beta_t (x - 75) + kappa_t (x - 75) 1{x < 75}` (eq. 5.2), background
#' mortality `exp(zeta_t)`, and Gamma-distributed frailty with mean 1 and
#' variance `sigma2` (Example 4.1). `exp(-sigma2 M)` is the mean frailty of the
#' survivors of a cohort, where `M(t, x)` cumulates the observed rates of the
#' cohort, background subtracted, from the first age of the data; years
#' before the first year take the rates of the first year (eqs. 4.7-4.8).
#' Deaths are Poisson (eq. 4.9). For a given `sigma2`, the indexes of every
#' year are the estimates of the EM algorithm of Appendix B; `sigma2`
#' maximises the profile pseudo-likelihood (eq. 4.14). `kappa` and `zeta`
#' follow random walks without drift (eq. 5.4). `alpha` and `beta` follow a
#' bivariate random walk with drift (Section 4.5) or, with `ecm`, the
#' error-correction model of Section 5.2.1 for two groups (eqs. 5.5-5.6): the
#' equilibrium differences are the differences at the jump-off year and the
#' common drift is the mean change of the group named by `ecm` (women in the
#' paper). The median forecast is then a random walk with this drift for
#' both groups, and in every path the differences between the groups stay
#' stationary. Rates are forecast by eqs. 4.15-4.17: the mean frailty follows
#' from the baseline cumulated along each cohort. Selection by frailty makes
#' the rates of improvement at high ages rise over time (Section 3.3).
#'
#' @section Spread (Section 5.4):
#' A target's rates are the trend times `exp(y_t' r(x))` (eq. 5.10), with
#' regressors `r_i(x) = min(1, max(0, x_i - x) / 20)` and knots `x_i` = 40,
#' 60, 80, 100, 120 (eq. 5.11). The spread parameters `y_t` of each year are
#' Poisson maximum likelihood estimates given the fitted trend (eq. 5.9); the
#' trend is estimated first and kept fixed, the two-stage estimation of
#' Jarner and Kryger (2011, concluding remarks). They
#' follow a VAR(1) model without mean, `y_t = A y_{t-1} + u_t` with normal
#' innovations (eq. 5.12), so the spread reverts towards zero. `spread_ar`
#' sets `A`; the paper recommends a diagonal `A` with `0 <= A_ii < 1` and
#' imposes `A = 0.99 I` for Denmark, a deviation half-life of about 69 years
#' (Section 5.4.1). Each target has a VAR(1) of its own, "a separate VAR
#' model" per target population (Section 5.4.1).
#'
#' @section Forecasts and uncertainty:
#' Draws include the innovations of the trend indexes and of the spread.
#' Trend and spread are independent (Jarner and Kryger 2011, eq. 5).
#' Parameters are fixed, as in Jarner and Kryger (2011, Section 2.3.3), who
#' leave parameter uncertainty aside. Forecasts start from the fitted rates of the last year;
#' the paper uses no jump-off correction (2022, Section 5.3). Predictive
#' draws (`mt_forecast(predictive = TRUE)`) add Poisson variation of the
#' deaths given the rates, the observation model of the papers.
#'
#' @section Ages:
#' Both papers use single ages 20-100 and model adult mortality only (2011,
#' Section 1.2). `saint()` needs consecutive single-year age groups from age
#' 20 without an open interval: the frailty term follows cohorts year by
#' year. Ages below 20 are refused. Ages above 100 are accepted; the
#' parametric trend extends to them (2022, Section 5.1.1). Drop an open
#' interval with `mt_subset(x, ages = c(20, 100))`; estimate single years
#' from age groups with [mt_ungroup()] or `step_ungroup()`.
#'
#' @section Corrections and package-specific choices:
#' Three formulas of Jarner and Jallbjørn (2022) are read as evident
#' misprints. In (B14) and (B16) the denominator `c F - G` is the rate `c F +
#' G` of (B7)-(B8). In the second line of (4.16) the baseline of year
#' `t_max` is that of year `t - 1`. The common drift of Section 5.2.1 divides
#' the sum of changes by `t_max - t_min + 1`; the package divides by the
#' number of changes, the empirical average that the text describes.
#'
#' Package-specific (reference unknown):
#'
#' * the fixed point of the EM algorithm is computed by Newton-Raphson
#'   sweeps (Fisher scoring) over the parameters of each year with the mean
#'   frailty held fixed, alternated with updates of the mean frailty; Section
#'   4.4 names Newton-Raphson as the alternative to EM, and the tests check
#'   that both give the same estimates;
#' * `sigma2` is searched from 0 by doubling (0.025, 0.05, ..., 3.2) and then
#'   by [stats::optimize()] around the best value; a warning (class
#'   `mt_warning_saint`) is given at the upper end;
#' * starting values, and a lower bound of `1e-10` on the background rate,
#'   with a warning when it binds (see `background`);
#' * on age ranges narrower than 20-100: `kappa` is set to 0 without ages on
#'   both sides of 75, and spread regressors that are 0 at all ages, or 1 at
#'   all ages like the level, are left out;
#' * the innovations of `kappa` and `zeta` are independent of each other and
#'   of `alpha` and `beta` (the paper gives their marginal distributions),
#'   and the spread innovations of different targets are independent;
#'   variances are maximum likelihood estimates;
#' * the adjustment coefficients and the innovation covariance of the
#'   error-correction model are Gaussian maximum likelihood estimates given
#'   the equilibrium differences and the drift; a warning is given when the
#'   model does not correct deviations (eigenvalues of eq. 5.8 outside (-1,
#'   1));
#' * `spread_ar = "diagonal"`: each `A_ii` is the least-squares coefficient
#'   of its component, bounded to `[0, 0.99]` (a warning names the strata
#'   where 0.99 binds); `"full"`: conditional maximum likelihood (the 2011
#'   paper uses Yule-Walker estimates of the same model), with a warning
#'   when `A` is not stationary.
#'
#' `background = FALSE` is an opt-in variant of this package: the SAINT
#' specification of Section 5.1 has a background, the class of models of
#' Section 4.2 allows it to be left out.
#'
#' Not implemented: the deterministic trend and the three regressors of
#' 2011, the error-correction model with free equilibrium and drift, and
#' frailty distributions other than Gamma (Appendix A).
#'
#' @param reference Reference strata: stratum labels, or a named list of
#'   values of stratum variables, for example `list(region = "Deutschland")`.
#' @param by Stratum variable(s) with a trend each, for example `"sex"`.
#'   `NULL`: one trend for all strata.
#' @param ecm `NULL`, or the group of `by` whose mean changes give the common
#'   drift of the error-correction model that couples the level and slope of
#'   two trends (Section 5.2.1), for example `"female"`. Needs `by` with two
#'   groups.
#' @param spread_ar `A` of the VAR(1) of the spread: `"diagonal"` (diagonal,
#'   estimated, `0 <= A_ii <= 0.99`), `"full"` (unrestricted, eq. 5.12), or a
#'   number `a` in `[0, 1)` for `A = a I` (the paper's choice for Denmark:
#'   0.99).
#' @param background Whether the trend has the background mortality
#'   `exp(zeta_t)` of eq. 5.3 (default). `FALSE` leaves it out, a model of
#'   the class of Section 4.2, which allows empty components. Use it where
#'   the data do not identify a background, for example rates smoothed over
#'   age; a warning says when the estimate of the background reaches its
#'   lower bound, which makes its random walk and the forecasts at young ages
#'   unstable.
#' @return A model specification (joint over strata).
#' @references
#' Jarner SF, Kryger EM (2011). Modelling adult mortality in small
#' populations: the SAINT model. *ASTIN Bulletin* 41(2): 377-418.
#' \doi{10.2143/AST.41.2.2136982}
#'
#' Jarner SF, Jallbjørn S (2022). The SAINT model: a decade later. *ASTIN
#' Bulletin* 52(2): 483-517. \doi{10.1017/asb.2021.37}
#' @seealso [li_lee()] and [coherent_fdm()], other models whose strata keep
#'   bounded distances to a common trend.
#' @examples
#' x <- mt_example_regions(years = 1995:2019, regions = 3, ages = "single")
#' x <- mt_aggregate_strata(mt_subset(x, ages = c(40, 100)), mt_hierarchy(region = "all", quintile = "all"))
#' x <- mt_subset(x, quintile == "all", drop = TRUE)
#' fit <- mt_fit(saint(reference = list(region = "all"), by = "sex"), x)
#' fc <- mt_forecast(fit, h = 10, n_draws = 100, seed = 1)
#' p <- mt_parameters(fit)
#' p[parameter == "sigma2"]
#' @export
saint <- function(reference, by = NULL, ecm = NULL, spread_ar = "diagonal", background = TRUE) {
  check_flag(background, "background")
  ok_ref <- (is.character(reference) && length(reference) > 0 && !anyNA(reference)) ||
    (is.list(reference) && length(reference) > 0 && !is.null(names(reference)) && all(nzchar(names(reference))))
  if (!ok_ref) {
    mt_abort(paste("`reference` must name the reference strata: stratum labels, or a named list (or named character",
                   "vector) of values of stratum variables, for example list(region = \"Deutschland\")."),
             class = "argument")
  }
  if (!is.null(by) && (!is.character(by) || !length(by) || anyNA(by))) {
    mt_abort("`by` must be NULL or names of stratum variables, for example \"sex\".", class = "argument")
  }
  if (!is.null(ecm)) {
    if (!is.character(ecm) || length(ecm) != 1L || is.na(ecm)) {
      mt_abort("`ecm` must be NULL or one group of `by`, for example \"female\".", class = "argument")
    }
    if (is.null(by)) mt_abort("`ecm` couples the trends of two groups and needs `by`.", class = "argument")
  }
  ok_ar <- (is.character(spread_ar) && length(spread_ar) == 1L && spread_ar %in% c("diagonal", "full")) ||
    (is.numeric(spread_ar) && length(spread_ar) == 1L && !is.na(spread_ar) && spread_ar >= 0 && spread_ar < 1)
  if (!ok_ar) {
    mt_abort("`spread_ar` must be \"diagonal\", \"full\" or a number in [0, 1).", class = "argument")
  }
  new_mt_spec("saint", list(reference = reference, by = by, ecm = ecm, spread_ar = spread_ar, background = background),
              joint = TRUE, min_periods = 8,
              label = paste0("SAINT", if (!is.null(by)) sprintf(" (trend by %s%s)", paste(by, collapse = ", "),
                                                                if (!is.null(ecm)) ", error correction" else ""),
                             if (!background) ", no background" else ""))
}

#' @export
fit_model.mt_spec_saint <- function(spec, D, E, ages, periods, strata = NULL, hierarchy = NULL, ...) {
  p <- spec$params
  x <- saint_check_ages(ages)
  ns <- dim(D)[3]
  st <- saint_strata(p, strata, ns)
  ng <- length(st$groups)
  # trend of each group from the pooled deaths and exposures of its reference
  # strata (Jarner and Kryger 2011, Section 3.1)
  trends <- lapply(seq_len(ng), function(g) {
    use <- st$group == g & st$reference
    saint_trend_fit(apply(D[, , use, drop = FALSE], 1:2, sum), apply(E[, , use, drop = FALSE], 1:2, sum), x,
                    background = !isFALSE(p$background))
  })
  names(trends) <- st$groups
  dyn <- saint_trend_dynamics(trends, p$ecm, st$groups)
  # spreads: every stratum except the single reference stratum of a group,
  # whose rates are the trend
  cols <- saint_spread_columns(x)
  R <- saint_regressors(x)[, cols, drop = FALSE]
  spreads <- lapply(seq_len(ns), function(s) {
    if (st$sole[s]) return(NULL)
    mu_ref <- trends[[st$group[s]]]$fitted
    y <- saint_spread_years(D[, , s], E[, , s], mu_ref, R)
    c(list(y = y), saint_var_fit(y, p$spread_ar))
  })
  fitted <- array(NA_real_, dim(D))
  for (s in seq_len(ns)) {
    mu_ref <- trends[[st$group[s]]]$fitted
    fitted[, , s] <- if (is.null(spreads[[s]])) mu_ref else mu_ref * exp(R %*% t(spreads[[s]]$y))
  }
  saint_report(trends, dyn, spreads, strata$label %||% as.character(seq_len(ns)), p$spread_ar)
  list(x = x, periods = periods, group = st$group, groups = st$groups, reference = st$reference, sole = st$sole,
       trends = trends, dynamics = dyn, R = R, spread_columns = cols, spreads = spreads, fitted = fitted)
}

#' @export
simulate_model.mt_spec_saint <- function(spec, fit, h, n_draws, ...) {
  na <- length(fit$x)
  ns <- length(fit$group)
  paths <- saint_index_paths(fit, h, n_draws)
  mu_ref <- lapply(seq_along(fit$trends), function(g) saint_trend_forecast(fit$trends[[g]], paths[[g]], h, n_draws))
  out <- array(NA_real_, c(na, h, ns, n_draws))
  for (s in seq_len(ns)) {
    base <- mu_ref[[fit$group[s]]]
    sp <- fit$spreads[[s]]
    if (is.null(sp)) {
      out[, , s, ] <- base
      next
    }
    y <- saint_var_paths(sp, h, n_draws)
    for (k in seq_len(h)) out[, k, s, ] <- base[, k, ] * exp(fit$R %*% y[, k, ])
  }
  out
}

# Checks ----------------------------------------------------------------------

# Integer ages of single-year closed age groups from age 20 (both papers use
# ages 20-100)
saint_check_ages <- function(ages) {
  lo <- ages$lower
  up <- ages$upper
  if (any(!is.finite(up)) || any(abs(up - lo - 1) > 1e-9) || any(abs(lo - round(lo)) > 1e-9) ||
      (length(lo) > 1 && any(abs(diff(lo) - 1) > 1e-9))) {
    mt_abort(paste("saint() needs consecutive single-year age groups without an open interval: the frailty term",
                   "cumulates the rates of each cohort year by year (Jarner and Jallbjorn 2022, eq. 4.7), and both",
                   "papers use ages 20-100. Drop the open interval with mt_subset(x, ages = c(20, 100)); estimate",
                   "single years from age groups with mt_ungroup() or step_ungroup()."), class = "model_data",
             call = NULL)
  }
  if (min(lo) < 20) {
    mt_abort(paste("saint() models adult mortality from age 20 (Jarner and Kryger 2011, Section 1.2; both papers use",
                   "ages 20-100); its baseline cannot describe the mortality of children and adolescents. Restrict",
                   "the ages with mt_subset(x, ages = c(20, 100))."), class = "model_data", call = NULL)
  }
  if (length(lo) < 10) mt_abort("saint() needs at least ten single ages.", class = "model_data", call = NULL)
  round(lo)
}

# Groups (one trend each), reference strata and the strata whose rates are
# the trend (the single reference stratum of a group)
saint_strata <- function(p, strata, ns) {
  if (is.null(strata) || nrow(strata) != ns) mt_abort("saint() needs the strata of the data.", class = "argument")
  vars <- setdiff(names(strata), c("stratum", "label"))
  ref <- p$reference
  if (!is.null(names(ref)) && all(nzchar(names(ref)))) {
    bad <- setdiff(names(ref), vars)
    if (length(bad)) {
      mt_abort(sprintf("`reference` names unknown stratum variable(s) %s; the strata have %s.",
                       paste(bad, collapse = ", "), paste(vars, collapse = ", ")), class = "argument", call = NULL)
    }
    is_ref <- rep(TRUE, ns)
    for (v in unique(names(ref))) {
      is_ref <- is_ref & as.character(strata[[v]]) %in% as.character(unlist(ref[names(ref) == v]))
    }
  } else {
    bad <- setdiff(ref, strata$label)
    if (length(bad)) {
      mt_abort(sprintf("Unknown reference stratum/strata %s. Available: %s.", paste(bad, collapse = ", "),
                       paste(utils::head(strata$label, 12), collapse = ", ")), class = "argument", call = NULL)
    }
    is_ref <- strata$label %in% ref
  }
  if (!any(is_ref)) mt_abort("No stratum matches `reference`.", class = "argument", call = NULL)
  group <- rep(1L, ns)
  groups <- "all"
  if (!is.null(p$by)) {
    bad <- setdiff(p$by, vars)
    if (length(bad)) {
      mt_abort(sprintf("`by` names unknown stratum variable(s) %s.", paste(bad, collapse = ", ")), class = "argument",
               call = NULL)
    }
    key <- do.call(paste, c(lapply(p$by, function(v) as.character(strata[[v]])), sep = "."))
    groups <- unique(key)
    group <- match(key, groups)
  }
  n_ref <- tabulate(group[is_ref], length(groups))
  if (any(n_ref == 0)) {
    mt_abort(sprintf("Group(s) %s of `by` have no reference stratum.", paste(groups[n_ref == 0], collapse = ", ")),
             class = "argument", call = NULL)
  }
  if (!is.null(p$ecm)) {
    if (length(groups) != 2L || !p$ecm %in% groups) {
      mt_abort(sprintf(paste("The error-correction model couples two trends: `by` must give two groups and `ecm`",
                             "must name one of them (groups: %s)."), paste(groups, collapse = ", ")),
               class = "argument", call = NULL)
    }
  }
  sole <- is_ref & n_ref[group] == 1L
  list(group = group, groups = groups, reference = is_ref, sole = sole)
}

# Trend -------------------------------------------------------------------------

# Design of the baseline log rate (eq. 5.2): log mu0(t, x) = alpha_t +
# beta_t (x - 75) + kappa_t (x - 75) 1{x < 75}; the excess slope kappa is
# identified only with ages on both sides of 75
saint_design <- function(x) {
  X <- cbind(alpha = 1, beta = x - 75)
  if (any(x < 75) && any(x > 75)) X <- cbind(X, kappa = (x - 75) * (x < 75))
  X
}

# Cumulated rates M(t, x) of eqs. 4.7-4.8 from the rates with the background
# subtracted, mt [age, period] (ages x_min, x_min + 1, ...): the sum along the
# cohort diagonal from age x_min to x - 1, with the rates of the first year
# for the years before it; M = 0 at age x_min
saint_cumulated <- function(mt) {
  na <- nrow(mt)
  nt <- ncol(mt)
  M <- matrix(0, na, nt)
  if (na > 1) {
    for (i in 2:na) {
      M[i, 1] <- M[i - 1, 1] + mt[i - 1, 1]
      if (nt > 1) M[i, -1] <- M[i - 1, -nt] + mt[i - 1, -nt]
    }
  }
  M
}

# Lower bound of the log background rate zeta (a background rate of about
# 1e-10; package-specific): where the data hold no background mortality, its
# estimate would go to minus infinity
saint_zeta_min <- -23

# Starting values of the trend (package-specific): the log-linear baseline
# fitted to all deaths, and a background of half the lowest observed rate of
# each year
saint_trend_init <- function(D, E, X) {
  m <- D / E
  theta <- matrix(c(log(sum(D) / sum(E)), 0.1, rep(0, ncol(X) - 2)), ncol(X), ncol(D))
  theta <- poisson_by_year(D, E, X, theta)
  zeta <- vapply(seq_len(ncol(D)), function(t) {
    pos <- m[, t][m[, t] > 0]
    log(0.5 * if (length(pos)) min(pos) else 1e-8)
  }, numeric(1))
  list(theta = theta, zeta = zeta)
}

# Trend for a fixed frailty variance sigma2: the EM algorithm of Jarner and
# Jallbjorn (2022, Appendix B) for the pseudo-likelihood model (4.9)-(4.10)
# with Gamma frailty (4.13), baseline (5.2) and background (5.3). Each
# iteration computes M(t, x) from the observed rates and the current
# background (4.7)-(4.8) and c = exp(-sigma2 M) (B12), splits the deaths into
# their expected senescent and background parts (B14, B16; the paper prints
# c F - G in their denominators, a misprint of the rate c F + G), and fits
# the baseline by Poisson maximum likelihood per year with offset c E (B13)
# and the background as deaths over exposures (B15). Iterations stop when no
# parameter changes by more than `tol`.
saint_trend_em <- function(D, E, x, sigma2, init = NULL, tol = 1e-9, max_iter = 10000) {
  X <- saint_design(x)
  na <- nrow(D)
  nt <- ncol(D)
  m <- D / E
  st <- init %||% saint_trend_init(D, E, X)
  theta <- st$theta
  zeta <- st$zeta
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    G <- matrix(exp(zeta), na, nt, byrow = TRUE)
    cc <- exp(-sigma2 * saint_cumulated(m - G))
    Fx <- exp(X %*% theta)
    mu <- cc * Fx + G
    theta_new <- poisson_by_year(D * cc * Fx / mu, cc * E, X, theta)
    zeta_new <- pmax(log(colSums(D * G / mu) / colSums(E)), saint_zeta_min)
    change <- max(abs(theta_new - theta), abs(zeta_new - zeta))
    theta <- theta_new
    zeta <- zeta_new
    if (!is.finite(change)) break
    if (change < tol) {
      converged <- TRUE
      break
    }
  }
  G <- matrix(exp(zeta), na, nt, byrow = TRUE)
  M <- saint_cumulated(m - G)
  Fx <- exp(X %*% theta)
  mu <- exp(-sigma2 * M) * Fx + G
  ll <- sum(ifelse(D > 0, D * log(mu * E), 0) - mu * E)
  if (!is.finite(ll)) ll <- -Inf
  list(theta = theta, zeta = zeta, sigma2 = sigma2, M = M, fitted = mu, baseline = Fx, loglik = ll,
       iterations = it, converged = converged, X = X)
}

# Trend without background mortality (background = FALSE; G empty in
# eq. 4.3): for a fixed sigma2, c = exp(-sigma2 M) follows from the observed
# rates alone and the baseline is the Poisson maximum likelihood estimate per
# year with offset c E (Section 4.4: the case mu_b = 0)
saint_trend_nobg <- function(D, E, x, sigma2, init = NULL) {
  X <- saint_design(x)
  st <- init %||% saint_trend_init(D, E, X)
  M <- saint_cumulated(D / E)
  cc <- exp(-sigma2 * M)
  theta <- poisson_by_year(D, cc * E, X, st$theta)
  rownames(theta) <- colnames(X)
  Fx <- exp(X %*% theta)
  mu <- cc * Fx
  ll <- sum(ifelse(D > 0, D * log(mu * E), 0) - mu * E)
  if (!is.finite(ll)) ll <- -Inf
  list(theta = theta, zeta = rep(-Inf, ncol(D)), sigma2 = sigma2, M = M, fitted = mu, baseline = Fx, loglik = ll,
       iterations = 1L, converged = all(is.finite(theta)), X = X)
}

# The same estimates as saint_trend_em(), computed faster: the fixed point of
# the EM algorithm maximises the likelihood of (4.9)-(4.10) in the baseline
# and background parameters with c = exp(-sigma2 M) held at the value of the
# fixed point. Here c is computed from the current background, all
# parameters of each year are fitted for this c by Fisher scoring with step
# halving (Newton-Raphson sweeps, the alternative to EM that Section 4.4
# names), and c is updated until no parameter changes by more than `tol`.
saint_trend_newton <- function(D, E, x, sigma2, init = NULL, tol = 1e-9, max_iter = 200) {
  X <- saint_design(x)
  na <- nrow(D)
  nt <- ncol(D)
  p <- ncol(X)
  m <- D / E
  st <- init %||% saint_trend_init(D, E, X)
  phi <- rbind(st$theta, pmax(st$zeta, saint_zeta_min))
  pairs <- which(upper.tri(diag(p + 1), diag = TRUE), arr.ind = TRUE)
  ll_year <- function(ph, cc, j = seq_len(nt)) {
    mu <- cc[, j, drop = FALSE] * exp(X %*% ph[1:p, , drop = FALSE]) +
      matrix(exp(ph[p + 1, ]), na, length(j), byrow = TRUE)
    colSums(ifelse(D[, j, drop = FALSE] > 0, D[, j, drop = FALSE] * log(mu), 0) - mu * E[, j, drop = FALSE])
  }
  rates <- function(ph) {
    G <- matrix(exp(ph[p + 1, ]), na, nt, byrow = TRUE)
    exp(-sigma2 * saint_cumulated(m - G)) * exp(X %*% ph[1:p, , drop = FALSE]) + G
  }
  converged <- FALSE
  mu_old <- rates(phi)
  for (outer in seq_len(max_iter)) {
    theta_old <- phi[1:p, , drop = FALSE]
    G <- matrix(exp(phi[p + 1, ]), na, nt, byrow = TRUE)
    cc <- exp(-sigma2 * saint_cumulated(m - G))
    cur <- ll_year(phi, cc)
    for (it in 1:50) {
      Fx <- exp(X %*% phi[1:p, , drop = FALSE])
      G <- matrix(exp(phi[p + 1, ]), na, nt, byrow = TRUE)
      mu <- cc * Fx + G
      r <- D / mu - E
      w <- E / mu
      sF <- cc * Fx
      g <- rbind(crossprod(X, r * sF), colSums(r * G))
      H <- array(0, c(p + 1, p + 1, nt))
      for (k in seq_len(nrow(pairs))) {
        i <- pairs[k, 1]
        j <- pairs[k, 2]
        Ji <- if (i <= p) sF * X[, i] else G
        Jj <- if (j <= p) sF * X[, j] else G
        v <- colSums(w * Ji * Jj)
        H[i, j, ] <- v
        H[j, i, ] <- v
      }
      step <- solve_sym_batch(H, g)
      proj <- function(s) {
        nw <- phi + s
        nw[p + 1, ] <- pmax(nw[p + 1, ], saint_zeta_min)
        nw
      }
      new <- proj(step)
      nl <- ll_year(new, cc)
      bad <- !is.finite(nl) | nl < cur - 1e-12 * abs(cur)
      half <- 0
      while (any(bad) && half < 30) {
        step[, bad] <- step[, bad] / 2
        new[, bad] <- proj(step)[, bad]
        nl[bad] <- ll_year(new[, bad, drop = FALSE], cc, which(bad))
        bad <- !is.finite(nl) | nl < cur - 1e-12 * abs(cur)
        half <- half + 1
      }
      moved <- max(abs(new - phi))
      gain <- max(nl - cur)
      phi <- new
      cur <- nl
      # stop at a negligible step, or a negligible gain of the log-likelihood
      # (flat directions, such as a background rate near zero)
      if (moved < tol / 10 || gain < 1e-10) break
    }
    mu_new <- rates(phi)
    change <- max(abs(phi[1:p, , drop = FALSE] - theta_old), abs(log(mu_new / mu_old)))
    mu_old <- mu_new
    if (!is.finite(change)) break
    if (change < tol) {
      converged <- TRUE
      break
    }
  }
  theta <- phi[1:p, , drop = FALSE]
  rownames(theta) <- colnames(X)
  zeta <- phi[p + 1, ]
  G <- matrix(exp(zeta), na, nt, byrow = TRUE)
  M <- saint_cumulated(m - G)
  Fx <- exp(X %*% theta)
  mu <- exp(-sigma2 * M) * Fx + G
  ll <- sum(ifelse(D > 0, D * log(mu * E), 0) - mu * E)
  if (!is.finite(ll)) ll <- -Inf
  list(theta = theta, zeta = zeta, sigma2 = sigma2, M = M, fitted = mu, baseline = Fx, loglik = ll,
       iterations = outer, converged = converged, X = X)
}

# Maximum likelihood of the trend: the profile pseudo-log-likelihood (4.14)
# over the frailty variance sigma2. The search (package-specific) evaluates
# sigma2 = 0, 0.025, 0.05, 0.1, ... (doubling up to `sigma2_max`) until the
# profile falls, then maximises it between the neighbours of the best value
# (stats::optimize()); each evaluation starts from the previous estimates.
# For each sigma2 the estimates of the EM algorithm are computed by
# saint_trend_newton() (saint_trend_nobg() without background).
saint_trend_fit <- function(D, E, x, sigma2_max = 3.2, background = TRUE) {
  X <- saint_design(x)
  state <- new.env(parent = emptyenv())
  state$init <- saint_trend_init(D, E, X)
  one <- function(s2, ...) {
    if (background) saint_trend_newton(D, E, x, s2, init = state$init, ...) else
      saint_trend_nobg(D, E, x, s2, init = state$init)
  }
  profile <- function(s2) {
    f <- one(s2, tol = 1e-8, max_iter = 50)
    if (all(is.finite(f$theta)) && all(!is.na(f$zeta))) state$init <- list(theta = f$theta, zeta = f$zeta)
    f$loglik
  }
  grid <- c(0, 0.025 * 2^(0:floor(log2(sigma2_max / 0.025))))
  ll <- rep(-Inf, length(grid))
  starts <- vector("list", length(grid))
  for (i in seq_along(grid)) {
    ll[i] <- profile(grid[i])
    starts[[i]] <- state$init
    if (i > 1 && ll[i] < max(ll[seq_len(i - 1)])) break
  }
  best <- which.max(ll)
  s2 <- grid[best]
  lo <- grid[max(best - 1, 1)]
  hi <- grid[min(best + 1, length(grid))]
  state$init <- starts[[best]]
  opt <- stats::optimize(profile, c(lo, hi), maximum = TRUE, tol = 1e-7)
  if (opt$objective > ll[best]) s2 <- opt$maximum else state$init <- starts[[best]]
  fit <- one(s2)
  fit$sigma2_bound <- s2 > 0.99 * grid[length(grid)]
  fit$zeta_bound <- background && any(fit$zeta <= saint_zeta_min + 1e-8)
  fit
}

# Time dynamics of the trend (Jarner and Jallbjorn 2022, Section 5.2): the
# excess slope kappa and the background zeta of each group follow random
# walks without drift (5.4); the level alpha and the slope beta follow a
# bivariate random walk with drift (Section 4.5) or, with `ecm`, the
# error-correction model (5.5)-(5.6) of the two groups. Variances are
# maximum likelihood estimates (sums of squares over the number of changes)
saint_trend_dynamics <- function(trends, ecm, groups) {
  rw0 <- function(v) sqrt(mean(diff(v)^2))
  out <- lapply(trends, function(tr) {
    th <- tr$theta
    d <- t(diff(t(th[1:2, , drop = FALSE])))
    drift <- rowMeans(d)
    res <- d - drift
    list(drift = drift, cov = tcrossprod(res) / ncol(res),
         sd_kappa = if (nrow(th) > 2) rw0(th[3, ]) else 0, sd_zeta = if (all(is.finite(tr$zeta))) rw0(tr$zeta) else 0)
  })
  ec <- NULL
  if (!is.null(ecm)) {
    lead <- match(ecm, groups)
    other <- 3L - lead
    tf <- trends[[lead]]$theta
    tm <- trends[[other]]$theta
    ec <- saint_ecm_fit(tf[1, ], tm[1, ], tf[2, ], tm[2, ])
    ec$lead <- lead
    ec$other <- other
  }
  list(groups = out, ecm = ec)
}

# Error-correction model (5.5)-(5.6) of Y = (alpha_f, alpha_m, beta_f,
# beta_m) with the restrictions of Section 5.2.1: the equilibrium differences
# C_D are the differences at the jump-off year, the common drift C_Delta the
# mean change of the leading group ("female"; the paper divides the sum of
# changes by t_max - t_min + 1, here by the number of changes, the
# empirical average its text describes). The adjustment coefficients a11,
# a21, a32, a42 and the covariance Omega of the innovations are Gaussian
# maximum likelihood estimates given C_D and C_Delta (iterated generalised
# least squares; package-specific, the paper does not state the estimator).
saint_ecm_fit <- function(af, am, bf, bm) {
  n_t <- length(af)
  Y <- cbind(af, am, bf, bm)
  n <- n_t - 1
  dY <- Y[-1, , drop = FALSE] - Y[-n_t, , drop = FALSE]
  CD <- c(af[n_t] - am[n_t], bf[n_t] - bm[n_t])
  ca <- (af[n_t] - af[1]) / n
  cb <- (bf[n_t] - bf[1]) / n
  Cd <- c(ca, ca, cb, cb)
  za <- (af - am)[-n_t] - CD[1]
  zb <- (bf - bm)[-n_t] - CD[2]
  Z <- cbind(za, za, zb, zb)
  Ys <- dY - matrix(Cd, n, 4, byrow = TRUE)
  ssz <- colSums(Z^2)
  a <- ifelse(ssz > 0, colSums(Ys * Z) / pmax(ssz, 1e-300), 0)
  for (it in 1:200) {
    e <- Ys - Z * matrix(a, n, 4, byrow = TRUE)
    Om <- crossprod(e) / n
    Oi <- tryCatch(solve(Om), error = function(err) NULL)
    if (is.null(Oi) || any(ssz == 0)) break
    a_new <- tryCatch(solve(Oi * crossprod(Z), colSums(Z * (Ys %*% Oi))), error = function(err) NULL)
    if (is.null(a_new)) break
    done <- max(abs(a_new - a)) < 1e-12
    a <- a_new
    if (done) break
  }
  a <- unname(a)
  e <- Ys - Z * matrix(a, n, 4, byrow = TRUE)
  A <- matrix(0, 4, 2)
  A[1, 1] <- a[1]
  A[2, 1] <- a[2]
  A[3, 2] <- a[3]
  A[4, 2] <- a[4]
  B <- cbind(c(1, -1, 0, 0), c(0, 0, 1, -1))
  list(A = A, B = B, CD = CD, C_delta = Cd, Omega = crossprod(e) / n, last = Y[n_t, ],
       lambda = c(alpha = 1 + a[1] - a[2], beta = 1 + a[3] - a[4]))
}

# Simulated paths of the trend indexes: per group a list of matrices
# [h, n_draws] alpha, beta, kappa (NULL when not identified) and zeta
saint_index_paths <- function(fit, h, n_draws) {
  dyn <- fit$dynamics
  rw0 <- function(last, sd) last + apply(matrix(stats::rnorm(h * n_draws, 0, sd), h), 2, cumsum)
  cs <- function(m) if (h == 1) m else apply(m, 2, cumsum)
  paths <- lapply(seq_along(fit$trends), function(g) {
    tr <- fit$trends[[g]]
    nt <- ncol(tr$theta)
    dg <- dyn$groups[[g]]
    ab <- NULL
    if (is.null(dyn$ecm)) {
      ab <- lapply(1:2, function(i) matrix(NA_real_, h, n_draws))
      inc <- draw_mvnorm(h * n_draws, dg$drift, dg$cov)
      for (i in 1:2) ab[[i]] <- tr$theta[i, nt] + matrix(cs(matrix(inc[i, ], h)), h)
    }
    list(alpha = ab[[1]], beta = ab[[2]],
         kappa = if (nrow(tr$theta) > 2) matrix(rw0(tr$theta[3, nt], dg$sd_kappa), h) else NULL,
         zeta = matrix(rw0(tr$zeta[nt], dg$sd_zeta), h))
  })
  ec <- dyn$ecm
  if (!is.null(ec)) {
    Yp <- matrix(ec$last, 4, n_draws)
    out <- array(NA_real_, c(4, h, n_draws))
    for (k in seq_len(h)) {
      eps <- draw_mvnorm(n_draws, rep(0, 4), ec$Omega)
      Yp <- Yp + ec$A %*% (crossprod(ec$B, Yp) - ec$CD) + ec$C_delta + eps
      out[, k, ] <- Yp
    }
    for (j in 1:2) {
      g <- if (j == 1) ec$lead else ec$other
      paths[[g]]$alpha <- matrix(out[j, , ], h)
      paths[[g]]$beta <- matrix(out[j + 2, , ], h)
    }
  }
  paths
}

# Forecast rates of the trend [age, h, draw] (eqs. 4.15-4.17): the baseline
# of the simulated indexes, mean frailty 1 / (1 + sigma2 M0) with M0 cumulated
# along each cohort from the jump-off year (4.16; the paper prints the
# baseline of year t_max in the second line of 4.16, a misprint of t - 1) and
# M0 = (exp(sigma2 M) - 1) / sigma2 in the jump-off year (Example 4.2), plus
# the background
saint_trend_forecast <- function(tr, path, h, n_draws) {
  na <- nrow(tr$M)
  nt <- ncol(tr$M)
  s2 <- tr$sigma2
  m0 <- if (s2 > 0) expm1(s2 * tr$M[, nt]) / s2 else tr$M[, nt]
  M0p <- matrix(m0, na, n_draws)
  Fp <- matrix(tr$baseline[, nt], na, n_draws)
  X <- tr$X
  out <- array(NA_real_, c(na, h, n_draws))
  for (k in seq_len(h)) {
    M0 <- rbind(0, M0p[-na, , drop = FALSE] + Fp[-na, , drop = FALSE])
    eta <- outer(X[, 1], path$alpha[k, ]) + outer(X[, 2], path$beta[k, ])
    if (!is.null(path$kappa)) eta <- eta + outer(X[, 3], path$kappa[k, ])
    Fk <- exp(eta)
    out[, k, ] <- Fk / (1 + s2 * M0) + matrix(exp(path$zeta[k, ]), na, n_draws, byrow = TRUE)
    M0p <- M0
    Fp <- Fk
  }
  out
}

# Spread -------------------------------------------------------------------------

# Knots x_i of the spread regressors r_i(x) = min(1, max(0, x_i - x) / 20)
# (Jarner and Jallbjorn 2022, eq. 5.11)
saint_knots <- c(40, 60, 80, 100, 120)

# Spread regressors (eq. 5.11) at the ages x, columns y1 ... y5
saint_regressors <- function(x) {
  R <- matrix(vapply(saint_knots, function(k) pmin(1, pmax(0, k - x) / 20), numeric(length(x))), length(x))
  colnames(R) <- paste0("y", seq_along(saint_knots))
  R
}

# Columns of the spread regressors identified on the ages x (all five on ages
# 20-100, the range of the paper): columns that are 0 at every age are left
# out, and of the columns that are 1 at every age only the one with the
# highest knot (the level) is kept
saint_spread_columns <- function(x) {
  R <- saint_regressors(x)
  keep <- colSums(abs(R)) > 0
  ones <- which(keep & apply(R == 1, 2, all))
  if (length(ones) > 1) keep[ones[-length(ones)]] <- FALSE
  idx <- unname(which(keep))
  q <- qr(R[, idx, drop = FALSE])
  sort(idx[q$pivot[seq_len(q$rank)]])
}

# Spread parameters y_t of each year: Poisson maximum likelihood of the
# deaths given the trend rates and exposures (eqs. 5.9-5.10); matrix
# [period, regressor]
saint_spread_years <- function(D, E, mu_ref, R) {
  y <- t(poisson_by_year(D, mu_ref * E, R, matrix(0, ncol(R), ncol(D))))
  colnames(y) <- colnames(R)
  y
}

# VAR(1) of the spread parameters without mean (eq. 5.12): A diagonal with
# 0 <= A_ii <= 0.99, each the least-squares coefficient of its component
# ("diagonal"); unrestricted by conditional maximum likelihood ("full",
# multivariate least squares); or fixed at a number times the identity.
# Omega is the covariance of the residuals (sums of squares over the number
# of changes)
saint_var_fit <- function(y, ar) {
  p <- ncol(y)
  n <- nrow(y) - 1
  Y1 <- y[-1, , drop = FALSE]
  Y0 <- y[-nrow(y), , drop = FALSE]
  bound <- rep(FALSE, p)
  if (is.numeric(ar)) {
    A <- diag(ar, p)
  } else if (ar == "diagonal") {
    ss <- colSums(Y0^2)
    raw <- ifelse(ss > 0, colSums(Y1 * Y0) / pmax(ss, 1e-300), 0)
    bound <- raw > 0.99
    A <- diag(pmin(pmax(raw, 0), 0.99), p)
  } else {
    At <- tryCatch(solve(crossprod(Y0), crossprod(Y0, Y1)), error = function(e) NULL)
    if (is.null(At)) {
      mt_abort(paste("The unrestricted VAR(1) of the spread cannot be estimated (singular cross-products); use",
                     "spread_ar = \"diagonal\"."), class = "saint", call = NULL)
    }
    A <- t(At)
  }
  e <- Y1 - Y0 %*% t(A)
  dimnames(A) <- list(colnames(y), colnames(y))
  Om <- crossprod(e) / n
  list(A = A, Omega = Om, bound = bound, radius = max(Mod(eigen(A, only.values = TRUE)$values)))
}

# Simulated spread parameters [regressor, h, draw] from the last estimate
# (2011, Section 2.3.3: y_{T+h} | y_T ~ N(A^h y_T, sum_i A^i Omega A^i'))
saint_var_paths <- function(sp, h, n_draws) {
  p <- ncol(sp$y)
  yp <- matrix(sp$y[nrow(sp$y), ], p, n_draws)
  out <- array(NA_real_, c(p, h, n_draws))
  for (k in seq_len(h)) {
    yp <- sp$A %*% yp + draw_mvnorm(n_draws, rep(0, p), sp$Omega)
    out[, k, ] <- yp
  }
  out
}

# Warnings of a fit (class mt_warning_saint)
saint_report <- function(trends, dyn, spreads, labels, ar) {
  which_trends <- function(flag) {
    paste(names(trends)[vapply(trends, function(tr) isTRUE(tr[[flag]]), NA)], collapse = ", ")
  }
  nc <- vapply(trends, function(tr) !isTRUE(tr$converged), logical(1))
  if (any(nc)) {
    mt_warn(sprintf("the estimation of the trend did not converge for %s.", paste(names(trends)[nc], collapse = ", ")),
            class = "saint")
  }
  if (nzchar(which_trends("sigma2_bound"))) {
    mt_warn(sprintf(paste("the frailty variance of the trend of %s is at the upper end of its search range; the trend",
                          "model may not suit the reference data."), which_trends("sigma2_bound")), class = "saint")
  }
  if (nzchar(which_trends("zeta_bound"))) {
    mt_warn(sprintf(paste("the background mortality of the trend of %s is at its lower bound (rate 1e-10) in some",
                          "years: the data do not identify it (for example rates smoothed over age, or no young adult",
                          "ages). Its random walk then has a large variance and forecasts at young ages are unstable;",
                          "consider background = FALSE."), which_trends("zeta_bound")), class = "saint")
  }
  ec <- dyn$ecm
  if (!is.null(ec) && any(abs(ec$lambda) >= 1)) {
    mt_warn(sprintf(paste("the error-correction model of the trends does not correct deviations from the",
                          "equilibrium (eigenvalues lambda_alpha = %.3f, lambda_beta = %.3f; Jarner and Jallbjorn",
                          "2022, eq. 5.8, need them below 1 in absolute value)."), ec$lambda[1], ec$lambda[2]),
            class = "saint")
  }
  if (identical(ar, "diagonal")) {
    b <- vapply(spreads, function(sp) any(sp$bound), logical(1))
    if (any(b)) {
      mt_warn(sprintf(paste("the autoregressive coefficient of a spread component is bounded at 0.99 in %d",
                            "stratum/strata (%s); their deviations from the trend persist for decades."), sum(b),
                      paste(utils::head(labels[b], 5), collapse = ", ")), class = "saint")
    }
  }
  if (identical(ar, "full")) {
    b <- vapply(spreads, function(sp) isTRUE(sp$radius >= 1), logical(1))
    if (any(b)) {
      mt_warn(sprintf(paste("the VAR(1) of the spread is not stationary in %d stratum/strata (%s); the",
                            "deviations from the trend are then unbounded. The papers restrict the model (Jarner and",
                            "Jallbjorn 2022, Section 5.4.1): spread_ar = \"diagonal\"."), sum(b),
                      paste(utils::head(labels[b], 5), collapse = ", ")), class = "saint")
    }
  }
  invisible()
}

# Numerics -------------------------------------------------------------------------

# Solve the symmetric positive definite systems H[, , t] b[, t] = g[, t] for
# all t at once (Cholesky by components; p is small)
solve_sym_batch <- function(H, g) {
  p <- dim(H)[1]
  n <- dim(H)[3]
  L <- array(0, c(p, p, n))
  for (j in seq_len(p)) {
    s <- H[j, j, ]
    if (j > 1) for (k in seq_len(j - 1)) s <- s - L[j, k, ]^2
    L[j, j, ] <- sqrt(pmax(s, 1e-300))
    if (j < p) {
      for (i in (j + 1):p) {
        s <- H[i, j, ]
        if (j > 1) for (k in seq_len(j - 1)) s <- s - L[i, k, ] * L[j, k, ]
        L[i, j, ] <- s / L[j, j, ]
      }
    }
  }
  z <- matrix(0, p, n)
  for (i in seq_len(p)) {
    s <- g[i, ]
    if (i > 1) for (k in seq_len(i - 1)) s <- s - L[i, k, ] * z[k, ]
    z[i, ] <- s / L[i, i, ]
  }
  b <- matrix(0, p, n)
  for (i in p:1) {
    s <- z[i, ]
    if (i < p) for (k in (i + 1):p) s <- s - L[k, i, ] * b[k, ]
    b[i, ] <- s / L[i, i, ]
  }
  b
}

# Poisson maximum likelihood of D[, t] ~ Poisson(W[, t] exp(X theta_t)) for
# every year t (columns), by Newton-Raphson with step halving from `theta`
# [p, nt]
poisson_by_year <- function(D, W, X, theta, tol = 1e-10, max_iter = 100) {
  p <- ncol(X)
  pairs <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  XX <- X[, pairs[, 1], drop = FALSE] * X[, pairs[, 2], drop = FALSE]
  ll <- function(th, j = seq_len(ncol(D))) {
    eta <- X %*% th
    colSums(D[, j, drop = FALSE] * eta - W[, j, drop = FALSE] * exp(eta))
  }
  cur <- ll(theta)
  for (it in seq_len(max_iter)) {
    mu <- W * exp(X %*% theta)
    g <- crossprod(X, D - mu)
    hv <- crossprod(XX, mu)
    H <- array(0, c(p, p, ncol(D)))
    for (k in seq_len(nrow(pairs))) {
      H[pairs[k, 1], pairs[k, 2], ] <- hv[k, ]
      H[pairs[k, 2], pairs[k, 1], ] <- hv[k, ]
    }
    step <- solve_sym_batch(H, g)
    new <- theta + step
    nl <- ll(new)
    bad <- !is.finite(nl) | nl < cur - 1e-12 * abs(cur)
    half <- 0
    while (any(bad) && half < 30) {
      step[, bad] <- step[, bad] / 2
      new[, bad] <- theta[, bad] + step[, bad]
      nl[bad] <- ll(new[, bad, drop = FALSE], which(bad))
      bad <- !is.finite(nl) | nl < cur - 1e-12 * abs(cur)
      half <- half + 1
    }
    theta <- new
    cur <- nl
    if (max(abs(step)) < tol) break
  }
  theta
}

# Rates of the trend model (4.10) with Gamma frailty (4.13) generated year by
# year, with M(t, x) cumulated from these rates as in (4.7)-(4.8); fitting
# deaths equal to these rates times the exposures returns the parameters
# (used by the tests). theta [p, nt] for the design saint_design(x).
saint_trend_rates <- function(theta, zeta, sigma2, x) {
  X <- saint_design(x)
  na <- length(x)
  nt <- length(zeta)
  Fx <- exp(X %*% theta)
  G <- exp(zeta)
  mu <- M <- matrix(0, na, nt)
  for (t in seq_len(nt)) {
    if (t == 1) {
      for (i in seq_len(na)) {
        if (i > 1) M[i, 1] <- M[i - 1, 1] + mu[i - 1, 1] - G[1]
        mu[i, 1] <- exp(-sigma2 * M[i, 1]) * Fx[i, 1] + G[1]
      }
    } else {
      if (na > 1) M[-1, t] <- M[-na, t - 1] + mu[-na, t - 1] - G[t - 1]
      mu[, t] <- exp(-sigma2 * M[, t]) * Fx[, t] + G[t]
    }
  }
  mu
}

# Parameters of a SAINT fit (mt_parameters()): per trend the indexes, the
# frailty variance and the time-series parameters; per stratum with a spread
# the spread parameters and their VAR(1)
saint_parameters <- function(fit) {
  f <- fit$fits
  lab <- fit$strata$label
  yrs <- as.character(fit$periods)
  row <- function(stratum, parameter, index, value) {
    data.table(stratum = stratum, parameter = parameter, index = as.character(index), value = as.numeric(value))
  }
  pairs <- function(nm) {
    k <- which(upper.tri(diag(length(nm)), diag = TRUE), arr.ind = TRUE)
    list(i = k[, 1], j = k[, 2], index = paste(nm[k[, 1]], nm[k[, 2]], sep = ","))
  }
  ec <- f$dynamics$ecm
  trend <- rbindlist(lapply(seq_along(f$trends), function(g) {
    tr <- f$trends[[g]]
    dg <- f$dynamics$groups[[g]]
    gl <- paste("trend", f$groups[g])
    th <- tr$theta
    bg <- all(is.finite(tr$zeta))
    rbind(row(gl, "alpha", yrs, th[1, ]), row(gl, "beta", yrs, th[2, ]),
          if (nrow(th) > 2) row(gl, "kappa", yrs, th[3, ]), if (bg) row(gl, "zeta", yrs, tr$zeta),
          row(gl, "sigma2", "", tr$sigma2), row(gl, "loglik", "", tr$loglik),
          if (is.null(ec)) {
            rbind(row(gl, "drift", c("alpha", "beta"), dg$drift),
                  row(gl, "cov", c("alpha,alpha", "alpha,beta", "beta,beta"), dg$cov[c(1, 3, 4)]))
          },
          if (nrow(th) > 2) row(gl, "sd_kappa", "", dg$sd_kappa), if (bg) row(gl, "sd_zeta", "", dg$sd_zeta))
  }))
  ecm <- NULL
  if (!is.null(ec)) {
    g <- f$groups[c(ec$lead, ec$other)]
    nm <- c(paste0("alpha_", g), paste0("beta_", g))
    pr <- pairs(nm)
    el <- "trend error correction"
    ecm <- rbind(row(el, c("a11", "a21", "a32", "a42"), "", ec$A[cbind(1:4, c(1, 1, 2, 2))]),
                 row(el, "C_D", c("alpha", "beta"), ec$CD),
                 row(el, "C_delta", c("alpha", "beta"), ec$C_delta[c(1, 3)]),
                 row(el, "lambda", c("alpha", "beta"), ec$lambda),
                 row(el, "Omega", pr$index, ec$Omega[cbind(pr$i, pr$j)]))
  }
  spread <- rbindlist(lapply(seq_along(f$spreads), function(s) {
    sp <- f$spreads[[s]]
    if (is.null(sp)) return(NULL)
    nm <- colnames(sp$y)
    pr <- pairs(nm)
    off <- sp$A - diag(diag(sp$A), nrow(sp$A))
    rbind(rbindlist(lapply(seq_along(nm), function(i) row(lab[s], nm[i], yrs, sp$y[, i]))),
          if (any(off != 0)) {
            k <- expand.grid(i = seq_along(nm), j = seq_along(nm))
            row(lab[s], "A", paste(nm[k$i], nm[k$j], sep = ","), sp$A[cbind(k$i, k$j)])
          } else {
            row(lab[s], "ar", nm, diag(sp$A))
          },
          row(lab[s], "Omega", pr$index, sp$Omega[cbind(pr$i, pr$j)]))
  }))
  rbind(trend, ecm, spread)
}
