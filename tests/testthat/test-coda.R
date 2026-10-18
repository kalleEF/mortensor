# Compositional data (CoDa) models: coda(), coda(coherent = TRUE), coda_causes()

h3_coda <- function() mt_hierarchy(sex = "both", region = "all", quintile = "all")
closure_cols <- function(M) sweep(M, 2, colSums(M), "/")
closure_rows <- function(M) M / rowSums(M)

# The published algorithm, transcribed from model.Oeppen(), predict.Oeppen()
# and get_dx_values() of MortalityForecast 0.8.1 (GPL-3, Marius D. Pascariu
# and Marie-Pier Bergeron-Boucher), with base R for the compositions package
# (acomp() closes rows, clr(), clrInv(), geometricmeanCol()). `dx`: [part,
# year]; `B`: the benchmark population's fit (coherent model). The index is
# forecast from the last year (the reference starts one step later, see
# ?coda), ARIMA(0,1,0) with drift: the ML drift is the mean change.
ref_oeppen <- function(dx, h, B = NULL) {
  clr <- function(M) log(M) - rowMeans(log(M))
  clr_inv <- function(M) closure_rows(exp(M))
  dx <- closure_rows(t(dx))
  ax <- exp(colMeans(log(dx)))
  ax <- ax / sum(ax)
  cdx <- sweep(dx, 2, ax, "/")
  if (!is.null(B)) cdx <- cdx / clr_inv(B$kt %*% t(B$bx))
  S <- svd(clr(closure_rows(cdx)))
  kt <- S$d[1] * S$u[, 1]
  bx <- S$v[, 1]
  bench <- function(k) if (is.null(B)) 1 else clr_inv(k %*% t(B$bx))
  fitted <- closure_rows(sweep(clr_inv(kt %*% t(bx)) * bench(B$kt), 2, ax, "*"))
  list(ax = ax, bx = bx, kt = kt, dx = dx, fitted = t(fitted), bench = bench)
}

ref_forecast <- function(r, kf, Bkf = NULL) {
  p <- closure_rows(sweep(closure_rows(exp(kf %*% t(r$bx))) * r$bench(Bkf), 2, r$ax, "*"))
  J <- r$dx[nrow(r$dx), ] / p[1, ]
  t(closure_rows(sweep(p, 2, J, "*"))[-1, , drop = FALSE])
}

test_that("CoDa models satisfy the model contract", {
  expect_mt_model(coda())
  expect_mt_model(coda(jump_off = "fitted"))
  expect_mt_model(coda(coherent = TRUE))
  x <- mt_example_causes(years = 1998:2012)
  expect_mt_model(coda_causes(), data = x)
})

test_that("CoDa constructors check their arguments", {
  expect_s3_class(coda(), "mt_spec_coda")
  expect_false(coda()$joint)
  expect_true(coda(common_by = "sex")$params$coherent)
  expect_true(coda(common_by = "sex")$joint)
  expect_match(coda(common_by = "sex")$label, "common by sex")
  expect_error(coda(common_by = "sex", coherent = FALSE), class = "mt_error_argument")
  expect_error(coda(coherent = "yes"), class = "mt_error_argument")
  expect_error(coda(jump_off = "actual"))
  expect_error(coda_causes(residual = "sometimes"), class = "mt_error_argument")
  expect_error(coda_causes(var = c("a", "b")), class = "mt_error_argument")
})

test_that("coda() reproduces the published algorithm on England & Wales men", {
  skip_if_not_installed("StMoMo")
  ew <- as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo"))
  fit <- mt_fit(coda(), ew)
  f <- fit$fits[[1]]
  ages <- fit$age
  comp <- closure_cols(exp(matrix(mortensor:::coda_forward(fit$observed[, , 1], ages), ncol = length(fit$periods))))
  r <- ref_oeppen(comp, 0)
  expect_equal(closure_cols(exp(f$alpha + outer(f$b, f$k))), r$fitted, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(f$fitted, mortensor:::coda_backward(log(r$fitted), ages), tolerance = 1e-10, ignore_attr = TRUE)
  # point forecast: draws without innovations and drift uncertainty
  h <- 15
  kt <- r$kt
  kf <- kt[length(kt)] + (0:h) * mean(diff(kt))
  fit$fits[[1]]$sigma_k <- 0
  fc <- mt_forecast(fit, h = h, n_draws = 2, seed = 1)
  expect_equal(fc$rates[, , 1, 1], mortensor:::coda_backward(log(ref_forecast(r, kf)), ages), tolerance = 1e-10,
               ignore_attr = TRUE)
  # the median of the draws is close to the point forecast
  fit$fits[[1]] <- f
  fc <- mt_forecast(fit, h = 5, n_draws = 1000, seed = 1)
  ratio <- apply(fc$rates[, , 1, ], 1:2, stats::median) / mortensor:::coda_backward(log(ref_forecast(r, kf[1:6])),
                                                                                     ages)
  expect_lt(stats::median(abs(ratio - 1)), 0.01)
})

test_that("the coherent CoDa model reproduces the published algorithm on French data", {
  skip_if_not_installed("demography")
  fr <- as_mt_data(mortensor:::pkg_data("fr.mort", "demography"))
  fr <- mt_window(fr, from = 1960, to = 2006)
  fr <- mt_regroup(fr, age = single_ages(0, 100, open = TRUE))
  fit <- mt_fit(coda(coherent = TRUE), fr)
  ages <- fit$age
  nt <- length(fit$periods)
  comp_of <- function(m) closure_cols(exp(matrix(mortensor:::coda_forward(m, ages), ncol = ncol(m))))
  comps <- lapply(1:2, function(s) comp_of(mortensor:::coda_positive_rates(fit$observed[, , s])))
  # benchmark: closed geometric mean of the two distributions
  bench <- closure_cols(sqrt(comps[[1]] * comps[[2]]))
  B <- ref_oeppen(bench, 0)
  h <- 10
  KT <- B$kt[nt]
  Kf <- KT + (0:h) * mean(diff(B$kt))
  # no innovations, no drift uncertainty: point forecasts
  fit$fits$commons[[1]]$sigma_k <- 0
  for (s in 1:2) fit$fits$strata[[s]]$ar$sigma <- 0
  fc <- mt_forecast(fit, h = h, n_draws = 2, seed = 1)
  for (s in 1:2) {
    r <- ref_oeppen(comps[[s]], 0, B = B)
    st <- fit$fits$strata[[s]]
    cm <- fit$fits$commons[[1]]
    ours <- closure_cols(exp(st$alpha + outer(cm$b, cm$k) + outer(st$b, st$k)))
    expect_equal(ours, r$fitted, tolerance = 1e-10, ignore_attr = TRUE)
    # deviations: AR(1) with mean by maximum likelihood, as forecast::Arima(order = c(1, 0, 0))
    a <- stats::arima(r$kt, order = c(1, 0, 0), method = "ML")
    expect_equal(abs(st$ar$phi), abs(stats::coef(a)[["ar1"]]), tolerance = 1e-6)
    mu <- stats::coef(a)[["intercept"]]
    kf <- c(r$kt[nt], mu + stats::coef(a)[["ar1"]]^(1:h) * (r$kt[nt] - mu))
    expect_equal(comp_of(matrix(fc$rates[, , s, 1], length(ages))), ref_forecast(r, kf, Kf), tolerance = 1e-8,
                 ignore_attr = TRUE)
  }
})

test_that("CoDa fits handle zero deaths, closed age ranges and open intervals", {
  m <- matrix(c(0, 0.01, 0.02, 0.04), 2)
  expect_warning(r <- mortensor:::coda_positive_rates(m), class = "mt_warning_sparse")
  expect_equal(r, matrix(c(0.005, 0.01, 0.02, 0.04), 2))
  x <- mt_example(deaths = "single", years = 1995:2012)
  fit <- mt_fit(coda(), x)
  f <- fit$fits[[1]]
  expect_equal(length(f$b), length(fit$age))
  expect_equal(sum(f$b), 0, tolerance = 1e-10)
  # the open interval follows its regression on k
  expect_equal(log(f$fitted[length(fit$age), ]), f$open$a + f$open$b * f$k, tolerance = 1e-10)
  # a closed age range adds the survivors as a part
  xc <- mt_set_grid(x, age = single_ages(30, 99, open = FALSE))
  xc <- mt_regroup(xc, age = single_ages(30, 99, open = FALSE))
  fc <- mt_fit(coda(), xc)
  expect_equal(length(fc$fits[[1]]$b), length(fc$age) + 1)
  expect_true(all(is.finite(fitted(fc)) & fitted(fc) > 0))
  p <- mortensor:::coda_parameters(fit)
  expect_true(all(c("a", "b", "k", "drift_k", "sd_k", "a_open", "b_open") %in% p$parameter))
  expect_equal(sum(p[stratum == "female" & parameter == "a"]$value), 1)
  expect_type(mortensor:::coda_assumptions(coda()), "character")
})

test_that("the observed jump-off starts at the observed deaths", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  fit <- mt_fit(coda(), x)
  fit$fits[[1]]$sigma_k <- 0
  fit$fits[[1]]$drift <- 0
  fc <- mt_forecast(fit, h = 2, n_draws = 1, seed = 1)
  obs <- fit$observed[, length(fit$periods), 1]
  expect_equal(fc$rates[-length(fit$age), 1, 1, 1], obs[-length(obs)], tolerance = 1e-8, ignore_attr = TRUE)
  expect_equal(fc$rates[length(fit$age), 1, 1, 1], obs[[length(obs)]], tolerance = 1e-8)
})

test_that("the coherent CoDa model pools bottom strata within groups", {
  x <- mt_example_regions(years = 2000:2014, regions = 2)
  xa <- mt_aggregate_strata(x, h3_coda())
  fit <- mt_fit(coda(common_by = "sex"), xa)
  cm <- fit$fits$commons
  expect_length(cm, 3)
  st <- fit$strata
  bottom <- which(st$sex != "both" & st$region != "all" & st$quintile != "all")
  expect_true(all(unlist(lapply(cm[unique(fit$fits$group[bottom])], `[[`, "members")) %in% bottom))
  fc <- mt_forecast(fit, h = 3, n_draws = 10, seed = 1)
  expect_true(all(is.finite(fc$rates) & fc$rates > 0))
  p <- mortensor:::coda_parameters(fit)
  expect_true(all(c("K", "B", "phi", "mean_k") %in% p$parameter))
  expect_error(mt_fit(coda(common_by = "nonsense"), xa), class = "mt_error_model_fit")
})

test_that("coda_causes() forecasts causes that add up to all causes in every draw", {
  x <- mt_example_causes(years = 1998:2012)
  fit <- mt_fit(coda_causes(), x)
  expect_s3_class(fit, "mt_fit_coda_causes")
  expect_false(fit$fits$residual)
  st <- fit$strata
  gap <- function(arr) {
    max(vapply(c("female", "male"), function(sx) {
      i <- st$stratum[st$sex == sx & st$cause == "all"]
      j <- st$stratum[st$sex == sx & st$cause != "all"]
      sums <- if (length(dim(arr)) == 4) apply(arr[, , j, , drop = FALSE], c(1, 2, 4), sum) else
        apply(arr[, , j, drop = FALSE], 1:2, sum)
      tot <- if (length(dim(arr)) == 4) arr[, , i, ] else arr[, , i]
      max(abs(sums - tot) / pmax(tot, 1e-12))
    }, 0))
  }
  expect_lt(gap(fitted(fit)), 1e-12)
  fc <- mt_forecast(fit, h = 4, n_draws = 30, seed = 1)
  expect_lt(gap(fc$rates), 1e-12)
  pr <- mt_forecast(fit, h = 4, n_draws = 30, seed = 1, predictive = TRUE)
  expect_true(pr$predictive)
  expect_lt(gap(pr$rates), 1e-12)
  expect_identical(pr$rates, mt_forecast(fit, h = 4, n_draws = 30, seed = 1, predictive = TRUE)$rates)
  expect_false(isTRUE(all.equal(pr$rates, fc$rates)))
  fc2 <- mt_forecast(fit, to = 2014, n_draws = 5, seed = 1)
  expect_equal(fc2$periods, 2013:2014)
  p <- mortensor:::coda_parameters(fit)
  expect_true(all(c("a", "b", "k") %in% p$parameter))
  bt <- mt_backtest(coda_causes(), x, origins = 2010, h = 2, n_draws = 20)
  expect_true(all(is.finite(bt$crps)))
})

test_that("coda_causes() adds a residual cause when the causes do not cover all deaths", {
  x <- mt_example_causes(years = 1998:2012)
  xs <- mt_subset(x, cause != "other")
  fit <- mt_fit(coda_causes(), xs)
  expect_true(fit$fits$residual)
  fc <- mt_forecast(fit, h = 3, n_draws = 10, seed = 1)
  st <- fit$strata
  i <- st$stratum[st$sex == "female" & st$cause == "all"]
  j <- st$stratum[st$sex == "female" & st$cause != "all"]
  expect_true(all(apply(fc$rates[, , j, , drop = FALSE], c(1, 2, 4), sum) < fc$rates[, , i, ]))
  expect_error(mt_fit(coda_causes(var = "nonsense"), x), class = "mt_error_model_fit")
  expect_error(mt_fit(coda_causes(total_level = "none"), x), class = "mt_error_model_fit")
})
