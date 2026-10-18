# Compositional data (CoDa) models: coda(), coda(coherent = TRUE), coda_causes()

h3_coda <- function() mt_hierarchy(sex = "both", region = "all", quintile = "all")
closure_cols <- function(M) sweep(M, 2, colSums(M), "/")

# French women and men 1960-2006 (HMD via demography), single ages 0-109
# and 110+, as in the data of Bergeron-Boucher et al. (2017)
fr_coda <- function() {
  fr <- as_mt_data(mortensor:::pkg_data("fr.mort", "demography"))
  fr <- mt_window(fr, from = 1960, to = 2006)
  mt_regroup(fr, age = single_ages(0, 110, open = TRUE))
}

test_that("CoDa models satisfy the model contract", {
  expect_mt_model(coda())
  expect_mt_model(coda(jump_off = "fitted"))
  expect_mt_model(coda(open = "regression"))
  expect_mt_model(coda(order = c(0, 1, 0)))
  expect_mt_model(coda(coherent = TRUE))
  expect_mt_model(coda(coherent = TRUE, order_dev = c(1, 0, 0), open = "regression"))
  x <- mt_example_causes(years = 1998:2012)
  expect_mt_model(coda_causes(), data = x)
  expect_mt_model(coda_causes(rank = 2, order = c(0, 2, 2)), data = x)
})

test_that("CoDa constructors check their arguments", {
  expect_s3_class(coda(), "mt_spec_coda")
  expect_false(coda()$joint)
  expect_identical(coda()$label, "CoDa")
  expect_equal(coda()$params$order, c(0, 1, 1))
  expect_equal(coda(coherent = TRUE)$params$order_dev, c(1, 1, 0))
  expect_match(coda(order = c(0, 1, 0))$label, "ARIMA\\(0,1,0\\) with drift")
  expect_match(coda(coherent = TRUE, order_dev = c(1, 0, 0))$label, "deviations ARIMA\\(1,0,0\\)")
  expect_true(coda(common_by = "sex")$params$coherent)
  expect_true(coda(common_by = "sex")$joint)
  expect_match(coda(common_by = "sex")$label, "common by sex")
  expect_error(coda(common_by = "sex", coherent = FALSE), class = "mt_error_argument")
  expect_error(coda(coherent = "yes"), class = "mt_error_argument")
  expect_error(coda(jump_off = "actual"))
  expect_error(coda(open = "none"))
  expect_error(coda(order = c(0, 3, 1)), class = "mt_error_argument")
  expect_error(coda(order = c(0, 1)), class = "mt_error_argument")
  expect_error(coda(order_dev = c(-1, 1, 0)), class = "mt_error_argument")
  expect_error(coda(drift = NA), class = "mt_error_argument")
  expect_error(coda_causes(residual = "sometimes"), class = "mt_error_argument")
  expect_error(coda_causes(var = c("a", "b")), class = "mt_error_argument")
  expect_error(coda_causes(rank = 0), class = "mt_error_argument")
  expect_match(coda_causes(rank = 2, order = c(0, 2, 2))$label, "rank 2, ARIMA\\(0,2,2\\)")
})

test_that("coda() reproduces the code of Bergeron-Boucher et al. (2017) on French women", {
  skip_if_not_installed("demography")
  # Reference values computed with the functions fit120() and CoDa() of
  # CoDaForecast.R, the supplementary code of Bergeron-Boucher et al. (2017,
  # Demographic Research 37(17)), evaluated from the authors' file (no licence
  # stated; not copied here) with forecast 8.21.1 and compositions 2.x, on
  # demography::fr.mort, women, 1960-2006, single ages 0-110+: fit120() on
  # the deaths and exposures, and CoDa(dx, t = 20) on the life-table deaths
  # of the parts below (ages 0-120), printed with signif(, 12). Parts 1, 41,
  # 81, 101 and 121 are ages 0, 40, 80, 100 and 120.
  ref_kannisto_1960 <- c(0.09274239225, 0.256822874, 0.5387981525, 0.7979544278, 0.9303197355)
  ref_kannisto_2006 <- c(0.03271236138, 0.1271036988, 0.3853540012, 0.7296890054, 0.9207787119)
  ref_ax <- c(0.00863262654749, 0.00144030162659, 0.03471334435, 0.00343439298251, 1.81196196789e-09)
  ref_bx <- c(-0.114609432179, -0.044536398492, -0.0335710738563, 0.126751530823, 0.195834062631)
  ref_kt <- c(-9.97561328784, 10.9905207688)
  ref_dx_h1 <- c(0.00306079502309, 0.000903897911362, 0.0232898996723, 0.0131308603612, 2.17031741396e-08)
  ref_dx_h20 <- c(0.000985043500515, 0.000504394303547, 0.0141652020898, 0.0281340845775, 8.00043751111e-08)
  fr <- fr_coda()
  D <- mt_array(fr, "deaths")
  E <- mt_array(fr, "exposure")
  layout <- mortensor:::coda_layout(mt_grid(fr)$age)
  expect_identical(layout$method, "kannisto")
  expect_length(layout$parts, 121)
  sc <- mortensor:::coda_schedule(D[, , 1], E[, , 1], layout, "women")
  # Kannisto model by Poisson maximum likelihood: the authors' Nelder-Mead
  # optimum differs by up to 0.15% (its log-likelihood is lower in every year)
  expect_equal(sc$m[c(81, 91, 101, 111, 121), 1], ref_kannisto_1960, tolerance = 2e-3)
  expect_equal(sc$m[c(81, 91, 101, 111, 121), 47], ref_kannisto_2006, tolerance = 2e-3)
  expect_equal(sc$m[1:80, ], D[1:80, , 1] / E[1:80, , 1], ignore_attr = TRUE)
  fit <- mortensor:::coda_core(mortensor:::coda_clr(sc$m, layout), 1, c(0, 1, 1), TRUE)
  pick <- c(1, 41, 81, 101, 121)
  expect_equal(exp(fit$alpha[pick]) / sum(exp(fit$alpha)), ref_ax, tolerance = 1e-9)
  expect_equal(fit$b[pick, 1], ref_bx, tolerance = 1e-9)
  expect_equal(fit$k[c(1, 47), 1], ref_kt, tolerance = 1e-9)
  # point forecast: ARIMA(0,1,1) with drift, eq. (20), observed jump-off (17)
  fit$jump_off <- "observed"
  pt <- mortensor:::coda_point(fit, 20)
  dx <- closure_cols(exp(pt$Z))
  expect_equal(dx[pick, 1], ref_dx_h1, tolerance = 1e-9)
  expect_equal(dx[pick, 20], ref_dx_h20, tolerance = 1e-9)
  # eq. (20): the corrected forecast is the last index plus the drift per year
  expect_equal(pt$k[, 1], fit$k[47, 1] + (1:20) * fit$ts[[1]]$const)
  # the same through the fit method of coda() (mt_fit() needs positive
  # exposures at every age of the grid, which ages 101-109 lack in 1960)
  f2 <- fit_model(coda(), D[, , 1], E[, , 1], mt_grid(fr)$age, 1960:2006)
  expect_equal(f2$b, fit$b)
  expect_equal(f2$k, fit$k)
})

test_that("the coherent CoDa model reproduces the code of Bergeron-Boucher et al. (2017)", {
  skip_if_not_installed("demography")
  # Reference values from CoDa.CH() of CoDaForecast.R (see the previous test),
  # women (1) and men (2) with the common factor of CoDa() fitted to the
  # life-table deaths of the mean of the rates of women and men (Appendix C):
  # deviation loadings, deviation index (first and last year) and the
  # forecast life-table deaths after 20 years at parts 1, 41, 81, 101, 121.
  ref_bx <- list(c(0.171414408323, 0.00103623127191, -0.0607034243649, -0.143259080825, 0.179627345153),
                 c(0.0335671616486, -0.0318254460927, -0.0628773714666, -0.0812786131215, 0.312746494013))
  ref_kt <- list(c(-1.09439267857, 1.41257149175), c(-1.29821931168, 0.518046987371))
  ref_dx_h20 <- list(c(0.000818753130233, 0.000463778684248, 0.0145212933125, 0.0268818865773, 6.60622479597e-08),
                     c(0.00145783705354, 0.00127238067568, 0.0243044231884, 0.0102564193476, 2.08428848694e-08))
  fr <- fr_coda()
  f <- fit_model(coda(coherent = TRUE), mt_array(fr, "deaths"), mt_array(fr, "exposure"), mt_grid(fr)$age,
                 1960:2006, strata = mt_strata(fr))
  pick <- c(1, 41, 81, 101, 121)
  expect_length(f$commons, 1)
  for (s in 1:2) {
    st <- f$strata[[s]]
    expect_equal(st$b[pick, 1], ref_bx[[s]], tolerance = 1e-9)
    expect_equal(st$k[c(1, 47), 1], ref_kt[[s]], tolerance = 1e-9)
    pt <- mortensor:::coda_point_coherent(f, s, 20)
    expect_equal(closure_cols(exp(pt$Z))[pick, 20], ref_dx_h20[[s]], tolerance = 1e-9)
    # deviations: ARIMA(1,1,0) without drift; the point forecast levels off
    expect_equal(st$ts[[1]]$order, c(1, 1, 0))
    expect_lt(abs(diff(pt$k[19:20, 1])), abs(diff(pt$k[1:2, 1])))
  }
  fit <- mt_fit(coda(coherent = TRUE), mt_regroup(fr, age = single_ages(0, 100, open = TRUE)))
  p <- mt_parameters(fit)
  expect_true(all(c("K", "B", "arima_K", "order_k", "arima_k", "sd_k", "kannisto_a") %in% p$parameter))
  expect_equal(p[parameter == "order_k" & stratum == "female"]$value, c(1, 1, 0))
})

test_that("ARIMA models of the index match stats::arima() and forecast::Arima()", {
  set.seed(3)
  k <- cumsum(c(0, -1 + stats::arima.sim(list(ma = 0.4), 30)))
  f <- mortensor:::index_arima_fit(k, c(0, 1, 1), TRUE)
  a <- stats::arima(k, order = c(0, 1, 1), xreg = cbind(drift = seq_along(k)), method = "ML")
  expect_equal(f$coef, stats::coef(a))
  mean_ours <- mortensor:::index_arima_paths(f, 5, 1, innovations = "none")
  mean_ref <- stats::predict(a, n.ahead = 5, newxreg = cbind(drift = length(k) + 1:5))$pred
  expect_equal(as.vector(mean_ours), as.vector(mean_ref), tolerance = 1e-6)
  # AR(1) of the changes without drift, and AR(1) with a mean
  for (ord in list(c(1, 1, 0), c(1, 0, 0))) {
    f <- mortensor:::index_arima_fit(k, ord, FALSE)
    a <- stats::arima(k, order = ord, method = "ML", include.mean = ord[2] == 0)
    expect_equal(f$coef, stats::coef(a))
    expect_equal(as.vector(mortensor:::index_arima_paths(f, 4, 1, innovations = "none")),
                 as.vector(stats::predict(a, n.ahead = 4)$pred), tolerance = 1e-8)
  }
  skip_if_not_installed("forecast")
  g <- forecast::Arima(k, order = c(0, 1, 1), include.drift = TRUE, method = "ML")
  expect_equal(unname(mortensor:::index_arima_fit(k, c(0, 1, 1), TRUE)$coef), unname(stats::coef(g)))
  # simulated paths start from the last value and spread with the horizon
  p <- with_seed_test(1, mortensor:::index_arima_paths(mortensor:::index_arima_fit(k, c(0, 1, 0), TRUE), 10, 2000))
  expect_equal(mean(p[10, ]), k[length(k)] + 10 * mean(diff(k)), tolerance = 0.05)
  expect_gt(stats::sd(p[10, ]), stats::sd(p[1, ]))
  # a series too short for the model falls back to the random walk with drift
  f <- mortensor:::index_arima_fit(1:4 + c(0, 0.1, -0.1, 0), c(0, 1, 1), TRUE)
  expect_true(f$fallback)
  expect_equal(f$const, mean(diff(1:4 + c(0, 0.1, -0.1, 0))))
})

test_that("zero deaths are replaced multiplicatively, keeping the total of each year", {
  D <- matrix(c(0, 4, 6, 10, 2, 0, 8, 10), 4)
  expect_warning(r <- mortensor:::coda_replace_zeros(D, "test"), class = "mt_warning_sparse")
  # delta = (2 / 2) / 20 = 0.05 of the deaths of each year
  expect_equal(colSums(r), colSums(D))
  expect_equal(r[1, 1], 0.05 * 20)
  expect_equal(r[-1, 1], D[-1, 1] * (1 - 0.05))
  expect_equal(r[2, 2], 0.05 * 20)
  expect_identical(mortensor:::coda_replace_zeros(D + 1, "test"), D + 1)
  expect_error(mortensor:::coda_replace_zeros(matrix(0, 2, 2), "test"), class = "mt_error_model_fit")
  expect_error(suppressWarnings(mortensor:::coda_replace_zeros(cbind(c(1, 2), 0), "test")),
               class = "mt_error_model_fit")
})

test_that("the Kannisto model is the Poisson maximum likelihood fit", {
  x <- 0:20
  mu <- stats::plogis(-2.5 + 0.11 * x)
  E <- 1e4 * exp(-0.12 * x)
  D <- round(E * mu)
  th <- mortensor:::kannisto_ml(D, E, x)
  m <- stats::plogis(th[1] + th[2] * x)
  # score equations of the Poisson likelihood
  expect_lt(max(abs(c(sum((1 - m) * (D - E * m)), sum(x * (1 - m) * (D - E * m))))), 1e-6)
  nll <- function(p) -sum(D * log(stats::plogis(p[1] + p[2] * x)) - E * stats::plogis(p[1] + p[2] * x))
  expect_lte(nll(th), stats::optim(c(-2, 0.1), nll)$value + 1e-8)
  expect_true(all(is.na(mortensor:::kannisto_ml(0 * D, E, x))))
  # a year without deaths at ages 80 and over takes the fit to all years
  Dm <- cbind(D, 0)
  expect_warning(k <- mortensor:::coda_kannisto(Dm, cbind(E, E), 80 + x, "test"), class = "mt_warning_sparse")
  expect_equal(k$par[, 2], mortensor:::kannisto_ml(rowSums(Dm), 2 * E, x))
})

test_that("the composition covers ages 80-120 where single ages allow it, otherwise the regression", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  ages <- mt_grid(x)$age
  lay <- mortensor:::coda_layout(ages)
  expect_identical(lay$method, "kannisto")
  expect_identical(lay$parts, c(as.character(30:79), as.character(80:120)))
  expect_identical(mortensor:::coda_layout(ages, "regression")$method, "regression")
  grouped <- age_groups(c(0, 1, seq(5, 90, 5)))
  expect_identical(mortensor:::coda_layout(grouped)$method, "regression")
  expect_error(mortensor:::coda_layout(grouped, "kannisto"), class = "mt_error_model_data")
  expect_identical(mortensor:::coda_layout(single_ages(30, 85, open = FALSE))$method, "closed")
  expect_null(mortensor:::coda_kannisto_ages(single_ages(0, 85)))
  # rates back from the composition: the open interval gets the deaths over
  # the person-years of the parts from its lower bound
  fit <- mt_fit(coda(), x)
  f <- fit$fits[[1]]
  expect_equal(length(f$b), 91)
  expect_equal(sum(f$b), 0, tolerance = 1e-10)
  expect_gt(f$k[nrow(f$k), 1], f$k[1, 1])
  P <- closure_cols(exp(f$alpha + f$b %*% t(f$k)))
  l <- apply(P, 2, function(p) rev(cumsum(rev(p))))
  L <- l - P / 2
  op <- 71:91
  expect_equal(f$fitted[71, ], colSums(P[op, ]) / colSums(L[op, ]), ignore_attr = TRUE)
  expect_equal(f$fitted[51, ], P[51, ] / L[51, ], ignore_attr = TRUE)
  # the composition of the observed schedule is recovered from its clr
  sc <- mortensor:::coda_schedule(fit$observed[, , 1] * 1e5, matrix(1e5, 71, 18), lay, "test")
  back <- mortensor:::coda_layout_rates(mortensor:::coda_clr(sc$m, lay), lay)
  expect_equal(back[1:70, ], sc$m[1:70, ], tolerance = 1e-10, ignore_attr = TRUE)
  # a closed age range without the extension adds the survivors as a part
  xc <- mt_set_grid(x, age = single_ages(30, 85, open = FALSE))
  xc <- mt_regroup(xc, age = single_ages(30, 85, open = FALSE))
  fc <- mt_fit(coda(), xc)
  expect_equal(length(fc$fits[[1]]$b), length(fc$age) + 1)
  expect_true(all(is.finite(fitted(fc)) & fitted(fc) > 0))
  # regression of the open interval on the index
  fr <- mt_fit(coda(open = "regression"), x)$fits[[1]]
  expect_equal(log(fr$fitted[71, ]), fr$open$a + fr$open$b * fr$k[, 1], tolerance = 1e-10)
  p <- mortensor:::coda_parameters(mt_fit(coda(open = "regression"), x))
  expect_true(all(c("a", "b", "k", "order_k", "arima_k", "sd_k", "a_open", "b_open") %in% p$parameter))
  expect_equal(sum(p[stratum == "female" & parameter == "a"]$value), 1)
  expect_type(mortensor:::coda_assumptions(coda()), "character")
  expect_type(mortensor:::coda_assumptions(coda(coherent = TRUE)), "character")
})

test_that("the observed jump-off starts at the observed deaths and draws are centred on the point forecast", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  for (open in c("regression", "kannisto")) {
    fit <- mt_fit(coda(open = open), x)
    f <- fit$fits[[1]]
    f$ts[[1]]$const <- 0
    f$ts[[1]]$theta <- 0
    pt <- mortensor:::coda_point(f, 2)
    expect_equal(pt$k[, 1], rep(f$k[nrow(f$k), 1], 2))
    m <- mortensor:::coda_layout_rates(pt$Z, f$layout, pt$open)
    obs <- fit$observed[, length(fit$periods), 1]
    if (open == "regression") {
      expect_equal(m[, 1], obs, tolerance = 1e-8, ignore_attr = TRUE)
    } else {
      expect_equal(m[1:50, 1], obs[1:50], tolerance = 1e-8, ignore_attr = TRUE)
    }
  }
  # the median of the draws in the first forecast year is the point forecast
  fit <- mt_fit(coda(), x)
  f <- fit$fits[[1]]
  fc <- mt_forecast(fit, h = 5, n_draws = 401, seed = 1)
  pt <- mortensor:::coda_point(f, 5)
  ptr <- mortensor:::coda_layout_rates(pt$Z, f$layout)
  med <- apply(fc$rates[, , 1, ], 1:2, stats::median)
  expect_lt(max(abs(med[, 1] / ptr[, 1] - 1)), 0.01)
  expect_lt(stats::median(abs(med[, 5] / ptr[, 5] - 1)), 0.02)
  # intervals widen with the horizon
  w <- apply(log(fc$rates[, , 1, ]), 1:2, stats::IQR)
  expect_gt(mean(w[, 5]), mean(w[, 1]))
})

test_that("the coherent CoDa model averages the rates of the bottom strata within groups", {
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
  expect_true(all(c("K", "B", "arima_K", "arima_k") %in% p$parameter))
  expect_error(mt_fit(coda(common_by = "nonsense"), xa), class = "mt_error_model_fit")
  # the common component is the CoDa model of the mean rates of the members
  x2 <- mt_example(deaths = "single", years = 1995:2012)
  f2 <- mt_fit(coda(coherent = TRUE), x2)
  D <- mt_array(x2, "deaths")
  E <- mt_array(x2, "exposure")
  lay <- f2$fits$layout
  sc <- lapply(1:2, function(s) mortensor:::coda_schedule(D[, , s], E[, , s], lay, "s"))
  core <- mortensor:::coda_core(mortensor:::coda_clr((sc[[1]]$m + sc[[2]]$m) / 2, lay), 1, c(0, 1, 1), TRUE)
  expect_equal(f2$fits$commons[[1]]$b, core$b)
  expect_equal(f2$fits$commons[[1]]$k, core$k)
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
  expect_true(all(c("a", "b", "k", "arima_k") %in% p$parameter))
  bt <- mt_backtest(coda_causes(), x, origins = 2010, h = 2, n_draws = 20)
  expect_true(all(is.finite(bt$crps)))
  # two components (Oeppen 2008, Section 7), each index with its own model
  f2 <- mt_fit(coda_causes(rank = 2, order = c(0, 2, 2)), x)
  g <- f2$fits$groups[[1]]
  expect_equal(dim(g$b), c(length(f2$fits$parts), 2))
  expect_length(g$ts, 2)
  expect_equal(g$ts[[1]]$order, c(0, 2, 2))
  expect_true(all(c("b1", "b2", "k1", "k2") %in% mortensor:::coda_parameters(f2)$parameter))
  expect_lt(gap(mt_forecast(f2, h = 3, n_draws = 10, seed = 1)$rates), 1e-12)
  # cause rates follow eq. (7) of Oeppen (2008): m[x, c] = m[x] d[x, c] / d[x]
  P <- closure_cols(exp(g$alpha + g$b %*% t(g$k)))
  R <- mortensor:::md_rates(log(P), f2$fits$layout, 3, mortensor:::coda_open_fitted(g$open, g$k))
  expect_equal(R[, , 1] / mortensor:::md_total(R), P[1:71, ] / (P[1:71, ] + P[72:142, ] + P[143:213, ]),
               tolerance = 1e-10, ignore_attr = TRUE)
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

test_that("coda() fits and forecasts a closed age range with either open-interval method", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  a <- single_ages(30, 90, open = FALSE)
  xc <- mt_regroup(mt_set_grid(x, age = a), age = a)
  for (open in c("regression", "auto")) {
    fc <- mt_forecast(mt_fit(coda(open = open), xc), h = 3, n_draws = 20, seed = 1)
    expect_true(all(is.finite(fc$rates)))
  }
})
