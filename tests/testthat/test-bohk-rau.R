short_mcmc <- c(chains = 2, iterations = 1700, burn_in = 200, thin = 3)

test_that("bohk_rau() fits and forecasts with both core models", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  set.seed(1)
  res <- expect_mt_model(bohk_rau(mcmc = short_mcmc), x, h = 4, n_draws = 30)
  f <- res$fit$fits[[1]]
  expect_identical(length(f$samples$sigma), 1000L)
  expect_identical(dim(f$romi), c(71L, 17L))
  # fitted rates are the smoothed rates; the forecast starts from the last of them
  expect_equal(stats::fitted(res$fit)[, "2012", 1], exp(f$jump_off), ignore_attr = TRUE)
  expect_lt(stats::median(abs(log(stats::fitted(res$fit) / res$fit$observed))), 0.1)
  set.seed(2)
  resl <- expect_mt_model(bohk_rau(core = "linear", draws = "samples", mcmc = short_mcmc), x, h = 4, n_draws = 30)
  expect_true(all(c("mu1", "omega2", "r") %in% names(resl$fit$fits[[1]]$samples)))
  p <- mt_parameters(resl$fit)
  expect_true(all(c("b1", "b2", "sigma", "mu1", "rhat") %in% p$parameter))
  expect_true(all(p[parameter == "b1", lower <= value & value <= upper]))
  expect_match(summary(res$fit)$assumptions, "Bohk-Ewald and Rau")
})

test_that("bohk_rau() follows equations 1 and 20 and applies the limits", {
  x <- mt_example(deaths = "single", years = 1995:2012)
  set.seed(3)
  fit <- mt_fit(bohk_rau(smooth = FALSE, limits = c(0.0099, 0.0101), mcmc = short_mcmc), x)
  D <- mt_array(x, "deaths")[, , 1]
  E <- mt_array(x, "exposure")[, , 1]
  m <- pmax(D, 0.5) / E
  f <- fit$fits[[1]]
  expect_equal(f$romi, 1 - m[, -1] / m[, -18], ignore_attr = TRUE)
  expect_equal(f$jump_off, log(m[, 18]), ignore_attr = TRUE)
  fc <- mt_forecast(fit, h = 5, n_draws = 20, seed = 1)
  r <- fc$rates[, , 1, ] / m[, 18]
  hh <- rep(1:5, each = 71)
  expect_true(all(r <= (1 - 0.0099)^hh + 1e-12 & r >= (1 - 0.0101)^hh - 1e-12))
})

test_that("bohk_rau() quantile draws are comonotone and match the sample paths in the median", {
  x <- mt_subset(mt_example(deaths = "single", years = 1990:2012), sex == "female")
  set.seed(4)
  fq <- mt_fit(bohk_rau(limits = NULL, mcmc = short_mcmc), x)
  fs <- fq
  fs$fits[[1]]$draws <- "samples"
  dq <- mt_forecast(fq, h = 3, n_draws = 400, seed = 2)$rates[, , 1, ]
  ds <- mt_forecast(fs, h = 3, n_draws = 400, seed = 2)$rates[, , 1, ]
  # the same draw is at the same rank at every age and horizon
  rk <- apply(dq[, 3, ], 1, rank)
  expect_true(all(rk == rk[, 1]))
  expect_true(all(apply(dq[, 1, ], 1, rank) == rk[, 1]))
  # medians agree; quantile draws are wider (comonotone over the years)
  mq <- apply(dq, 1:2, stats::median)
  ms <- apply(ds, 1:2, stats::median)
  expect_lt(max(abs(log(mq / ms))), 0.01)
  wq <- apply(log(dq[, 3, ]), 1, stats::IQR)
  ws <- apply(log(ds[, 3, ]), 1, stats::IQR)
  expect_gt(stats::median(wq / ws), 1)
})

test_that("the samplers of bohk_rau() recover the parameters of simulated improvements", {
  set.seed(3)
  na <- 12
  tt <- 1:30
  b1 <- log(0.02) + stats::rnorm(na, 0, 0.2)
  b2 <- -0.15 + stats::rnorm(na, 0, 0.03)
  R <- exp(b1 + outer(b2, log(tt))) + matrix(stats::rnorm(na * 30, 0, 0.001), na)
  mc <- c(chains = 2, iterations = 2200, burn_in = 200, thin = 4)
  s <- mortensor:::br_mcmc_exp(R, tt, mc, 0.001)
  expect_lt(abs(stats::median(s$sigma) / 0.001 - 1), 0.1)
  fit <- list(core = "exponential", samples = s)
  curve <- apply(mortensor:::br_curve(fit, 31), 1, stats::median)
  expect_lt(max(abs(curve / exp(b1 + b2 * log(31)) - 1)), 0.05)
  expect_lt(max(s$rhat), 1.1)
  a1 <- 0.015 + stats::rnorm(na, 0, 0.003)
  a2 <- 2e-4 + stats::rnorm(na, 0, 5e-5)
  R2 <- a1 + outer(a2, tt) + matrix(stats::rnorm(na * 30, 0, 0.002), na)
  s2 <- mortensor:::br_mcmc_lin(R2, tt, mc)
  expect_lt(abs(stats::median(s2$sigma) / 0.002 - 1), 0.1)
  expect_lt(abs(stats::median(s2$mu1) - mean(a1)), 0.003)
  expect_lt(abs(stats::median(s2$mu2) - mean(a2)), 1e-4)
  expect_lt(max(s2$rhat), 1.1)
})

test_that("the P-spline surface of bohk_rau() equals the PCLM engine at the same smoothing parameters", {
  x <- mt_example(deaths = "single", years = 2000:2012)
  x <- mt_regroup(mt_set_grid(x, age = single_ages(60, 89, open = FALSE)), age = single_ages(60, 89, open = FALSE))
  g <- mt_grid(x)
  D <- mt_array(x, "deaths")[, , 1]
  E <- mt_array(x, "exposure")[, , 1]
  m <- mortensor:::br_smooth(D, E, g$age, g$period$lower)
  xa <- mortensor:::cell_ages(g$age$lower, g$age$upper)
  xp <- g$period$lower + 0.5
  ba <- mortensor:::pclm_basis(xa, min(xa), max(xa), kr = 5)
  bp <- mortensor:::pclm_basis(xp, min(xp), max(xp), kr = 5)
  pen <- list(kronecker(diag(ncol(bp$B)), ba$P), kronecker(bp$P, diag(ncol(ba$B))))
  ref <- mortensor:::pclm_fit(as.vector(D), Matrix::Diagonal(length(D)), as.vector(E), kronecker(bp$B, ba$B), pen,
                              lambda = attr(m, "lambda"))
  expect_lt(max(abs(ref$gamma / as.vector(m) - 1)), 1e-5)
})

test_that("bohk_rau() smooths age 0 over time only", {
  skip_if_not_installed("StMoMo")
  x <- mt_subset(as_mt_data(mortensor:::pkg_data("EWMaleData", "StMoMo")), periods = c(1990, 2011))
  g <- mt_grid(x)
  D <- mt_array(x, "deaths")[, , 1]
  E <- mt_array(x, "exposure")[, , 1]
  m <- mortensor:::br_smooth_rates(D, E, g$age, g$period$lower)
  surface <- mortensor:::br_smooth(D, E, g$age, g$period$lower)
  err <- function(r, a) stats::median(abs(log(r[a, ] / (D[a, ] / E[a, ]))))
  expect_lt(err(m, 1), 0.03)
  # one surface over ages 0-100 cannot follow the fall of the rates after
  # age 0 and misses the rates at ages 1-4
  expect_gt(err(surface, 2:5), 2 * err(m, 2:5))
  expect_equal(m[-1, ], mortensor:::br_smooth(D[-1, ], E[-1, ], list(lower = g$age$lower[-1],
                                                                       upper = g$age$upper[-1]), g$period$lower),
               ignore_attr = TRUE)
})

test_that("bohk_rau() checks its arguments", {
  expect_error(bohk_rau(limits = c(0.03, 0.01)), class = "mt_error_argument")
  expect_error(bohk_rau(limits = 0.01), class = "mt_error_argument")
  expect_error(bohk_rau(mcmc = c(chain = 2)), class = "mt_error_argument")
  expect_error(bohk_rau(mcmc = c(iterations = 100, burn_in = 100)), class = "mt_error_argument")
  expect_error(bohk_rau(smooth = NA), class = "mt_error_argument")
  expect_identical(bohk_rau(mcmc = c(chains = 3))$params$mcmc,
                   c(chains = 3, iterations = 5200, burn_in = 200, thin = 5))
  expect_identical(bohk_rau()$params$limits, c(0.005, 0.035))
  expect_match(bohk_rau(core = "linear")$label, "linear")
})
