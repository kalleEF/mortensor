test_that("life table identities hold", {
  m <- 3e-4 + 2e-5 * exp(0.1 * (30:100))
  lt <- life_table(m, single_ages(30, 100))
  expect_equal(lt$lx[-1], lt$lx[-nrow(lt)] - lt$dx[-nrow(lt)])
  expect_equal(sum(lt$dx), lt$lx[1])
  expect_equal(lt$Tx[1], sum(lt$Lx))
  expect_equal(lt$ex, lt$Tx / lt$lx)
  expect_equal(tail(lt$qx, 1), 1)
})

test_that("D10: life expectancy uses interval widths (single and abridged ages agree)", {
  ages <- 0:109
  m <- 2e-5 * exp(0.095 * ages) + 2e-4
  e_single <- life_table(c(m[1:100], m[101]), single_ages(0, 100))$ex[1]
  grp <- c(0, 1, seq(5, 100, 5))
  ax <- age_groups(grp)
  # abridged rates from the single-year schedule with stationary weights
  lt <- life_table(c(m[1:100], m[101]), single_ages(0, 100))
  gid <- findInterval(0:100, grp)
  mg <- tapply(lt$dx, gid, sum) / tapply(lt$Lx, gid, sum)
  e_abr <- life_table(as.numeric(mg), ax)$ex[1]
  expect_lt(abs(e_single - e_abr), 0.1)
})

test_that("Gompertz closure of the open interval matches the exact integral", {
  b <- 0.1
  a <- 1e-5
  m <- a * exp(b * (30:99 + 0.5))
  open_rate <- 0.5
  lt_g <- life_table(c(m, open_rate), single_ages(30, 100), open = "gompertz")
  # exact remaining life expectancy at 100 under mu(x) = a exp(b x)
  surv <- function(t) exp(-a / b * (exp(b * (100 + t)) - exp(b * 100)))
  e100 <- stats::integrate(surv, 0, 40)$value
  expect_lt(abs(tail(lt_g$ex, 1) - e100), 0.05)
  expect_error(life_table(1:3 / 10, single_ages(0, 5)), class = "mt_error_lifetable")
})

test_that("life_expectancy works on arrays and checks `at`", {
  x <- mt_example(years = 2001:2003, deaths = "single")
  m <- mt_array(x, "deaths") / mt_array(x, "exposure")
  e <- life_expectancy(m, mt_grid(x)$age, at = 65)
  expect_equal(dim(e), c(3, 2))
  expect_true(all(e > 15 & e < 30))
  expect_error(life_expectancy(m, mt_grid(x)$age, at = 65.5), class = "mt_error_lifetable")
})

test_that("expected deaths use future exposures", {
  x <- mt_example(years = 2001:2012, deaths = "single")
  fc <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 10, seed = 1)
  mu <- mt_deaths(fc, noise = FALSE)
  expect_equal(mu[, 1, , 1], fc$rates[, 1, , 1] * fc$exposure_last[, 1, ])
  E2 <- fc$exposure_last[, 1, ] * 2
  expect_equal(mt_deaths(fc, exposure = E2, noise = FALSE), 2 * mu)
  expect_error(suppressWarnings(mt_deaths(fc, exposure = matrix(1, 3, 2))), class = "mt_error_argument")
})

test_that("HMD tables are converted", {
  d <- data.frame(Year = 2000, Age = 0:2, Female = c(10, 1, 1), Male = c(12, 2, 1), Total = c(22, 3, 2),
                  OpenInterval = c(FALSE, FALSE, TRUE))
  e <- transform(d, Female = 1000, Male = 1000, Total = 2000)
  x <- mt_from_hmd(d, e)
  expect_equal(mt_strata(x)$sex, c("female", "male"))
  expect_equal(tail(mt_grid(x)$age$labels, 1), "2+")
  expect_error(mt_read_hmd("SWE", username = "", password = ""), class = "mt_error_adapter")
})

test_that("probability of dying, age-standardised rates and YLL are computed per path", {
  x <- mt_example(deaths = "single", years = 2001:2015)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 40, seed = 1)
  p <- prob_dying(fc, from = 30, to = 70)
  expect_equal(nrow(p), 6)
  expect_true(all(p$q0.5 > 0 & p$q0.5 < 1))
  # equals 1 - l70 / l30 of the median path's life table within sampling error
  m <- apply(fc$rates[, 1, 1, ], 1, median)
  lt <- life_table(m, fc$age)
  expect_equal(p$q0.5[1], 1 - lt$lx[41] / lt$lx[1], tolerance = 0.05)
  a <- asr(fc)
  expect_true(all(a$mean > 0))
  esp <- standard_population()
  expect_equal(sum(esp$weight), 1e5)
  w <- mortensor:::standard_weights(fc$age, esp)
  expect_equal(sum(w), sum(esp$weight[esp$age_lo >= 30]))
  y <- yll(fc)
  expect_true(all(y$mean > 0))
  expect_error(prob_dying(fc, from = 70, to = 30), class = "mt_error_argument")
})

test_that("forecasts are exported as summaries and draws", {
  x <- mt_example(deaths = "single", years = 2005:2015)
  fc <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 5, seed = 1)
  f1 <- mt_export(fc, tempdir())
  f2 <- mt_export(fc, tempdir(), what = "draws")
  on.exit(unlink(c(f1, f2)))
  s <- data.table::fread(f1)
  d <- data.table::fread(f2)
  expect_equal(nrow(s), 71 * 2 * 2)
  expect_equal(nrow(d), 71 * 2 * 2 * 5)
  expect_match(basename(f1), "origin2015_summary.csv")
  expect_error(mt_export(fc, file.path(tempdir(), "missing-dir")), class = "mt_error_argument")
})

test_that("standard weights cover open grid intervals and paths can be returned", {
  esp <- standard_population("esp2013")
  w <- standard_weights(age_groups(c(seq(50, 85, 5), 90)), esp)
  expect_equal(w[length(w)], sum(esp$weight[esp$age_lo >= 90]))
  expect_equal(sum(w), sum(esp$weight[esp$age_lo >= 50]))
  w1 <- standard_weights(single_ages(0, 100), esp)
  expect_equal(sum(w1), sum(esp$weight))
  x <- mt_example(deaths = "single", years = 2005:2019)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 20, seed = 1)
  a <- asr(fc, paths = TRUE)
  expect_identical(dim(a), c(3L, 2L, 20L))
  expect_equal(apply(a, 1:2, stats::median)[, 1], asr(fc)[stratum == 1]$q0.5, ignore_attr = TRUE)
  expect_identical(dim(prob_dying(fc, 30, 70, paths = TRUE)), c(3L, 2L, 20L))
})
