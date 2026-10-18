# Rotated Lee-Carter (Li, Lee and Gerland 2013), lee_carter(rotate = TRUE)

# Lee-Carter-type data with ages 0-100+ and a known age pattern of decline
# that is steep at young ages (as in low-mortality populations): women with
# e0 from about 80 to 86, men from about 72 to 77 over 1981-2010. Large
# exposures, so that the fit is close to the truth.
rot_data <- function(years = 1981:2010, seed = 1) {
  set.seed(seed)
  ages <- 0:100
  braw <- ifelse(ages < 15, 2 - ages / 15, ifelse(ages < 65, 0.8, 1.2 - 0.8 * pmax(ages - 80, 0) / 20))
  b <- braw / sum(braw)
  rows <- list()
  for (sx in c("female", "male")) {
    lvl <- if (sx == "female") 0.55 else 1.6
    a <- log(5e-3 * (ages == 0) + 2e-4 * exp(-0.3 * ages) + 1.5e-4 + lvl * 2.5e-5 * exp(0.1 * ages))
    k <- -2 * (seq_along(years) - mean(seq_along(years))) + stats::rnorm(length(years), 0, 0.5)
    m <- exp(a + outer(b, k))
    E <- matrix(1e6, length(ages), length(years))
    D <- matrix(stats::rpois(length(m), m * E), length(ages))
    rows[[sx]] <- data.table::data.table(sex = sx, year = rep(years, each = length(ages)), age_lo = ages,
                                         age_hi = ifelse(ages == 100, Inf, ages + 1), deaths = as.vector(D),
                                         exposure = as.vector(E))
  }
  long <- data.table::rbindlist(rows)
  mt_data(deaths = long[, c("sex", "year", "age_lo", "age_hi", "deaths")],
          exposure = long[, c("sex", "year", "age_lo", "age_hi", "exposure")], by = "sex",
          age = c("age_lo", "age_hi"), omega = 110)
}

test_that("rotation weights and the ultimate pattern agree with MortCast", {
  # Expected values computed with MortCast 2.8-0 (GPL >= 2; Sevcikova, Li and
  # Gerland), by sourcing its R/LC.R: `b` below with names 0, 1, 5, ..., 100
  # (5-year age groups), ultimate.bx(b), and rotate.leecarter(b, ultimate.bx(b),
  # e0) for e0 = 75, 80, 85, 90, 95, 101.9, 102 and 110, of which the columns
  # for 75, 85, 101.9 and 110 are kept; printed with dput(control = "digits17").
  b_raw <- c(0.09, 0.07, 0.06, 0.055, 0.04, 0.035, 0.035, 0.034, 0.033, 0.035, 0.038, 0.042, 0.046, 0.05, 0.052,
             0.053, 0.05, 0.045, 0.038, 0.03, 0.022, 0.015)
  b <- b_raw / sum(b_raw)
  age <- age_groups(c(0, 1, seq(5, 100, 5)))
  bu_mortcast <- c(rep(0.05033881897386254, 15), 0.051306873184898349, 0.048402710551790899, 0.043562439496611809,
                   0.036786060019361085, 0.029041626331074537, 0.021297192642787996, 0.014520813165537268)
  B_mortcast <- matrix(c(
    0.092975206611570244, 0.072314049586776868, 0.061983471074380167, 0.056818181818181823, 0.041322314049586778,
    0.036157024793388434, 0.036157024793388434, 0.035123966942148761, 0.034090909090909095, 0.036157024793388434,
    0.039256198347107439, 0.043388429752066117, 0.047520661157024795, 0.051652892561983473, 0.053719008264462811,
    0.054752066115702477, 0.051652892561983473, 0.046487603305785122, 0.039256198347107439, 0.030991735537190084,
    0.022727272727272728, 0.015495867768595042,
    0.078075316384689591, 0.064634493649133762, 0.057914082281355861, 0.054553876597466903, 0.044473259545800038,
    0.041113053861911081, 0.041113053861911081, 0.040441012725133292, 0.039768971588355503, 0.041113053861911081,
    0.043129177272244454, 0.045817341819355623, 0.048505506366466784, 0.051193670913577946, 0.05253775318713353,
    0.053548094594578405, 0.050517070372243775, 0.045465363335019396, 0.038392973482905268, 0.030310242223346266,
    0.022227510963787264, 0.015155121111673133,
    0.050339905757156404, 0.050339379113096473, 0.050339115791066501, 0.050338984130051521, 0.05033858914700657,
    0.050338457485991583, 0.050338457485991583, 0.050338431153788589, 0.050338404821585594, 0.050338457485991583,
    0.050338536482600574, 0.050338641811412566, 0.05033874714022455, 0.050338852469036535, 0.050338905133442531,
    0.051306961001393335, 0.048402793397540889, 0.0435625140577868, 0.036786122982131075, 0.02904167603852453,
    0.021297229094917991, 0.014520838019262265,
    bu_mortcast), 22)
  bu <- mortensor:::lc_ultimate_b(b, age)
  expect_equal(bu, bu_mortcast, tolerance = 1e-12)
  expect_equal(sum(bu), 1)
  w <- mortensor:::rotation_weights(c(75, 85, 101.9, 110))
  expect_equal(outer(b, 1 - w) + outer(bu, w), B_mortcast, tolerance = 1e-12)
  # weights: MortCast's (wt) at 85, 90, 95; 0 below 80, 1 from 102
  expect_equal(mortensor:::rotation_weights(c(70, 80, 85, 90, 95, 102, 120)),
               c(0, 0, 0.34946417959909837, 0.65486073394528510, 0.87767898956725554, 1, 1), tolerance = 1e-14)
  expect_equal(mortensor:::rotation_weights(90, 75, 100, 1), 0.5 * (1 + sin(pi / 2 * (2 * 15 / 25 - 1))))
})

test_that("the ultimate pattern generalises to single years and keeps groups without deaths flat", {
  age <- single_ages(0, 100)
  b <- c(seq(0.02, 0.012, length.out = 15), seq(0.008, 0.01, length.out = 50), seq(0.012, 0.004, length.out = 36))
  b <- b / sum(b)
  bu <- mortensor:::lc_ultimate_b(b, age)
  lo <- age$lower
  m <- mean(bu[lo >= 15 & lo < 65])
  expect_equal(bu[lo < 65], rep(m, 65))
  expect_equal(bu[lo == 65], m)
  expect_equal(bu[lo >= 65] / b[lo >= 65], rep(bu[lo == 65] / b[lo == 65], 36))
  expect_equal(sum(bu), 1)
  # an age group without deaths (b = 0) keeps no trend
  fixed <- lo == 3
  b0 <- b
  b0[fixed] <- 0
  bu0 <- mortensor:::lc_ultimate_b(b0, age, fixed)
  expect_equal(bu0[fixed], 0)
  expect_equal(sum(bu0), sum(b0))
  # no positive pattern at 65: no ultimate pattern
  b1 <- b
  b1[lo == 65] <- -0.001
  expect_null(mortensor:::lc_ultimate_b(b1, age))
})

test_that("rotation settings are checked and shown in the label", {
  expect_match(lee_carter(rotate = TRUE)$label, "rotated")
  expect_identical(lee_carter()$label, "Lee-Carter")
  expect_false(lee_carter()$params$rotate)
  s <- lee_carter(rotate = TRUE, rotation = c(e0_end = 100))
  expect_equal(s$params$rotation, c(e0_start = 80, e0_end = 100, power = 0.5))
  expect_match(s$label, "e0 80-100")
  expect_equal(lee_carter(rotation = list(power = 1))$params$rotation[["power"]], 1)
  expect_error(lee_carter(rotate = "yes"), class = "mt_error_argument")
  expect_error(lee_carter(rotate = TRUE, rotation = c(80, 102, 0.5)), class = "mt_error_argument")
  expect_error(lee_carter(rotate = TRUE, rotation = c(start = 80)), class = "mt_error_argument")
  expect_error(lee_carter(rotate = TRUE, rotation = c(e0_start = 90, e0_end = 85)), class = "mt_error_argument")
  expect_error(lee_carter(rotate = TRUE, rotation = c(power = 0)), class = "mt_error_argument")
})

test_that("rotation needs ages from 0 with age groups starting at 15 and 65", {
  expect_error(mortensor:::rotation_check(single_ages(30, 100)), class = "mt_error_model_data")
  expect_error(mortensor:::rotation_check(age_groups(seq(0, 100, 10))), class = "mt_error_model_data")
  expect_silent(mortensor:::rotation_check(age_groups(c(0, 1, seq(5, 100, 5)))))
  expect_warning(mortensor:::rotation_check(single_ages(0, 100, open = FALSE)), class = "mt_warning_model_data")
  x <- mt_example(deaths = "single", years = 2001:2012)
  expect_error(suppressWarnings(mt_fit(lee_carter(rotate = TRUE), x)), class = "mt_error_model_fit",
               regexp = "ages from 0")
})

test_that("rotated draws keep the life expectancy of the unrotated draws path by path", {
  x <- rot_data()
  f0 <- mt_fit(lee_carter(), x)
  f1 <- mt_fit(lee_carter(rotate = TRUE), x)
  expect_equal(fitted(f1), fitted(f0))
  expect_length(f1$fits[[1]]$b_ultimate, 101)
  fc0 <- mt_forecast(f0, h = 30, n_draws = 20, seed = 3)
  fc1 <- mt_forecast(f1, h = 30, n_draws = 20, seed = 3)
  e0 <- life_expectancy(fc0, paths = TRUE)
  e1 <- life_expectancy(fc1, paths = TRUE)
  expect_lt(max(abs(e1 - e0)), 1e-6)
  expect_true(all(e0[, "female", ] > 80))
  # women: decline slower at ages 0-14 and faster at 15-64 than without rotation
  lo <- fc0$age$lower
  dec <- function(fc, rows) mean(log(fc$rates[rows, 1, "female", ] / fc$rates[rows, 30, "female", ]))
  expect_lt(dec(fc1, lo < 15), dec(fc0, lo < 15))
  expect_gt(dec(fc1, lo >= 15 & lo < 65), dec(fc0, lo >= 15 & lo < 65))
  expect_true(all(is.finite(fc1$rates) & fc1$rates > 0))
})

test_that("rotation leaves paths below e0_start unchanged, per stratum", {
  x <- rot_data()
  f0 <- mt_fit(lee_carter(), x)
  f1 <- mt_fit(lee_carter(rotate = TRUE), x)
  fc0 <- mt_forecast(f0, h = 3, n_draws = 10, seed = 5)
  fc1 <- mt_forecast(f1, h = 3, n_draws = 10, seed = 5)
  e0 <- life_expectancy(fc0, paths = TRUE)
  # men stay below 80 and are not rotated; women are above 80 and are
  expect_true(all(e0[, "male", ] < 80))
  expect_identical(fc1$rates[, , "male", ], fc0$rates[, , "male", ])
  expect_false(isTRUE(all.equal(fc1$rates[, , "female", ], fc0$rates[, , "female", ])))
  # a rotation that starts above every path changes nothing
  f2 <- mt_fit(lee_carter(rotate = TRUE, rotation = c(e0_start = 150, e0_end = 160)), x)
  expect_identical(mt_forecast(f2, h = 3, n_draws = 10, seed = 5)$rates, fc0$rates)
})

test_that("from e0_end on, the age pattern of decline is the ultimate pattern", {
  x <- rot_data()
  for (jump_off in c("fitted", "observed")) {
    fit <- mt_fit(lee_carter(rotate = TRUE, rotation = c(e0_start = 50, e0_end = 60), jump_off = jump_off), x)
    fc <- mt_forecast(fit, h = 5, n_draws = 5, seed = 1)
    for (s in 1:2) {
      f <- fit$fits[[s]]
      bu <- f$b_ultimate
      J <- if (jump_off == "observed") f$last_log_rate else log(f$fitted[, ncol(f$fitted)])
      for (j in 1:5) {
        # log rates move from the jump-off along bu only
        d <- unname(log(fc$rates[, , s, j]) - J)
        kap <- colSums(d) / sum(bu)
        expect_equal(unname(d), outer(bu, kap), tolerance = 1e-8)
      }
    }
  }
})

test_that("rotation works with the observed jump-off, the bootstrap, SVD and a chosen period", {
  x <- rot_data()
  check_e0 <- function(spec, h = 10, n = 12) {
    unrot <- spec
    unrot$params$rotate <- FALSE
    # the same bootstrap resamples in both fits
    f0 <- with_seed_test(1, mt_fit(unrot, x))
    f1 <- with_seed_test(1, mt_fit(spec, x))
    e0 <- life_expectancy(mt_forecast(f0, h = h, n_draws = n, seed = 2), paths = TRUE)
    e1 <- life_expectancy(mt_forecast(f1, h = h, n_draws = n, seed = 2), paths = TRUE)
    expect_lt(max(abs(e1 - e0)), 1e-6)
    f1
  }
  check_e0(lee_carter(rotate = TRUE, jump_off = "observed"))
  check_e0(lee_carter(rotate = TRUE, method = "svd", adjust = "e0"))
  check_e0(lee_carter(rotate = TRUE, choose_period = TRUE, min_years = 15))
  fb <- check_e0(lee_carter(rotate = TRUE, bootstrap = TRUE, n_boot = 3), n = 6)
  expect_true(all(vapply(fb$fits[[1]]$boot, function(s) length(s$b_ultimate) == 101, NA)))
  expect_mt_model(lee_carter(rotate = TRUE), data = x, h = 3, n_draws = 6)
})

test_that("the index solver matches target life expectancies", {
  age <- single_ages(0, 100)
  lo <- age$lower
  J <- log(5e-3 * (lo == 0) + 1e-4 + 2e-5 * exp(0.1 * lo))
  B <- matrix(c(rep(1 / 101, 101), seq(2, 1, length.out = 101) / sum(seq(2, 1, length.out = 101))), 101)
  target <- c(85, 95)
  kap <- mortensor:::lc_solve_e0(J, B, target, c(0, 0), age)
  e0 <- life_expectancy(exp(J + B * rep(kap, each = 101)), age)
  expect_lt(max(abs(e0 - target)), 1e-9)
})
