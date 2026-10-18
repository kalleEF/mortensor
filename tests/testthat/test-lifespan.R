# Life tables of France (HMD via demography::fr.mort): single ages 0-109 and
# 110+, or deaths and exposures summed to the age groups starting at `breaks`
fr_life_tables <- function() {
  fr <- mortensor:::pkg_data("fr.mort", "demography")
  grouped <- function(s, yr, breaks) {
    m <- fr$rate[[s]][, as.character(yr)]
    p <- fr$pop[[s]][, as.character(yr)]
    g <- findInterval(0:110, breaks)
    D <- tapply(ifelse(is.na(m), 0, m * p), g, sum)
    E <- tapply(p, g, sum)
    list(m = as.vector(D / E), age = age_groups(breaks))
  }
  list(female_2000_single = list(m = fr$rate$female[, "2000"], age = single_ages(0, 110)),
       male_1987_single = list(m = fr$rate$male[, "1987"], age = single_ages(0, 110)),
       female_1950_5y = grouped("female", 1950, seq(0, 100, 5)),
       female_2000_abridged = grouped("female", 2000, c(0, 1, seq(5, 85, 5))))
}

test_that("lifespan variation agrees with LifeIneq", {
  skip_if_not_installed("demography")
  # Reference values at ages 0 and 65 from LifeIneq 00.04.09 (R/indices.R,
  # sourced from github.com/alysonvanraalte/LifeIneq): ineq_edag(),
  # ineq_sd(), ineq_cov(), ineq_H(), ineq_gini(distribution_type = "aad" and
  # "rl") and ineq_iqr() on the columns of mortensor's life table (dx, lx
  # with radix 1, ex, ax with the open interval's ax = Lx / dx), with
  # set.seed(1) before the calls (ineq_iqr() perturbs tied lx at random).
  ref <- list(
    female_2000_single = rbind(
      edagger = c(9.663024709, 7.008585465),
      sd = c(14.18425557, 8.211834891),
      cv = c(0.1712647262, 0.09521290908),
      entropy = c(0.116674102, 0.3298611941),
      gini = c(0.08453618794, 0.0537459495),
      gini_rl = c(0.08453618794, 0.2181679582),
      iqr = c(13.41011018, 11.22595719)),
    male_1987_single = rbind(
      edagger = c(12.4459491, 7.202768889),
      sd = c(17.07315081, 7.945863166),
      cv = c(0.2369938974, 0.09934097606),
      entropy = c(0.1727633064, 0.480640979),
      gini = c(0.1237094046, 0.0568842088),
      gini_rl = c(0.1237094046, 0.3036167286),
      iqr = c(19.32465572, 12.05458686)),
    female_1950_5y = rbind(
      edagger = c(14.64681576, 6.98491648),
      sd = c(22.05606492, 7.485704412),
      cv = c(0.3193239563, 0.09399239359),
      entropy = c(0.2120541073, 0.4770597137),
      gini = c(0.1577954255, 0.05295861841),
      gini_rl = c(0.1577954255, 0.2880634363),
      iqr = c(19.20655623, 11.21954679)),
    female_2000_abridged = rbind(
      edagger = c(11.35334869, 8.84333507),
      sd = c(13.95217267, 8.481925776),
      cv = c(0.1681096598, 0.09811803231),
      entropy = c(0.1367964426, 0.412350776),
      gini = c(0.07637345279, 0.04400371641),
      gini_rl = c(0.07637345279, 0.1773722725),
      iqr = c(16.58540941, 14.49182808)))
  lts <- fr_life_tables()
  for (nm in names(ref)) {
    lt <- lts[[nm]]
    r <- ref[[nm]]
    for (k in 1:2) {
      at <- c(0, 65)[k]
      ours <- function(ms, ...) lifespan_variation(lt$m, ms, lt$age, at = at, ...)
      rv <- function(ms) unname(r[ms, k])
      info <- paste(nm, "at", at)
      expect_equal(ours("edagger"), rv("edagger"), tolerance = 1e-8, info = info)
      expect_equal(ours("entropy"), rv("entropy"), tolerance = 1e-8, info = info)
      expect_equal(ours("gini"), rv("gini"), tolerance = 1e-8, info = info)
      expect_equal(ours("gini", distribution = "remaining_life"), rv("gini_rl"), tolerance = 1e-8, info = info)
      # LifeIneq's ineq_var() offsets the ages at death by the first ages of
      # the table, which differs above age 0 for the abridged table
      if (nm != "female_2000_abridged" || at == 0) {
        expect_equal(ours("sd"), rv("sd"), tolerance = 1e-8, info = info)
        expect_equal(ours("var"), rv("sd")^2, tolerance = 1e-8, info = info)
        expect_equal(ours("cv"), rv("cv"), tolerance = 1e-8, info = info)
      }
      # quartiles: linear interpolation of the survival curve here, a monotone
      # spline in LifeIneq (which extrapolates beyond the open interval's
      # start, so the abridged table with 85+ is not compared)
      if (grepl("single", nm)) expect_lt(abs(ours("iqr") - rv("iqr")), 0.03)
      if (grepl("5y", nm)) expect_lt(abs(ours("iqr") - rv("iqr")), 0.4)
    }
  }
})

test_that("lifespan variation agrees with sampled ages at death", {
  skip_if_not_installed("demography")
  lt <- fr_life_tables()$female_2000_single
  r <- mortensor:::lt_columns(matrix(lt$m, ncol = 1), lt$age)
  na <- length(lt$m)
  # ages at death: interval by the life-table deaths, uniform within closed
  # intervals, exponential (mean Lx / dx) in the open interval
  n <- 2e5
  set.seed(2)
  j <- sample.int(na, n, replace = TRUE, prob = r$dx[, 1])
  closed <- j < na
  y <- lt$age$lower[j] + ifelse(closed, stats::runif(n), stats::rexp(n, r$dx[na, 1] / r$Lx[na, 1]))
  ys <- sort(y)
  gini <- 2 * sum(ys * seq_len(n)) / (n * sum(ys)) - (n + 1) / n
  q <- stats::quantile(y, c(0.25, 0.75), names = FALSE)
  expect_lt(abs(lifespan_variation(lt$m, "sd", lt$age) - stats::sd(y)), 0.1)
  expect_lt(abs(lifespan_variation(lt$m, "gini", lt$age) - gini), 0.001)
  expect_lt(abs(lifespan_variation(lt$m, "iqr", lt$age) - diff(q)), 0.1)
})

test_that("lifespan variation of a constant hazard", {
  # exponential lifespans with mean 20: e-dagger = e, H = 1, sd = e, Gini of
  # the remaining life 1/2, interquartile range log(3) / m
  m <- rep(0.05, 301)
  a <- single_ages(0, 300)
  lv <- function(ms) lifespan_variation(m, ms, a, distribution = "remaining_life")
  expect_equal(lv("edagger"), 20, tolerance = 1e-6)
  expect_equal(lv("entropy"), 1, tolerance = 1e-6)
  expect_lt(abs(lv("sd") - 20), 0.01)
  expect_lt(abs(lv("gini") - 0.5), 0.001)
  expect_lt(abs(lv("iqr") - log(3) / 0.05), 0.01)
})

test_that("lifespan measures are consistent with each other and the life table", {
  m <- 5e-5 * exp(0.095 * (0:100)) + 2e-4
  a <- single_ages(0, 100)
  lt <- life_table(m, a)
  e65 <- lt$ex[lt$age_lo == 65]
  lv <- function(ms, ...) lifespan_variation(m, ms, a, at = 65, ...)
  expect_equal(lv("entropy"), lv("edagger") / e65)
  expect_equal(lv("var"), lv("sd")^2)
  expect_equal(lv("cv"), lv("sd") / (65 + e65))
  expect_equal(lv("cv", distribution = "remaining_life"), lv("sd") / e65)
  expect_equal(lv("gini", distribution = "remaining_life"), lv("gini") * (65 + e65) / e65)
  # the Gini coefficient equals its double-sum definition, also when the
  # ages at death are not increasing (user ax)
  ax <- rep(0.5, 101)
  ax[1:2] <- c(0.9, 0.05)
  r <- mortensor:::lt_columns(matrix(m, ncol = 1), a, ax = ax)
  p <- r$dx[, 1]
  y <- a$lower + c(r$ax[-101, 1], r$Lx[101, 1] / r$dx[101, 1])
  mu <- sum(p * y)
  g <- sum(abs(outer(y, y, "-")) * outer(p, p)) / (2 * mu)
  expect_equal(lifespan_variation(m, "gini", a, ax = ax), g)
  # lifespans vary less at higher ages
  expect_gt(lifespan_variation(m, "sd", a, at = 30), lv("sd"))
})

test_that("lifespan_variation works on arrays, draws and data", {
  x <- mt_example(deaths = "single", years = 2015:2019)
  g <- mt_grid(x)
  m <- mt_array(x, "deaths") / mt_array(x, "exposure")
  v <- lifespan_variation(m, "edagger", g$age, at = 65)
  expect_equal(dim(v), c(5, 2))
  expect_equal(dimnames(v), dimnames(m)[-1])
  expect_equal(v[2, 1], lifespan_variation(m[, 2, 1], "edagger", g$age, at = 65))
  obs <- lifespan_variation(x, measure = c("edagger", "gini"), at = 65, probs = c(0.1, 0.9))
  expect_equal(nrow(obs), 2 * 5 * 2)
  expect_true(all(c("sex", "period", "mean", "q0.1", "q0.9", "measure", "at") %in% names(obs)))
  expect_true(all(obs$q0.1 <= obs$q0.9))
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 30, seed = 1)
  p <- lifespan_variation(fc, "sd", at = 30, paths = TRUE)
  expect_equal(dim(p), c(3, 2, 30))
  s <- lifespan_variation(fc, "sd", at = 30)
  expect_equal(s$mean, as.vector(apply(p, c(1, 2), mean)))
  expect_equal(s$at, rep(30, 6))
  # path by path: the same as the values of the paths' rates
  expect_equal(p[1, 2, 5], lifespan_variation(fc$rates[, 1, 2, 5], "sd", fc$age, at = 30))
})

test_that("lifespan_variation checks its arguments", {
  m <- 5e-5 * exp(0.095 * (0:100))
  a <- single_ages(0, 100)
  expect_error(lifespan_variation(m, "variance", a), class = "mt_error_argument")
  expect_error(lifespan_variation(m, c("sd", "gini"), a), class = "mt_error_argument")
  expect_error(lifespan_variation(m, "sd", a, at = 65.5), class = "mt_error_lifetable")
  expect_error(lifespan_variation(m, "sd", single_ages(0, 100, open = FALSE)), class = "mt_error_lifetable")
  expect_error(lifespan_variation(m[-1], "sd", a), class = "mt_error_lifetable")
  x <- mt_example(deaths = "single", years = 2015:2019)
  fc <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 10, seed = 1)
  expect_error(lifespan_variation(fc, c("sd", "gini"), paths = TRUE), class = "mt_error_argument")
  # a missing rate gives missing values, not an error
  mm <- m
  mm[50] <- NA
  expect_true(is.na(lifespan_variation(mm, "gini", a)))
})

test_that("lifespan variation is a quantity of comparisons, contrasts and inequality indices", {
  x <- mt_example(deaths = "single", years = 2005:2019)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 3, n_draws = 50, seed = 1)
  lv <- lifespan_variation(fc, measure = "gini", at = 30)
  tab <- mt_compare_forecasts(list(lc = fc), quantity = "lifespan_variation", measure = "gini", at = 30,
                              wide = FALSE, level = 0.95)
  expect_equal(tab$med, lv$q0.5, tolerance = 1e-12)
  expect_equal(tab$lo, lv$q0.025, tolerance = 1e-12)
  paths <- lifespan_variation(fc, measure = "edagger", paths = TRUE)
  ct <- mt_contrast(fc, var = "sex", ref = "female", quantity = "lifespan_variation")
  d <- paths[, "male", ] - paths[, "female", ]
  expect_equal(ct$mean, unname(rowMeans(d)), tolerance = 1e-12)
  expect_error(mt_compare_forecasts(list(fc), quantity = "lifespan_variation", measure = c("sd", "gini")),
               class = "mt_error_argument")
  expect_error(mt_contrast(fc, var = "sex", ref = "female", quantity = "lifespan_variation", measure = "x"),
               class = "mt_error_argument")
  skip_if_not_installed("ggplot2")
  p <- mt_plot_compare(list(lc = fc), quantity = "lifespan_variation", measure = "sd", data = x)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$y, "Standard deviation")
})
