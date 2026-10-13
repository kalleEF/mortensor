# Behaviour with defective or unusual input data (see dev/stress-tests.R)

test_that("age and period labels of many sources are parsed", {
  p <- parse_age_labels(c("Under 1", "1 to 4", "5 - 9", "85 years and over", "90 & over", "0-4 years", "Age 85+",
                          "1 bis 4 Jahre", "15 bis unter 20 Jahre", "1 and under 5", "95 Jahre und mehr"))
  expect_equal(p$age_lo, c(0, 1, 5, 85, 90, 0, 85, 1, 15, 1, 95))
  expect_equal(p$age_hi, c(1, 5, 10, Inf, Inf, 5, Inf, 5, 20, 5, Inf))
  q <- parse_period_labels(c("2001", "2001-2005", "2001/02", "2001-02"))
  expect_equal(q$period_lo, rep(2001, 4))
  expect_equal(q$period_hi, c(2002, 2006, 2002, 2002))
  expect_error(parse_period_labels("spring"), class = "mt_error_axis")
})

test_that("totals, unknown ages, separators and label mismatches are handled with warnings", {
  d <- data.frame(year = 2020, sex = "f", age = c("0-49", "50+", "All ages", "Unknown"),
                  deaths = c("1,000", "9,000", "10,100", "100"), exposure = c(1e6, 2e5, 1.2e6, NA))
  w <- character()
  x <- withCallingHandlers(mt_data(deaths = d, exposure = d[1:2, ], by = "sex"),
                           warning = function(cnd) {
                             w <<- c(w, class(cnd)[1])
                             invokeRestart("muffleWarning")
                           })
  expect_true(all(c("mt_warning_totals", "mt_warning_unknown_age") %in% w))
  expect_equal(sum(mt_obs(x)$value), 10100)
  expect_equal(mt_obs(x)$value, c(1010, 9090))
  e <- data.frame(year = 2020, sex = "F", age = c("0-49", "50+"), exposure = "1.000,5")
  expect_warning(suppressWarnings(mt_data(deaths = d[1:2, ], exposure = e, by = "sex", decimal = ","),
                                  classes = "mt_warning_markers"), class = "mt_warning_strata")
  y <- suppressWarnings(mt_data(deaths = d[1:2, ], exposure = e, by = "sex", decimal = ","))
  expect_equal(mt_obs(y, "exposure")$value, c(1000.5, 1000.5))
  expect_warning(suppressWarnings(mt_data(deaths = d[1:2, ], exposure = e, by = "sex"),
                                  classes = c("mt_warning_strata", "mt_warning_missing")),
                 class = "mt_warning_markers")
})

test_that("errors name the problem: duplicates, overlaps and missing cells", {
  d <- data.frame(year = 2020, region = c("A", "B"), age = "90+", deaths = 5, exposure = 100)
  expect_error(mt_data(deaths = d, exposure = d), "columns not used: region", class = "mt_error_obs")
  o <- data.frame(year = 2020, age = c("85+", "85-89", "90+"), deaths = 1:3, exposure = 100)
  expect_error(mt_data(deaths = o, exposure = o), "Overlapping age groups", class = "mt_error_obs")
  g <- expand.grid(year = 2001:2008, age = c("80-84", "85+"), stringsAsFactors = FALSE)
  g$deaths <- 10
  g$exposure <- 100
  x <- mt_data(deaths = g[g$year != 2005, ], exposure = g)
  expect_error(mt_fit(rwd(), x), "No deaths in 2005", class = "mt_error_model_data")
  z <- g
  z$exposure[1] <- 0
  expect_error(mt_fit(rwd(), mt_data(deaths = g, exposure = z)), "zero or negative", class = "mt_error_model_data")
  h <- g
  h$deaths[1] <- 1000
  expect_warning(mt_fit(rwd(), mt_data(deaths = h, exposure = g)), class = "mt_warning_data_quality")
})

test_that("a level that is the total of the others is detected", {
  d <- expand.grid(year = 2001:2006, age = c("60-64", "65+"), sex = c("f", "m"), stringsAsFactors = FALSE)
  d$deaths <- 10 + seq_len(nrow(d))
  t <- aggregate(deaths ~ year + age, d, sum)
  t$sex <- "persons"
  d <- rbind(d, t[, names(d)])
  d$exposure <- 1000
  x <- mt_data(deaths = d, exposure = d, by = "sex")
  v <- mt_validate(x)
  expect_true("total among levels" %in% v$check)
  expect_false("total among levels" %in% mt_validate(x, hierarchy = mt_hierarchy(sex = "persons"))$check)
})

test_that("groups of ages that change over time are regrouped to common groups", {
  d <- data.frame(year = rep(2001:2002, c(2, 3)), age = c("80-84", "85+", "80-84", "85-89", "90+"),
                  deaths = c(10, 20, 10, 12, 9))
  p <- data.frame(year = rep(2001:2002, each = 2), age = c("80-84", "85+"), exposure = 100)
  x <- mt_data(deaths = d, exposure = p)
  expect_error(mt_fit(rwd(), x), "common", class = "mt_error_model_data")
  y <- mt_regroup(x, age = "common")
  expect_identical(mt_grid(y)$age$labels, c("80-84", "85+"))
  expect_equal(mt_obs(y)$value, c(10, 20, 10, 21))
})

test_that("Lee-Carter copes with age groups without deaths and reports failures once", {
  set.seed(2)
  E <- matrix(200, 6, 15)
  D <- matrix(stats::rpois(90, E * c(0, 0.001, 0.003, 0.01, 0.03, 0.1)), 6)
  D[1, ] <- 0
  ages <- mt_axis(c(0, 20, 40, 60, 70, 80), c(20, 40, 60, 70, 80, Inf))
  expect_warning(est <- fit_model(lee_carter(), D, E, ages, 2001:2015), class = "mt_warning_sparse")
  expect_true(all(is.finite(est$a)))
  expect_equal(est$b[1], 0)
})

test_that("implausible forecast rates are reported", {
  arr <- array(0.01, c(2, 3, 1, 10))
  arr[1, 3, 1, 1:4] <- 2
  expect_warning(check_draws_plausible(arr, mt_axis(c(0, 50), c(50, Inf)), data.table(label = "s"), "Model"),
                 "40.0% of the draws", class = "mt_warning_implausible")
})

test_that("crossed hierarchies of regions and causes aggregate exposures once over causes", {
  g <- expand.grid(year = 2015:2020, age = "80+", region = c("r1", "r2"), cause = c("c1", "c2", "c3"),
                   stringsAsFactors = FALSE)
  set.seed(3)
  g$deaths <- stats::rpois(nrow(g), 50)
  p <- expand.grid(year = 2015:2020, age = "80+", region = c("r1", "r2"), stringsAsFactors = FALSE)
  p$exposure <- ifelse(p$region == "r1", 1000, 3000)
  x <- mt_data(deaths = g, exposure = p, by = c("region", "cause"))
  expect_error(mt_hierarchy(region = "all", additive = c("deaths", "rates")), class = "mt_error_hierarchy")
  h <- mt_hierarchy(region = "all", cause = list(total = "all", additive = "rates"))
  xa <- mt_aggregate_strata(x, h)
  st <- mt_strata(xa)
  e <- mt_obs(xa, "exposure")
  expect_equal(unique(e[stratum == st[region == "all" & cause == "all"]$stratum]$value), 4000)
  expect_equal(unique(e[stratum == st[region == "r1" & cause == "all"]$stratum]$value), 1000)
  expect_true(all(mt_validate(xa, hierarchy = h)$severity == "info"))
  fc <- mt_forecast(mt_fit(reconciled(rwd(), h, method = "wls"), xa), h = 2, n_draws = 20, seed = 1)
  a <- fc$rates[1, 2, , ]
  s <- fc$strata
  expect_equal(a[s$region == "all" & s$cause == "all", ], colSums(a[s$region == "all" & s$cause != "all", ]),
               tolerance = 1e-8)
  dep <- a[s$region != "all" & s$cause == "all", ] * c(1000, 3000)
  expect_equal(a[s$region == "all" & s$cause == "all", ] * 4000, colSums(dep), tolerance = 1e-8)
})

test_that("suppressed cells need totals or bounds and respect bounds", {
  d <- data.frame(region = c("A", "B", "C", "all"), year = 2020, age = "90+", deaths = c(10, NA, NA, 40),
                  exposure = c(100, 300, 100, 500))
  x <- suppressWarnings(mt_data(deaths = d, exposure = d, by = "region"))
  expect_error(mt_fill_suppressed(x, mt_hierarchy(region = "Total")), class = "mt_error_fill")
  expect_error(mt_fill_suppressed(x), class = "mt_error_fill")
  expect_warning(y <- mt_fill_suppressed(x, mt_hierarchy(region = "all"), lower = 1, upper = 2),
                 class = "mt_warning_fill")
  expect_true(all(mt_obs(y)[provenance == "estimated"]$value <= 2))
  d2 <- d[1:3, ]
  x2 <- suppressWarnings(mt_data(deaths = d2, exposure = d2, by = "region"))
  z <- mt_fill_suppressed(x2, lower = 1, upper = 4)
  v <- mt_obs(z)[provenance == "estimated"]$value
  expect_equal(v, c(4, 4))
  z2 <- mt_fill_suppressed(x2, lower = 1, upper = 40)
  expect_equal(mt_obs(z2)[provenance == "estimated"]$value, c(30, 10))
})

test_that("smoothing leaves age 0 alone and respects group widths", {
  d <- data.frame(year = 2019, age = c("0", "1-4", "5-9", "10-14", "15-19", "20-24"),
                  deaths = c(40, 4, 2, 3, 8, 12), exposure = c(1e4, 4e4, 5e4, 5e4, 5e4, 5e4))
  x <- mt_data(deaths = d, exposure = d)
  s <- mt_smooth(x, method = "whittaker")
  expect_equal(mt_obs(s)$value[1], 40)
  expect_identical(mt_obs(s)$provenance[1], "observed")
})

test_that("future exposures are matched by names, not by position", {
  x <- mt_example(deaths = "single", years = 2001:2015)
  fc <- mt_forecast(mt_fit(rwd(), x), h = 2, n_draws = 10, seed = 1)
  E <- mt_array(x, "exposure")[, "2015", , drop = FALSE]
  d1 <- mt_deaths(fc, exposure = E, noise = FALSE)
  rev_E <- E[, , rev(dimnames(E)[[3]]), drop = FALSE]
  d2 <- mt_deaths(fc, exposure = rev_E, noise = FALSE)
  expect_equal(d1, d2)
  expect_warning(mt_deaths(fc, exposure = unname(E), noise = FALSE), class = "mt_warning_exposure_order")
  d3 <- mt_deaths(fc, exposure = mt_subset(x, periods = c(2015, 2015)), noise = FALSE)
  expect_equal(d1, d3)
  bad <- E[, , 1, drop = FALSE]
  expect_error(mt_deaths(fc, exposure = bad), class = "mt_error_argument")
})

test_that("observed indicators, contrasts and inequality indices", {
  x <- mt_example_regions(years = 2015:2019, regions = 2)
  e <- life_expectancy(x)
  expect_true(all(e$q0.025 <= e$q0.5 & e$q0.5 <= e$q0.975))
  a <- asr(x, standard = "esp1976")
  expect_true(all(a$standard == "esp1976" & a$per == 1e5))
  ct <- mt_contrast(x, var = "quintile", ref = "Q1")
  expect_setequal(unique(ct$level), paste0("Q", 2:5))
  expect_true(all(c("sex", "region", "p_positive") %in% names(ct)))
  iq <- mt_inequality(x, "quintile", order = paste0("Q", 5:1))
  expect_setequal(unique(iq$index), c("sii", "rii"))
  expect_error(mt_contrast(x, var = "quintile", ref = "Q9"), class = "mt_error_argument")
  expect_equal(sum(standard_population("who2000")$weight), 100035)
  expect_equal(sum(standard_population("esp1976")$weight), 1e5)
})

test_that("relabelling, window steps, exports and conversions work", {
  d <- data.frame(Sex = rep(1:2, each = 6), year = 2001:2006, age = "90+", deaths = 5:16, exposure = 100)
  x <- mt_data(deaths = d, exposure = d, by = "Sex")
  y <- mt_relabel(x, "Sex", c("1" = "male", "2" = "female"))
  expect_setequal(mt_strata(y)$Sex, c("male", "female"))
  w <- mt_prep_apply(mt_prep(step_window(3)), y)
  expect_identical(mt_grid(w)$period$labels, c("2004", "2005", "2006"))
  fc <- mt_forecast(mt_fit(rwd(), y), h = 2, n_draws = 5, seed = 1)
  f <- mt_export(asr(fc), tempdir(), file = "asr-test.csv", digits = 3)
  expect_true(file.exists(f))
  expect_warning(mt_forecast(mt_fit(rwd(), y), h = 3, to = 2008, n_draws = 5), class = "mt_warning_argument")
  skip_if_not_installed("demography")
  xe <- mt_example(deaths = "single", years = 2001:2005)
  dd <- as_demogdata(xe)
  expect_s3_class(dd, "demogdata")
  expect_setequal(names(dd$rate), mt_strata(xe)$label)
  skip_if_not_installed("StMoMo")
  sd <- as_stmomodata(xe, "male")
  expect_s3_class(sd, "StMoMoData")
  expect_equal(sum(sd$Dxt), sum(mt_array(xe, "deaths")[, , "male"]))
})
