test_that("GBE exports are parsed with carried row labels and statistical symbols", {
  g <- mt_read_gbe(test_path("fixtures", "gbe-synthetic.csv"))
  expect_setequal(unique(g$region), c("Deutschland", "Bremen"))
  expect_identical(sort(unique(g$year)), c(2019L, 2020L))
  expect_equal(nrow(g), 2 * 2 * 2 * 2 * 3)
  expect_true(all(g[symbol == "-"]$value == 0))
  expect_true(all(is.na(g[symbol %in% c(".", "X")]$value)))
  expect_equal(g[region == "Bremen" & year == 2020 & sex == "Weiblich" & grepl("^Alle", cause) &
                   grepl("^85", age)]$value, 4)
})

test_that("suppressed cells are raked to the published totals", {
  d <- data.frame(region = c("A", "B", "C", "all"), year = 2020, age = "90+", deaths = c(10, NA, NA, 40),
                  exposure = c(100, 300, 100, 500))
  x <- suppressWarnings(mt_data(deaths = d, exposure = d, by = "region"))
  expect_warning(mt_validate(x), NA)
  y <- mt_fill_suppressed(x, mt_hierarchy(region = "all"))
  o <- mt_obs(y)
  expect_equal(o[label == "B"]$value, 22.5)
  expect_equal(o[label == "C"]$value, 7.5)
  expect_identical(o[label == "B"]$provenance, "estimated")
  # two crossed variables: the margins identify the cells here
  tr <- round(outer(c(A = 50, B = 30, C = 20), c(c1 = 0.5, c2 = 0.3, c3 = 0.2)) * 10)
  g <- expand.grid(region = c("A", "B", "C", "all"), cause = c("c1", "c2", "c3", "all"), stringsAsFactors = FALSE)
  full <- rbind(cbind(tr, all = rowSums(tr)), all = c(colSums(tr), sum(tr)))
  g$deaths <- full[cbind(match(g$region, rownames(full)), match(g$cause, colnames(full)))]
  g$exposure <- 1000
  g$year <- 2020
  g$age <- "90+"
  obs <- g
  obs$deaths[obs$region == "B" & obs$cause %in% c("c2", "c3")] <- NA
  obs$deaths[obs$region == "C" & obs$cause == "c3"] <- NA
  x2 <- suppressWarnings(mt_data(deaths = obs, exposure = obs, by = c("region", "cause")))
  y2 <- mt_fill_suppressed(x2, mt_hierarchy(region = "all", cause = "all"))
  o2 <- mt_obs(y2)[provenance == "estimated"]
  expect_equal(o2[label == "B.c2"]$value, 90, tolerance = 1e-6)
  expect_equal(o2[label == "C.c3"]$value, 40, tolerance = 1e-6)
  expect_lt(attr(y2, "fill")$max_gap, 1e-6)
  expect_error(mt_fit(rwd(), x), class = "mt_error_model_data")
})

test_that("31 December populations give exposures and the German data read consistently", {
  p <- data.frame(year = 2018:2020, age = "90+", population = c(100, 110, 130))
  x <- mt_as_exposure(mt_add_population(mt_data(), p), type = "dec31")
  e <- mt_obs(x, "exposure")
  expect_equal(e$period_lo, c(2019, 2020))
  expect_equal(e$value, c(105, 120))
  dd <- as.data.frame(de_deaths)
  dd <- dd[dd$year %in% 2017:2019, ]
  expect_warning(y <- mt_data(deaths = dd[dd$cause == "all", setdiff(names(dd), "cause")],
                              population = de_population[de_population$year %in% 2016:2019, ],
                              by = c("region", "sex")), class = "mt_warning_missing")
  y <- mt_as_exposure(y, type = "dec31")
  ag <- unique(mt_obs(y, "deaths")[, list(age_lo)])[order(age_lo)]
  y <- suppressWarnings(mt_regroup(y, age = ag$age_lo))
  y <- mt_fill_suppressed(y, mt_hierarchy(region = "Deutschland"))
  y <- mt_subset(y, region != "Deutschland")
  expect_identical(nrow(mt_strata(y)), 32L)
  expect_identical(mt_grid(y)$age$labels[c(1, 2, 18)], c("0", "1-14", "90+"))
  pub <- sum(dd$deaths[dd$region == "Deutschland" & dd$cause == "all"])
  expect_equal(sum(mt_obs(y)$value), pub, tolerance = 1e-6)
  expect_true(all(is.finite(mt_obs(y, "exposure")$value)))
  z <- suppressWarnings(mt_data(deaths = dd[dd$year >= 2018, ],
                                population = de_population[de_population$year >= 2017, ],
                                by = c("region", "sex", "cause")), classes = "mt_warning_missing")
  z <- mt_as_exposure(z, type = "dec31")
  z <- mt_subset(z, region %in% c("Berlin", "Bremen", "Deutschland"), ages = c(80, Inf))
  expect_warning(z <- mt_group_levels(z, "cause", de_cause_groups, other = "rest", keep = "all",
                                      additive = "rates"), class = "mt_warning_group")
  expect_setequal(unique(mt_strata(z)$cause), c("all", names(de_cause_groups), "rest"))
})

test_that("text values with suppression markers are read as missing or zero", {
  d <- data.frame(year = 2020, age = c("80-84", "85-89", "90+"), deaths = c("12", " . ", "-"),
                  exposure = c("1000", "800", "500"))
  expect_warning(x <- mt_data(deaths = d, exposure = d), class = "mt_warning_missing")
  o <- mt_obs(x)
  expect_identical(o$value, c(12, 0))
  expect_true(is.numeric(mt_obs(x, "exposure")$value))
  d$deaths[1] <- "<5"
  expect_warning(suppressWarnings(mt_data(deaths = d, exposure = d), classes = "mt_warning_missing"),
                 class = "mt_warning_markers")
  x2 <- suppressWarnings(mt_data(deaths = d, exposure = d, missing = c("<5", "."), zero = c("-", "0*")),
                         classes = "mt_warning_missing")
  expect_identical(nrow(mt_obs(x2)), 1L)
})

test_that("tables without a stratum variable are shared over its levels", {
  d <- expand.grid(year = 2019:2020, age = "90+", sex = c("f", "m"), cause = c("a", "b"), stringsAsFactors = FALSE)
  d$deaths <- seq_len(nrow(d))
  p <- expand.grid(year = 2019:2020, age = "90+", sex = c("f", "m"), stringsAsFactors = FALSE)
  p$exposure <- c(100, 110, 90, 95)
  x <- mt_data(deaths = d, exposure = p, by = c("sex", "cause"))
  e <- mt_obs(x, "exposure")
  expect_identical(nrow(e), 8L)
  st <- mt_strata(x)
  expect_equal(e[stratum == st[sex == "f" & cause == "b"]$stratum]$value, c(100, 110))
})

test_that("levels are grouped with summed or shared exposures", {
  d <- expand.grid(year = 2020, age = "90+", region = c("r1", "r2", "r3"), stringsAsFactors = FALSE)
  d$deaths <- c(10, 20, 30)
  d$exposure <- c(100, 200, 300)
  x <- mt_data(deaths = d, exposure = d, by = "region")
  y <- mt_group_levels(x, "region", list(north = c("r1", "r2")))
  expect_setequal(mt_strata(y)$region, c("north", "r3"))
  st <- mt_strata(y)
  expect_equal(mt_obs(y)[stratum == st[region == "north"]$stratum]$value, 30)
  expect_equal(mt_obs(y, "exposure")[stratum == st[region == "north"]$stratum]$value, 300)
  z <- mt_group_levels(x, "region", list(a = "r1"), other = "b")
  expect_setequal(mt_strata(z)$region, c("a", "b"))
  expect_error(mt_group_levels(x, "region", list(a = "r9")), class = "mt_error_argument")
  expect_error(mt_group_levels(x, "region", list(a = "r1", b = "r1")), class = "mt_error_argument")
  c2 <- expand.grid(year = 2020, age = "90+", cause = c("all", "c1", "c2", "c3"), stringsAsFactors = FALSE)
  c2$deaths <- c(60, 10, 20, 30)
  c2$exposure <- 1000
  xc <- mt_data(deaths = c2, exposure = c2, by = "cause")
  yc <- mt_group_levels(xc, "cause", list(c12 = c("c1", "c2")), other = "rest", keep = "all", additive = "rates")
  st <- mt_strata(yc)
  expect_equal(mt_obs(yc, "exposure")[stratum == st[cause == "c12"]$stratum]$value, 1000)
  expect_equal(mt_obs(yc)[stratum == st[cause == "rest"]$stratum]$value, 30)
  expect_equal(mt_obs(yc)[stratum == st[cause == "all"]$stratum]$value, 60)
})

test_that("mt_subset restricts strata, ages and periods", {
  x <- mt_example(years = 2001:2019)
  y <- mt_subset(x, sex == "female", ages = c(60, Inf), periods = c(2010, 2019))
  expect_true(all(mt_strata(y)$sex == "female"))
  expect_gte(min(mt_obs(y)$age_lo), 60)
  expect_identical(range(mt_obs(y)$period_lo), c(2010, 2019))
  expect_error(mt_subset(x, sex == "none"), class = "mt_error_argument")
  expect_error(mt_subset(x, ages = 60), class = "mt_error_argument")
  z <- mt_subset(x, sex == "male", drop = TRUE)
  expect_identical(z$by, character())
  expect_identical(mt_strata(z)$label, "all")
})
