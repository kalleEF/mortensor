# Stress tests: simulated data with realistic defects (dev/stress-tests.md).
# Each scenario builds a data set with one defect, runs the usual pipeline
# and records whether mortensor stops with an error, warns, or continues
# silently, together with the message. Run from the package root:
#   Rscript dev/stress-tests.R
suppressMessages(devtools::load_all(quiet = TRUE))
library(data.table)

# A clean data set: deaths and populations by sex and region, 5-year age
# groups 0-4, ..., 85+, years 2001-2020, Gompertz-like rates improving by
# 1.5% a year, Poisson deaths
sim_clean <- function(regions = c("North", "South", "Island"), years = 2001:2020, pop_scale = c(2e5, 1e5, 5e3),
                      seed = 1) {
  set.seed(seed)
  lo <- seq(0, 85, 5)
  g <- CJ(region = regions, sex = c("female", "male"), year = years, age_lo = lo)
  g[, age := ifelse(age_lo == 85, "85+", sprintf("%d-%d", age_lo, age_lo + 4))]
  g[, population := round(pop_scale[match(region, regions)] / 18 * exp(-age_lo / 60))]
  g[, rate := exp(-9.5 + 0.085 * (age_lo + 2.5) + 0.25 * (sex == "male") - 0.015 * (year - 2001))]
  g[age_lo == 0, rate := rate * 8]
  g[, deaths := rpois(.N, population * pmin(rate, 0.9))]
  g[, c("age_lo", "rate") := NULL]
  as.data.frame(g)
}

run_case <- function(name, defect, f) {
  warns <- character()
  res <- withCallingHandlers(tryCatch(f(), error = function(e) e),
                             warning = function(w) {
                               warns <<- c(warns, paste0("[", class(w)[1], "] ", conditionMessage(w)))
                               invokeRestart("muffleWarning")
                             })
  err <- if (inherits(res, "error")) paste0("[", class(res)[1], "] ", conditionMessage(res)) else ""
  outcome <- if (nzchar(err)) "error" else if (length(warns)) "warning" else "silent"
  detail <- if (!inherits(res, "error") && is.character(res)) res else ""
  data.table(case = name, defect = defect, outcome = outcome, error = err,
             warnings = paste(unique(warns), collapse = " | "), detail = detail)
}

pipeline <- function(d, p = d, by = c("region", "sex"), ...) {
  x <- mt_data(deaths = d, population = p, by = by, ...)
  x <- mt_as_exposure(x)
  fit <- mt_fit(lee_carter(), x)
  fc <- mt_forecast(fit, h = 5, n_draws = 50, seed = 1)
  e0 <- life_expectancy(fc, at = 0)
  sprintf("ok: %d strata, e0 2025 median %.1f-%.1f", nrow(mt_strata(x)), min(e0[period == max(period)]$q0.5),
          max(e0[period == max(period)]$q0.5))
}

clean <- sim_clean()
cases <- list()
add <- function(name, defect, f) cases[[length(cases) + 1]] <<- run_case(name, defect, f)

add("clean", "none (reference)", function() pipeline(clean))

# Coding of variables -------------------------------------------------------
add("total_rows_age", "rows with age 'All ages' (table totals) in the deaths and populations", function() {
  t <- as.data.table(clean)[, .(deaths = sum(deaths), population = sum(population), age = "All ages"),
                           by = .(region, sex, year)]
  pipeline(rbind(as.data.table(clean), t, use.names = TRUE))
})
add("total_rows_sex", "rows with sex 'Persons' (both sexes) next to female and male", function() {
  t <- as.data.table(clean)[, .(deaths = sum(deaths), population = sum(population), sex = "Persons"),
                           by = .(region, year, age)]
  pipeline(rbind(as.data.table(clean), t, use.names = TRUE))
})
add("labels_mismatch", "sex coded 'female'/'male' in deaths, 'F'/'M' in populations", function() {
  p <- copy(as.data.table(clean))[, sex := ifelse(sex == "female", "F", "M")]
  pipeline(clean, p)
})
add("labels_case_space", "region 'North ' (trailing space) and 'south' (case) in the populations", function() {
  p <- copy(as.data.table(clean))[region == "North", region := "North "][region == "South", region := "south"]
  pipeline(clean, p)
})
add("na_stratum", "sex missing (NA) in some death rows", function() {
  d <- copy(as.data.table(clean))
  d[sample(.N, 20), sex := NA]
  pipeline(d)
})
add("year_factor", "year and age stored as factors (read.csv with stringsAsFactors)", function() {
  d <- copy(as.data.table(clean))[, `:=`(year = factor(year), age = factor(age))]
  pipeline(d)
})
add("year_character_fy", "years as financial years '2001/02'", function() {
  d <- copy(as.data.table(clean))[, year := sprintf("%d/%02d", year, (year + 1) %% 100)]
  pipeline(d)
})
add("age_labels_variants", "age labels '0 to 4', '5 - 9', '85 years and over', 'Under 1'", function() {
  l <- c("Under 1", "1 to 4", "5 - 9", "10 to 14", "85 years and over", "85 and older", "90 & over", "00-04",
         "0-4 years", "85+ years", "Age 85+", "95 Jahre und mehr", "unter 1 Jahr", "1 bis 4 Jahre")
  r <- lapply(l, function(z) tryCatch(paste(parse_age_labels(z)[, c("age_lo", "age_hi")], collapse = ","),
                                      error = function(e) "FAILED"))
  paste(sprintf("%s -> %s", l, unlist(r)), collapse = "; ")
})
add("age_numeric_open", "single ages 0-100 where 100 means 100+ (open, as in many tables), no age_open", function() {
  d <- CJ(sex = c("female", "male"), year = 2001:2020, age = 0:100)
  d[, population := round(1e5 * exp(-age / 40))]
  d[age == 100, population := 2000]
  d[, deaths := rpois(.N, population * pmin(exp(-9.5 + 0.09 * age), 0.6))]
  pipeline(d, by = "sex")
})
add("thousands_separator", "populations as text with thousands separators ('12,345')", function() {
  p <- copy(as.data.table(clean))[, population := format(population, big.mark = ",", trim = TRUE)]
  pipeline(clean, p)
})
add("decimal_comma", "populations as text with decimal commas ('1234,5'), as in German exports", function() {
  p <- copy(as.data.table(clean))[, population := sub("\\.", ",", sprintf("%.1f", population + 0.5))]
  pipeline(clean, p)
})

# Structure of the tables ---------------------------------------------------
add("duplicate_rows", "a block of death rows appears twice (appended file)", function() {
  d <- as.data.table(clean)
  pipeline(rbind(d, d[year == 2010]))
})
add("missing_year", "year 2010 missing in all tables", function() {
  pipeline(as.data.table(clean)[year != 2010])
})
add("missing_year_one_stratum", "year 2010 missing for one region in the deaths only", function() {
  d <- as.data.table(clean)
  pipeline(d[!(year == 2010 & region == "Island")], d)
})
add("population_missing_stratum", "populations missing for one region", function() {
  d <- as.data.table(clean)
  pipeline(d, d[region != "Island"])
})
add("different_open_interval", "deaths to 90+ (85-89, 90+), populations to 85+", function() {
  d <- as.data.table(clean)
  top <- d[age == "85+"]
  d2 <- rbind(d[age != "85+"], copy(top)[, `:=`(age = "85-89", deaths = round(deaths * 0.6))],
              copy(top)[, `:=`(age = "90+", deaths = deaths - round(deaths * 0.6))])
  pipeline(d2, d)
})
add("grouping_change_over_time", "deaths 85+ until 2010, 85-89 and 90+ from 2011 (populations 85+)", function() {
  d <- as.data.table(clean)
  late <- d[age == "85+" & year > 2010]
  d2 <- rbind(d[!(age == "85+" & year > 2010)], copy(late)[, `:=`(age = "85-89", deaths = round(deaths * 0.6))],
              copy(late)[, `:=`(age = "90+", deaths = deaths - round(deaths * 0.6))])
  pipeline(d2, d)
})
add("overlapping_groups", "deaths contain both '85+' and '85-89', '90+' for the same year", function() {
  d <- as.data.table(clean)
  top <- d[age == "85+" & year == 2015]
  d2 <- rbind(d, copy(top)[, `:=`(age = "85-89", deaths = round(deaths * 0.6))],
              copy(top)[, `:=`(age = "90+", deaths = deaths - round(deaths * 0.6))])
  pipeline(d2, d)
})
add("short_series", "only 4 years of data", function() pipeline(as.data.table(clean)[year >= 2017]))

# Values ----------------------------------------------------------------------
add("zero_population", "population 0 in a cell with deaths > 0 (coding error)", function() {
  p <- copy(as.data.table(clean))
  p[region == "Island" & age == "85+" & year == 2005, population := 0]
  pipeline(clean, p)
})
add("deaths_exceed_population", "deaths larger than the population in one cell", function() {
  d <- copy(as.data.table(clean))
  d[region == "North" & age == "40-44" & year == 2005, deaths := population * 2]
  pipeline(d)
})
add("negative_deaths", "negative death count (correction record)", function() {
  d <- copy(as.data.table(clean))
  d[region == "South" & age == "50-54" & year == 2012, deaths := -3]
  pipeline(d)
})
add("suppressed_lt5", "counts below 5 published as '<5'", function() {
  d <- copy(as.data.table(clean))[, deaths := ifelse(deaths < 5, "<5", as.character(deaths))]
  pipeline(d)
})
add("suppressed_lt5_declared", "counts below 5 published as '<5', declared with missing = '<5'", function() {
  d <- copy(as.data.table(clean))[, deaths := ifelse(deaths < 5, "<5", as.character(deaths))]
  pipeline(d, missing = "<5")
})
add("suppressed_no_totals", "suppressed cells ('.') without published totals, then mt_fill_suppressed()", function() {
  d <- copy(as.data.table(clean))[, deaths := ifelse(deaths < 3, ".", as.character(deaths))]
  x <- mt_as_exposure(mt_data(deaths = d, population = d, by = c("region", "sex")))
  x <- mt_fill_suppressed(x, mt_hierarchy(region = "all"))
  "filled"
})
add("sparse_small_region", "region of 5,000 people: many zero death cells (Lee-Carter)", function() {
  x <- mt_as_exposure(mt_data(deaths = clean, population = clean, by = c("region", "sex")))
  x <- mt_subset(x, region == "Island")
  o <- mt_obs(x)
  fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 50, seed = 1)
  r <- fc$rates
  sprintf("ok: %.0f%% zero cells; forecast rates finite: %s; max/min ratio of median rates at age 0-4: %.1f",
          100 * mean(o$value == 0), all(is.finite(r)), max(apply(r[1, , 1, ], 1, median)) /
            min(apply(r[1, , 1, ], 1, median)))
})
add("sparse_rwd_apc", "same sparse region with rwd() and li_lee()", function() {
  x <- mt_as_exposure(mt_data(deaths = clean, population = clean, by = c("region", "sex")))
  out <- vapply(list(rwd = rwd(), li_lee = li_lee(), tucker3 = tucker3()), function(sp) {
    r <- tryCatch({
      fc <- mt_forecast(mt_fit(sp, x), h = 5, n_draws = 50, seed = 1)
      sprintf("finite=%s", all(is.finite(fc$rates)))
    }, error = function(e) paste("ERROR", conditionMessage(e)))
    r
  }, "")
  paste(names(out), out, sep = ": ", collapse = "; ")
})
add("all_zero_age_group", "an age group with zero deaths in every year (small region, ages 5-9)", function() {
  d <- copy(as.data.table(clean))[region == "Island"][age %in% c("5-9", "10-14"), deaths := 0]
  pipeline(d, by = c("region", "sex"))
})
add("covid_shock", "a mortality shock (+20%) in the last two years", function() {
  d <- copy(as.data.table(clean))[year >= 2019, deaths := round(deaths * 1.2)]
  pipeline(d)
})

# Scale ------------------------------------------------------------------------
add("many_strata", "400 districts x 2 sexes x 18 ages x 20 years (288,000 rows): reading and grouping time", function() {
  d <- sim_clean(regions = sprintf("D%03d", 1:400), pop_scale = rep(5e4, 400))
  t0 <- Sys.time()
  x <- mt_data(deaths = d, population = d, by = c("region", "sex"))
  t1 <- Sys.time()
  x <- mt_as_exposure(x)
  t2 <- Sys.time()
  y <- mt_group_levels(x, "region", function(r) paste0("G", (as.integer(substr(r, 2, 4)) - 1) %/% 40 + 1))
  t3 <- Sys.time()
  h <- mt_hierarchy(region = "all", sex = "both")
  z <- mt_aggregate_strata(y, h)
  t4 <- Sys.time()
  s <- function(a, b) round(as.numeric(b - a, units = "secs"), 1)
  sprintf("mt_data %ss, mt_as_exposure %ss, mt_group_levels %ss, mt_aggregate_strata %ss", s(t0, t1), s(t1, t2),
          s(t2, t3), s(t3, t4))
})

res <- rbindlist(cases)
saveRDS(res, Sys.getenv("STRESS_OUT", "dev/stress-results.rds"))
options(width = 200)
for (i in seq_len(nrow(res))) {
  r <- res[i]
  cat(sprintf("\n## %s (%s)\n%s\n", r$case, r$outcome, r$defect))
  if (nzchar(r$error)) cat("  error:", r$error, "\n")
  if (nzchar(r$warnings)) cat("  warnings:", r$warnings, "\n")
  if (nzchar(r$detail)) cat("  result:", r$detail, "\n")
}
