# Benchmarks of the forecasting pipeline, shown in vignettes/benchmarks.Rmd:
# 1. interval calibration on simulated data with known truth,
# 2. rolling-origin backtests on four real data sets,
# 3. the behaviour on input tables with realistic defects.
# The results are stored in vignettes/benchmarks-results.rds; the vignette
# reports them as they are. There are no thresholds.
#
# Run from the package root:
#   Rscript dev/benchmarks.R
# BENCH_QUICK=true: short run (fewer replicates, origins and draws).
# BENCH_PARTS=sim,real,defects: the parts to run (default all); parts not run
# keep their stored results. BENCH_CORES: cores for the backtests (default 2).
# BENCH_OUT: output file.

suppressMessages(devtools::load_all(quiet = TRUE))
suppressMessages(library(data.table))
quick <- identical(Sys.getenv("BENCH_QUICK"), "true")
parts <- strsplit(Sys.getenv("BENCH_PARTS", "sim,real,defects"), ",")[[1]]
cores <- as.integer(Sys.getenv("BENCH_CORES", "2"))
out_file <- Sys.getenv("BENCH_OUT", "vignettes/benchmarks-results.rds")
log_step <- function(...) message(format(Sys.time(), "%H:%M:%S "), sprintf(...))
res <- if (file.exists(out_file)) readRDS(out_file) else list()

# Scores of a backtest: overall and by horizon, for rates and life expectancy
score_tables <- function(bt) {
  le <- "target" %in% names(bt) && "life_expectancy" %in% bt$target
  list(rate = mt_score(bt, by = "model"),
       rate_h = mt_score(bt, by = c("model", "horizon")),
       le = if (le) mt_score(bt, by = "model", target = "life_expectancy"),
       le_h = if (le) mt_score(bt, by = c("model", "horizon"), target = "life_expectancy"))
}

# ---- 1. Known truth -----------------------------------------------------------
# simulate_lc_dgp(): the Lee-Carter model with Poisson deaths, 10 strata, 5-year
# age groups 30-94, a random walk with drift per stratum. sd_obs = 0.03 adds
# independent normal variation to the log rates of every cell (extra-Poisson
# variation, which no model below except the negative binomial one assumes).
if ("sim" %in% parts) {
  n_rep <- if (quick) 3 else 20
  n_train <- 30
  h <- 10
  n_draws <- if (quick) 200 else 500
  generators <- list("Lee-Carter, Poisson deaths" = 0, "Lee-Carter, extra variation" = 0.03)
  models <- list("Lee-Carter" = lee_carter(), "Lee-Carter (negative binomial)" = lee_carter(family = "nb"),
                 "random walk with drift" = rwd(),
                 "ARIMA(0,1,1) with drift" = arima_series(order = c(0, 1, 1), drift = TRUE),
                 "Lee-Carter + random walk" = ensemble(lc = lee_carter(), rw = rwd()))
  t0 <- Sys.time()
  sim <- rbindlist(lapply(names(generators), function(g) {
    rbindlist(lapply(seq_len(n_rep), function(r) {
      log_step("simulation: %s, replicate %d", g, r)
      x <- simulate_lc_dgp(n_t = n_train + h, seed = 1000 + r, sd_obs = generators[[g]])
      bt <- mt_compare(models, x, origins = 2000 + n_train, h = h, n_draws = n_draws, seed = r, cores = cores)
      cbind(generator = g, rep = r, as.data.table(bt))
    }))
  }))
  res$sim <- list(
    settings = list(n_rep = n_rep, n_train = n_train, h = h, n_draws = n_draws, sd_obs = unlist(generators)),
    overall = sim[, mt_score(.SD, by = "model"), by = generator],
    by_horizon = sim[, mt_score(.SD, by = c("model", "horizon")), by = generator],
    minutes = as.numeric(Sys.time() - t0, units = "mins"))
}

# ---- 2. Real data -------------------------------------------------------------
if ("real" %in% parts) {
  n_draws <- if (quick) 100 else 300
  single <- list("Lee-Carter" = lee_carter(), "Lee-Carter (negative binomial)" = lee_carter(family = "nb"),
                 "random walk with drift" = rwd(), "Lee-Carter + random walk" = ensemble(lc = lee_carter(), rw = rwd()),
                 "Lee-Carter, recalibrated" = recalibrated(lee_carter()))
  joint <- c(single, list("Li-Lee" = li_lee(common_by = "sex")))

  # England & Wales men (HMD via StMoMo), single ages 0-100, 1961-2011
  ew <- as_mt_data(StMoMo::EWMaleData)
  # France women (HMD via demography), single ages 0-100, 1950-2006
  fr <- mt_subset(as_mt_data(demography::fr.mort), sex == "female", ages = c(0, 100), periods = c(1950, 2006),
                  drop = TRUE)
  # England (ONS), sex x deprivation quintile, 5-year groups 0-4, ..., 90+, 2001-2019
  ons <- mt_data(deaths = ons_data, population = ons_data, by = c("sex", "imd_decile"), age = "age_group",
                 value = c(deaths = "all_causes", population = "population"))
  ons <- mt_as_exposure(ons, type = "midyear")
  ons <- mt_group_levels(ons, "imd_decile", function(d) paste0("Q", ceiling(as.numeric(d) / 2)), name = "quintile")
  # Germany (Destatis), 16 Laender x sex, all causes, age groups 0, 1-14, 15-19, ..., 90+, 1998-2019;
  # suppressed counts estimated from the published totals, as in vignette("laender-germany")
  de <- suppressWarnings(mt_data(deaths = de_deaths, population = de_population, by = c("region", "sex", "cause")))
  de <- mt_as_exposure(de, type = "dec31")
  de <- mt_subset(de, periods = c(1998, 2019))
  de <- mt_regroup(de, age = "deaths")
  de <- mt_fill_suppressed(de, mt_hierarchy(region = "Deutschland", cause = list(total = "all", additive = "rates")))
  de <- mt_subset(de, cause == "all", region != "Deutschland", drop = TRUE)

  sets <- list(
    ew = list(label = "England & Wales men", data = ew, models = single,
              origins = if (quick) c(1996, 2001) else 1996:2001, h = 10),
    france = list(label = "France women", data = fr, models = single,
                  origins = if (quick) c(1991, 1996) else 1991:1996, h = 10),
    ons = list(label = "England by sex and deprivation quintile", data = ons, models = joint,
               origins = if (quick) c(2010, 2014) else 2010:2014, h = 5),
    germany = list(label = "German Laender by sex", data = de, models = joint,
                   origins = if (quick) c(2010, 2014) else 2010:2014, h = 5))
  res$real <- lapply(sets, function(s) {
    log_step("backtest: %s", s$label)
    t0 <- Sys.time()
    bt <- mt_compare(s$models, s$data, origins = s$origins, h = s$h, n_draws = n_draws,
                     quantities = "life_expectancy", cores = cores)
    c(list(label = s$label, origins = s$origins, h = s$h, n_draws = n_draws,
           ages = range(mt_grid(s$data)$age$lower), periods = range(mt_grid(s$data)$period$lower),
           n_strata = nrow(mt_strata(s$data))),
      score_tables(bt), list(minutes = as.numeric(Sys.time() - t0, units = "mins")))
  })
}

# ---- 3. Defective input tables --------------------------------------------------
# A clean data set (two regions of 200,000 and 100,000 people, two sexes,
# 5-year age groups 0-4, ..., 85+, 2001-2020, Poisson deaths) and one
# realistic defect per scenario. Each scenario runs the usual pipeline and
# records whether mortensor stops with an error, warns or continues silently.
# A third region of 5,000 people is used in the scenarios on small populations.
if ("defects" %in% parts) {
  sim_clean <- function(regions = c("North", "South"), years = 2001:2020, pop_scale = c(2e5, 1e5), seed = 1) {
    with_seed(seed, {
      lo <- seq(0, 85, 5)
      g <- CJ(region = regions, sex = c("female", "male"), year = years, age_lo = lo)
      g[, age := ifelse(age_lo == 85, "85+", sprintf("%d-%d", age_lo, age_lo + 4))]
      g[, population := round(pop_scale[match(region, regions)] / 18 * exp(-age_lo / 60))]
      g[, rate := exp(-9.5 + 0.085 * (age_lo + 2.5) + 0.25 * (sex == "male") - 0.015 * (year - 2001))]
      g[age_lo == 0, rate := rate * 8]
      g[, deaths := stats::rpois(.N, population * pmin(rate, 0.9))]
      g[, c("age_lo", "rate") := NULL]
      g[]
    })
  }
  run_case <- function(group, defect, f) {
    warns <- character()
    out <- withCallingHandlers(tryCatch(f(), error = function(e) e), warning = function(w) {
      warns <<- c(warns, sprintf("[%s] %s", class(w)[1], conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
    err <- if (inherits(out, "error")) sprintf("[%s] %s", class(out)[1], conditionMessage(out)) else ""
    data.table(group = group, defect = defect,
               outcome = if (nzchar(err)) "error" else if (length(warns)) "warning" else
                 if (is.character(out) && startsWith(out, "mt_validate")) "mt_validate() finding" else "no message",
               error = err, warnings = paste(unique(warns), collapse = " | "),
               result = if (is.character(out)) out else "")
  }
  # the workflow of vignette("own-data"): read, exposures, mt_validate(), fit, forecast
  pipeline <- function(d, p = d, by = c("region", "sex"), ...) {
    x <- mt_as_exposure(mt_data(deaths = d, population = p, by = by, ...), type = "midyear")
    v <- mt_validate(x)
    found <- if (nrow(v)) sprintf("mt_validate(): %s; ", paste(sprintf("%s (%s)", v$check, v$severity), collapse = ", "))
    fc <- mt_forecast(mt_fit(lee_carter(), x), h = 5, n_draws = 50, seed = 1)
    e0 <- life_expectancy(fc, at = 0)[period == max(period)]
    sprintf("%s%d strata; e0 in the last year %.1f-%.1f", found %||% "", nrow(mt_strata(x)), min(e0$q0.5),
            max(e0$q0.5))
  }
  clean <- sim_clean()
  small <- sim_clean(regions = c("North", "South", "Island"), pop_scale = c(2e5, 1e5, 5e3))
  share_above_1 <- function(fc) mean(apply(fc$rates, length(dim(fc$rates)), function(r) any(r > 1)))
  top_split <- function(d) rbind(copy(d)[, `:=`(age = "85-89", deaths = round(deaths * 0.6))],
                                 copy(d)[, `:=`(age = "90+", deaths = deaths - round(deaths * 0.6))])
  cases <- list(
    list("none", "clean data (reference)", function() pipeline(clean)),
    list("coding", "rows with age 'All ages' (table totals)", function() {
      t <- clean[, .(deaths = sum(deaths), population = sum(population), age = "All ages"), by = .(region, sex, year)]
      pipeline(rbind(clean, t, use.names = TRUE))
    }),
    list("coding", "rows with sex 'Persons' next to female and male", function() {
      t <- clean[, .(deaths = sum(deaths), population = sum(population), sex = "Persons"), by = .(region, year, age)]
      pipeline(rbind(clean, t, use.names = TRUE))
    }),
    list("coding", "sex 'female'/'male' in the deaths, 'F'/'M' in the populations", function() {
      pipeline(clean, copy(clean)[, sex := ifelse(sex == "female", "F", "M")])
    }),
    list("coding", "region 'North ' (trailing space) and 'south' (case) in the populations", function() {
      pipeline(clean, copy(clean)[region == "North", region := "North "][region == "South", region := "south"])
    }),
    list("coding", "sex missing (NA) in 20 death rows", function() {
      d <- copy(clean)
      rows <- with_seed(2, sample(nrow(d), 20))
      d[rows, sex := NA]
      pipeline(d)
    }),
    list("coding", "year and age stored as factors", function() {
      pipeline(copy(clean)[, `:=`(year = factor(year), age = factor(age))])
    }),
    list("coding", "years as financial years '2001/02'", function() {
      pipeline(copy(clean)[, year := sprintf("%d/%02d", year, (year + 1) %% 100)])
    }),
    list("coding", "single ages 0-90 where 90 means 90+, without age_open", function() {
      d <- CJ(sex = c("female", "male"), year = 2001:2020, age = 0:90)
      d[, population := round(1e5 * exp(-age / 40))]
      d[age == 90, population := round(sum(1e5 * exp(-(90:105) / 40)))]
      mu <- d$population * pmin(exp(-9.5 + 0.09 * d$age), 0.6)
      d[, deaths := with_seed(3, stats::rpois(length(mu), mu))]
      pipeline(d, by = "sex")
    }),
    list("coding", "populations as text with thousands separators ('12,345')", function() {
      pipeline(clean, copy(clean)[, population := format(population, big.mark = ",", trim = TRUE)])
    }),
    list("coding", "populations as text with decimal commas ('1234,5')", function() {
      pipeline(clean, copy(clean)[, population := sub("\\.", ",", sprintf("%.1f", population + 0.5))])
    }),
    list("structure", "a block of death rows appears twice", function() pipeline(rbind(clean, clean[year == 2010]))),
    list("structure", "year 2010 missing in all tables", function() pipeline(clean[year != 2010])),
    list("structure", "year 2010 missing for one region in the deaths", function() {
      pipeline(clean[!(year == 2010 & region == "South")], clean)
    }),
    list("structure", "populations missing for one region", function() pipeline(clean, clean[region != "South"])),
    list("structure", "deaths to 90+ (85-89, 90+), populations to 85+", function() {
      pipeline(rbind(clean[age != "85+"], top_split(clean[age == "85+"])), clean)
    }),
    list("structure", "deaths 85+ until 2010, 85-89 and 90+ from 2011", function() {
      pipeline(rbind(clean[!(age == "85+" & year > 2010)], top_split(clean[age == "85+" & year > 2010])), clean)
    }),
    list("structure", "both '85+' and '85-89', '90+' in the deaths of one year", function() {
      pipeline(rbind(clean, top_split(clean[age == "85+" & year == 2015])), clean)
    }),
    list("structure", "only 4 years of data", function() pipeline(clean[year >= 2017])),
    list("values", "population 0 in a cell with deaths", function() {
      pipeline(clean, copy(clean)[region == "South" & age == "85+" & year == 2005, population := 0])
    }),
    list("values", "deaths larger than the population in one cell", function() {
      pipeline(copy(clean)[region == "North" & age == "40-44" & year == 2005, deaths := population * 2])
    }),
    list("values", "a negative death count", function() {
      pipeline(copy(clean)[region == "South" & age == "50-54" & year == 2012, deaths := -3])
    }),
    list("values", "counts below 5 published as '<5'", function() {
      pipeline(copy(clean)[, deaths := ifelse(deaths < 5, "<5", as.character(deaths))])
    }),
    list("values", "counts below 5 published as '<5', declared with missing = '<5'", function() {
      pipeline(copy(clean)[, deaths := ifelse(deaths < 5, "<5", as.character(deaths))], missing = "<5")
    }),
    list("values", "suppressed cells ('.') without published totals, then mt_fill_suppressed()", function() {
      d <- copy(clean)[, deaths := ifelse(deaths < 3, ".", as.character(deaths))]
      x <- mt_as_exposure(mt_data(deaths = d, population = d, by = c("region", "sex")), type = "midyear")
      mt_fill_suppressed(x, mt_hierarchy(region = "all"))
      "filled"
    }),
    list("values", "an age group with zero deaths in every year", function() {
      pipeline(copy(clean)[region == "South" & age == "5-9", deaths := 0])
    }),
    list("values", "a mortality shock (+20%) in the last two years", function() {
      pipeline(copy(clean)[year >= 2019, deaths := round(deaths * 1.2)])
    }),
    list("small population", "a region of 5,000 people, Lee-Carter for each region and sex", function() {
      x <- mt_as_exposure(mt_data(deaths = small, population = small, by = c("region", "sex")), type = "midyear")
      fc <- mt_forecast(mt_fit(lee_carter(), mt_subset(x, region == "Island")), h = 5, n_draws = 100, seed = 1)
      sprintf("%.0f%% of its cells without deaths; draws with a rate above 1: %.0f%%",
              100 * mean(mt_obs(mt_subset(x, region == "Island"))$value == 0), 100 * share_above_1(fc))
    }),
    list("small population", "the same region with rwd(), li_lee(common_by = \"sex\") and tucker3()", function() {
      x <- mt_as_exposure(mt_data(deaths = small, population = small, by = c("region", "sex")), type = "midyear")
      out <- vapply(list(rwd = rwd(), li_lee = li_lee(common_by = "sex"), tucker3 = tucker3()), function(sp) {
        fc <- mt_forecast(mt_fit(sp, x), h = 5, n_draws = 100, seed = 1)
        isl <- grep("Island", dimnames(fc$rates)[[3]])
        sprintf("%.0f%%", 100 * mean(apply(fc$rates[, , isl, , drop = FALSE], 4, function(r) any(r > 1))))
      }, "")
      paste("draws with a rate above 1 in the small region:", paste(names(out), out, sep = " ", collapse = ", "))
    }),
    list("scale", "400 districts x 2 sexes x 18 ages x 20 years (288,000 rows)", function() {
      d <- sim_clean(regions = sprintf("D%03d", 1:400), pop_scale = rep(5e4, 400))
      t0 <- Sys.time()
      x <- mt_as_exposure(mt_data(deaths = d, population = d, by = c("region", "sex")), type = "midyear")
      t1 <- Sys.time()
      y <- mt_group_levels(x, "region", function(r) paste0("G", (as.integer(substr(r, 2, 4)) - 1) %/% 40 + 1))
      z <- mt_aggregate_strata(y, mt_hierarchy(region = "all", sex = "both"))
      t2 <- Sys.time()
      sprintf("reading %.1f s; grouping and totals %.1f s", as.numeric(t1 - t0, units = "secs"),
              as.numeric(t2 - t1, units = "secs"))
    }))
  log_step("defective data: %d scenarios", length(cases))
  res$defects <- rbindlist(lapply(cases, function(cs) run_case(cs[[1]], cs[[2]], cs[[3]])))
}

res$settings <- list(date = Sys.Date(), quick = quick,
                     versions = vapply(c("mortensor", "StMoMo", "demography", "forecast"),
                                       function(p) as.character(utils::packageVersion(p)), ""),
                     R = R.version.string)
saveRDS(res, out_file, compress = "xz")
log_step("saved %s", out_file)
