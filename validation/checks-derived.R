# Proposed validation checks for lifespan variation (lifespan_variation())
# and excess mortality (mt_excess()), in the format of checks.R: each check
# returns a data frame with columns `metric` and `value`. To use them, source
# this file after checks.R in run_validation.R
# (source(file.path(vdir, "checks-derived.R"), local = TRUE)) or move the
# checks into checks.R, and add the thresholds below to thresholds.csv.
#
# Proposed lines for thresholds.csv (check,metric,comparator,threshold,note):
#
# lifespan_reference,max_rel_diff_lifeineq,<=,1e-8,"vs LifeIneq 00.04.09; 3e-10 at 2026-10 (10-digit constants)"
# lifespan_reference,max_iqr_diff_lifeineq,<=,0.03,"single ages; spline in LifeIneq; 0.017 at 2026-10"
# lifespan_reference,mc_sd_diff,<=,0.08,"1e6 sampled ages at death; 0.036 at 2026-10 (Monte Carlo error)"
# lifespan_reference,mc_gini_diff,<=,0.0005,"0.00012 at 2026-10"
# lifespan_reference,mc_iqr_diff,<=,0.05,"0.015 at 2026-10"
# excess_simulation,cov95_null,>=,0.9,"interval covers 0 without a shock; 0.975 at 2026-10 (20 replicates)"
# excess_simulation,cov95_null,<=,0.995,""
# excess_simulation,cov95_shock,>=,0.9,"interval covers the injected deaths; 0.983 at 2026-10"
# excess_simulation,rel_bias_shock,<=,0.03,"mean relative error of the median excess; 0.007 at 2026-10"
# excess_simulation,seconds,<=,120,"runtime budget; 28 s at 2026-10"
# excess_france,excess_2003,>=,8000,"15,600 at 2026-10 (published for August 2003: about 15,000)"
# excess_france,excess_2003,<=,25000,""
# excess_france,excess_2003_lower,>=,0,"lower 95% limit; 6,200 at 2026-10"
# excess_france,excess_2004,<=,0,"deficit after the heatwave; -29,700 at 2026-10"
# excess_france,p_score_1918_female,>=,20,"43% at 2026-10"
# excess_france,p_score_1918_male,>=,50,"99% at 2026-10 (male excess 1914-1917 points to war deaths)"

if (!exists("checks")) checks <- list()

# France (HMD via demography::fr.mort) life tables: single ages 0-109 and
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
       female_1950_5y = grouped("female", 1950, seq(0, 100, 5)))
}

checks$lifespan_reference <- list(
  suite = "quick",
  description = paste("lifespan_variation() on France life tables (fr.mort: women 2000 and men 1987 in single ages to",
                      "110+, women 1950 in 5-year groups to 100+) at ages 0 and 65: largest relative difference",
                      "from LifeIneq 00.04.09 (e-dagger, entropy, Gini, standard deviation; reference values as in",
                      "tests/testthat/test-lifespan.R), largest difference of the interquartile range (single ages),",
                      "and the difference of the standard deviation, Gini coefficient and interquartile range from",
                      "1e6 ages at death sampled from the life table of women 2000 (uniform within intervals,",
                      "exponential in the open interval)."),
  run = function() {
    ref <- list(
      female_2000_single = rbind(edagger = c(9.663024709, 7.008585465), sd = c(14.18425557, 8.211834891),
                                 entropy = c(0.116674102, 0.3298611941), gini = c(0.08453618794, 0.0537459495),
                                 iqr = c(13.41011018, 11.22595719)),
      male_1987_single = rbind(edagger = c(12.4459491, 7.202768889), sd = c(17.07315081, 7.945863166),
                               entropy = c(0.1727633064, 0.480640979), gini = c(0.1237094046, 0.0568842088),
                               iqr = c(19.32465572, 12.05458686)),
      female_1950_5y = rbind(edagger = c(14.64681576, 6.98491648), sd = c(22.05606492, 7.485704412),
                             entropy = c(0.2120541073, 0.4770597137), gini = c(0.1577954255, 0.05295861841),
                             iqr = c(19.20655623, 11.21954679)))
    lts <- fr_life_tables()
    rel <- 0
    iqr <- 0
    for (nm in names(ref)) {
      lt <- lts[[nm]]
      for (k in 1:2) {
        for (ms in c("edagger", "sd", "entropy", "gini")) {
          v <- lifespan_variation(lt$m, ms, lt$age, at = c(0, 65)[k])
          rel <- max(rel, abs(v / ref[[nm]][ms, k] - 1))
        }
        if (grepl("single", nm)) {
          iqr <- max(iqr, abs(lifespan_variation(lt$m, "iqr", lt$age, at = c(0, 65)[k]) - ref[[nm]]["iqr", k]))
        }
      }
    }
    lt <- lts$female_2000_single
    r <- mortensor:::lt_columns(matrix(lt$m, ncol = 1), lt$age)
    na <- length(lt$m)
    n <- 1e6
    set.seed(1)
    j <- sample.int(na, n, replace = TRUE, prob = r$dx[, 1])
    y <- lt$age$lower[j] + ifelse(j < na, stats::runif(n), stats::rexp(n, r$dx[na, 1] / r$Lx[na, 1]))
    ys <- sort(y)
    gini <- 2 * sum(ys * seq_len(n)) / (n * sum(ys)) - (n + 1) / n
    data.frame(metric = c("max_rel_diff_lifeineq", "max_iqr_diff_lifeineq", "mc_sd_diff", "mc_gini_diff",
                          "mc_iqr_diff"),
               value = c(rel, iqr, abs(lifespan_variation(lt$m, "sd", lt$age) - stats::sd(y)),
                         abs(lifespan_variation(lt$m, "gini", lt$age) - gini),
                         abs(lifespan_variation(lt$m, "iqr", lt$age) - diff(stats::quantile(y, c(0.25, 0.75))))))
  })

# One replicate of the excess simulation: mt_example() with deaths at ages
# 70+ in 2016 and 2017 raised by 10%, Lee-Carter baseline 1990-2014, excess
# by sex and age group (30-69, 70+) and in total for 2015-2019, with the
# injected deaths as truth
excess_replicate <- function(r, n_draws = 500) {
  x <- mt_example(years = 1990:2019, deaths = "single", seed = r)
  d <- x$obs$deaths
  hit <- d$period_lo %in% c(2016, 2017) & d$age_lo >= 70
  add <- round(1.1 * d$value[hit]) - d$value[hit]
  inj <- data.table::data.table(sex = x$strata$sex[d$stratum[hit]], period = d$period_lo[hit], add = add)
  x$obs$deaths$value[hit] <- d$value[hit] + add
  ex <- mt_excess(lee_carter(), x, baseline = c(1990, 2014), ages = c(30, 70), n_draws = n_draws, seed = r)
  tot <- mt_excess(lee_carter(), x, baseline = c(1990, 2014), by = character(), n_draws = n_draws, seed = r)
  tr <- inj[, list(true = sum(add)), by = c("sex", "period")]
  ex <- merge(ex, tr, by = c("sex", "period"), all.x = TRUE)
  ex$true[ex$age == "30-69" | is.na(ex$true)] <- 0
  tt <- inj[, list(true = sum(add)), by = "period"]
  tot <- merge(tot, tt, by = "period", all.x = TRUE)
  tot$true[is.na(tot$true)] <- 0
  rbind(data.frame(level = "sex_age", true = ex$true, med = ex$excess_q0.5, lo = ex$excess_q0.025,
                   hi = ex$excess_q0.975),
        data.frame(level = "total", true = tot$true, med = tot$excess_q0.5, lo = tot$excess_q0.025,
                   hi = tot$excess_q0.975))
}

checks$excess_simulation <- list(
  suite = "quick",
  description = paste("mt_excess() on simulated data with known truth (mt_example, single ages, 1990-2019, 20",
                      "replicates): deaths at ages 70+ in 2016 and 2017 raised by 10%, Lee-Carter baseline",
                      "1990-2014, excess 2015-2019 by sex and age group (30-69, 70+) and in total; 95% coverage of",
                      "zero where there is no shock and of the injected deaths where there is one, and the mean",
                      "relative error of the median excess."),
  run = function() {
    t0 <- proc.time()[["elapsed"]]
    res <- do.call(rbind, lapply(1:20, excess_replicate))
    secs <- proc.time()[["elapsed"]] - t0
    covered <- res$lo <= res$true & res$true <= res$hi
    null <- res$true == 0
    shock <- res$true > 0
    data.frame(metric = c("cov95_null", "cov95_shock", "cov95_total_null", "rel_bias_shock", "seconds"),
               value = c(mean(covered[null & res$level == "sex_age"]), mean(covered[shock]),
                         mean(covered[null & res$level == "total"]),
                         abs(mean((res$med[shock] - res$true[shock]) / res$true[shock])), secs))
  })

checks$excess_france <- list(
  suite = "quick",
  description = paste("mt_excess() on France (HMD via demography::fr.mort), ages 0-99 and 100+: excess deaths in",
                      "2003 (August heatwave) and 2004 against a Lee-Carter baseline 1990-2002 (both sexes, median",
                      "and 95% interval), the loss of life expectancy at birth of women in 2003, and the P-scores of",
                      "1918 against a 1900-1913 baseline (ages 0-94 and 95+)."),
  run = function() {
    fr <- as_mt_data(mortensor:::pkg_data("fr.mort", "demography"))
    x <- mt_regroup(mt_subset(fr, periods = c(1990, 2004)), age = 0:100)
    tot <- mt_excess(lee_carter(), x, baseline = c(1990, 2002), by = character(), n_draws = 1000, seed = 1)
    le <- mt_excess(lee_carter(), x, baseline = c(1990, 2002), periods = 2003, quantity = "life_expectancy",
                    n_draws = 1000, seed = 1)
    y <- mt_regroup(mt_subset(fr, periods = c(1900, 1918)), age = 0:95)
    e18 <- mt_excess(lee_carter(), y, baseline = c(1900, 1913), periods = 1918, n_draws = 1000, seed = 1)
    data.frame(metric = c("excess_2003", "excess_2003_lower", "excess_2003_upper", "excess_2004",
                          "le_excess_2003_female", "p_score_1918_female", "p_score_1918_male"),
               value = c(tot$excess_q0.5[tot$period == 2003], tot$excess_q0.025[tot$period == 2003],
                         tot$excess_q0.975[tot$period == 2003], tot$excess_q0.5[tot$period == 2004],
                         le$excess_q0.5[le$sex == "female"], e18$p_score_q0.5[e18$sex == "female"],
                         e18$p_score_q0.5[e18$sex == "male"]))
  })
