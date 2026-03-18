# Create synthetic example mortality dataset for the mortensor package
#
# Generates realistic mortality rates following a Gompertz-Makeham structure
# with sex differentials and deprivation gradients.
# Run this script from the package root: source("data-raw/create_example_data.R")
#
# Uses devtools::load_all() so the package does not need to be installed first.

devtools::load_all(".")

set.seed(42)

# --- Parameters ---
ages  <- seq(30, 90, by = 5)       # 13 five-year age groups (30-34 to 90+)
years <- 2001:2020                   # 20 years of data
sexes <- c("Male", "Female")
quintiles <- paste0("Q", 1:5)       # 5 deprivation quintiles (Q1 = least, Q5 = most)

n_a <- length(ages)
n_t <- length(years)
groups <- expand.grid(sex = sexes, quintile = quintiles,
                      stringsAsFactors = FALSE)
n_g <- nrow(groups)

# --- Gompertz-Makeham baseline ---
# log(mx) = alpha + beta * age + gamma * year
# Typical UK-like parameters
alpha_base <- -10.5   # intercept (log scale)
beta_age   <-  0.085  # age slope (Gompertz)
gamma_year <- -0.015  # period trend (mortality improvement)

# --- Build 3D arrays ---
mx_arr  <- array(NA_real_, dim = c(n_a, n_t, n_g))
pop_arr <- array(NA_real_, dim = c(n_a, n_t, n_g))

for (g in seq_len(n_g)) {
  sex <- groups$sex[g]
  quint <- groups$quintile[g]
  q_num <- as.integer(sub("Q", "", quint))

  # Sex differential: males have ~30% higher mortality (log scale shift)
  sex_shift <- ifelse(sex == "Male", 0.3, 0.0)

  # Deprivation gradient: ~8% increase per quintile step (log scale)
  dep_shift <- (q_num - 1) * 0.08

  for (t in seq_len(n_t)) {
    for (a in seq_len(n_a)) {
      log_mx <- alpha_base + beta_age * ages[a] +
                gamma_year * (years[t] - 2010) +
                sex_shift + dep_shift +
                rnorm(1, 0, 0.04)  # random noise
      mx_arr[a, t, g] <- exp(log_mx)
    }
    # Population: declining with age, slight growth over time
    pop_arr[, t, g] <- round(50000 * exp(-0.03 * (ages - 30)) *
                               (1 + 0.005 * (years[t] - 2001)) *
                               runif(n_a, 0.9, 1.1))
  }
}

# Add COVID-like shock in 2020 (last year)
covid_idx <- which(years == 2020)
mx_arr[, covid_idx, ] <- mx_arr[, covid_idx, ] * exp(rnorm(n_a * n_g, 0.15, 0.05))

# --- Create mort_data object ---
example_mortality <- mort_data(
  mx     = mx_arr,
  ages   = ages,
  years  = years,
  groups = groups,   # data.frame with sex, quintile columns
  pop    = pop_arr,
  cause  = "All causes"
)

# --- Save ---
usethis::use_data(example_mortality, overwrite = TRUE)

cat("Created example_mortality dataset:\n")
cat(sprintf("  Ages: %d-%d+ (%d groups)\n", min(ages), max(ages), n_a))
cat(sprintf("  Years: %d-%d (%d years)\n", min(years), max(years), n_t))
cat(sprintf("  Groups: %d (%s x %s)\n", n_g,
            paste(unique(groups$sex), collapse = "/"),
            paste(unique(groups$quintile), collapse = "/")))
