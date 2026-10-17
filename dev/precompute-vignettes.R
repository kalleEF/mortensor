# Precompute vignettes whose chunks take several minutes (MCMC, backtests). Run from the
# package root; the vignettes use the development version of the package:
#   Rscript dev/precompute-vignettes.R
# Each knit takes 5-20 minutes. The figures are written next to the vignettes.
suppressMessages(devtools::load_all(".", quiet = TRUE))
old <- setwd("vignettes")
on.exit(setwd(old))
for (v in c("bayesian-apc", "model-choice", "inequality-england", "laender-germany")) {
  knitr::knit(paste0(v, ".Rmd.orig"), output = paste0(v, ".Rmd"))
}
