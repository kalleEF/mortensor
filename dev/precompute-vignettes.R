# Precompute vignettes whose chunks take several minutes (MCMC, backtests). Run from the
# package root after installing the current version:
#   R CMD INSTALL . && Rscript dev/precompute-vignettes.R
old <- setwd("vignettes")
on.exit(setwd(old))
knitr::knit("bayesian-apc.Rmd.orig", output = "bayesian-apc.Rmd")
knitr::knit("regions.Rmd.orig", output = "regions.Rmd")
knitr::knit("deprivation.Rmd.orig", output = "deprivation.Rmd")
