#!/usr/bin/env Rscript

library(devtools)
library(roxygen2)
library(testthat)

setwd('/ua/pliu/repe/pram/')

  document()
# install(quick=F, reload=T, build_vignettes=T, threads=4)
# test()
  check(cran=F)

# build_vignettes()

# test( filter = 'runPRAM' )
