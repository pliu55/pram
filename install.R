#!/usr/bin/env Rscript

library(devtools)
library(roxygen2)
library(testthat)

document()
install(quick=F, reload=T, build_vignettes=T, threads=4)
test()
