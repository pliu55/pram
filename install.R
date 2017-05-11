#!/bin/env Rscript

library(devtools)
library(roxygen2)

setwd('/ua/pliu/repe/pram/')
document()
install(quick=T, threads=4)


library(testthat)
test()
