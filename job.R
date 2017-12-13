#!/bin/env Rscript

library(devtools)
library(roxygen2)
library(testthat)

setwd('/ua/pliu/repe/pram/')
document()
install(quick=T, reload=F, threads=4)


#test( filter = 'buildModel' )
test( )
