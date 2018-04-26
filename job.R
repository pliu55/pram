#!/bin/env Rscript

library(devtools)
library(roxygen2)
library(testthat)

setwd('/ua/pliu/repe/pram/')
document()
install(quick=T, reload=F, threads=4)

## deal with prepIgBam's selAlnInGRanges: readGAlignment, read in chunk
## debug by ~/repe/gata/70_prepIgBam/

# test( filter = 'screenModel' )
# test( filter = 'prepIgBam' )
# test( filter = 'defIgRanges' )
# test( filter = 'buildModel' )
# test( filter = 'evalModel' )
# test( filter = 'selModel' )
 test()
