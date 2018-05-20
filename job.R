#!/usr/bin/env Rscript

library(devtools)
library(roxygen2)
library(testthat)
# library(BiocStyle)

setwd('/ua/pliu/repe/pram/')
# build_vignettes()
  document()
  install(quick=T, reload=F, build_vignettes=F, threads=4)

## deal with prepIgBam's selAlnInGRanges: readGAlignment, read in chunk
## debug by ~/repe/gata/70_prepIgBam/

# test( filter = 'screenModel' )
# test( filter = 'prepIgBam' )
# test( filter = 'defIgRanges' )
# test( filter = 'buildModel' )
# test( filter = 'evalModel' )
# test( filter = 'selModel' )
# test( filter = 'GTF' )
  test()
