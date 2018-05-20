#' Predict intergenic transcript models from RNA-seq and screen them by ChIP-seq
#'
#' @inheritParams defIgRanges
##                to inherit `in_gtf` and `genome`
#'
#' @inheritParams buildModel
##                to inherit `in_bamv`, `out_gtf`, and `cufflinks`
#'
#' @inheritParams screenModel
##                to inherit `in_bedv`, `training_gtf`, and `training_tpms`
#'
#' @return  NULL
#'
#' @export
#'
setGeneric( 'runPRAM',
           function(in_gtf, in_bamv, out_gtf, cufflinks, genome, in_bedv,
                    training_gtf, training_tpms) standardGeneric('runPRAM'))

setMethod( 'runPRAM',
c('character', 'vector', 'character', 'character', 'character', 'vector',
  'character', 'character'),
function(in_gtf, in_bamv, out_gtf, cufflinks='', genome=NULL, in_bedv,
         training_gtf, training_tpms) {

})
