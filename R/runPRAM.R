#' @title Predict intergenic transcript models from RNA-seq
## and screen them by ChIP-seq
#'
#' @inheritParams defIgRanges
##                to inherit `in_gtf`
#'
#' @inheritParams buildModel
##                to inherit `in_bamv`, `out_gtf`
#'
## @inheritParams screenModel
##                to inherit `in_bedv`, `training_gtf`, and `training_tpms`
#'
#' @return  None
#'
#' @importFrom  tools  file_path_sans_ext
#'
#' @export
#'
#' @examples
#'
#' in_gtf = system.file('extdata/demo/in.gtf', package='pram')
#'
#' in_bamv = c( system.file('extdata/demo/SZP.bam', package='pram'),
#'              system.file('extdata/demo/TLC.bam', package='pram') )
#'
#' pred_out_gtf = tempfile(fileext='.gtf')
#'
#' \donttest{
#' runPRAM(in_gtf, in_bamv, pred_out_gtf)
#' }
#'
#'
## in_bedv = c( system.file('extdata/demo/H3K79me2.bed.gz', package='pram'),
##              system.file('extdata/demo/POLR2.bed.gz',    package='pram') )
##
## training_tpms = c( system.file('extdata/demo/AED1.isoforms.results',
##                                package='pram'),
##                    system.file('extdata/demo/AED2.isoforms.results',
##                                package='pram') )
##
## training_gtf = system.file('extdata/demo/training.gtf', package='pram')
##
## screen_out_gtf = tempfile(fileext='.gtf')
##
## \donttest{
## runPRAM(in_gtf, in_bamv, screen_out_gtf, in_bedv, training_tpms,
##         training_gtf)
## }
##
##
#runPRAM <- function(in_gtf, in_bamv, out_gtf, in_bedv, training_tpms,
#                    training_gtf) {
runPRAM <- function(in_gtf, in_bamv, out_gtf) {
    foutbam = finbam = NULL
    chromgrs = getMaxChromGRangesFromBams(in_bamv)

    iggrs = defIgRanges(in_gtf, chromgrs)

    bamdt = data.table(finbam = in_bamv)
    bamdt[, foutbam := tempfile(pattern=file_path_sans_ext(basename(finbam)),
                                fileext='.ig.bam')]
    apply(bamdt, 1, function(x) prepIgBam(x[['finbam']], iggrs, x[['foutbam']]))

    fgtf_all_mdl = tempfile(pattern='pram_all_mdl.', fileext='.gtf')
    fgtf_sel_mdl = tempfile(pattern='pram_sel_mdl.', fileext='.gtf')

    buildModel(bamdt$foutbam, fgtf_all_mdl)

    selModel(fgtf_all_mdl, fgtf_sel_mdl, min_n_exon=2, min_tr_len=200,
             info_keys = c('transcript_id'))

    file.copy(fgtf_sel_mdl, out_gtf, overwrite=TRUE)

  # if ( missing(in_bedv) | missing(training_tpms) | missing(training_gtf) ) {
  #     file.copy(fgtf_sel_mdl, out_gtf, overwrite=T)
  # } else {
  #     screenModel(in_bedv, training_tpms, training_gtf, fgtf_sel_mdl, out_gtf)
  # }
}


#' @importFrom  GenomicRanges  makeGRangesFromDataFrame
#'
setGeneric( 'getMaxChromGRangesFromBams',
            function(fbams) standardGeneric('getMaxChromGRangesFromBams'))

setMethod('getMaxChromGRangesFromBams', 'vector',
function(fbams) {
    seqlength = NULL
    dt = rbindlist(lapply(fbams, getIdxStatsFromBam))
    maxdt = dt[, list(end = max(seqlength)), by=seqnames]
    maxdt[, start := 1]
    grs = makeGRangesFromDataFrame(maxdt, keep.extra.columns=FALSE)

    return(grs)
})


#' @importFrom  Rsamtools   idxstatsBam  indexBam
#' @importFrom  data.table  data.table
#'
setGeneric( 'getIdxStatsFromBam',
            function(fbam) standardGeneric('getIdxStatsFromBam'))

setMethod('getIdxStatsFromBam', 'character',
function(fbam) {
    fbai = paste0(fbam, '.bai')
    if ( ! file.exists(fbai) ) indexBam(fbam)
    dt = data.table(idxstatsBam(fbam))
    dt[, fbam := fbam]

    return(dt)
})
