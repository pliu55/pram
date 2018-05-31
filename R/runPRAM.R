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
#' @importFrom  tools  file_path_sans_ext
#'
#' @export
#'
setGeneric( 'runPRAM',
           function(in_gtf, in_bamv, out_gtf, in_bedv,
                    training_tpms, training_gtf) standardGeneric('runPRAM'))

setMethod( 'runPRAM',
c('character', 'vector', 'character', 'vector', 'vector', 'character'),
function(in_gtf, in_bamv, out_gtf, in_bedv, training_tpms, training_gtf) {

    chromgrs = getMaxChromGRangesFromBams(in_bamv)

    iggrs = defIgRanges(in_gtf, chromgrs)

    tmpdir = paste0(tempdir(), '/runPRAM/')
    bamdt = data.table(finbam = in_bamv)
    bamdt[, foutbam := paste0(tmpdir, file_path_sans_ext(basename(finbam)),
                              '.ig.bam')]
    apply(bamdt, 1, function(x) prepIgBam(x[['finbam']], iggrs, x[['foutbam']]))

    fgtf_all_mdl = paste0(tmpdir, 'all_mdl.gtf')
    buildModel(fbamv_ig, fgtf_all_mdl)

    selModel(fgtf_all_mdl, out_gtf)

    if ( ( ! missing(in_bedv)) & (! missing(training_gtf)) &
         ( ! missing(training_tpms)) ) {
        screenModel(in_bedv, training_tpms, training_gtf, out_gtf, out_gtf)
    }
})



#' @importFrom  Rsamtools      idxstatsBam
#' @importFrom  GenomicRanges  makeGRangesFromDataFrame
#'
setGeneric( 'getMaxChromGRangesFromBams',
            function(fbams) standardGeneric('getMaxChromGRangesFromBams'))

setMethod('getMaxChromGRangesFromBams', 'vector',
function(fbams) {
    dt = rbindlist(lapply(fbams,
                          function(x) data.table(idxstatsBam(x))[, fbam := x]))
    maxdt = dt[, list(end = max(seqlength)), by=seqnames]
    maxdt[, start := 1]
    grs = makeGRangesFromDataFrame(maxdt, keep.extra.columns=F)

    return(grs)
})
