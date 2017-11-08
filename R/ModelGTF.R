ModelGTF = setClass(
    'ModelGTF',
    slots = list( model_method = 'character' ),
    contains = 'GTF'
)


setMethod(
    'initialize',
    signature('ModelGTF'),
    function(.Object) {
        .Object = callNextMethod(.Object)
        .Object@model_method = character()
        return(.Object)
    }
)


#' @export
#'
setMethod(
    'initFromGTFFile',
    signature('ModelGTF', 'character', 'vector'),
    function(obj, fgtf, infokeys, ...) {
        obj = callNextMethod(obj, fgtf, infokeys, ...)
        ecl = list(...)
        if ( ! is.null(ecl$model_method) ) {
            renameGnTrID(obj@grangedt, ecl$model_method)
        }

        return(obj)
    }
)


renameGnTrID <- function(exondt, runid) {
    exondt[, ori_lab := ifelse(strand == '+', 'plus',
                               ifelse(strand == '-', 'minus', NA) )]
    exondt[, ign := .GRP, by=gene_id]
    exondt[, itr := .GRP, by=transcript_id]
    exondt[, gene_id := paste0(runid, '_', ori_lab, '_', chrom, '.', ign)]
    exondt[, transcript_id := paste0(gene_id, '.', itr)]
    exondt[, c('ori_lab', 'ign', 'itr') := NULL]
}
