#' @export
ModelGTF = setClass( 'ModelGTF', contains = 'GTF' )

setMethod('initialize', 'ModelGTF',
    function(.Object, fgtf, infokeys, origin, model_method) {
        .Object = callNextMethod(.Object, fgtf, infokeys, origin=origin)
        renameGnTrID(.Object@exon, model_method)
        return(.Object)
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
