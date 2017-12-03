#' Evaluate transcript model
#'
#' Evaluate trascript model's precision and recall on exon nucleotide,
#' splice junctions, and transcript's splice pattern
#'
#' @param mdltr a Transcript object for models
#' @param tgttr a Transcript object for targets
#'
#' @import data.table
#'
#' @return a data table of precision, recall, number of true positive,
#'         false negative, false positive for all three evaluated features
#'
#' @export
#'
evalModel <- function(mdltr, tgttr) {
    mdlexondt = getExon(mdltr)
    tgtexondt = getExon(tgttr)

    exonnucdt = evalMdlExonNuc(mdlexondt, tgtexondt)

    mdl_ol_tgtdt = findMdlTrOLTgtTr(mdlexondt, tgtexondt)

    mdljncdt = getJnc(mdltr)
    tgtjncdt = getJnc(tgttr)

    setnames(mdljncdt, 'trid', 'mdlid')
    mdljncdt[, mdljncid := .I]
    tgtjncdt[, tgtjncid := .I]

    tgtjnc_in_mdldt = findTargetJncInModel(tgtjncdt, mdljncdt)
    mdljnc_in_tgtdt = findModelJncInTarget(mdljncdt, tgtjncdt)

    jncdt = evalMdlJnc(tgtjnc_in_mdldt, mdljnc_in_tgtdt, mdl_ol_tgtdt, tgtjncdt)
    outdt = rbind(exonnucdt, jncdt)

    return(outdt)
}


#' Evaluate transcript model's splice junction
#'
#' @return a data table of precision, recall, number of true positive,
#'         false negative, false positive at each junction and all junctions in
#'         a transcript
#'
evalMdlJnc <- function(tgtjnc_in_mdldt, mdljnc_in_tgtdt, mdl_ol_tgtdt,
                       tgtjncdt) {
    tpfndt = calTpFn4Jnc(tgtjnc_in_mdldt, tgtjncdt)
    fpdt = calFp4Jnc(mdljnc_in_tgtdt, mdl_ol_tgtdt, tgtjncdt)

    dt = merge(tpfndt, fpdt, by='feat', all=T)
    dt[, `:=`( precision = ntp/(ntp + nfp),
               recall    = ntp/(ntp + nfn) )]

    return(dt)
}

#' find if model's splice junction exists in a target transcript
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors mcols
#'
findModelJncInTarget <- function(mdljncdt, tgtjncdt) {
    mdlgrs = makeGRangesFromDataFrame(mdljncdt, keep.extra.columns=T)
    tgtgrs = makeGRangesFromDataFrame(tgtjncdt, keep.extra.columns=T)

    ol = findOverlaps(mdlgrs, tgtgrs, type='equal', ignore.strand=F)
    oldt = data.table(as.data.frame(ol))
    oldt[, `:=`( mdljncid = mcols(mdlgrs)$mdljncid[queryHits],
                 trid     = mcols(tgtgrs)$trid[subjectHits] )]
    dt = merge(mdljncdt, oldt[, list(mdljncid, trid)], by='mdljncid', all=T)
    dt[, mdljncid := NULL]
    return(dt)
}


#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#'
findTargetJncInModel <- function(tgtjncdt, mdljncdt) {
    mdlgrs = makeGRangesFromDataFrame(mdljncdt, keep.extra.columns=T)
    tgtgrs = makeGRangesFromDataFrame(tgtjncdt, keep.extra.columns=T)

    ol = findOverlaps(tgtgrs, mdlgrs, type='equal', ignore.strand=F)
    oldt = data.table(as.data.frame(ol))
    oldt[, `:=`( tgtjncid = mcols(tgtgrs)$tgtjncid[queryHits],
                 mdlid    = mcols(mdlgrs)$mdlid[subjectHits] )]
    dt = merge(tgtjncdt, oldt[, list(tgtjncid, mdlid)], by='tgtjncid', all=T)
    dt[, tgtjncid := NULL]
    return(dt)
}


#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges findOverlaps
#'
findMdlTrOLTgtTr <- function(mdlexondt, tgtexondt) {
    mdldt = getTrFromExon(mdlexondt, id_col_name='trid')
    trdt  = getTrFromExon(tgtexondt, id_col_name='trid')

    setnames(mdldt, 'trid', 'mdlid')

    mdlgrs = makeGRangesFromDataFrame(mdldt, keep.extra.columns=T)
    trgrs  = makeGRangesFromDataFrame(trdt, keep.extra.columns=T)

    ol = findOverlaps(mdlgrs, trgrs, type='any', ignore.strand=F)
    oldt = data.table(as.data.frame(ol))
    oldt[, `:=`( mdlid = mcols(mdlgrs)$mdlid[queryHits],
                 trid  = mcols(trgrs)$trid[subjectHits] )]
    outdt = merge(mdldt[, list(mdlid, nexon)], oldt[, list(mdlid, trid)],
                  by='mdlid', all.x=T)

    return(outdt)
}


calTpFn4Jnc <- function(tgtjnc_in_mdldt, tgtjncdt) {
    tgtdt = unique(tgtjncdt, by=c('trid', 'njnc'))
    dt = subset(tgtjnc_in_mdldt, (trid %in% tgtdt[, trid]) & (! is.na(mdlid)))

    ## at individual junction level
    ## - TP: number of target junctions exists in any model
    ## - FN: number of target junctions does not exist in any model
    indidt = dt[, list(nprd = length(unique(ijnc))), by=trid]
    indi_prddt = merge(tgtdt, indidt, by='trid', all.x=T)
    indi_prddt[, nprd := ifelse(is.na(nprd), 0, nprd)]
    indi_prddt[, n_not_prd := njnc - nprd]
    indi_tp = sum(indi_prddt[, nprd])
    indi_fn = sum(indi_prddt[, n_not_prd])

    ## at transcript level
    ## - TP: number of targets w/ all junctions exist in a model
    ## - FN: number of targets w/ a junction does not exist in any model
    mdl_njncdt = dt[, list(nmdljnc = length(unique(ijnc))),
                    by=list(trid, njnc, mdlid)]
    trdt = subset(mdl_njncdt, njnc == nmdljnc)
    tr_prddt = data.table(trid = tgtdt[, trid])
    tr_prddt[, nprd := ifelse(trid %in% trdt[, trid], 1, 0)]
    tr_prddt[, n_not_prd := ifelse(nprd == 1, 0, 1)]
    tr_tp = sum(tr_prddt[, nprd])
    tr_fn = sum(tr_prddt[, n_not_prd])

    outdt = data.table( feat = c( 'indi_jnc', 'tr_jnc' ),
                        ntp  = c( indi_tp, tr_tp ),
                        nfn  = c( indi_fn, tr_fn ) )

    return(outdt)
}


calFp4Jnc <- function(mdljnc_in_tgtdt, mdl_ol_tgtdt, tgtjncdt) {
    tgtdt = unique(tgtjncdt, by=c('trid', 'njnc'))

    indt = subset(mdljnc_in_tgtdt, trid %in% tgtdt[, trid])
    oldt = subset(mdl_ol_tgtdt, trid %in% tgtdt[, trid])

    ## at individual junction level
    ## FP: number of model junctions not existing in model's overlapping target
    indidt = indt[, list(nprd = length(unique(ijnc))), by=mdlid]
    indi_prddt = merge(oldt, indidt, by='mdlid', all.x=T)
    indi_prddt[, nprd := ifelse(is.na(nprd), 0, nprd)]
    indi_prddt[, n_not_prd := ifelse(nexon == 1, 0, nexon - 1 - nprd)]
    if ( nrow(subset(indi_prddt, n_not_prd < 0)) > 0 ) {
        cat('error: n_not_prd < 0:', in_method, "\n")
    }
    indi_fp = sum(indi_prddt[, n_not_prd])


    ## at transcript level
    ## FP: number of models w/ a junction not existing in its overlapping target
    ##     or a single-exon model
    tr_fp = nrow(subset(indi_prddt, (nexon == 1) | (n_not_prd > 0) ))

    outdt = data.table( feat = c( 'indi_jnc', 'tr_jnc' ),
                        nfp  = c( indi_fp, tr_fp ) )

    return(outdt)
}


#' Evaluate transcript model's exon nucleotide
#'
#' @inheritParams evalModel
#'
#' @return a data.table of precision, recall, number of true positive,
#          false negative, and false positive
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges intersect
#' @importFrom GenomicRanges setdiff
#' @importFrom IRanges width
#'
evalMdlExonNuc <- function(mdldt, tgtdt) {
    mdlgrs = makeGRangesFromDataFrame(mdldt, keep.extra.columns=T)
    tgtgrs = makeGRangesFromDataFrame(tgtdt, keep.extra.columns=T)

    olgrs = intersect(mdlgrs, tgtgrs, ignore.strand=F)

    # - TP: number of target nucleotides exists in any model
    # - FN: number of target nucleotides does not exist in any model
    ntp = sum(width(olgrs))
    nfn = sum(width(tgtgrs)) - ntp

    # - FP: number of model nuc not existing in model's OL target
    # - every model overlaps with a target
    not_ol_mdlgrs = setdiff(mdlgrs, olgrs)
    nfp = sum(width(not_ol_mdlgrs))

    precision = ntp/(ntp + nfp)
    recall    = ntp/(ntp + nfn)

    outdt = data.table( feat      = 'exon_nuc',
                        ntp       = ntp,
                        nfn       = nfn,
                        nfp       = nfp,
                        precision = precision,
                        recall    = recall )
    return(outdt)
}
