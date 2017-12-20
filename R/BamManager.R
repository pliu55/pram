#' @import methods
#'
BamManager = setClass('BamManager',
                      slots = list( fparentbams = 'vector',
                                    foutbam     = 'character',
                                    chrom       = 'character',
                                    ori         = 'character',
                                    strand      = 'character',
                                    grs         = 'GRanges' )
)


setGeneric('fparentbams', function(x) standardGeneric('fparentbams'))
setGeneric('foutbam',     function(x) standardGeneric('foutbam'))
setGeneric('chrom',       function(x) standardGeneric('chrom'))
setGeneric('ori',         function(x) standardGeneric('ori'))
## strand was defined in BioGenerics, cannot set here, just setMethod is fine
#setGeneric('strand',      function(x) standardGeneric('strand'))
setGeneric('grs',         function(x) standardGeneric('grs'))

setGeneric('fparentbams<-', function(x, value) standardGeneric('fparentbams<-'))
setGeneric('foutbam<-',     function(x, value) standardGeneric('foutbam<-'))
setGeneric('chrom<-',       function(x, value) standardGeneric('chrom<-'))
setGeneric('ori<-',         function(x, value) standardGeneric('ori<-'))
setGeneric('strand<-',      function(x, value) standardGeneric('strand<-'))
setGeneric('grs<-',         function(x, value) standardGeneric('grs<-'))

setMethod('fparentbams', 'BamManager', function(x) x@fparentbams)
setMethod('foutbam',     'BamManager', function(x) x@foutbam)
setMethod('chrom',       'BamManager', function(x) x@chrom)
setMethod('ori',         'BamManager', function(x) x@ori)
setMethod('strand',      'BamManager', function(x) x@strand)
setMethod('grs',         'BamManager', function(x) x@grs)

setReplaceMethod('fparentbams', 'BamManager',
                 function(x, value) {x@fparentbams=value; x})
setReplaceMethod('foutbam', 'BamManager',
                 function(x, value) {x@foutbam=value; x})
setReplaceMethod('chrom', 'BamManager', function(x, value) {x@chrom=value; x})
setReplaceMethod('ori',   'BamManager', function(x, value) {x@ori=value; x})
setReplaceMethod('strand', 'BamManager', function(x, value) {x@strand=value; x})
setReplaceMethod('grs', 'BamManager', function(x, value) {x@grs=value; x})


setMethod('show', 'BamManager',
    function(object) {
        cat('fparentbams:', fparentbams(object), "\n")
        cat('foutbam:',     foutbam(object),     "\n")
        cat('chrom:',       chrom(object),       "\n")
        cat('ori:',         ori(object),         "\n")
        cat('strand:',      strand(object),      "\n")
        print(grs(object))
    }
)


#' @importFrom  BiocGenerics  strand
#' @importFrom  GenomeInfoDb  seqnames
#' @importFrom  tools         file_path_sans_ext
#'
setMethod('initialize', 'BamManager',
    function(.Object, fparentbams, foutbam, chrom=NULL, ori=NULL, allgrs=NULL){
        fparentbams(.Object) = c(fparentbams)
        foutbam(.Object)     = foutbam
        chrom(.Object)       = chrom
        ori(.Object)         = ori
        strand(.Object)      = convertOri2Strand(ori)
        if ( ! is.null(allgrs) ) {
            grs(.Object) = allgrs[ seqnames(allgrs) == chrom &
                                   strand(allgrs) == ori ]
        }

        return(.Object)
    }
)


#' @importFrom  data.table  setkey data.table
#' @importFrom  Rsamtools   indexBam filterBam BamFile ScanBamParam
#' @importFrom  IRanges     filterRules
#'
setGeneric('filterParentBamByGRanges',
           function(bam, prm) standardGeneric('filterParentBamByGRanges'))
setMethod('filterParentBamByGRanges', c('BamManager', 'Param'),
    function(bam, prm) {
        finbam = fparentbams(bam)
        if ( length(finbam) > 1 ) {
            msg = paste0('BamManager::filterParentBamByGRanges: ',
                         "length(fparentbams) > 1\n")
            stop(msg)
        }
        grs = grs(bam)

        finbam_index = paste0(finbam, '.bai')
        if ( ! file.exists(finbam_index) ) indexBam(finbam)

        mate2flag = fr1ststrand2mate2flag(prm)[[strand(bam)]]
        flag_1stmate = mate2flag[['1stmate']]
        flag_2ndmate = mate2flag[['2ndmate']]

        alns_1st = selAlnInGRanges(grs, finbam, flag_1stmate)
        alns_2nd = selAlnInGRanges(grs, finbam, flag_2ndmate)
        alns = c(alns_1st, alns_2nd)

        sel_alndt = selAlnByMateMaxNDup(alns, maxunindupaln(prm),
                                        maxmulndupaln(prm))

        filter_func = function(x) {
            dt = data.table(as.data.frame(x))
            dt[, qname_HI := paste0(qname, '_', HI)]
            ## splice reads may be selected twice, need to only select them once
            dt[, i := seq_len(.N), by=list(qname, flag, HI)]
            dt[, to_select := ifelse((qname_HI %in% sel_alndt[, rdid]) & (i==1),
                                     T, F)]
            return(dt[, to_select])
        }

        filter_rules = FilterRules(list(tmp=filter_func))

        inbam = BamFile(finbam, yieldSize=maxyieldsize(prm))

        filterBam(inbam, foutbam(bam), filter=filter_rules,
                  param=ScanBamParam(what=c('qname', 'flag'), tag=c('HI'),
                                     which=grs))
    }
)


#' @importFrom  data.table         data.table
#' @importFrom  GenomicRanges      start
#' @importFrom  S4Vectors          mcols
#' @importFrom  GenomicAlignments  cigar
#'
selAlnByMateMaxNDup <- function(alns, max_uni_ndup, max_mul_ndup) {
    alndt = data.table( qname = mcols(alns)$qname,
                        HI    = mcols(alns)$HI,
                        cigar = cigar(alns),
                        start = start(alns),
                        flag  = mcols(alns)$flag )

    alndt[, `:=`( is_mul = ifelse(bitwAnd(flag, 0x100) > 0, T, F),
                  is_rd1 = ifelse(bitwAnd(flag, 0x40)  > 0, T, F),
                  is_rd2 = ifelse(bitwAnd(flag, 0x80)  > 0, T, F))]

    rd1dt = subset(alndt, is_rd1)
    rd2dt = subset(alndt, is_rd2)

    vnames = c('cigar', 'start', 'is_mul')
    setnames(rd1dt, vnames, paste0(vnames, '1'))
    setnames(rd2dt, vnames, paste0(vnames, '2'))
    rd1dt[, `:=`( is_rd1 = NULL, is_rd2 = NULL, flag = NULL )]
    rd2dt[, `:=`( is_rd1 = NULL, is_rd2 = NULL, flag = NULL )]

    all_matedt = merge(rd1dt, rd2dt, by=c('qname', 'HI'), all=T)

    matedt = subset(all_matedt, (! is.na(cigar1)) & (! is.na(cigar2)))
    matedt[, rdid := paste0(qname, '_', HI)]

    uni_matedt = subset(matedt, (! is_mul1) & (! is_mul2))
    mul_matedt = subset(matedt, is_mul1 | is_mul2)

    ## uni_matedt's HI are all equal to 1
    ## lean to keep alignments of multi-reads that map to fewer loci
    setkey(mul_matedt, HI)

    uni_matedt[, dupi := seq_len(.N), by=list(cigar1, start1, cigar2, start2)]
    mul_matedt[, dupi := seq_len(.N), by=list(cigar1, start1, cigar2, start2)]

    sel_uni_matedt = subset(uni_matedt, dupi <= max_uni_ndup)
    sel_mul_matedt = subset(mul_matedt, dupi <= max_mul_ndup)

    alndt[, rdid := paste0(qname, '_', HI)]
    sel_alndt = subset(alndt, rdid %in% c( sel_uni_matedt[, rdid],
                                           sel_mul_matedt[, rdid] ))

    return(sel_alndt)
}


#' @importFrom  Rsamtools          ScanBamParam
#' @importFrom  GenomicAlignments  readGAlignments
#'
selAlnInGRanges <- function(iggrs, finbam, flag_mate) {
    # read bam file for paired-end reads that overlap with intergenic regions
    # speed up the reading of bam file
    bamprm = ScanBamParam(flag=flag_mate, tag=c('HI'), what=c('flag', 'qname'),
                          which=iggrs)
    all_aln = readGAlignments(finbam, use.names=T, param=bamprm)

    out_aln = all_aln
    if ( length(all_aln) > 0 ) {
        # 1. get reads that fall outside intergenic regions
        outside_aln = subset(all_aln, ( start(all_aln) < start(iggrs)) |
                                      ( end(all_aln)   > end(iggrs))   )

        out_aln = subset(all_aln, ! names(all_aln) %in% names(outside_aln))
    }

    return(out_aln)
}
