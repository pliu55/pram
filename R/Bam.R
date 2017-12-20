#' @import methods
#'
Bam = setClass('Bam',
               slots = list( fparentbam = 'character',
                             foutbam    = 'character',
                             chrom      = 'character',
                             ori        = 'character',
                             strand     = 'character',
                             grs        = 'GRanges',
                             prm        = 'Param' )
)


setGeneric('fparentbam',   function(x) standardGeneric('fparentbam'))
setGeneric('foutbam',      function(x) standardGeneric('foutbam'))
setGeneric('chrom',        function(x) standardGeneric('chrom'))
setGeneric('ori',          function(x) standardGeneric('ori'))
## strand was defined in BioGenerics, cannot set here, just setMethod is fine
#setGeneric('strand',       function(x) standardGeneric('strand'))
setGeneric('grs',          function(x) standardGeneric('grs'))
setGeneric('prm',          function(x) standardGeneric('prm'))

setGeneric('fparentbam<-', function(x, value) standardGeneric('fparentbam<-'))
setGeneric('foutbam<-',    function(x, value) standardGeneric('foutbam<-'))
setGeneric('chrom<-',      function(x, value) standardGeneric('chrom<-'))
setGeneric('ori<-',        function(x, value) standardGeneric('ori<-'))
setGeneric('strand<-',     function(x, value) standardGeneric('strand<-'))
setGeneric('grs<-',        function(x, value) standardGeneric('grs<-'))
setGeneric('prm<-',        function(x, value) standardGeneric('prm<-'))

setMethod('fparentbam', 'Bam', function(x) x@fparentbam)
setMethod('foutbam',    'Bam', function(x) x@foutbam)
setMethod('chrom',      'Bam', function(x) x@chrom)
setMethod('ori',        'Bam', function(x) x@ori)
setMethod('strand',     'Bam', function(x) x@strand)
setMethod('grs',        'Bam', function(x) x@grs)
setMethod('prm',        'Bam', function(x) x@prm)

setReplaceMethod('fparentbam', 'Bam',
                 function(x, value) {x@fparentbam=value; x})
setReplaceMethod('foutbam', 'Bam', function(x, value) {x@foutbam=value; x})
setReplaceMethod('chrom',   'Bam', function(x, value) {x@chrom=value;   x})
setReplaceMethod('ori',     'Bam', function(x, value) {x@ori=value;     x})
setReplaceMethod('strand',  'Bam', function(x, value) {x@strand=value;  x})
setReplaceMethod('grs',     'Bam', function(x, value) {x@grs=value;     x})
setReplaceMethod('prm',     'Bam', function(x, value) {x@prm=value;     x})


setMethod('show', 'Bam',
    function(object) {
        cat('fparentbam:', fparentbam(object), "\n")
        cat('foutbam:',    foutbam(object), "\n")
        cat('chrom:',      chrom(object),   "\n")
        cat('ori:',        ori(object),     "\n")
        cat('strand:',     strand(object),  "\n")
        print(grs(object))
    }
)


#' @importFrom  BiocGenerics  strand
#' @importFrom  GenomeInfoDb  seqnames
#' @importFrom  tools         file_path_sans_ext
#'
setMethod('initialize', 'Bam',
    function(.Object, fparentbam, foutbam, chrom, ori, allgrs, prm){
        fparentbam(.Object) = fparentbam
        foutbam(.Object)    = foutbam
        chrom(.Object)      = chrom
        ori(.Object)        = ori
        strand(.Object)     = convertOri2Strand(ori)
        grs(.Object)        = allgrs[ seqnames(allgrs) == chrom &
                                      strand(allgrs) == ori ]
        prm(.Object)        = prm

        return(.Object)
    }
)


#' @importFrom  data.table    setkey data.table
#' @importFrom  Rsamtools     indexBam filterBam BamFile ScanBamParam
#' @importFrom  IRanges       filterRules
#'
setGeneric('filterParentBamByGRanges',
           function(bam) standardGeneric('filterParentBamByGRanges'))
setMethod('filterParentBamByGRanges', 'Bam',
    function(bam) {
        grs = grs(bam)
        prm = prm(bam)
        finbam = fparentbam(bam)

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


initBam <- function(fparentbams, iggrs, prm) {
    dt = data.table( chrom = as.character(seqnames(iggrs)),
                     ori   = as.character(strand(iggrs)) )
    setkey(dt, NULL)
    unidt = unique(dt, by=c('chrom', 'ori'))
    indidt = rbindlist(lapply(fparentbams,
                              function(x) copy(unidt)[, fparentbam := x]))

    indidt[, strand := convertOri2Strand(ori)]
    bams = apply(indidt, 1, initBamFromParentChromOri, iggrs, prm)

    return(bams)
}


initBamFromParentChromOri <- function(vec, iggrs, prm) {
    fparentbam = vec[['fparentbam']]
    chrom      = vec[['chrom']]
    ori        = vec[['ori']]
    strand     = vec[['strand']]
    foutbam    = paste0(tmpdir(prm), file_path_sans_ext(basename(fparentbam)),
                        '.', chrom, '.', strand, '.bam')

    bam = Bam(fparentbam, foutbam, chrom, ori, iggrs, prm)
}


#' filter bam files for a given intergenic regions
#'
#' @importFrom  parallel      mcmapply
#'
#' @return a data.table contains input and filtered bam files by chromosome and
#'         strand of a given intergenic genomic ranges
#'
#filterBam4IG <- function(cmbdt, iggrs, prm) {
#   nthr = nthreads(prm)
#   if ( nthr == 1 ) {
#       mapply(filterBam4IGByChrOri, cmbdt$finbam, cmbdt$chrom, cmbdt$ori,
#              MoreArgs=list(cmbdt=cmbdt, in_iggrs=iggrs, prm=prm))
#   } else {
#       mcmapply(filterBam4IGByChrOri, cmbdt$finbam, cmbdt$chrom, cmbdt$ori,
#                MoreArgs=list(cmbdt=cmbdt, in_iggrs=iggrs, prm=prm),
#                mc.cores=nthr)
#   }
#}


#' @importFrom  data.table    setkey data.table
#' @importFrom  GenomeInfoDb  seqnames
#' @importFrom  BiocGenerics  strand
#' @importFrom  Rsamtools     indexBam filterBam BamFile ScanBamParam
#' @importFrom  IRanges       filterRules
#'
#filterBam4IGByChrOri <- function(in_finbam, in_chrom, in_ori, cmbdt, in_iggrs,
#                                prm) {
#   selcmbdt = cmbdt[ finbam == in_finbam & chrom == in_chrom & ori == in_ori]
#   strand = selcmbdt$strand
#   figbam = selcmbdt$figbam

#   iggrs = in_iggrs[ seqnames(in_iggrs) == in_chrom &
#                     strand(in_iggrs) == in_ori ]

#   finbam_index = paste0(in_finbam, '.bai')
#   if ( ! file.exists(finbam_index) ) indexBam(in_finbam)

#   strand2mate2flag = fr1ststrand2mate2flag(prm)
#   flag_1stmate = strand2mate2flag[[strand]][['1stmate']]
#   flag_2ndmate = strand2mate2flag[[strand]][['2ndmate']]

#   alns_1st = selAlnInGRanges(iggrs, in_finbam, flag_1stmate)
#   alns_2nd = selAlnInGRanges(iggrs, in_finbam, flag_2ndmate)
#   alns = c(alns_1st, alns_2nd)

#   max_uni_ndup = maxunindupaln(prm)
#   max_mul_ndup = maxmulndupaln(prm)

#   sel_alndt = selAlnByMateMaxNDup(alns, max_uni_ndup, max_mul_ndup)

#   filter_func = function(x) {
#       dt = data.table(as.data.frame(x))
#       dt[, qname_HI := paste0(qname, '_', HI)]
#       ## splice reads may be selected twice, need to only select them once
#       dt[, i := seq_len(.N), by=list(qname, flag, HI)]
#       dt[, to_select := ifelse((qname_HI %in% sel_alndt[, rdid]) & (i==1),
#                                T, F)]
#       return(dt[, to_select])
#   }

#   filter_rules = FilterRules(list(tmp=filter_func))

#   max_yield_size = maxyieldsize(prm)
#   inbam = BamFile(in_finbam, yieldSize=max_yield_size)

#   filterBam(inbam, figbam, filter=filter_rules,
#             param=ScanBamParam(what=c('qname', 'flag'), tag=c('HI'),
#                                which=iggrs))
#}


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
