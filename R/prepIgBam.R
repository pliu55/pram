#' @title Extract alignments in intergenic regions from BAM files
#'
#' @param  finbam   Full name of an input RNA-seq BAM file.
#'                  Currently, PRAM only supports strand-specific paired-end
#'                  data with the first mate on the right-most of transcript
#'                  coordinate, i.e., 'fr-firststrand' by Cufflinks's
#'                  definition
#'
#' @param  iggrs  A GenomicRanges object defining intergenic regions
#'
#' @param  foutbam  Full name of an output BAM file to save all alignment fell
#'                  into intergenic regions
#'
#' @param  max_uni_n_dup_aln  Maximum number of uniquely mapped fragments
#'                            to report per each alignment.
#'                            Default: 10
#'
#' @param  max_mul_n_dup_aln  Maximum number of multi-mapping fragments
#'                            to report per each alignment.
#'                            Default: 10
#'
#' @return None
#'
#' @export
#'
#' @examples
#'
#' finbam = system.file('extdata/bam/CMPRep2.sortedByCoord.raw.bam',
#'                      package='pram')
#'
#' iggrs = GenomicRanges::GRanges('chr10:77236000-77247000:+')
#'
#' foutbam = tempfile(fileext='.bam')
#'
#' prepIgBam(finbam, iggrs, foutbam)
#'
#'
prepIgBam <- function(finbam, iggrs, foutbam,
                      max_uni_n_dup_aln=10, max_mul_n_dup_aln=10) {

    prm = new('Param')
    maxunindupaln(prm) = max_uni_n_dup_aln
    maxmulndupaln(prm) = max_mul_n_dup_aln

    ori_iggrs = addOri(iggrs)

    extractBam(finbam, ori_iggrs, foutbam, prm)
}


#' @importFrom  BiocGenerics  strand strand<-
#'
addOri <- function(iggrs) {
    w_ori_iggrs  = iggrs[ ( strand(iggrs) == '+' ) |
                          ( strand(iggrs) == '-' ) ]
    no_ori_iggrs = iggrs[ ( strand(iggrs) != '+' ) &
                          ( strand(iggrs) != '-' ) ]

    plus_iggrs  = copy(no_ori_iggrs)
    minus_iggrs = copy(no_ori_iggrs)

    strand(plus_iggrs)  = '+'
    strand(minus_iggrs) = '-'

    outgrs = c( w_ori_iggrs, plus_iggrs, minus_iggrs )

    return(outgrs)
}


#' @importFrom  data.table  setkey data.table
#' @importFrom  Rsamtools   indexBam filterBam BamFile ScanBamParam
#' @importFrom  S4Vectors   FilterRules
#'
extractBam <- function(finbam, iggrs, foutbam, prm) {
    qname_HI = qname = HI = i = flag = to_select = rdid = NULL
    finbam_index = paste0(finbam, '.bai')
    if ( ! file.exists(finbam_index) ) indexBam(finbam)

    fr1ststrand2mate2flag = fr1ststrand2mate2flag(prm)
    plus_flag_1st  = fr1ststrand2mate2flag[['plus']][['1stmate']]
    plus_flag_2nd  = fr1ststrand2mate2flag[['plus']][['2ndmate']]

    minus_flag_1st = fr1ststrand2mate2flag[['minus']][['1stmate']]
    minus_flag_2nd = fr1ststrand2mate2flag[['minus']][['2ndmate']]

    plus_iggrs  = iggrs[ strand(iggrs) == '+' ]
    minus_iggrs = iggrs[ strand(iggrs) == '-' ]

    plus_alns_1st  = selAlnInGRanges(plus_iggrs, finbam, plus_flag_1st)
    plus_alns_2nd  = selAlnInGRanges(plus_iggrs, finbam, plus_flag_2nd)

    minus_alns_1st = selAlnInGRanges(minus_iggrs, finbam, minus_flag_1st)
    minus_alns_2nd = selAlnInGRanges(minus_iggrs, finbam, minus_flag_2nd)

    alns = c(plus_alns_1st, plus_alns_2nd, minus_alns_1st, minus_alns_2nd)

    sel_alndt = selAlnByMateMaxNDup(alns, maxunindupaln(prm),
                                    maxmulndupaln(prm))

    filter_func = function(x) {
        dt = data.table(as.data.frame(x))
        dt[, qname_HI := paste0(qname, '_', HI)]
        ## splice reads may be selected twice, need to only select them once
        dt[, i := seq_len(.N), by=list(qname, flag, HI)]
        dt[, to_select := ifelse((qname_HI %in% sel_alndt[, rdid]) & (i==1),
                                 TRUE, FALSE)]
        return(dt[, to_select])
    }

    filter_rules = FilterRules(list(tmp=filter_func))

    inbam = BamFile(finbam, yieldSize=maxyieldsize(prm))

    filterBam(inbam, foutbam, filter=filter_rules,
              param=ScanBamParam(what=c('qname', 'flag'), tag=c('HI'),
                                 which=iggrs))
}


#' @importFrom  data.table         data.table
#' @importFrom  GenomicRanges      start
#' @importFrom  S4Vectors          mcols
#' @importFrom  GenomicAlignments  cigar
#'
selAlnByMateMaxNDup <- function(alns, max_uni_ndup, max_mul_ndup) {
    flag = is_rd1 = is_rd2 = cigar1 = cigar2 = rdid = qname = HI = is_mul1 = is_mul2 = dupi = start1 = start2 = NULL
    alndt = data.table( qname = mcols(alns)$qname,
                        HI    = mcols(alns)$HI,
                        cigar = cigar(alns),
                        start = start(alns),
                        flag  = mcols(alns)$flag )

    alndt[, `:=`( is_mul = ifelse(bitwAnd(flag, 0x100) > 0, TRUE, FALSE),
                  is_rd1 = ifelse(bitwAnd(flag, 0x40)  > 0, TRUE, FALSE),
                  is_rd2 = ifelse(bitwAnd(flag, 0x80)  > 0, TRUE, FALSE))]

    rd1dt = subset(alndt, is_rd1)
    rd2dt = subset(alndt, is_rd2)

    vnames = c('cigar', 'start', 'is_mul')
    setnames(rd1dt, vnames, paste0(vnames, '1'))
    setnames(rd2dt, vnames, paste0(vnames, '2'))
    rd1dt[, `:=`( is_rd1 = NULL, is_rd2 = NULL, flag = NULL )]
    rd2dt[, `:=`( is_rd1 = NULL, is_rd2 = NULL, flag = NULL )]

    all_matedt = merge(rd1dt, rd2dt, by=c('qname', 'HI'), all=TRUE)

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


#' @importFrom  IRanges            subsetByOverlaps
#' @importFrom  Rsamtools          ScanBamParam
#' @importFrom  GenomicAlignments  readGAlignments
#'
selAlnInGRanges <- function(iggrs, finbam, flag_mate) {
    bamprm = ScanBamParam(flag=flag_mate, tag=c('HI'), what=c('flag', 'qname'),
                          which=iggrs)
    alns = readGAlignments(finbam, use.names=FALSE, param=bamprm)
    olalns = subsetByOverlaps(alns, iggrs, type='within', ignore.strand=TRUE)

    return(olalns)
}
