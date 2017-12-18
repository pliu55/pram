#' @importFrom  Rsamtools  mergeBam indexBam
#'
poolBamByChromOri <- function(in_chrom, in_ori, indidt) {
    seldt = indidt[ chrom == in_chrom & ori == in_ori]
    figbams = seldt$figbam
    fplbam  = unique(seldt$fplbam)

    mergeBam(figbams, fplbam, overwrite=T)
    indexBam(fplbam)
}


#' filter bam files for a given intergenic regions
#'
#' @importFrom  GenomeInfoDb  seqnames
#' @importFrom  BiocGenerics  strand
#' @importFrom  data.table    rbindlist data.table
#' @importFrom  tools         file_path_sans_ext
#' @importFrom  parallel      mcmapply
#'
#' @return a data.table contains input and filtered bam files by chromosome and
#'         strand of a given intergenic genomic ranges
#'
filterBam4IG <- function(fbams, iggrs, prm) {
    dt = data.table( chrom = as.character(seqnames(iggrs)),
                     ori   = as.character(strand(iggrs)) )
    setkey(dt, NULL)
    unidt = unique(dt, by=c('chrom', 'ori'))
    cmbdt = rbindlist(lapply(fbams, function(x) copy(dt)[, finbam := x]))
    cmbdt[, strand := convertOri2Strand(ori)]

    cmbdt[, figbam := paste0(tmpdir(prm),
                             file_path_sans_ext(basename(finbam)),
                             '.', chrom, '.', strand, '.bam')]

    nthr = nthreads(prm)
    if ( nthr == 1 ) {
        mapply(filterBam4IGByChrOri, cmbdt$finbam, cmbdt$chrom, cmbdt$ori,
               MoreArgs=list(cmbdt=cmbdt, in_iggrs=iggrs, prm=prm))
    } else {
        mcmapply(filterBam4IGByChrOri, cmbdt$finbam, cmbdt$chrom, cmbdt$ori,
                 MoreArgs=list(cmbdt=cmbdt, in_iggrs=iggrs, prm=prm),
                 mc.cores=nthr)
    }

    return(cmbdt)
}


#' @importFrom  data.table    setkey data.table
#' @importFrom  GenomeInfoDb  seqnames
#' @importFrom  BiocGenerics  strand
#' @importFrom  Rsamtools     indexBam filterBam BamFile ScanBamParam
#' @importFrom  IRanges       filterRules
#'
filterBam4IGByChrOri <- function(in_finbam, in_chrom, in_ori, cmbdt, in_iggrs,
                                 prm) {
    selcmbdt = cmbdt[ finbam == in_finbam & chrom == in_chrom & ori == in_ori]
    strand = selcmbdt$strand
    figbam = selcmbdt$figbam

    iggrs = in_iggrs[ seqnames(in_iggrs) == in_chrom &
                      strand(in_iggrs) == in_ori ]

    finbam_index = paste0(in_finbam, '.bai')
    if ( ! file.exists(finbam_index) ) indexBam(in_finbam)

    strand2mate2flag = fr1ststrand2mate2flag(prm)
    flag_1stmate = strand2mate2flag[[strand]][['1stmate']]
    flag_2ndmate = strand2mate2flag[[strand]][['2ndmate']]

    alns_1st = selAlnInGRanges(iggrs, in_finbam, flag_1stmate)
    alns_2nd = selAlnInGRanges(iggrs, in_finbam, flag_2ndmate)
    alns = c(alns_1st, alns_2nd)

    max_uni_ndup = maxunindupaln(prm)
    max_mul_ndup = maxmulndupaln(prm)

    sel_alndt = selAlnByMateMaxNDup(alns, max_uni_ndup, max_mul_ndup)

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

    max_yield_size = maxyieldsize(prm)
    inbam = BamFile(in_finbam, yieldSize=max_yield_size)

    filterBam(inbam, figbam, filter=filter_rules,
              param=ScanBamParam(what=c('qname', 'flag'), tag=c('HI'),
                                 which=iggrs))
}


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
