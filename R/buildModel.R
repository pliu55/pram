#' Build transcript models from aligned RNA-seq data
#'
#' @param  fbams  character vector of a list of RNA-seq BAM files.  Currently,
#'                PRAM only supports strand-specific paired-end data with the
#'                first mate on the right-most of transcript coordinate, i.e.,
#'                'fr-firststrand' by Cufflinks's definition
#' @param  iggrs  a GenomicRanges object defining intergenic regions
#' @param  outdir  a character string defining the full name of a directory for
#'                 saving output files. PRAM will a folder named 'pram_tmp/'
#                  under this directory to save temporary files.
#' @param  method  a character string defining PRAM's model building method.
#'                 Must be one of 'pooling+cufflinks', 'pooling+stringtie',
#'                 'cufflinks+cuffmerge', 'stringtie+merging', or
#'                 'cufflinks+taco'.  Default: 'pooling+cufflinks'
#' @param  nthr  an integer defining the number of threads to-be-used.
#'               Default: 1
#'
#' @param  cufflinks  Cufflinks binary with full path. Default: ''
#'
#' @param  stringtie  StringTie binary with full path. Default: ''
#'
#' @return  NULL
#'
#' @export
#'
buildModel <- function(fbams, iggrs, outdir, method='pooling+cufflinks',
                       nthr=1, cufflinks='', stringtie='') {

    method2func = list( 'pooling+cufflinks'   = modelByPoolingCufflinks,
                        'pooling+stringtie'   = modelByPoolingStringTie,
                        'cufflinks+cuffmerge' = modelByCufflinksCuffmerge,
                        'stringtie+merging'   = modelByStringtieMerge,
                        'cufflinks+taco'      = modelByCufflinksTACO )

    method2label = list( 'pooling+cufflinks'   = 'plcf',
                         'pooling+stringtie'   = 'plst',
                         'cufflinks+cuffmerge' = 'cfmg',
                         'stringtie+merging'   = 'stmg',
                         'cufflinks+taco'      = 'cftc' )

    lo_method = tolower(method)
    if ( ! lo_method %in% names(method2func) ) {
        msg = paste0('method= ', method, ' is not implemented. Must be one of ',
                     modeling_methods, "\n")
        stop(msg)
    }

    tmp_dir = paste0(outdir, '/pram_tmp/')
    if ( ! file.exists(tmp_dir) ) dir.create(tmp_dir, recursive=T)

    prm = new('Param')
    outdir(prm) = paste0(outdir, '/')
    tmpdir(prm) = tmp_dir
    nthreads(prm) = nthr
    cufflinks(prm) = cufflinks
    stringtie(prm) = stringtie


    func  = method2func[[lo_method]]
    label = method2label[[lo_method]]
    func(fbams, iggrs, label, prm)
}


#' @importFrom  parallel       mcmapply
#' @importFrom  GenomicRanges  GRanges
#'
modelByPoolingCufflinks <- function(fbams, iggrs, label, prm) {
    origin = label
    info_keys = c('gene_id', 'transcript_id')

    indidt = filterBam4IG(fbams, iggrs, prm)

    indidt[, `:=`( fplbam = paste0(tmpdir(prm),
                                   origin, '.', chrom, '.', strand, '.bam'),
                   outdir = paste0(tmpdir(prm),
                                   origin, '_', chrom, '_', strand, '/') )]

    pldt = unique(indidt[, .(chrom, ori, strand, fplbam, outdir)],
                  by=c('chrom', 'ori'))

    nthr = nthreads(prm)
    gtfs = NULL
    if ( nthr == 1 ) {
        mapply(poolBamByChromOri, pldt$chrom, pldt$ori,
               MoreArgs=list(indidt=indidt))

        gtfs = mapply(buildCufflinksModelByChromOri, pldt$chrom, pldt$ori,
                      SIMPLIFY=F,  MoreArgs=list(pldt=pldt, origin=origin,
                      info_keys=info_keys, prm=prm))
    } else {
        mcmapply(poolBamByChromOri, pldt$chrom, pldt$ori,
                 MoreArgs=list(indidt=indidt), mc.cores=nthr)

        gtfs = mcmapply(buildCufflinksModelByChromOri, pldt$chrom, pldt$ori,
                        SIMPLIFY=F, MoreArgs=list(pldt=pldt, origin=origin,
                        info_keys=info_keys, prm=prm), mc.cores=nthr)
    }

    gtf = new('GTF')
    fgtf(gtf)     = paste0(outdir(prm), origin, '.gtf')
    origin(gtf)   = origin
    infokeys(gtf) = info_keys
    grangedt(gtf) = rbindlist(lapply(gtfs, function(x) grangedt(x)))
    writeGTF(gtf, fgtf(gtf), to_append=F)
    cat('File writtne:', fgtf(gtf), "\n")

    fout_grs = paste0(outdir(prm), origin, '.rda')
    grs = makeGRangesFromDataFrame(grangedt(gtf), keep.extra.columns=T)
    save(grs, file=fout_grs)
    cat('R obj saved:', fout_grs, "\n")

    return(indidt)
}


modelByPoolingStringTie <- function(bamdt, prm) {}
modelByCufflinksCuffmerge <- function(bamdt, prm) {}
modelByStringtieMerge <- function(bamdt, prm) {}
modelByCufflinksTACO <- function(bamdt, prm) {}


buildCufflinksModelByChromOri <- function(in_chrom, in_ori, pldt, origin,
                                          info_keys, prm) {
    seldt = pldt[chrom == in_chrom & ori == in_ori]
    fplbam = seldt$fplbam
    strand = seldt$strand
    outdir = seldt$outdir

    label = paste0(origin, '.', strand, '.', in_chrom)

    modelByCufflinks(outdir, label, fplbam, prm)

    fgtf = paste0(outdir, 'transcripts.gtf')
    gtf = new('GTF')
    gtf = initFromGTFFile(gtf, fgtf, info_keys, origin=origin)

    return(gtf)
}


modelByCufflinks <- function(outdir, label, finbam, prm) {
    if ( ! file.exists(outdir) ) dir.create(outdir, recursive=T)
    setwd(outdir)

    fout = paste0(outdir, 'run.out')
    ferr = paste0(outdir, 'run.err')

    ## not use '--frag-bias-correct' or '--multi-read-correct'
    args = c( cufflinks(prm),
              '-o', outdir,
              '-p 1',
              '--library-type', libtype(prm),
              '--min-isoform-fraction',    minisoformfraction(prm),
              '--max-multiread-fraction',  maxmultireadfraction(prm),
              '--min-frags-per-transfrag', minfragspertransfrag(prm),
              '--label', label,
              '--quiet',
              '--no-update-check', finbam)

    cat('nohup', args, "\n")
    system2('nohup', args=args, stdout=fout, stderr=ferr)
}


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
