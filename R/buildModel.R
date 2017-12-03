#' Build transcript models from aligned RNA-seq data
#'
#' @param  fbams  character vector of a list of RNA-seq BAM files.  Currently,
#'                PRAM only supports strand-specific paired-end data with the
#'                first mate on the right-most of transcript coordinate, i.e.,
#'                'fr-firststrand' by Cufflinks's definition
#' @param  iggrs  a GenomicRanges object defining intergenic regions
#' @param  outdir  a character string defining the full name of a directory for
#'                 saving output files
#' @param  method  a character string defining PRAM's model building method.
#'                 Must be one of 'pooling+cufflinks', 'pooling+stringtie',
#'                 'cufflinks+cuffmerge', 'stringtie+merging', or
#'                 'cufflinks+taco'.  Default: 'pooling+cufflinks'
#' @param  nthr  an integer defining the number of threads to-be-used.
#'               Default: 1
#'
#' @return  NULL
#'
#' @export
#'
buildModel <- function(fbams, iggrs, outdir, method='pooling+cufflinks',
                       nthr=1) {

    tmpdir = paste0(tempdir(), '/')

    bamdt = filterBam4IG(fbams, iggrs, nthr, tmpdir)

    browser()
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
filterBam4IG <- function(fbams, iggrs, nthr, tmpdir) {
    dt = data.table( chrom = as.character(seqnames(iggrs)),
                     ori   = as.character(strand(iggrs)) )
    setkey(dt, NULL)
    unidt = unique(dt, by=c('chrom', 'ori'))
    cmbdt = rbindlist(lapply(fbams, function(x) copy(dt)[, finbam := x]))
    cmbdt[, strand := ifelse(ori == '+', 'plus',
                             ifelse( ori == '-', 'minus', NA))]

    cmbdt[, foutbam := paste0(tmpdir, file_path_sans_ext(basename(finbam)),
                              '.', chrom, '.', strand, '.bam')]

   #mapply(filterBam4IGByChrOri, cmbdt$finbam, cmbdt$chrom, cmbdt$ori,
   #       MoreArgs=list(cmbdt=cmbdt, iggrs=iggrs))

    mcmapply(filterBam4IGByChrOri, cmbdt$finbam, cmbdt$chrom, cmbdt$ori,
             MoreArgs=list(cmbdt=cmbdt, iggrs=iggrs), SIMPLIFY=F, mc.cores=nthr)

    return(cmbdt)
}


#' @importFrom  GenomeInfoDb  seqnames
#' @importFrom  BiocGenerics  strand
#' @importFrom  Rsamtools     indexBam filterBam BamFile ScanBamParam
#'
filterBam4IGByChrOri <- function(in_finbam, in_chrom, in_ori, cmbdt, iggrs) {
    foutbam = cmbdt[ finbam == in_finbam & chrom == in_chrom & ori == in_ori,
                     foutbam ]

    finbam_index = paste0(in_finbam, '.bai')
    if ( ! file.exists(finbam_index) ) indexBam(in_finbam)

    prm = new('Param')
    inbam = BamFile(in_finbam, yieldSize=getMaxYieldSize(prm))
    grs = iggrs[ seqnames(iggrs) == in_chrom & strand(iggrs) == in_ori ]
    filterBam(inbam, foutbam, param=ScanBamParam(which=grs))
}
