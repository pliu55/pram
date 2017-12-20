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
#' @param  cufflinks  Cufflinks executable file. Required by method
#'                    'pooling+cufflinks' and 'cufflinks+cuffmerge'. Default: ''
#'
#' @param  cuffmerge  Cuffmerge executable file. Required by method
#'                    'cufflinks+cuffmerge'. Default: ''
#'
#' @param  stringtie  StringTie executable file. Required by method
#'                    'pooling+stringtie' and 'stringtie+merging'. Default: ''
#'
#' @param  taco       TACO executable file. Required by method
#'                    'cufflinks+taco'. Default: ''
#'
#' @return  NULL
#'
#' @importFrom  GenomeInfoDb  seqnames
#' @importFrom  BiocGenerics  strand
#' @importFrom  tools         file_path_sans_ext
#' @importFrom  data.table    rbindlist data.table
#'
#' @export
#'
buildModel <- function(fbams, iggrs, outdir, method='pooling+cufflinks',
                       nthr=1, cufflinks='', cuffmerge='', stringtie='',
                       taco='') {

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

    prm = new('Param')

    lo_method = tolower(method)
    checkMethodBin(lo_method, cufflinks, cuffmerge, stringtie, taco, prm)


    tmp_dir = paste0(outdir, '/pram_tmp/')
    if ( ! file.exists(tmp_dir) ) dir.create(tmp_dir, recursive=T)

    outdir(prm) = paste0(outdir, '/')
    tmpdir(prm) = tmp_dir
    nthreads(prm) = nthr
    chromoridt(prm) = getUniChromOriDt(iggrs)


    bams = initBy1ParentChromOri(fbams, iggrs, prm)

    if ( nthr == 1 ) {
        lapply(bams, filterParentBamByGRanges, prm)
    } else {
        mclapply(bams, filterParentBamByGRanges, prm, mc.cores=nthr)
    }

    func  = method2func[[lo_method]]
    label = method2label[[lo_method]]

    gtfs = func(bams, label, prm)

    outputModel(gtfs, label, prm)
}


checkMethodBin <- function(lo_method, cufflinks, cuffmerge, stringtie, taco,
                           prm) {
    if ( lo_method  == 'pooling+cufflinks' ) {
        checkCufflinksBin(cufflinks, prm)
    } else if ( lo_method %in% c('pooling+stringtie', 'stringtie+merging') ) {
        checkStringTieBin(stringtie, prm)
    } else if ( lo_method == 'cufflinks+cuffmerge' ) {
        checkCuffmergeBin(cuffmerge, prm)
        checkCufflinksBin(cufflinks, prm)
    } else if ( lo_method == 'cufflinks+taco' ) {
        checkTacoBin(taco, prm)
        checkCufflinksBin(cufflinks, prm)
    } else if ( ! lo_method %in% names(method2func) ) {
        msg = paste0('method= ', method, ' is not implemented. Must be one of ',
                     modeling_methods, "\n")
        stop(msg)
    }
}


outputModel <- function(gtfs, label, prm) {
    gtf = new('GTF')
    fgtf(gtf)     = paste0(outdir(prm), label, '.gtf')
    origin(gtf)   = label
    infokeys(gtf) = gtfinfokeys(prm)
    grangedt(gtf) = rbindlist(lapply(gtfs, function(x) grangedt(x)))
    writeGTF(gtf, fgtf(gtf), to_append=F)
    cat('File writtne:', fgtf(gtf), "\n")

    fout_grs = paste0(outdir(prm), label, '.rda')
    grs = makeGRangesFromDataFrame(grangedt(gtf), keep.extra.columns=T)
    save(grs, file=fout_grs)
    cat('R obj saved:', fout_grs, "\n")
}


modelByPoolingCufflinks <- function(bams, label, prm) {
    gtfs = modelByPoolingBams('cufflinks', bams, label, prm)
    return(gtfs)
}


modelByPoolingStringTie <- function(bams, label, prm) {
    gtfs = modelByPoolingBams('stringtie', bams, label, prm)
    return(gtfs)
}


modelByCufflinksCuffmerge <- function(bamdt, prm) {
}


modelByStringtieMerge <- function(bamdt, prm) {
}


modelByCufflinksTACO <- function(bamdt, prm) {
}


#' @importFrom  parallel       mcmapply
#' @importFrom  GenomicRanges  GRanges
#'
modelByPoolingBams <- function(method, bams, label, prm) {
    dt = chromoridt(prm)

    nthr = nthreads(prm)
    gtfs = NULL
    if ( nthr == 1 ) {
        plbams = mapply(poolBamByChromOri, dt$chrom, dt$ori,
                        MoreArgs=list(bams=bams, label=label, prm=prm))
        gtfs = lapply(plbams, modelByChromOriBam, method, prm)
    } else {
    }

    return(gtfs)
}


#' @importFrom  Rsamtools  mergeBam indexBam
#'
poolBamByChromOri <- function(chrom, ori, bams, label, prm) {
    fparentbaml = lapply(bams, function(x) {
                         if ( chrom(x) == chrom & ori(x) == ori )
                            return(foutbam(x))})

    fparentbams = unlist(fparentbaml)

    fplbam = paste0(tmpdir(prm), label, '.', chrom, '.', convertOri2Strand(ori),
                    '.bam')

    mergeBam(fparentbams, fplbam, overwrite=T)
    indexBam(fplbam)

    plbam = BamManager(fparentbams, fplbam, chrom, ori)

    return(plbam)
}


#' @importFrom  tools  file_path_sans_ext
#'
initBy1ParentChromOri <- function(fparentbams, iggrs, prm) {
    chromoridt = chromoridt(prm)

    bams = c()
    for ( fparentbam in fparentbams ) {
        for ( i in 1:nrow(chromoridt) ) {
            chrom = chromoridt[i]$chrom
            ori   = chromoridt[i]$ori
            strand = convertOri2Strand(ori)
            foutbam = paste0(tmpdir(prm),
                             file_path_sans_ext(basename(fparentbam)), '.',
                             chrom, '.', strand, '.bam')

            bams = c(bams, BamManager(fparentbam, foutbam, chrom, ori, iggrs))
        }
    }

    return(bams)
}


#' @importFrom  tools file_path_sans_ext
#'
modelByChromOriBam <- function(bam, method, prm) {
    finbam = foutbam(bam)
    chrom  = chrom(bam)
    strand = strand(bam)
    outdir = outdir(prm)

    label = file_path_sans_ext(basename(finbam))

    if ( ! file.exists(outdir) ) dir.create(outdir, recursive=T)
    setwd(outdir)

    fout = paste0(outdir, 'run.out')
    ferr = paste0(outdir, 'run.err')
    fout_gtf = paste0(outdir, 'transcripts.gtf')

    method2func = list( 'cufflinks' = getCufflinksArgs,
                        'stringtie' = getStringTieArgs )

    func = method2func[[method]]
    args = func(outdir, label, finbam, fout_gtf, prm)

    system2('nohup', args=args, stdout=fout, stderr=ferr)

    gtf = new('GTF')
    gtf = initFromGTFFile(gtf, fout_gtf, gtfinfokeys(prm), origin=label)

    return(gtf)
}


getCufflinksArgs <- function(outdir, label, finbam, fout_gtf, prm) {
    ## not use '--frag-bias-correct' or '--multi-read-correct'
    args = c( cufflinks(prm),
              '-o', outdir,
              '-p 1',
              '--library-type', cufflinkslibtype(prm),
              '--min-isoform-fraction',    minisoformfraction(prm),
              '--max-multiread-fraction',  maxmultireadfraction(prm),
              '--min-frags-per-transfrag', minfragspertransfrag(prm),
              '--label', label,
              '--quiet',
              '--no-update-check', finbam)

    return(args)
}


getStringTieArgs <- function(outdir, label, finbam, fout_gtf, prm) {
    args = c( stringtie(prm),
              finbam,
              '-o', fout_gtf,
              stringtielibtype(prm),
              '-l', label,
              '-f', minisoformfraction(prm),
              '-M', maxmultireadfraction(prm),
              '-c', minfragspertransfrag(prm),
              '-p 1' )

    return(args)
}
