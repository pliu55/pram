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


    func  = method2func[[lo_method]]
    label = method2label[[lo_method]]
    func(fbams, iggrs, label, prm)
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

    ## build models in a separate function so that can be tested from pooled bam
    browser()


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

   #cat('nohup', args, "\n")
    system2('nohup', args=args, stdout=fout, stderr=ferr)
}
