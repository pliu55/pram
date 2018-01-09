#' Build transcript models from aligned RNA-seq data
#'
#' @param  fbams  A character vector of input RNA-seq BAM file(s).  Currently,
#'                PRAM only supports strand-specific paired-end data with the
#'                first mate on the right-most of transcript coordinate, i.e.,
#'                'fr-firststrand' by Cufflinks's definition
#'
#' @param  outdir  A character string defining the full name of a directory for
#'                 saving output files. PRAM will a folder named 'pram_tmp/'
#'                 under this directory to save temporary files.
#'
#' @param  mode  A character string defining PRAM's model building mode.
#'               Must be one of 'pooling+cufflinks', 'pooling+stringtie',
#'               'cufflinks+cuffmerge', 'stringtie+merging', 'cufflinks+taco',
#'               'cufflinks', or 'stringtie'. Default: 'pooling+cufflinks'
#'
#' @param  nthreads  An integer defining the number of threads to-be-used.
#'                   Default: 1
#'
#' @param  cufflinks  Cufflinks executable file.  Required by mode 'cufflinks',
#'                    'pooling+cufflinks' and 'cufflinks+cuffmerge'.
#'                    For mode 'cufflinks+cuffmerge',  executable files of
#'                    Cuffmerge, Cuffcompare, and gtf_to_sam from the Cufflinks
#'                    suite are assumed to be under the same folder as
#'                    Cufflinks.  Default: ''
#'
#' @param  stringtie  StringTie executable file.  Required by mode 'stringtie',
#'                    'pooling+stringtie' and 'stringtie+merging'. Default: ''
#'
#' @param  taco       TACO executable file. Required by mode
#'                    'cufflinks+taco'. Default: ''
#'
#' @return  NULL
#'
#' @export
#'
buildModel <- function(fbams, outdir, mode='pooling+cufflinks', nthreads=1,
                       cufflinks='', stringtie='', taco='') {

    prm = new('Param')
    fuserbams(prm) = fbams
    outdir(prm)    = ifelse(grepl('/$', outdir, perl=T), outdir,
                            paste0(outdir, '/'))
    mode(prm)      = tolower(mode)
    nthreads(prm)  = nthreads
    cufflinks(prm) = cufflinks
    stringtie(prm) = stringtie
    taco(prm)      = taco

    checkModeBin(prm)

    tmpdir = paste0(outdir, 'pram_tmp/')
    if ( ! file.exists(tmpdir) ) dir.create(tmpdir, recursive=T)
    tmpdir(prm) = tmpdir

    if ( mode %in% c('pooling+cufflinks', 'pooling+stringtie') ) {
        prm = def1StepManager(prm)
    } else if ( mode %in% c('cufflinks+cuffmerge', 'stringtie+merging',
                            'cufflinks+taco') ) {
        prm = def2StepManager(prm)
    } else if ( mode %in% c( 'cufflinks', 'stringtie' ) ) {
        prm = defCSManager(prm)
    }

    splitUserBams(prm)

    mode2func = list( 'pooling+cufflinks'   = modelByPoolingCufflinks,
                      'pooling+stringtie'   = modelByPoolingStringTie,
                      'cufflinks+cuffmerge' = modelByCufflinksCuffmerge,
                      'stringtie+merging'   = modelByStringtieMerge,
                      'cufflinks+taco'      = modelByCufflinksTACO,
                      'cufflinks'           = modelByCufflinks,
                      'stringtie'           = modelByStringTie )

    func = mode2func[[mode(prm)]]
    func(prm)

    outputModel(prm)
}


#' @importFrom  parallel  mclapply
#'
splitUserBams <- function(prm) {
    nthr = nthreads(prm)
    managerdt = managerdt(prm)

    if ( nthr == 1 ) {
        lapply(1:nrow(managerdt), splitUserBamByChromOri, prm)
    } else if ( nthr > 1 ) {
        mclapply(1:nrow(managerdt), splitUserBamByChromOri, prm, mc.cores=nthr)
    }
}


splitUserBamByChromOri <- function(i, prm) {
    managerdt = managerdt(prm)
    fuserbam     = managerdt[i]$fuserbam
    fchromoribam = managerdt[i]$fchromoribam
    chrom        = managerdt[i]$chrom
    strand       = managerdt[i]$strand

    filterBamByChromOri(fuserbam, fchromoribam, chrom, strand, prm)
}


#' @importFrom  Rsamtools  filterBam BamFile ScanBamParam
#' @importFrom  S4Vectors  FilterRules
#' @importFrom  GenomicRanges  GRanges
#'
filterBamByChromOri <- function(fuserbam, fchromoribam, chrom, strand, prm) {
    fr1ststrand2mate2flag = fr1ststrand2mate2flag(prm)
    flag1st = fr1ststrand2mate2flag[[strand]][['1stmate']]

    grs = GRanges(paste0(chrom, ':1-', maxchromlen(prm)))
    mate1st = scanBam(fuserbam, param=ScanBamParam(flag=flag1st, tag=c('HI'),
                                                   what=c('qname'), which=grs))
    seldt = data.table()
    if ( length(mate1st[[1]]$qname) > 0 ) {
        seldt = data.table( qname = mate1st[[1]]$qname,
                            HI    = mate1st[[1]]$tag$HI )
        seldt[, qname_HI := paste0(qname, '_', HI)]
    }

    ## scanBam first for qname flag, and HI, save them and filterBam
    filter_func = function(x) {
        dt = data.table(as.data.frame(x))
        dt[, qname_HI := paste0(qname, '_', HI)]
        dt[, to_select := ifelse(qname_HI %in% seldt$qname_HI, T, F)]
        return(dt$to_select)
    }

    filter_rules = FilterRules(list(tmp=filter_func))
    inbam = BamFile(fuserbam, yieldSize=maxyieldsize(prm))
    filterBam(inbam, fchromoribam, filter=filter_rules,
              param=ScanBamParam(what=c('qname'), tag=c('HI'), which=grs))
}


#' @importFrom  Rsamtools  indexBam
#'
genBamChromOri <- function(fbam)  {
    fbai = paste0(fbam, '.bai')
    if ( ! file.exists(fbai) ) indexBam(fbam)

    idxdt = data.table(idxstatsBam(fbam))
    idxdt[, fuserbam := fbam]
    setnames(idxdt, 'seqnames', 'chrom')

    plusdt  = idxdt[mapped > 0, .(fuserbam, chrom)]
    minusdt = idxdt[mapped > 0, .(fuserbam, chrom)]

    plusdt[,  strand := 'plus']
    minusdt[, strand := 'minus']

    dt = rbind(plusdt, minusdt)
    dt[, ori := convertStrand2Ori(strand)]

    return(dt)
}


#' @importFrom  tools  file_path_sans_ext
#'
defCSManager <- function(prm) {
    fuserbams = fuserbams(prm)
    nthr = nthreads(prm)

    dt = data.table()
    if ( nthr == 1 ) {
        dt = rbindlist(lapply(fuserbams, genBamChromOri))
    } else if ( nthr > 1 ) {
        dt = rbindlist(mclapply(fuserbams, genBamChromOri, mc.cores=nthr))
    }

    dt[, `:=`( label = mode2label(prm)[[mode(prm)]],
               bamid = file_path_sans_ext(basename(fuserbam)) )]

    dt[, tag := paste0(bamid, '_', chrom, '_', strand) ]

    dt[, `:=`( fchromoribam = paste0(tmpdir(prm), tag, '.bam'),
               fmdlbam      = paste0(tmpdir(prm), tag, '.bam'),
               mdldir       = paste0(tmpdir(prm), tag, '/'),
               fmdlgtf      = paste0(tmpdir(prm), tag, '/transcripts.gtf'),
               foutgtf      = paste0(outdir(prm), bamid, '_', label, '.gtf') )]

    all_chromoridt = dt[, .(chrom, ori)]
    chromoridt(prm) = unique(all_chromoridt, by=c('chrom', 'ori'))
    managerdt(prm) = dt

    return(prm)
}


#' @importFrom  tools  file_path_sans_ext
#'
def1StepManager <- function(prm) {
    fuserbams = fuserbams(prm)
    nthr = nthreads(prm)

    dt = data.table()
    if ( nthr == 1 ) {
        dt = rbindlist(lapply(fuserbams, genBamChromOri))
    } else if ( nthr > 1 ) {
        dt = rbindlist(mclapply(fuserbams, genBamChromOri, mc.cores=nthr))
    }

    dt[, label := mode2label(prm)[[mode(prm)]] ]
    dt[, tag := paste0(label, '_', chrom, '_', strand) ]

    dt[, `:=`( fchromoribam = paste0(tmpdir(prm),
                                     file_path_sans_ext(basename(fuserbam)),
                                     '.', tag, '.bam'),
               fmdlbam = paste0(tmpdir(prm), tag, '.bam'),
               mdldir  = paste0(tmpdir(prm), tag, '/'),
               fmdlgtf = paste0(tmpdir(prm), tag, '/transcripts.gtf'),
               foutgtf = paste0(outdir(prm), label, '.gtf') )]

    all_chromoridt = dt[, .(chrom, ori)]
    chromoridt(prm) = unique(all_chromoridt, by=c('chrom', 'ori'))
    managerdt(prm) = dt

    return(prm)
}


#' @importFrom  tools  file_path_sans_ext
#'
def2StepManager <- function(prm) {
    fuserbams = fuserbams(prm)
    nthr = nthreads(prm)

    dt = data.table()
    if ( nthr == 1 ) {
        dt = rbindlist(lapply(fuserbams, genBamChromOri))
    } else if ( nthr > 1 ) {
        dt = rbindlist(mclapply(fuserbams, genBamChromOri, mc.cores=nthr))
    }

    mode = mode(prm)
    dt[, label := mode2label(prm)[[mode(prm)]] ]
    dt[, `:=`( mdlpref = paste0(tmpdir(prm),
                                file_path_sans_ext(basename(fuserbam)),
                                '_', label, '_', chrom, '_', strand),
               mrgpref = paste0(tmpdir(prm), label, '_', chrom, '_', strand),
               mrgbase = ifelse(mode %in% c('cufflinks+cuffmerge',
                                            'stringtie+merging'), 'merged',
                                ifelse(mode == 'cufflinks+taco', 'assembly',
                                       NA)) )]

    dt[, `:=`( fchromoribam = paste0(mdlpref, '.bam'),
               fmdlbam      = paste0(mdlpref, '.bam'),
               mdldir       = paste0(mdlpref, '/'),
               fmdlgtf      = paste0(mdlpref, '/transcripts.gtf'),
               fmrglist = paste0(mrgpref, '_gtfs.list'),
               fmrg_out = paste0(mrgpref, '_run.out'),
               fmrg_err = paste0(mrgpref, '_run.err'),
               mrgdir   = paste0(mrgpref, '/'),
               fmrggtf  = paste0(mrgpref, '/', mrgbase, '.gtf'),
               foutgtf = paste0(outdir(prm), label, '.gtf'))]

    all_chromoridt = dt[, .(chrom, ori)]
    setkey(all_chromoridt, NULL)
    chromoridt(prm) = unique(all_chromoridt, by=c('chrom', 'ori'))
    managerdt(prm) = dt

    return(prm)
}


outputModel <- function(prm) {
    ## for 1-step or 2-step builder, there is only one foutgtf
    ## for cs builder, there may be multiple foutgtf

    dt = managerdt(prm)
    nthr = nthreads(prm)
    info_keys = gtfinfokeys(prm)
    foutgtfs = unique(dt$foutgtf)
    lapply(foutgtfs, outputModelByGTF, dt, info_keys, nthr)
}


outputModelByGTF <- function(in_foutgtf, managerdt, info_keys, nthr) {
    dt = managerdt[ foutgtf == in_foutgtf ]

    label = unique(dt$label)
    fmdlgtfs = unique(dt$fmdlgtf)

    grdt = data.table()
    if ( nthr == 1 ) {
        grdt = rbindlist(lapply(fmdlgtfs, getGRangeDT, info_keys, label))
    } else if ( nthr > 1 ) {
        grdt = rbindlist(mclapply(fmdlgtfs, getGRangeDT, info_keys, label,
                                  mc.cores=nthr))
    }

    gtf = new('GTF')
    fgtf(gtf)     = in_foutgtf
    origin(gtf)   = label
    infokeys(gtf) = info_keys
    grangedt(gtf) = grdt

    writeGTF(gtf, in_foutgtf, append=F)
}


modelByPoolingCufflinks <- function(prm) {
    modelByPoolingBams('cufflinks', prm)
}


modelByPoolingStringTie <- function(prm) {
    modelByPoolingBams('stringtie', prm)
}


modelByCufflinksCuffmerge <- function(prm) {
    modelByMethod('cufflinks', prm)
    mergeModels('cuffmerge', prm)
}


modelByStringtieMerge <- function(prm) {
    modelByMethod('stringtie', prm)
    mergeModels('stringtiemerge', prm)
}


modelByCufflinksTACO <- function(prm) {
    modelByMethod('cufflinks', prm)
    mergeModels('taco', prm)
}


modelByCufflinks <- function(prm) {
    modelByMethod('cufflinks', prm)
}


modelByStringTie <- function(prm) {
    modelByMethod('stringtie', prm)
}


mergeModels <- function(method, prm) {
    dt = chromoridt(prm)
    nthr = nthreads(prm)

    if ( nthr == 1 ) {
        mapply(mergeModelsByChromOri, dt$chrom, dt$ori,
               MoreArgs=list(method=method, prm=prm))
    } else {
        mcmapply(mergeModelsByChromOri, dt$chrom, dt$ori,
                 MoreArgs=list(method=method, prm=prm), mc.cores=nthr)
    }
}


mergeModelsByChromOri <- function(in_chrom, in_ori, method, prm) {
    dt = managerdt(prm)[ chrom == in_chrom & ori == in_ori ]
    fmdlgtfs = unique(dt$fmdlgtf)
    fmrggtf  = unique(dt$fmrggtf)
    mrgdir   = unique(dt$mrgdir)
    label    = unique(dt$label)
    fmrglist = unique(dt$fmrglist)
    fmrg_out = unique(dt$fmrg_out)
    fmrg_err = unique(dt$fmrg_err)

    if ( file.exists(mrgdir) ) unlink(mrgdir, recursive=T, force=T)

    ## taco does not work on an existing directory
    if ( method %in% c( 'cuffmerge', 'stringtiemerge' ) ) {
        dir.create(mrgdir, recursive=T)
        setwd(mrgdir)
    }

    sel = sapply(fmdlgtfs, function(x) length(grep('^#', readLines(x),
                                                   perl=T, invert=T)) > 0 )

    fsel_mdlgtfs = fmdlgtfs[sel]

   #gtf = new('GTF')
    if ( length(fsel_mdlgtfs) > 0 ) {
        write(paste0(fsel_mdlgtfs, collapse="\n"), fmrglist)

        method2func = list( 'cuffmerge'      = getCuffmergeArgs,
                            'stringtiemerge' = getStringTieMergeArgs,
                            'taco'           = getTacoArgs )

        func = method2func[[method]]
        args = func(fmrglist, mrgdir, fmrggtf, prm)

        if ( method == 'cuffmerge' ) {
            Sys.setenv(PATH=paste0(dirname(cufflinks(prm)), ':',
                                   Sys.getenv('PATH')))
        }
        system2('nohup', args=args, stdout=fmrg_out, stderr=fmrg_err)

       #gtf = initFromGTFFile(gtf, fmrggtf, gtfinfokeys(prm), origin=label)
    } else {
        write('no model was build from bam', fmrglist)
    }

   #return(gtf)
}


#' @importFrom  parallel       mcmapply
#'
modelByPoolingBams <- function(method, prm) {
    alldt = managerdt(prm)[, .(chrom, ori, fmdlbam)]
    setkey(alldt, NULL)
    dt = unique(alldt, by=c('chrom', 'ori', 'fmdlbam'))

    nthr = nthreads(prm)
    gtfs = NULL
    if ( nthr == 1 ) {
        mapply(poolBamByChromOri, dt$chrom, dt$ori, MoreArgs=list(prm=prm))
       #gtfl = lapply(dt$fmdlbam, modelByChromOriBam, method, prm)
        lapply(dt$fmdlbam, modelByChromOriBam, method, prm)
    } else if ( nthr > 1 ) {
        mcmapply(poolBamByChromOri, dt$chrom, dt$ori, MoreArgs=list(prm=prm),
                 mc.cores=nthr)
       #gtfl = mclapply(dt$fmdlbam, modelByChromOriBam, method, prm,
       #                mc.cores=nthr)
        mclapply(dt$fmdlbam, modelByChromOriBam, method, prm, mc.cores=nthr)
    }

   #return(gtfl)
}


#' @importFrom  Rsamtools  mergeBam indexBam
#'
poolBamByChromOri <- function(in_chrom, in_ori, prm) {
    managerdt = managerdt(prm)
    dt = managerdt[ chrom == in_chrom & ori == in_ori ]

    fchromoribams = dt$fchromoribam
    fmdlbam = unique(dt$fmdlbam)

    mergeBam(fchromoribams, fmdlbam, overwrite=T)
    indexBam(fmdlbam)
}


#' @importFrom  parallel  mcmapply
#'
modelByMethod <- function(method, prm) {
    dt = managerdt(prm)
    nthr = nthreads(prm)
    if ( nthr == 1 ) {
        lapply(dt$fmdlbam, modelByChromOriBam, method, prm)
    } else if ( nthr > 1 ) {
        mclapply(dt$fmdlbam, modelByChromOriBam, method, prm, mc.cores=nthr)
    }
}


#modelByChromOriBam <- function(in_chrom, in_ori, in_fmdlbam, method, prm) {
modelByChromOriBam <- function(in_fmdlbam, method, prm) {
    managerdt = managerdt(prm)
   #dt = managerdt[ chrom == in_chrom & ori == in_ori & fmdlbam == in_fmdlbam ]
    dt = managerdt[ fmdlbam == in_fmdlbam ]
    mdldir  = unique(dt$mdldir)
    fmdlbam = unique(dt$fmdlbam)
    fmdlgtf = unique(dt$fmdlgtf)
    label   = unique(dt$label)

    if ( file.exists(mdldir) ) unlink(mdldir, recursive=T, force=T)
    dir.create(mdldir, recursive=T)
    setwd(mdldir)

    fout = paste0(mdldir, 'run.out')
    ferr = paste0(mdldir, 'run.err')

    method2func = list( 'cufflinks' = getCufflinksArgs,
                        'stringtie' = getStringTieArgs )

    func = method2func[[method]]
    args = func(mdldir, label, fmdlbam, fmdlgtf, prm)

    system2('nohup', args=args, stdout=fout, stderr=ferr)

   #gtf = new('GTF')
   #gtf = initFromGTFFile(gtf, fmdlgtf, gtfinfokeys(prm), origin=label)

   #return(gtf)
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


getCuffmergeArgs <- function(fin_gtfs, outdir, fout_gtf, prm) {
    args = c( paste0(dirname(cufflinks(prm)), '/cuffmerge'),
              '-o', outdir,
             #'--ref-sequence', prm$fgnmfa,
              '--min-isoform-fraction', minisoformfraction(prm),
              '--num-threads 1',
              fin_gtfs )

    return(args)
}


getStringTieMergeArgs <- function(fin_gtfs, outdir, fout_gtf, prm) {
    args = c( stringtie(prm),
              '--merge',
              '-o', fout_gtf,
              '-F', mintrfpkmtoinclude(prm),
              '-T', mintrtpmtoinclude(prm),
              '-f', minisoformfraction(prm),
              fin_gtfs )

    return(args)
}


getTacoArgs <- function(fin_gtfs, outdir, fout_gtf, prm) {
    args = c( taco(prm),
              '--output-dir', outdir,
              '--num-processes 1',
              '--filter-min-expr', mintrtpmtoinclude(prm),
              '--isoform-frac',    minisoformfraction(prm),
              '--no-assemble-unstranded',
              fin_gtfs )

    return(args)
}


checkModeBin <- function(prm) {
    mode = mode(prm)
    if ( mode %in% c('pooling+cufflinks', 'cufflinks') ) {
        checkCufflinksBin(prm)
    } else if ( mode %in% c('pooling+stringtie', 'stringtie+merging',
                            'stringtie') ) {
        checkStringTieBin(prm)
    } else if ( mode == 'cufflinks+cuffmerge' ) {
        checkCuffmergeRequiredBins(prm)
    } else if ( mode == 'cufflinks+taco' ) {
        checkTacoBin(prm)
        checkCufflinksBin(prm)
    } else {
        msg = paste0('mode= ', mode, ' is not implemented. Must be one of ',
                     paste0(names(mode2label(prm)), collapse=', '), "\n")
        stop(msg)
    }
}
