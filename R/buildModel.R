#' @title Build transcript models from aligned RNA-seq data
#'
#' @param  in_bamv  A character vector of input BAM file(s). If mode 'cf'
#'                or 'st' is used, only one input RNA-seq BAM file is allowed.
#'                Currently, PRAM only supports strand-specific paired-end data
#'                with the first mate on the right-most of transcript
#'                coordinate, i.e., 'fr-firststrand' by Cufflinks's definition.
#'
#' @param  out_gtf  An output GTF file of predicted transcript models
#'
#' @param  method  A character string defining PRAM's model building method.
#'               Current available methods are:
#'               \itemize{
#'                   \item plcf: pooling   + cufflinks
#'                   \item plst: pooling   + stringtie
#'                   \item cfmg: cufflinks + cuffmerge
#'                   \item stmg: stringtie + merging
#'                   \item cftc: cufflinks + taco
#'                   \item cf:   cufflinks
#'                   \item st:   stringtie
#'               }
#'               Default: 'plcf'
#'
#' @param  nthreads  An integer defining the number of threads to-be-used.
#'                   Default: 1
#'
#' @param  tmpdir  A character string defining the full name of a folder for
#'                 saving temporary files. If not tmpdir is give, PRAM will
#'                 use R's tempdir().
#'
#' @param  keep_tmpdir  Whether to keep temporary files afterwards.
#'                      Default: False
#'
#' @param  cufflinks  Cufflinks executable.  Required by mode 'plcf',
#'                    'cfmg', and 'cf'.  For mode 'cfmg', executable files of
#'                    Cuffmerge, Cuffcompare, and gtf_to_sam from the Cufflinks
#'                    suite are assumed to be under the same folder as
#'                    Cufflinks.
#'                    All the executables are available to download for
#'                    Linux \url{http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz}
#'                    and MacOS \url{http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.OSX_x86_64.tar.gz}
#'                    Default: ''
#'
#' @param  stringtie  StringTie executable file.  Required by mode
#'                    'plst', 'stmg', and 'st'.
#'                    Default: ''
#'
#' @param  taco       TACO executable file. Required by mode 'cftc'.
#'                    Default: ''
#'
#' @param  cufflinks_ref_fa  Genome reference fasta file for Cufflinks. If
#'                           supplied, will be used for cufflinks's
#'                           '--frag-bias-correct' and cuffmerge's
#'                           '--ref-sequence' options.
#'                           Default: ''
#'
#' @return  NULL
#'
#' @export
#'
#' @examples
#'
#' fbams = c( system.file('extdata/bam/CMPRep1.sortedByCoord.clean.bam',
#'                        package='pram'),
#'            system.file('extdata/bam/CMPRep2.sortedByCoord.clean.bam',
#'                        package='pram') )
#'
#' foutgtf = tempfile(fileext='.gtf')
#'
#' buildModel(fbams, foutgtf, method='plst')
#'
#'
buildModel <- function(in_bamv, out_gtf, method='plcf', nthreads=1, tmpdir=NULL,
                       keep_tmpdir=FALSE, cufflinks='', stringtie='', taco='',
                       cufflinks_ref_fa='') {
    finbamv = in_bamv
    foutgtf = out_gtf
    mode    = method

    prm = new('Param')
    fuserbams(prm) = finbamv
    foutgtf(prm)   = foutgtf
    mode(prm)      = tolower(mode)
    nthreads(prm)  = nthreads
    cufflinks(prm) = cufflinks
    stringtie(prm) = stringtie
    taco(prm)      = taco
    cufflinksreffa(prm) = cufflinks_ref_fa

    tmpdir(prm) = createTmpdir(tmpdir, mode)

    prm = checkArgs(prm)

    if ( mode %in% c('plcf', 'plst') ) {
        prm = def1StepManager(prm)
    } else if ( mode %in% c('cfmg', 'stmg', 'cftc') ) {
        prm = def2StepManager(prm)
    } else if ( mode %in% c( 'cf', 'st' ) ) {
        prm = defCSManager(prm)
    }

    splitUserBams(prm)

    mode2func = list( 'plcf' = modelByPoolingCufflinks,
                      'plst' = modelByPoolingStringTie,
                      'cfmg' = modelByCufflinksCuffmerge,
                      'stmg' = modelByStringtieMerge,
                      'cftc' = modelByCufflinksTACO,
                      'cf'   = modelByCufflinks,
                      'st'   = modelByStringTie )

    func = mode2func[[mode(prm)]]
    func(prm)

    ## model will be removed if:
    ## - strand is not '+' or '-'
    ## - strand does not match the one where bam was derived
    outputCorrectStrandModel(prm)

    if ( ! keep_tmpdir ) {
        unlink(tmpdir(prm), recursive=T, force=T)
    }
}


createTmpdir <- function(tmpdir, mode) {
    if ( is.null(tmpdir) ) tmpdir = tempdir()

    sub_tmpdir = paste0(tmpdir, '/tmp_pram_', mode, '/')
    while ( file.exists(sub_tmpdir) ) {
        sub_tmpdir = paste0(tmpdir, '/tmp_pram_', mode, '_',
                            sample.int(999999, size=1), '/')
    }

    dir.create(sub_tmpdir, recursive=T)

    return(sub_tmpdir)
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


#' @importFrom  data.table     data.table
#' @importFrom  Rsamtools      filterBam BamFile ScanBamParam scanBam indexBam
#' @importFrom  S4Vectors      FilterRules
#' @importFrom  GenomicRanges  GRanges
#'
filterBamByChromOri <- function(fuserbam, fchromoribam, chrom, strand, prm) {
    fr1ststrand2mate2flag = fr1ststrand2mate2flag(prm)
    flag1st = fr1ststrand2mate2flag[[strand]][['1stmate']]

    grs = GRanges(paste0(chrom, ':1-', maxchromlen(prm)))
    fuserbai = paste0(fuserbam, '.bai')
    if ( ! file.exists(fuserbai) ) indexBam(fuserbam)

    ## Rsamtools version < 1.32.2 will give warnings on MacOS when reading tag
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


#' @importFrom  data.table  data.table
#' @importFrom  Rsamtools   indexBam idxstatsBam
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

    dt[, `:=`( mode  = mode(prm),
               bamid = file_path_sans_ext(basename(fuserbam)) )]

    dt[, tag := paste0(bamid, '_', chrom, '_', strand) ]

    dt[, `:=`( fchromoribam = paste0(tmpdir(prm), tag, '.bam'),
               fmdlbam      = paste0(tmpdir(prm), tag, '.bam'),
               mdldir       = paste0(tmpdir(prm), tag, '/'),
               fmdlgtf      = paste0(tmpdir(prm), tag, '/transcripts.gtf'),
               foutgtf      = paste0(tmpdir(prm), tag, '/transcripts.gtf') )]

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

    dt[, mode := mode(prm) ]
    dt[, tag := paste0(mode, '_', chrom, '_', strand) ]

    dt[, `:=`( fchromoribam = paste0(tmpdir(prm),
                                     file_path_sans_ext(basename(fuserbam)),
                                     '.', tag, '.bam'),
               fmdlbam = paste0(tmpdir(prm), tag, '.bam'),
               mdldir  = paste0(tmpdir(prm), tag, '/'),
               fmdlgtf = paste0(tmpdir(prm), tag, '/transcripts.gtf'),
               foutgtf = paste0(tmpdir(prm), tag, '/transcripts.gtf') )]

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

    dt[, mode := mode(prm) ]
    dt[, tag := paste0(mode, '_', chrom, '_', strand)]
    dt[, `:=`( mdlpref = paste0(tmpdir(prm),
                                file_path_sans_ext(basename(fuserbam)),
                                '_', tag),
               mrgpref = paste0(tmpdir(prm), tag),
               mrgbase = ifelse(mode %in% c('cfmg', 'stmg'), 'merged',
                                ifelse(mode == 'cftc', 'assembly', NA)) )]

    dt[, `:=`( fchromoribam = paste0(mdlpref, '.bam'),
               fmdlbam      = paste0(mdlpref, '.bam'),
               mdldir       = paste0(mdlpref, '/'),
               fmdlgtf      = paste0(mdlpref, '/transcripts.gtf'),
               fmrglist = paste0(mrgpref, '_gtfs.list'),
               fmrg_out = paste0(mrgpref, '_run.out'),
               fmrg_err = paste0(mrgpref, '_run.err'),
               mrgdir   = paste0(mrgpref, '/'),
               fmrggtf  = paste0(mrgpref, '/', mrgbase, '.gtf'),
               foutgtf  = paste0(mrgpref, '/pram_out.gtf') )]

    all_chromoridt = dt[, .(chrom, ori)]
    setkey(all_chromoridt, NULL)
    chromoridt(prm) = unique(all_chromoridt, by=c('chrom', 'ori'))
    managerdt(prm) = dt

    return(prm)
}


outputCorrectStrandModel <- function(prm) {
    dt        = managerdt(prm)
    nthr      = nthreads(prm)
    info_keys = gtfinfokeys(prm)
    foutgtf   = foutgtf(prm)
    mode      = mode(prm)

    stranddt = unique(dt[, .(ori, foutgtf)], by=c('ori', 'foutgtf'))

    grdt = data.table()
    if ( nthr == 1 ) {
        grdt = rbindlist(lapply(stranddt$foutgtf, getCorrectStrandExon,
                                stranddt, info_keys, mode))
    } else if ( nthr > 1 ) {
        grdt = rbindlist(mclapply(stranddt$foutgtf, getCorrectStrandExon,
                                  stranddt, info_keys, mode, mc.cores=nthr))
    }

    gtf = new('GTF')
    fgtf(gtf)     = foutgtf
    origin(gtf)   = mode
    infokeys(gtf) = info_keys
    grangedt(gtf) = grdt

    writeGTF(gtf, foutgtf, append=F)
}


getCorrectStrandExon <- function(fgtf, stranddt, info_keys, mode) {
    ori = stranddt[ foutgtf == fgtf ]$ori
    outdt = data.table()
    if ( file.exists(fgtf) ) {
        gtf = new('GTF')
        gtf = initFromGTFFile(gtf, fgtf, info_keys, origin=mode)
        grdt = grangedt(gtf)
        if ( nrow(grdt) > 0 ) {
            outdt = grdt[ (feature == 'exon') & (strand == ori) ]
        }
    }

    return(outdt)
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
    foutgtf  = unique(dt$foutgtf)
    mrgdir   = unique(dt$mrgdir)
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

        ## rename trid and geneid by chrom, strand, and id
        renameGTFTrGeneID(fmrggtf, foutgtf, prm)

    } else {
        write('no model was build from bam', fmrglist)
    }
}


renameGTFTrGeneID <- function(fingtf, foutgtf, prm) {
    gtf = new('GTF')
    info_keys = gtfinfokeys(prm)
    mode = mode(prm)
    gtf = initFromGTFFile(gtf, fingtf, info_keys, origin=mode)
    grdt = grangedt(gtf)[ feature == 'exon' ]

    grdt[, strand_label := ifelse(strand == '+', 'plus',
                                  ifelse(strand == '-', 'minus', NA))]
    grdt[, runid := paste0(mode, '_', chrom, '_', strand_label)]

    if ( mode == 'cfmg' ) {
        grdt[, `:=`( itr = tstrsplit(transcript_id, '_', fixed=T)[[2]],
                     ign = tstrsplit(gene_id,       '_', fixed=T)[[2]] )]
    } else if ( mode == 'stmg' ) {
        grdt[, `:=`( itr = tstrsplit(transcript_id, '.', fixed=T)[[3]],
                     ign = tstrsplit(gene_id,       '.', fixed=T)[[2]] )]
    } else if ( mode == 'cftc' ) {
        grdt[, `:=`( itr = tstrsplit(transcript_id, 'TU', fixed=T)[[2]],
                     ign = tstrsplit(gene_id,       'G',  fixed=T)[[2]] )]
    }

    grdt[, `:=`( trid = paste0(runid, '.', as.integer(ign), '.',
                               as.integer(itr)),
                 gnid = paste0(runid, '.', as.integer(ign)) )]

    grdt[, c('gene_id', 'transcript_id', 'strand_label', 'runid', 'itr',
             'ign') := NULL]
    setnames(grdt, c('trid', 'gnid'), c('transcript_id', 'gene_id'))

    outgtf = new('GTF')
    outgtf = initFromDataTable(outgtf, grdt, info_keys, origin=mode)
    writeGTF(outgtf, foutgtf, append=F)
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
        lapply(dt$fmdlbam, modelByChromOriBam, method, prm)
    } else if ( nthr > 1 ) {
        mcmapply(poolBamByChromOri, dt$chrom, dt$ori, MoreArgs=list(prm=prm),
                 mc.cores=nthr)
        mclapply(dt$fmdlbam, modelByChromOriBam, method, prm, mc.cores=nthr)
    }
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


modelByChromOriBam <- function(in_fmdlbam, method, prm) {
    managerdt = managerdt(prm)
    dt = managerdt[ fmdlbam == in_fmdlbam ]
    mdldir  = unique(dt$mdldir)
    fmdlbam = unique(dt$fmdlbam)
    fmdlgtf = unique(dt$fmdlgtf)
    tag     = unique(dt$tag)

    if ( file.exists(mdldir) ) unlink(mdldir, recursive=T, force=T)
    dir.create(mdldir, recursive=T)
    setwd(mdldir)

    fout = paste0(mdldir, 'run.out')
    ferr = paste0(mdldir, 'run.err')

    method2func = list( 'cufflinks' = getCufflinksArgs,
                        'stringtie' = getStringTieArgs )

    func = method2func[[method]]
    args = func(mdldir, tag, fmdlbam, fmdlgtf, prm)

    system2('nohup', args=args, stdout=fout, stderr=ferr)
}


getCufflinksArgs <- function(outdir, tag, finbam, fout_gtf, prm) {
    args = c( cufflinks(prm),
              '-o', outdir,
              '-p 1' )

    if ( cufflinksreffa(prm) != '' ) {
        args = c( args,
                  '--frag-bias-correct', cufflinksreffa(prm) )
    }

    args = c( args,
              '--multi-read-correct',

              '--library-type', cufflinkslibtype(prm),
              '--min-isoform-fraction',    minisoformfraction(prm),
              '--max-multiread-fraction',  maxmultireadfraction(prm),
              '--min-frags-per-transfrag', minfragspertransfrag(prm),
              '--label', tag,
              '--quiet',
              '--no-update-check', finbam)

    return(args)
}


getStringTieArgs <- function(outdir, tag, finbam, fout_gtf, prm) {
    args = c( stringtie(prm),
              finbam,
              '-o', fout_gtf,
              stringtielibtype(prm),
              '-l', tag,
              '-f', minisoformfraction(prm),
              '-M', maxmultireadfraction(prm),
              '-c', minfragspertransfrag(prm),
              '-p 1' )

    return(args)
}


getCuffmergeArgs <- function(fin_gtfs, outdir, fout_gtf, prm) {
    args = c( paste0(dirname(cufflinks(prm)), '/cuffmerge'),
              '-o', outdir )

    if ( cufflinksreffa(prm) != '' ) {
        args = c( args,
                  '--ref-sequence', cufflinksreffa(prm) )

    }

    args = c( args,
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


checkArgs <- function(prm) {
    mode = mode(prm)
    if ( mode %in% c('plcf', 'cf') ) {
        prm = checkCufflinksBin(prm)
    } else if ( mode %in% c('plst', 'stmg', 'st') ) {
        prm = checkStringTieBin(prm)
    } else if ( mode == 'cfmg' ) {
        prm = checkCuffmergeRequiredBins(prm)
    } else if ( mode == 'cftc' ) {
        prm = checkTacoBin(prm)
        prm = checkCufflinksBin(prm)
    } else {
        msg = paste0('mode "', mode, '" is not implemented. Must be one of ',
                     "[plcf, plst, cfmg, stmg, cftc, cf, st]\n")
        stop(msg)
    }

    for ( fbam in fuserbams(prm) ) {
        if ( ! file.exists(fbam) ) {
            stop(paste0('Cannot find input BAM: ', fbam, "\n"))
        }
    }

    if ( (mode %in% c( 'cf', 'st' )) & (length(fuserbams(prm)) > 1) ) {
        msg = paste0('length(finbamv) > 1. Only one input bam file is allowed ',
                     'for mode "', mode, '"\n')
        stop(msg)
    }
    return(prm)
}
