library(data.table)
suppressMessages(library(GenomicAlignments))

main <- function() {
    context('buildModel')

    prm = new('Param')
    outdir = paste0(tempdir(), '/')

    testFilterBamByChromOri('chr10', 'plus',  outdir, prm)
    testFilterBamByChromOri('chr12', 'minus', outdir, prm)

    if ( ( grepl('biostat.wisc.edu', Sys.info()[['nodename']], fixed=TRUE) &
           ( Sys.info()[['user']] == 'pliu' ) ) |
         ( ( Sys.info()[['sysname']] == 'Darwin' ) & 
           ( Sys.info()[['user']] == 'peng' ) &
           ( file.exists('/ua/pliu/repe/pram') ) ) ) {
        testBuild(prm, outdir)
    }
}


testBuild <- function(prm, outdir) {
    fbams = c( system.file('extdata/bam/CMPRep1.sortedByCoord.clean.bam',
                           package='pram'),
               system.file('extdata/bam/CMPRep2.sortedByCoord.clean.bam',
                           package='pram') )

    if ( getOS() == 'LINUX' ) {
        cufflinks(prm) = '/ua/pliu/local/cufflinks-2.2.1/cufflinks'
        stringtie(prm) = '/ua/pliu/local/stringtie-1.3.3/stringtie'
        taco(prm)      = '/ua/pliu/local/taco-0.7.0/taco_run'
    } else if ( getOS() == 'OSX' ) {
        cufflinks(prm) = '/ua/pliu/local/osx/cufflinks-2.1.1/cufflinks'
        stringtie(prm) = '/ua/pliu/local/osx/stringtie-1.3.3b/stringtie'
        taco(prm)      = '/ua/pliu/local/osx/taco-v0.7.0/taco_run'
    }

    ## test Cufflinks, StringTie, or TACO and define them in `prm`
    prm = checkCufflinksBin(prm)
    testBinCF(prm)

    prm = checkCuffmergeRequiredBins(prm)
    testBinCFMG(prm)

    prm = checkStringTieBin(prm)
    testBinST(prm)

    prm = checkTacoBin(prm)
    testBinTC(prm)

    #nthr = 4
    nthr = 1
    fout_cf_gtfs = paste0(outdir, 'CMPRep', 1:2, '.sortedByCoord.clean_cf.gtf')
    fout_st_gtfs = paste0(outdir, 'CMPRep', 1:2, '.sortedByCoord.clean_st.gtf')


    cufflinks = cufflinks(prm)
    stringtie = stringtie(prm)
    taco      = taco(prm)

    for ( i in 1:length(fbams) ) {
        testBuildByCF(fbams[i], fout_cf_gtfs[i], nthr, cufflinks)
    }

    testBuildByPLCF(fbams, outdir, nthr, cufflinks)
    testBuildByCFMG(fbams, outdir, nthr, cufflinks)

    testBuildByCFTC(fbams, outdir, nthr, cufflinks, taco)
    for ( i in 1:length(fbams) ) {
        testBuildByST(fbams[i], fout_st_gtfs[i], nthr, stringtie)
    }

    testBuildByPLST(fbams, outdir, nthr, stringtie)
    testBuildBySTMG(fbams, outdir, nthr, stringtie)
}


testFilterBamByChromOri <- function(chrom, strand, outdir, prm) {
    fin = system.file('extdata/bam/CMPRep1.sortedByCoord.clean.bam',
                      package='pram')

    fcmp = system.file( paste0('extdata/bam/CMPRep1.sortedByCoord.clean.',
                               chrom, '.', strand, '.bam'), package='pram')

    fout = paste0(outdir, chrom, '.', strand, '.bam')

    filterBamByChromOri(fin, fout, chrom, strand, prm)

    grs = GRanges(paste0(chrom, ':1-', maxchromlen(prm)))
    bamprm = ScanBamParam(what=c('flag', 'qname'), which=grs)

    fcmpbai = paste0(fcmp, '.bai')
    if ( ! file.exists(fcmpbai)) indexBam(fcmp)

    cmp_alns = readGAlignments(fcmp, use.names=FALSE, param=bamprm)
    out_alns = readGAlignments(fout, use.names=FALSE, param=bamprm)

    test_that( paste0('buildModel::testFilterBamByChromOri: ', chrom, strand),
               expect_identical(cmp_alns, out_alns) )
}


testBuildByCFTC <- function(fbams, outdir, nthr, cufflinks, taco) {
    foutgtf = paste0(outdir, 'build_cftc.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbams, foutgtf, method='cftc', nthreads=nthr,
               cufflinks=cufflinks, taco=taco)
    test_that(paste0('buildModel::testBuildByCFTC: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildBySTMG <- function(fbams, outdir, nthr, stringtie) {
    foutgtf = paste0(outdir, 'build_stmg.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbams, foutgtf, method='stmg', nthreads=nthr,
               stringtie=stringtie)
    test_that(paste0('buildModel::testBuildBySTMG: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByCFMG <- function(fbams, outdir, nthr, cufflinks) {
    foutgtf = paste0(outdir, 'build_cfmg.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbams, foutgtf, method='cfmg', nthreads=nthr,
               cufflinks=cufflinks)
    test_that(paste0('buildModel::testBuildByCFMG: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByCF <- function(fbam, foutgtf, nthr, cufflinks) {
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbam, foutgtf, method='cf', nthreads=nthr, cufflinks=cufflinks)
    test_that(paste0('buildModel::testBuildByCF: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByST <- function(fbam, foutgtf, nthr, stringtie) {
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbam, foutgtf, method='st', nthreads=nthr, stringtie=stringtie)
    test_that(paste0('buildModel::testBuildByST: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByPLCF <- function(fbams, outdir, nthr, cufflinks) {
    foutgtf = paste0(outdir, 'build_plcf.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbams, foutgtf, method='plcf', nthreads=nthr,
               cufflinks=cufflinks )
    test_that(paste0('buildModel::testBuildByPLCF: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByPLST <- function(fbams, outdir, nthr, stringtie) {
    foutgtf = paste0(outdir, 'build_plst.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbams, foutgtf, method='plst', nthreads=nthr,
               stringtie=stringtie)
    test_that(paste0('buildModel::testBuildByPLST: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBinCF <- function(prm) {
    test_that('buildModel::testBinCF',
              expect_true( file.exists(cufflinks(prm))))
}


testBinCFMG <- function(prm) {
    cfdir = dirname(cufflinks(prm))
    cuffmerge   = paste0(cfdir, '/cuffmerge')
    cuffcompare = paste0(cfdir, '/cuffcompare')
    gtftosam    = paste0(cfdir, '/gtf_to_sam')

    test_that('buildModel::testBinCFMG',
              expect_true( file.exists(cufflinks(prm))))
    test_that('buildModel::testBinCFMG',
              expect_true( file.exists(cuffmerge)))
    test_that('buildModel::testBinCFMG',
              expect_true( file.exists(cuffcompare)))
    test_that('buildModel::testBinCFMG',
              expect_true( file.exists(gtftosam)))
}


testBinST <- function(prm) {
    test_that('buildModel::testBinST',
              expect_true( file.exists(stringtie(prm))))
}


testBinTC <- function(prm) {
    test_that('buildModel::testBinTC',
              expect_true( file.exists(taco(prm))))
}


main()
