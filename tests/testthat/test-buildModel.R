library(data.table)
suppressMessages(library(GenomicAlignments))

main <- function() {
    context('buildModel')

    prm = new('Param')
    os = getOS()
    outdir = paste0(tempdir(), '/')

    testFilterBamByChromOri('chr10', 'plus',  outdir, prm)
    testFilterBamByChromOri('chr12', 'minus', outdir, prm)


    fbams = c( system.file('extdata/bam/CMPRep1.sortedByCoord.clean.bam',
                           package='pram'),
               system.file('extdata/bam/CMPRep2.sortedByCoord.clean.bam',
                           package='pram') )


    testBinPLCF(fbams, outdir, 'plcf', prm, os)
    testBinPLST(fbams, outdir, 'plst', prm, os)
    testBinCFMG(fbams, outdir, 'cfmg', prm, os)
    testBinCFTC(fbams, outdir, 'cftc', prm, os)


    nthr = 4
    fout_cf_gtfs = paste0(outdir, 'CMPRep', 1:2, '.sortedByCoord.clean_cf.gtf')
    fout_st_gtfs = paste0(outdir, 'CMPRep', 1:2, '.sortedByCoord.clean_st.gtf')


    cufflinks = '/ua/pliu/local/cufflinks-2.2.1/cufflinks'
    stringtie = '/ua/pliu/local/stringtie-1.3.3/stringtie'
    taco      = '/ua/pliu/local/taco-0.7.0/taco_run'
    fgnmfa    = '/tier2/deweylab/pliu/genome/cufflinks_mm10_male/genome.fa'
    if ( file.exists(cufflinks) ) {
        for ( i in 1:length(fbams) ) {
            testBuildByCF(fbams[i], fout_cf_gtfs[i], nthr, cufflinks)
        }

        testBuildByPLCF(fbams, outdir, nthr, cufflinks, fgnmfa)
        testBuildByCFMG(fbams, outdir, nthr, cufflinks, fgnmfa)

        if ( file.exists(taco) ) {
            testBuildByCFTC(fbams, outdir, nthr, cufflinks, taco, fgnmfa)
        }
    }

    if ( file.exists(stringtie) ) {
        for ( i in 1:length(fbams) ) {
            testBuildByST(fbams[i], fout_st_gtfs[i], nthr, stringtie)
        }

        testBuildByPLST(fbams, outdir, nthr, stringtie)
        testBuildBySTMG(fbams, outdir, nthr, stringtie)
    }
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

    cmp_alns = readGAlignments(fcmp, use.names=F, param=bamprm)
    out_alns = readGAlignments(fout, use.names=F, param=bamprm)

    test_that( paste0('buildModel::testFilterBamByChromOri: ', chrom, strand),
               expect_identical(cmp_alns, out_alns) )
}


testBuildByCFTC <- function(fbams, outdir, nthr, cufflinks, taco, fgnmfa) {
    foutgtf = paste0(outdir, 'build_cftc.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbams, foutgtf, mode='cftc', nthreads=nthr,
               cufflinks=cufflinks, taco=taco, cufflinks_ref_fa=fgnmfa)
    test_that(paste0('buildModel::testBuildByCFTC: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildBySTMG <- function(fbams, outdir, nthr, stringtie) {
    foutgtf = paste0(outdir, 'build_stmg.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbams, foutgtf, mode='stmg', nthreads=nthr, stringtie=stringtie)
    test_that(paste0('buildModel::testBuildBySTMG: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByCFMG <- function(fbams, outdir, nthr, cufflinks, fgnmfa) {
    foutgtf = paste0(outdir, 'build_cfmg.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbams, foutgtf, mode='cfmg', nthreads=nthr, cufflinks=cufflinks,
               cufflinks_ref_fa=fgnmfa)
    test_that(paste0('buildModel::testBuildByCFMG: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByCF <- function(fbam, foutgtf, nthr, cufflinks) {
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbam, foutgtf, mode='cf', nthreads=nthr, cufflinks=cufflinks)
    test_that(paste0('buildModel::testBuildByCF: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByST <- function(fbam, foutgtf, nthr, stringtie) {
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbam, foutgtf, mode='st', nthreads=nthr, stringtie=stringtie)
    test_that(paste0('buildModel::testBuildByST: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByPLCF <- function(fbams, outdir, nthr, cufflinks, fgnmfa) {
    foutgtf = paste0(outdir, 'build_plcf.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbams, foutgtf, mode='plcf', nthreads=nthr, cufflinks=cufflinks,
               cufflinks_ref_fa=fgnmfa)
    test_that(paste0('buildModel::testBuildByPLCF: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByPLST <- function(fbams, outdir, nthr, stringtie) {
    foutgtf = paste0(outdir, 'build_plst.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    buildModel(fbams, foutgtf, mode='plst', nthreads=nthr, stringtie=stringtie)
    test_that(paste0('buildModel::testBuildByPLST: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBinPLCF <- function(fbams, outdir, mode, prm, os) {
    foutgtf = paste0(outdir, 'bin_plcf.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    url = os2cufflinks_url(prm)[[os]]
    test_that('buildModel::testBinPLCF',
              expect_error(
                  buildModel(fbams, foutgtf, mode),
                  regexp=paste0('cufflinks not found: \n',
                                'It can be downloaded at ', url, "\n"),
                  ignore.case=T))
}


testBinPLST <- function(fbams, outdir, mode, prm, os) {
    foutgtf = paste0(outdir, 'bin_plst.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    url = os2stringtie_url(prm)[[os]]
    test_that('buildModel::testBinPLST',
              expect_error(
                  buildModel(fbams, foutgtf, mode),
                  regexp=paste0('StringTie not found: \n',
                                'It can be downloaded at ', url, "\n"),
                  ignore.case=T))
}


testBinCFMG <- function(fbams, outdir, mode, prm, os) {
    foutgtf = paste0(outdir, 'bin_cfmg.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    url = os2cufflinks_url(prm)[[os]]
    test_that('buildModel::testBinCFMG',
              expect_error(
                  buildModel(fbams, foutgtf, mode),
                  regexp=paste0('Cufflinks suite can be downloaded at ', url,
                                "\n"),
                  ignore.case=T))
}


testBinCFTC <- function(fbams, outdir, mode, prm, os) {
    foutgtf = paste0(outdir, 'bin_cftc.gtf')
    if ( file.exists(foutgtf) ) file.remove(foutgtf)
    url = os2taco_url(prm)[[os]]
    test_that('buildModel::testBinCFTC',
              expect_error(
                  buildModel(fbams, foutgtf, mode),
                  regexp=paste0('TACO not found: \n',
                                'It can be downloaded at ', url, "\n"),
                  ignore.case=T))
}


main()
