library(data.table)
suppressMessages(library(GenomicAlignments))

main <- function() {
    context('buildModel')

    prm = new('Param')

    testFilterBamByChromOri('chr10', 'plus',  prm)
    testFilterBamByChromOri('chr12', 'minus', prm)


    fbams = c( system.file('extdata/bam/CMPRep1.sortedByCoord.clean.bam',
                           package='pram'),
               system.file('extdata/bam/CMPRep2.sortedByCoord.clean.bam',
                           package='pram') )

    testBin(fbams, prm)


    nthr = 4
    outdir = paste0(tempdir(), '/')
    fout_cf_gtfs = paste0(outdir, 'CMPRep', 1:2, '.sortedByCoord.clean_cf.gtf')
    fout_st_gtfs = paste0(outdir, 'CMPRep', 1:2, '.sortedByCoord.clean_st.gtf')

    cufflinks = '/ua/pliu/local/cufflinks-2.2.1/cufflinks'
    stringtie = '/ua/pliu/local/stringtie-1.3.3/stringtie'
    taco      = '/ua/pliu/local/taco-0.7.0/taco_run'
    if ( file.exists(cufflinks) ) {
        testBuildByCF(fbams, outdir, nthr, cufflinks, fout_cf_gtfs)

        testBuildByPLCF(fbams, outdir, nthr, cufflinks)
        testBuildByCFMG(fbams, outdir, nthr, cufflinks)

        if ( file.exists(taco) ) {
            testBuildByCFTC(fbams, outdir, nthr, cufflinks, taco)
        }
    }

    if ( file.exists(stringtie) ) {
        testBuildByST(fbams, outdir, nthr, stringtie, fout_st_gtfs)

        testBuildByPLST(fbams, outdir, nthr, stringtie)
        testBuildBySTMG(fbams, outdir, nthr, stringtie)
    }
}


testBin <- function(fbams,prm) {
    os = getOS()
    nthr = 1
    testBinPLCF(fbams, 'pooling+cufflinks',   nthr, prm, os)
    testBinPLST(fbams, 'pooling+stringtie',   nthr, prm, os)
    testBinCFMG(fbams, 'cufflinks+cuffmerge', nthr, prm, os)
    testBinCFTC(fbams, 'cufflinks+taco',      nthr, prm, os)
}


testFilterBamByChromOri <- function(chrom, strand, prm) {
    outdir = paste0(tempdir(), '/')
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


testBuildByCFTC <- function(fbams, outdir, nthr, cufflinks, taco) {
    buildModel(fbams, outdir, mode='cufflinks+taco', nthreads=nthr,
               cufflinks=cufflinks, taco=taco)
    foutgtf = paste0(outdir, 'cftc.gtf')
    test_that(paste0('buildModel::testBuildByCFTC: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildBySTMG <- function(fbams, outdir, nthr, stringtie) {
    buildModel(fbams, outdir, mode='stringtie+merging', nthreads=nthr,
               stringtie=stringtie)
    foutgtf = paste0(outdir, 'stmg.gtf')
    test_that(paste0('buildModel::testBuildBySTMG: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByCFMG <- function(fbams, outdir, nthr, cufflinks) {
    buildModel(fbams, outdir, mode='cufflinks+cuffmerge', nthreads=nthr,
               cufflinks=cufflinks)
    foutgtf = paste0(outdir, 'cfmg.gtf')
    test_that(paste0('buildModel::testBuildByCFMG: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByCF <- function(fbams, outdir, nthr, cufflinks, foutgtfs) {
    buildModel(fbams, outdir, mode='cufflinks', nthreads=nthr,
               cufflinks=cufflinks)

    lapply(foutgtfs,
           function(foutgtf) {
               test_that(paste0('buildModel::testBuildByCF: ', foutgtf),
                         expect_true( file.exists(foutgtf) ))
           })
}


testBuildByST <- function(fbams, outdir, nthr, stringtie, foutgtfs) {
    buildModel(fbams, outdir, mode='stringtie', nthreads=nthr,
               stringtie=stringtie)

    lapply(foutgtfs,
           function(foutgtf) {
               test_that(paste0('buildModel::testBuildByST: ', foutgtf),
                         expect_true( file.exists(foutgtf) ))
           })
}


testBuildByPLCF <- function(fbams, outdir, nthr, cufflinks) {
    buildModel(fbams, outdir, mode='pooling+cufflinks', nthreads=nthr,
               cufflinks=cufflinks)
    foutgtf = paste0(outdir, 'plcf.gtf')
    test_that(paste0('buildModel::testBuildByPLCF: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBuildByPLST <- function(fbams, outdir, nthr, stringtie) {
    buildModel(fbams, outdir, mode='pooling+stringtie', nthreads=nthr,
               stringtie=stringtie)
    foutgtf = paste0(outdir, 'plst.gtf')
    test_that(paste0('buildModel::testBuildByPLST: ', foutgtf),
              expect_true( file.exists(foutgtf) ))
}


testBinPLCF <- function(fbams, mode, nthr, prm, os) {
    outdir = paste0(tempdir(), '/')
    url = os2cufflinks_url(prm)[[os]]
    test_that('buildModel::testBinPLCF',
              expect_error(
                  buildModel(fbams, outdir, mode, nthr),
                  regexp=paste0('cufflinks not found: \n',
                                'It can be downloaded at ', url, "\n"),
                  ignore.case=T))
}


testBinPLST <- function(fbams, mode, nthr, prm, os) {
    outdir = paste0(tempdir(), '/')
    url = os2stringtie_url(prm)[[os]]
    test_that('buildModel::testBinPLST',
              expect_error(
                  buildModel(fbams, outdir, mode, nthr),
                  regexp=paste0('StringTie not found: \n',
                                'It can be downloaded at ', url, "\n"),
                  ignore.case=T))
}


testBinCFMG <- function(fbams, mode, nthr, prm, os) {
    outdir = paste0(tempdir(), '/')
    url = os2cufflinks_url(prm)[[os]]
    test_that('buildModel::testBinCFMG',
              expect_error(
                  buildModel(fbams, outdir, mode, nthr),
                  regexp=paste0('Cufflinks suite can be downloaded at ', url,
                                "\n"),
                  ignore.case=T))
}


testBinCFTC <- function(fbams, mode, nthr, prm, os) {
    outdir = paste0(tempdir(), '/')
    url = os2taco_url(prm)[[os]]
    test_that('buildModel::testBinCFTC',
              expect_error(
                  buildModel(fbams, outdir, mode, nthr),
                  regexp=paste0('TACO not found: \n',
                                'It can be downloaded at ', url, "\n"),
                  ignore.case=T))
}


main()
