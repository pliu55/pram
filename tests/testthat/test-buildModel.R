library(data.table)
suppressMessages(library(GenomicAlignments))

main <- function() {
    context('buildModel')

    fbams = c( system.file('extdata/bam/CMPRep1.sortedByCoord.raw.bam',
                           package='pram'),
               system.file('extdata/bam/CMPRep2.sortedByCoord.raw.bam',
                           package='pram') )

    seqinfo = Seqinfo( c('chr10', 'chr12') )
    iggrs = c( GRanges( 'chr10:77236000-77247000:+', seqinfo = seqinfo ),
               GRanges( 'chr12:32095000-32125000:-', seqinfo = seqinfo ) )

    nthr = 1

    outdir = tempdir()

    prm = new('Param')
    os = getOS()

    cufflinks = '/ua/pliu/local/cufflinks-2.2.1/cufflinks'

    buildModel(fbams, iggrs, outdir, method='pooling+cufflinks',
               cufflinks=cufflinks)
    browser()


    testBinPLCF(fbams, iggrs, outdir, 'pooling+cufflinks',   nthr, prm, os)
    testBinPLST(fbams, iggrs, outdir, 'pooling+stringtie',   nthr, prm, os)
    testBinCFMG(fbams, iggrs, outdir, 'cufflinks+cuffmerge', nthr, prm, os)
    testBinCFTC(fbams, iggrs, outdir, 'cufflinks+taco',      nthr, prm, os)
}


testBinPLCF <- function(fbams, iggrs, outdir, method, nthr, prm, os) {
    url = os2cufflinks_url(prm)[[os]]
    test_that('buildModel::testBinPLCF',
              expect_error(
                  buildModel(fbams, iggrs, outdir, method, nthr),
                  regexp=paste0('Cufflinks not found: \n',
                                'It can be downloaded at ', url, "\n"),
                  ignore.case=T))
}


testBinPLST <- function(fbams, iggrs, outdir, method, nthr, prm, os) {
    url = os2stringtie_url(prm)[[os]]
    test_that('buildModel::testBinPLST',
              expect_error(
                  buildModel(fbams, iggrs, outdir, method, nthr),
                  regexp=paste0('StringTie not found: \n',
                                'It can be downloaded at ', url, "\n"),
                  ignore.case=T))
}


testBinCFMG <- function(fbams, iggrs, outdir, method, nthr, prm, os) {
    url = os2cufflinks_url(prm)[[os]]
    test_that('buildModel::testBinCFMG',
              expect_error(
                  buildModel(fbams, iggrs, outdir, method, nthr),
                  regexp=paste0('Cuffmerge not found: \n',
                                'It can be downloaded at ', url, "\n"),
                  ignore.case=T))
}


testBinCFTC <- function(fbams, iggrs, outdir, method, nthr, prm, os) {
    url = os2taco_url(prm)[[os]]
    test_that('buildModel::testBinCFTC',
              expect_error(
                  buildModel(fbams, iggrs, outdir, method, nthr),
                  regexp=paste0('TACO not found: \n',
                                'It can be downloaded at ', url, "\n"),
                  ignore.case=T))
}


main()
