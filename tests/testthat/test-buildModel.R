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

    testFilterBam4IG(fbams, iggrs, outdir, nthr)

    prm = new('Param')
    os = getOS()
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


testFilterBam4IG <- function(frawbams, iggrs, outdir, nthr) {
    tmpdir = paste0(outdir, 'pram_tmp/')
    if ( ! file.exists(tmpdir) ) dir.create(tmpdir, recursive=T)

    prm = new('Param')
    outdir(prm) = outdir
    tmpdir(prm) = tmpdir
    nthreads(prm) = nthr

    bamdt = filterBam4IG(frawbams, iggrs, prm)
    bamdt[, rnaseqid := tstrsplit(basename(finbam), '.', fixed=T)[[1]]]

    cleandt = data.table( rbind(
        c('chr10', '+', 'CMPRep1',
          system.file('extdata/bam/CMPRep1.sortedByCoord.clean.chr10.plus.bam',
                      package='pram')),

        c('chr10', '+', 'CMPRep2',
          system.file('extdata/bam/CMPRep2.sortedByCoord.clean.chr10.plus.bam',
                      package='pram')),

        c('chr12', '-', 'CMPRep1',
          system.file('extdata/bam/CMPRep1.sortedByCoord.clean.chr12.minus.bam',
                      package='pram')),

        c('chr12', '-', 'CMPRep2',
          system.file('extdata/bam/CMPRep2.sortedByCoord.clean.chr12.minus.bam',
                      package='pram')) ))

    setnames(cleandt, c('chrom', 'ori', 'rnaseqid', 'fcleanbam'))
    mrgdt = merge(bamdt, cleandt, by=c('chrom', 'ori', 'rnaseqid'), all=T)

    apply(mrgdt, 1, testFilterBam4IGByFile, iggrs)
}


testFilterBam4IGByFile <- function(vec, iggrs) {
   #bamprm = ScanBamParam(tag=c('HI'), what=c('flag', 'qname'), which=iggrs)
    bamprm = ScanBamParam(what=c('flag', 'qname'), which=iggrs)

    ig_alns    = readGAlignments(vec[['figbam']],    use.names=F, param=bamprm)
    clean_alns = readGAlignments(vec[['fcleanbam']], use.names=F, param=bamprm)

    test_that( paste0('buildModel::testFilterBam4IGByFile::', vec[['rnaseqid']],
                      '.', vec[['chrom']], '.', vec[['strand']]),
               expect_identical(ig_alns, clean_alns)
             )
}


main()
