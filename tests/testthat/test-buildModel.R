suppressMessages(library(GenomicRanges))

main <- function() {
    context('buildModel')

    test()
}


test <- function() {
    fbam1 = system.file('extdata/bam/CMPRep1.sortedByCoord.bam', package='pram')
    fbam2 = system.file('extdata/bam/CMPRep2.sortedByCoord.bam', package='pram')

    fbams = c(fbam1, fbam2)
    seqinfo = Seqinfo( c('chr10', 'chr12') )
    tgtgrs = c( GRanges( 'chr10:77236000-77247000:+', seqinfo = seqinfo ),
                GRanges( 'chr12:32095000-32125000:-', seqinfo = seqinfo ) )
    outdir = tempdir()
    nthr = 2

    buildModel(fbams, tgtgrs, outdir, nthr)
}


#   test_that(paste0('evalModel::testBenchmarkMethod::', method), {
#       expect_equal(indi_jnc_tp, evaldt[feat=='indi_jnc', ntp])
#       expect_equal(indi_jnc_fn, evaldt[feat=='indi_jnc', nfn])


main()
