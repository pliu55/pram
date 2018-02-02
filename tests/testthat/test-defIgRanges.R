suppressMessages(library(GenomicRanges))

main <- function() {
    context('defIgRanges')

    fgtf = system.file('extdata/gtf/defIgRanges_in.gtf', package='pram')
    fmm9 = system.file('extdata/chromsize/hg38.tsv.gz',  package='pram')
    genome = 'hg38'
    chroms = c('chr10', 'chr11')

    gn_iggrs = c(GRanges( 'chr10:100001-133797422', seqinfo = Seqinfo(chroms)),
                 GRanges( 'chr11:30001-39999',      seqinfo = Seqinfo(chroms)),
                 GRanges( 'chr11:100001-135086622', seqinfo = Seqinfo(chroms)))


    testGenome(fgtf, genome, gn_iggrs, chroms)

    testChromSizeFile(fgtf, fmm9, gn_iggrs, chroms)
}


testGenome <- function(fgtf, genome, gn_iggrs, chroms) {
    iggrs = defIgRanges(fgtf, genome=genome, chroms=chroms)
    seqlevels(iggrs) = seqlevels(gn_iggrs)
    test_that(paste0('genome=', genome), expect_identical(iggrs, gn_iggrs))
}


testChromSizeFile <- function(fgtf, fmm9, gn_iggrs, chroms) {
    iggrs = defIgRanges(fgtf, fchromsize=fmm9, chroms=chroms)
    seqlevels(iggrs) = seqlevels(gn_iggrs)
    test_that(paste0('fchromsize=', fmm9), expect_identical(iggrs, gn_iggrs))
}

main()
