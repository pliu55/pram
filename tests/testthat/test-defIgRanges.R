suppressMessages(library(GenomicRanges))

main <- function() {
    context('defIgRanges')

    fgtf = system.file('extdata/gtf/defIgRanges_in.gtf', package='pram')
    fmm9 = system.file('extdata/chromsize/hg38.tsv.gz',  package='pram')
    genome = 'hg38'

    seqinfo = Seqinfo( c('chr10', 'chr11') )
    gn_iggrs = c( GRanges( 'chr10:100001-133797422', seqinfo = seqinfo),
                  GRanges( 'chr11:30001-39999',      seqinfo = seqinfo),
                  GRanges( 'chr11:100001-135086622', seqinfo = seqinfo) )

    testGenome(fgtf, genome, gn_iggrs)

    testChromSizeFile(fgtf, fmm9, gn_iggrs)
}


testGenome <- function(fgtf, genome, gn_iggrs) {
    iggrs = defIgRanges(fgtf, genome=genome)
    int_grs = intersect(iggrs, gn_iggrs)
    seqlevels(int_grs) = seqlevels(gn_iggrs)
    test_that(paste0('genome=', genome), expect_identical(int_grs, gn_iggrs))
}


testChromSizeFile <- function(fgtf, fmm9, gn_iggrs) {
    iggrs = defIgRanges(fgtf, fchromsize=fmm9)
    int_grs = intersect(iggrs, gn_iggrs)
    seqlevels(int_grs) = seqlevels(gn_iggrs)
    test_that(paste0('fchromsize=', fmm9), expect_identical(int_grs, gn_iggrs))
}

main()
