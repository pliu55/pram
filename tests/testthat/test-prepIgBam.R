suppressMessages(library(GenomicAlignments))

main <- function() {
    context('prepIgBam')

    rnaseqids = c('CMPRep1', 'CMPRep2')
    seqinfo = Seqinfo( c('chr10', 'chr12') )
    iggrs = c( GRanges( 'chr10:77236000-77247000:+', seqinfo = seqinfo ),
               GRanges( 'chr12:32095000-32125000:-', seqinfo = seqinfo ),
               GRanges( 'chr12:77038000-77054000:-', seqinfo = seqinfo ) )

    outdir = paste0(tempdir(), '/')

    lapply(rnaseqids, testExtractIgBam, iggrs, outdir)
}


testExtractIgBam <- function(rnaseqid, iggrs, outdir) {
    finbam  = system.file(paste0('extdata/bam/', rnaseqid,
                                 '.sortedByCoord.raw.bam'), package='pram')
    fcmpbam = system.file(paste0('extdata/bam/', rnaseqid,
                                 '.sortedByCoord.clean.bam'), package='pram')
    foutbam = paste0(outdir, rnaseqid, '.ig.bam')

    prepIgBam(finbam, iggrs, foutbam)

    bamprm = ScanBamParam(what=c('flag', 'qname'), which=iggrs)
    cmp_alns = readGAlignments(fcmpbam, use.names=FALSE, param=bamprm)
    out_alns = readGAlignments(foutbam, use.names=FALSE, param=bamprm)

    test_that( paste0(foutbam, ' vs ', fcmpbam),
               expect_identical(cmp_alns, out_alns) )
}


main()
