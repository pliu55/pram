library(data.table)
library(tools)
suppressMessages(library(GenomicAlignments))

main <- function() {
    context('prepIgBam')

    fbams = c( system.file('extdata/bam/CMPRep1.sortedByCoord.raw.bam',
                           package='pram'),
               system.file('extdata/bam/CMPRep2.sortedByCoord.raw.bam',
                           package='pram') )

    seqinfo = Seqinfo( c('chr10', 'chr12') )
    iggrs = c( GRanges( 'chr10:77236000-77247000:+', seqinfo = seqinfo ),
               GRanges( 'chr12:32095000-32125000:-', seqinfo = seqinfo ),
               GRanges( 'chr12:77038000-77054000:-', seqinfo = seqinfo ) )

    outdir = paste0(tempdir(), '/')
    nthr = 2

    prepIgBam(fbams, iggrs, outdir, nthreads=nthr)
    lapply(fbams, testPrepIgBam, iggrs, outdir)
}


testPrepIgBam <- function(finbam, iggrs, outdir) {
    foutbam = paste0(outdir, file_path_sans_ext(basename(finbam)), '.ig.bam')
    fcleanbam = gsub('raw', 'clean', finbam)

    bamprm = ScanBamParam(what=c('flag', 'qname'), which=iggrs)
    out_alns   = readGAlignments(foutbam,   use.names=F, param=bamprm)
    clean_alns = readGAlignments(fcleanbam, use.names=F, param=bamprm)

    test_that( paste0('BamManager::testPrepIgBam:', foutbam, ' vs ', fcleanbam),
               expect_identical(out_alns, clean_alns) )
}

main()
