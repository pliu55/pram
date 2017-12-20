library(data.table)
suppressMessages(library(GenomicAlignments))

main <- function() {
    context('BamManager')

    fbams = c( system.file('extdata/bam/CMPRep1.sortedByCoord.raw.bam',
                           package='pram'),
               system.file('extdata/bam/CMPRep2.sortedByCoord.raw.bam',
                           package='pram') )

    seqinfo = Seqinfo( c('chr10', 'chr12') )
    iggrs = c( GRanges( 'chr10:77236000-77247000:+', seqinfo = seqinfo ),
               GRanges( 'chr12:32095000-32125000:-', seqinfo = seqinfo ) )

    prm = new('Param')
    outdir(prm)   = paste0(tempdir(), '/')
    tmpdir(prm)   = paste0(outdir(prm), 'pram_tmp/')
    nthreads(prm) = 1
    chromoridt(prm) = getUniChromOriDt(iggrs)

    if ( ! file.exists(tmpdir(prm)) ) dir.create(tmpdir(prm), recursive=T)

    testFilterParentBamByGRanges(fbams, iggrs, prm)
}


testFilterParentBamByGRanges <- function(fparentbams, iggrs, prm) {

    bams = initBy1ParentChromOri(fparentbams, iggrs, prm)

    lapply(bams, filterParentBamByGRanges, prm)

    lapply(bams, testFilterParentBamByGRangesByFile, iggrs)
}


testFilterParentBamByGRangesByFile <- function(bam, iggrs) {
   foutbam = foutbam(bam)
   fcleanbam = system.file(paste0('extdata/bam/',
                                  gsub('raw', 'clean', basename(foutbam))),
                           package='pram')

   bamprm = ScanBamParam(what=c('flag', 'qname'), which=iggrs)

   ig_alns    = readGAlignments(foutbam,   use.names=F, param=bamprm)
   clean_alns = readGAlignments(fcleanbam, use.names=F, param=bamprm)

   test_that( paste0('BamManager::testFilterParentBamByGRangesByFile::',
                     foutbam, ' vs ', fcleanbam),
              expect_identical(ig_alns, clean_alns) )
}


main()
