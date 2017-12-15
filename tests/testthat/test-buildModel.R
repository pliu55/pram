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

    ## to-be-changed later
    outdir = '/tier2/deweylab/scratch/pliu/repe/pram/'
    cufflinks = '/ua/pliu/local/cufflinks-2.2.1/cufflinks'
    stringtie = '/ua/pliu/local/stringtie-1.3.3/stringtie'

    testFilterBam4IG(fbams, iggrs, outdir, nthr)

    testBuildModel(fbams, iggrs, outdir, 'pooling+cufflinks', nthr,
                   cufflinks=cufflinks, stringtie='')
}


testBuildModel <- function(fbams, iggrs, outdir, method, nthr, cufflinks,
                           stringtie) {
    outdt = buildModel(fbams, iggrs, outdir, method, nthr, cufflinks, stringtie)
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
