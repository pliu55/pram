library(data.table)

main <- function() {
    context('trainModelFilter')

    fbeds = c( system.file('extdata/bed/GM12878POLR2AphosphoS5Rep1.bed',
                           package='pram'),
               system.file('extdata/bed/GM12878POLR2AphosphoS5Rep2.bed',
                           package='pram') )

    ftpms = c( system.file('extdata/rsem/GM12878Rep1.isoforms.results',
                           package='pram'),
               system.file('extdata/rsem/GM12878Rep2.isoforms.results',
                           package='pram') )

    fgtf = system.file('extdata/gtf/rf_tr_chr6.gtf', package='pram')

   #outdir = paste0(tempdir(), '/')
    outdir = '/tier2/deweylab/scratch/pliu/repe/pram/'

    rf = trainModelFilter(fbeds, ftpms, fgtf, outdir, expr_min_tpm=1,
                          cv_n_folds=10, nthreads=1)
}

main()
