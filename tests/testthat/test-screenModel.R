library(data.table)

main <- function() {
#   context('screenModel')

#   fbeds = c( system.file('extdata/bed/GM12878POLR2AphosphoS5Rep1.bed',
#                          package='pram'),
#              system.file('extdata/bed/GM12878POLR2AphosphoS5Rep2.bed',
#                          package='pram') )

#   ftpms = c( system.file('extdata/rsem/GM12878Rep1_training.isoforms.results',
#                          package='pram'),
#              system.file('extdata/rsem/GM12878Rep2_training.isoforms.results',
#                          package='pram') )

#   fgtf_training = system.file('extdata/gtf/rf_training.gtf', package='pram')
#   fgtf_testing  = system.file('extdata/gtf/rf_testing.gtf',  package='pram')

#   tmpdir = '/tier2/deweylab/scratch/pliu/repe/pram/'
#   fgtf_out = paste0(tmpdir, 'rf.gtf')
#   if ( file.exists(fgtf_out) ) file.remove(fgtf_out)

#  #outl = trainModelClassifier(fbeds, ftpms, fgtf_training, outdir,
#  #                            expr_min_tpm=1, cv_n_folds=10, nthreads=1)

#   screenModel(fbeds, ftpms, fgtf_training, fgtf_testing, fgtf_out, tmpdir,
#               expr_min_tpm=1, cv_n_folds=10, nthreads=1)

#   test_that(paste0('screenModel:', fgtf_out),
#             expect_true( file.exists(fgtf_out) ) )
}

main()
