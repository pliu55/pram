library(data.table)

main <- function() {
    context('runPRAM')

    pred_out_gtf   = tempfile(pattern='runPRAM_pred.',   fileext='.gtf')
    screen_out_gtf = tempfile(pattern='runPRAM_screen.', fileext='.gtf')

    in_gtf = system.file('extdata/demo/in.gtf', package='pram')

    in_bamv = c( system.file('extdata/demo/SZP.bam', package='pram'),
                 system.file('extdata/demo/TLC.bam', package='pram') )

  # in_bedv = c( system.file('extdata/demo/H3K79me2.bed.gz', package='pram'),
  #              system.file('extdata/demo/POLR2.bed.gz',    package='pram') )

  # training_tpms = c( system.file('extdata/demo/AED1.isoforms.results',
  #                                package='pram'),
  #                    system.file('extdata/demo/AED2.isoforms.results',
  #                                package='pram') )

  # training_gtf = system.file('extdata/demo/training.gtf', package='pram')

    if ( ( grepl('biostat.wisc.edu', Sys.info()[['nodename']], fixed=TRUE) &
           ( Sys.info()[['user']] == 'pliu' ) ) |
         ( ( Sys.info()[['sysname']] == 'Darwin' ) & 
           ( Sys.info()[['user']] == 'peng' ) ) ) {
        testPred(in_gtf, in_bamv, pred_out_gtf)
    }

  # testPredScreen(in_gtf, in_bamv, screen_out_gtf, in_bedv, training_tpms,
  #                training_gtf)
}


testPred <- function(in_gtf, in_bamv, out_gtf) {
    runPRAM(in_gtf, in_bamv, out_gtf)
    test_that('runPRAM::testPred',
              expect_true( file.exists(out_gtf)))
}


# testPredScreen <- function(in_gtf, in_bamv, out_gtf, in_bedv, training_tpms,
#                          training_gtf) {
#   runPRAM(in_gtf, in_bamv, out_gtf, in_bedv, training_tpms, training_gtf)
#   test_that('runPRAM::testPredScreen',
#             expect_true( file.exists(out_gtf)))
# }


main()
