library(data.table)

main <- function() {
    context('runPRAM')

    pred_out_gtf   = tempfile(pattern='runPRAM_pred.',   fileext='.gtf')
    screen_out_gtf = tempfile(pattern='runPRAM_screen.', fileext='.gtf')

    in_gtf = system.file('extdata/demo/in.gtf', package='pram')

    in_bamv = c( system.file('extdata/demo/SZP.bam', package='pram'),
                 system.file('extdata/demo/TLC.bam', package='pram') )

    if ( ( grepl('biostat.wisc.edu', Sys.info()[['nodename']], fixed=TRUE) &
           ( Sys.info()[['user']] == 'pliu' ) ) |
         ( ( Sys.info()[['sysname']] == 'Darwin' ) & 
           ( Sys.info()[['user']] == 'peng' ) &
           ( file.exists('/ua/pliu/repe/pram') ) ) ) {
        testPred(in_gtf, in_bamv, pred_out_gtf)
    }
}


testPred <- function(in_gtf, in_bamv, out_gtf) {
    st = ''
    if ( getOS() == 'LINUX' ) {
        st = '/ua/pliu/local/stringtie-1.3.3/stringtie'
    } else if ( getOS() == 'OSX' ) {
        st = '/ua/pliu/local/osx/stringtie-1.3.3b/stringtie'
    }
    runPRAM(in_gtf, in_bamv, out_gtf, method='plst', stringtie=st)
    test_that('runPRAM::testPred',
              expect_true( file.exists(out_gtf)))
}


main()
