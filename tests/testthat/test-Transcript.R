main <- function() {
    context('Transcript')

    ftgt = system.file('extdata/benchmark/tgt.tsv.gz', package='pram')
    tgtdt = data.table(read.table(ftgt, header=TRUE, sep="\t"))
    tgttr = Transcript(tgtdt)

    test_that('Transcript', {
          expect_true(validObject(tgttr))
    } )
}

main()
