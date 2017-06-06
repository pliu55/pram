main <- function() {
    context('Transcript')

    ftgt = system.file('extdata/benchmark/tgt.tsv', package='pram')
    tgtdt = fread(ftgt, header=T, sep="\t")
    tgttr = Transcript(tgtdt)

    test_that('Transcript', {
          expect_true(validObject(tgttr))
    } )
}

main()
