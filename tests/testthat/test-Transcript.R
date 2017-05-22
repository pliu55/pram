main <- function() {
    context('Transcript')

    data(tgtexondt)

    tr = Transcript(tgtexondt)

    test_that('Transcript', {
          expect_true(validObject(tr))
    } )
}

main()
