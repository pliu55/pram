main <- function() {
    context('ModelGTF')

    testSTMG()
    testSTPL()

    testCFMG()
    testCFPL()

    testTACO()
}

testSTMG <- function() {
    nlines = 17
    fgtf = system.file('extdata/stmg_minus_chr18.gtf', package='pram')
    origin = 'StringTieMerge'
    in_infokeys = c('gene_id', 'transcript_id')

    gtf = ModelGTF(fgtf, in_infokeys, origin, 'stmg')
    exondt = getExon(gtf)

    test_that( 'StringTieMerge GTF', {
        expect_equal(getOrigin(gtf), origin)
        expect_equal(nrow(exondt), nlines)
        expect_true(all(grepl('^stmg_minus_chr18', exondt[, gene_id], perl=T)))
    })
}


testSTPL <- function() {
    nlines = 18
    fgtf = system.file('extdata/stpl_minus_chr18.gtf', package='pram')
    origin = 'StringTiePool'
    in_infokeys = c('gene_id', 'transcript_id')

    gtf = ModelGTF(fgtf, in_infokeys, origin, 'stpl')
    exondt = getExon(gtf)

    test_that( 'StringTiePool GTF', {
        expect_equal(getOrigin(gtf), origin)
        expect_equal(nrow(exondt), nlines)
        expect_true(all(grepl('^stpl_minus_chr18', exondt[, gene_id], perl=T)))
    })
}


testCFMG <- function() {
    nlines = 17
    fgtf = system.file('extdata/cfmg_minus_chr18.gtf', package='pram')
    origin = 'CufflinksMerge'
    in_infokeys = c('gene_id', 'transcript_id')

    gtf = ModelGTF(fgtf, in_infokeys, origin=origin, 'cfmg')
    exondt = getExon(gtf)

    test_that( 'CufflinksMerge GTF', {
        expect_equal(getOrigin(gtf), origin)
        expect_equal(nrow(exondt), nlines)
        expect_true(all(grepl('^cfmg_minus_chr18', exondt[, gene_id], perl=T)))
    })
}


testCFPL <- function() {
    nlines = 19
    fgtf = system.file('extdata/cfpl_minus_chr18.gtf', package='pram')
    origin = 'CufflinksPool'
    model_method = 'cfpl'
    in_infokeys = c('gene_id', 'transcript_id')

    gtf = ModelGTF(fgtf, in_infokeys, origin, model_method)
    exondt = getExon(gtf)

    test_that( 'CufflinksPool GTF', {
        expect_equal(getOrigin(gtf), origin)
        expect_equal(nrow(exondt), nlines)
        expect_true(all(grepl('^cfpl_minus_chr18', exondt[, gene_id], perl=T)))
    })
}


testTACO <- function() {
    nlines = 19
    fgtf = system.file('extdata/taco_minus_chr18.gtf', package='pram')
    origin = 'TACO'
    model_method = 'taco'
    in_infokeys = c('gene_id', 'transcript_id')

    gtf = ModelGTF(fgtf, in_infokeys, origin, model_method)
    exondt = getExon(gtf)

    test_that( 'TACO GTF', {
        expect_equal(getOrigin(gtf), origin)
        expect_equal(nrow(exondt), nlines)
        expect_true(all(grepl('^taco_minus_chr18', exondt[, gene_id], perl=T)))
    })
}


main()
