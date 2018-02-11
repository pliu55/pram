main <- function() {
    context('selModel')

    fin_gtf  = system.file('extdata/gtf/selModel_in.gtf',  package='pram')
    fcmp_gtf = system.file('extdata/gtf/selModel_out.gtf', package='pram')
    fout_gtf = paste0(tempdir(), '/pram_selModel_out.gtf')

    selModel(fin_gtf, fout_gtf, min_n_exon=2, min_tr_len=200,
             info_keys=c('transcript_id', 'gene_id'))

    cmp_lines = readLines(fcmp_gtf)
    out_lines = readLines(fout_gtf)

    test_that(paste0(fout_gtf, ' vs ', fcmp_gtf),
              expect_identical(cmp_lines, out_lines))
}

main()
