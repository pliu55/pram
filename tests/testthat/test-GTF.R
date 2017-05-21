context('GTF')

nlines = 25
in_fgtf = system.file('extdata/gnc_v24_minus_chr18.gtf', package='pram')
in_origin = 'GENCODE'
in_infokeys = c('gene_id', 'transcript_id')

gtf = GTF(in_fgtf, in_infokeys, origin=in_origin)

out_fgtf = 'tmp.gtf'
writeGTF(gtf, out_fgtf, F)
new_gtf = GTF(out_fgtf, in_infokeys, origin=in_origin)
new_exondt = getExon(new_gtf)
lines = readLines(out_fgtf)

file.remove(out_fgtf)

exondt = getExon(gtf)

test_that('GTF', {
          expect_equal(nrow(exondt), nlines)
          expect_equal(names(exondt),
                       c('chrom', 'start', 'end', 'strand', in_infokeys))
          expect_identical(new_exondt, exondt)
          expect_equal(length(lines), nlines)
          }
)
