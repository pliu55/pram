#' Select transcript models
#'
#' @param  fin_gtf  Character of an input GTF file that contains
#'                  transcript models. Required to have 'transcript_id' in the
#'                  attribute column (column 9)
#'
#' @param  fout_gtf  Character of an output GTF file that contains selected
#'                   transcript models
#'
#' @param  min_n_exon  Minimium number of exons a transcript model required to
#'                     have
#'                     Default: 2
#'
#' @param  min_tr_len  Minimium length (bp) of exon(s) and intron(s) a
#'                     transcript model required to have
#'                     Default: 200
#'
#' @param  info_keys  A vector of characters defining the attributes in input
#'                    GTF file's column 9 to be saved in the output GTF file.
#'                    'transcript_id' will always be saved.
#'                    Default: c( 'transcript_id' )
#'
#' @export
#'
setGeneric('selModel',
           function(fin_gtf, fout_gtf, min_n_exon, min_tr_len, info_keys)
           standardGeneric('selModel'))

setMethod(
'selModel',
c('character', 'character', 'numeric', 'numeric', 'vector'),
function(fin_gtf, fout_gtf, min_n_exon=2, min_tr_len=200,
         info_keys = c('transcript_id') ) {
    in_gtf  = new('GTF')
    out_gtf = new('GTF')
    out_infokeys = unique(c('transcript_id', info_keys))

    in_gtf = initFromGTFFile(in_gtf, fin_gtf, infokeys=out_infokeys)
    grdt = grangedt(in_gtf)

    exondt = grdt[ feature == 'exon' ]
    dt = exondt[, list( n_exon = .N,
                        tr_len = max(end) - min(start) ), by=transcript_id]

    sel_trids = dt[ ( n_exon >= min_n_exon  ) &
                    ( tr_len >= min_tr_len ) ]$transcript_id
    sel_grdt = grdt[ transcript_id %in% sel_trids ]

    origin(out_gtf)   = origin(in_gtf)
    infokeys(out_gtf) = out_infokeys
    grangedt(out_gtf) = sel_grdt

    writeGTF(out_gtf, fout_gtf, append=F)
})
