#' @export
Transcript = setClass('Transcript',
    slots = list(
        source = 'character',  ## where this tr is defined, e.g. GENCODE
        exon   = 'data.table', ## exons with trid and geneid
        jnc    = 'data.table', ## splic junctions with trid and geneid
        tr     = 'data.table'  ## transcript with trid and geneid
    ),

    validity = function(object) {
        exon_required_cols = c('chrom', 'start', 'end', 'strand', 'trid')
        jnc_required_cols  = c(exon_required_cols, 'njnc', 'ijnc')
        tr_required_cols   = c(exon_required_cols, 'nexon')

        valid = T
        msg = NULL
        for ( required_col in exon_required_cols ) {
            if ( ! required_col %in% names(object@exon ) ) {
                valid = F
                msg = c(msg, 'Transcript::exon ', required_col, ' is missing')
            }
        }

        for ( required_col in jnc_required_cols ) {
            if ( ! required_col %in% names(object@jnc ) ) {
                valid = F
                msg = c(msg, 'Transcript::jnc ', required_col, ' is missing')
            }
        }

        for ( required_col in tr_required_cols ) {
            if ( ! required_col %in% names(object@tr ) ) {
                valid = F
                msg = c(msg, 'Transcript::tr ', required_col, ' is missing')
            }
        }

        if ( valid ) return(valid) else return(msg)
    }
)


setMethod('show', 'Transcript',
    function(object) {
        cat('source:', object@source, "\n")
        cat("exon:\n")
        print(object@exon)
    }
)

setMethod('initialize', 'Transcript',
    function(.Object, exon, ...) {
        .Object@exon = exon
        .Object@jnc  = getTrJncFromExon(exon)
        .Object@tr   = getTrFromExon(exon)

        validObject(.Object)
        return(.Object)
    }
)


#' obtain junction info from exon
#'
#' @param exondt a data.table of exons with columns: chrom, start, end, strand,
#'                and transcript ID
#'
#' @importFrom data.table data.table
#' @importFrom data.table setcolorder
#' @importFrom IRanges IRanges
#' @importFrom IRanges gaps
#' @importFrom S4Vectors split
#' @importFrom BiocGenerics unlist
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#'
#' @return a data table of junctions with columns: chrom, start, end, strand,
#'         transcript ID, number of junctions, and junction's index in
#'         transcript
#'
getTrJncFromExon <- function(exondt) {
    ## gaps() cannot be applied to GRangesList, have to use IRangesList instead
    exonirs = IRanges( start = exondt[, start],
                       end   = exondt[, end],
                       names = exondt[, trid] )
    exonirsl = split(exonirs, names(exonirs))
    jncirs = unlist(gaps(exonirsl))
    jncdt = data.table( trid  = names(jncirs),
                        start = start(jncirs),
                        end   = end(jncirs) )
    jncdt[, `:=`( ijnc = seq_along(.I),
                  njnc = .N ), by=trid]
    strdt = data.table( chrom  = exondt[, chrom],
                        strand = exondt[, strand],
                        trid   = exondt[, trid] )
    uni_strdt = unique(strdt, by=c('strand', 'trid'))
    outdt = merge(jncdt, uni_strdt, by='trid', all.x=T)
    setcolorder(outdt, c('chrom', 'start', 'end', 'strand', 'trid', 'njnc',
                         'ijnc'))
    return(outdt)
}


#' get transcript genomic ranges from exons
#'
#' @param exondt a data table of exon with columns: chrom, start, end, strand,
#'               and transcript ID
#' @param id_col_name the column name of transcript ID, default is 'trid'
#'
#' @return a data table of transcript genomic ranges with columns: chrom, start,
#'         end, strand, number of exons, and transcript ID
#'
getTrFromExon <- function(exondt, id_col_name='trid') {
    trdt = exondt[, list( chrom  = last(chrom),
                          start  = min(start),
                          end    = max(end),
                          strand = last(strand),
                          nexon  = .N ), by=id_col_name]
    return(trdt)
}

