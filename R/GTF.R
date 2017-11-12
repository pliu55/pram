#' @import methods
#' @import data.table
#'
GTF = setClass( 'GTF',
                slots = list( fgtf     = 'character',
                              origin   = 'character',
                              infokeys = 'vector',
                              grangedt = 'data.table' )
)


setGeneric('getFgtf',     function(x) standardGeneric('getFgtf'))
setGeneric('getOrigin',   function(x) standardGeneric('getOrigin'))
setGeneric('getInfokeys', function(x) standardGeneric('getInfokeys'))
setGeneric('getGrangedt', function(x) standardGeneric('getGrangedt'))
setGeneric('writeGTF',
           function(x, fout, to_append) standardGeneric('writeGTF'))
setGeneric('initFromGTFFile',
           function(obj, fgtf, infokeys, ...) standardGeneric('initFromGTFFile'))
setGeneric('initFromGRanges',
           function(obj, grs) standardGeneric('initFromGRanges'))

## 2nd argument to be named as 'value'
setGeneric('fgtf<-', function(x, value) standardGeneric('fgtf<-'))


setMethod('getFgtf',     'GTF', function(x) x@fgtf)
setMethod('getOrigin',   'GTF', function(x) x@origin)
setMethod('getInfokeys', 'GTF', function(x) x@infokeys)
setMethod('getGrangedt', 'GTF', function(x) x@grangedt)
setReplaceMethod('fgtf', 'GTF', function(x, value) {x@fgtf = value; x})


setMethod('show', 'GTF',
    function(object) {
        cat('fgtf:',     getFgtf(object),     "\n")
        cat('origin:',   getOrigin(object),   "\n")
        cat('infokeys:', getInfokeys(object), "\n")
        cat("granges:\n")
        print(getGrangedt(object))
    }
)


setMethod(
    'initialize',
    signature('GTF'),
    function(.Object) {
        .Object@fgtf     = character()
        .Object@origin   = character()
        .Object@infokeys = vector()
        .Object@grangedt = data.table()
        return(.Object)
    }
)

#' initialize a GTF object from a GTF file with a vector of info keys to parse
#'
#' @param obj a GTF object to be initialized
#' @param fgtf file name with full path to a GTF file
#' @param infokeys a vector of characters to define to-be-extracted entries in GTF file's column 9
#'
#' @import data.table
#'
#' @export
#'
setMethod(
    'initFromGTFFile',
    signature('GTF', 'character', 'vector'),
    function(obj, fgtf, infokeys, ...) {
        obj = GTF()
        obj@fgtf     = fgtf
        obj@infokeys = infokeys

        ecl = list(...)
        obj@origin = ifelse(is.null(ecl$origin), 'UNKNOWN', ecl$origin)

        lines = readLines(fgtf)
        nskip = 0
        for ( line in lines ) {
            if ( ! grepl('^#', line, perl=T) ) break
            nskip = nskip + 1
        }

        dt = fread(fgtf, header=F, sep="\t", colClasses=c('character',
                   'NULL', 'character', rep('integer', 2), 'NULL',
                   'character', 'NULL', 'character'), skip=nskip)
        setnames(dt, 1:6,
                 c('chrom', 'feature', 'start', 'end', 'strand', 'info'))
       #dt = subset(indt, ftr == 'exon')

        if ( length(infokeys) > 0 ) {
            for ( infokey in infokeys ) {
                dt[, eval(infokey) := gsub(paste0('.*', infokey, ' ([^;]+);.*'),
                                           '\\1', info, perl=T) ]
                dt[, eval(infokey) := gsub('"', '', get(infokey), fixed=T)]
            }
        }
        dt[, info := NULL ]

        obj@grangedt = dt
        return(obj)
    }
)


#' construct a GTF object from a GenomicRanges object
#'
#' @param obj a GTF object to be initialized
#' @param grs a GenomicRanges object to define ranges in GTF file
#'
#' @import data.table
#' @importFrom GenomicRanges as.data.frame
#'
#' @export
#'
setMethod(
    'initFromGRanges',
    c('GTF', 'GRanges'),
    function(obj, grs) {
        obj = GTF()
        dt = data.table(as.data.frame(grs))
        setnames(dt, 'seqnames', 'chrom')
        dt[, width := NULL]
        obj@grangedt = dt

        return(obj)
    }
)


#' write a GTF object to a GTF file
#'
#' @param x a GTF object
#' @param fout a character object of GTF file name
#' @param to_append a boolean to inidicate if to append to a GTF file or not
#'
#' @export
#'
setMethod('writeGTF',
    c('GTF', 'character', 'logical'),
    function(x, fout, to_append) {
        outdt = copy(getGrangedt(x))
        if ( ! ('feature' %in% names(outdt)) ) {
            outdt[, feature := 'unknown']
        }
        ori_field = ifelse(length(getOrigin(x)) == 0, 'UNKNOWN', getOrigin(x))
        outdt[, `:=`( source = ori_field,
                      score  = '.',
                      frame  = '.'     ) ]

        info_keys = getInfokeys(x)
        for ( infokey in info_keys ) {
            outdt[, eval(infokey) := paste0(infokey, ' "', get(infokey), '"')]
        }
        outdt[, irow := .I]
        outdt[, attr := paste0(paste0(.SD, collapse='; '), ';'),
              by=irow, .SDcols=info_keys]
        outdt = outdt[, list(chrom, source, feature, start, end, score, strand,
                             frame, attr)]
        write.table(outdt, fout, quote=F, sep="\t", col.names=F, row.names=F,
                    append=to_append)
    }
)



readGTF4Transcripts <- function(fin) {
  dt <- readGTF4ExonGeneIDTrID(fin)
  trdt <- dt[, list( chrom  = last(chrom),
                     start  = min(start),
                     end    = max(end),
                     strand = last(strand),
                     geneid = last(geneid) ), by=trid ]
  return(trdt)
}


readGTF4Genes <- function(fin) {
  dt <- readGTF4ExonGeneIDTrID(fin)
  genedt <- dt[, list( chrom  = last(chrom),
                       start  = min(start),
                       end    = max(end),
                       strand = last(strand),
                       trids  = paste0(unique(trid), collapse=',')
                     ), by=geneid]
  return(genedt)
}
