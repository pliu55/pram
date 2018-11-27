#' @import methods
#' @import data.table
#'
GTF = setClass( 'GTF',
                slots = list( fgtf     = 'character',
                              origin   = 'character',
                              infokeys = 'vector',
                              grangedt = 'data.table' )
)


setGeneric('fgtf',     function(x) standardGeneric('fgtf'))
setGeneric('origin',   function(x) standardGeneric('origin'))
setGeneric('infokeys', function(x) standardGeneric('infokeys'))
setGeneric('grangedt', function(x) standardGeneric('grangedt'))
setGeneric('writeGTF',
           function(x, fout, append) standardGeneric('writeGTF'))
setGeneric('initFromGTFFile',
           function(obj, fgtf, infokeys, ...) standardGeneric('initFromGTFFile'))
setGeneric('initFromGRanges',
           function(obj, grs) standardGeneric('initFromGRanges'))
setGeneric('initFromDataTable', function(obj, dt, infokeys, ...)
                                    standardGeneric('initFromDataTable'))

## 2nd argument to be named as 'value'
setGeneric('fgtf<-',     function(x, value) standardGeneric('fgtf<-'))
setGeneric('origin<-',   function(x, value) standardGeneric('origin<-'))
setGeneric('infokeys<-', function(x, value) standardGeneric('infokeys<-'))
setGeneric('grangedt<-', function(x, value) standardGeneric('grangedt<-'))


setMethod('fgtf',     'GTF', function(x) x@fgtf)
setMethod('origin',   'GTF', function(x) x@origin)
setMethod('infokeys', 'GTF', function(x) x@infokeys)
setMethod('grangedt', 'GTF', function(x) x@grangedt)

setReplaceMethod('fgtf',     'GTF', function(x, value) {x@fgtf = value; x})
setReplaceMethod('origin',   'GTF', function(x, value) {x@origin = value; x})
setReplaceMethod('infokeys', 'GTF', function(x, value) {x@infokeys = value; x})
setReplaceMethod('grangedt', 'GTF', function(x, value) {x@grangedt = value; x})


setMethod('show', 'GTF',
    function(object) {
        cat('fgtf:',     fgtf(object),     "\n")
        cat('origin:',   origin(object),   "\n")
        cat('infokeys:', infokeys(object), "\n")
        cat("grangedt:\n")
        print(grangedt(object))
    }
)


setMethod(
    'initialize',
    'GTF',
    function(.Object) {
        .Object@fgtf     = character()
        .Object@origin   = character()
        .Object@infokeys = vector()
        .Object@grangedt = data.table()
        return(.Object)
    }
)

#  initialize a GTF object from a GTF file with a vector of info keys to parse
#
#  @param obj a GTF object to be initialized
#  @param fgtf file name with full path to a GTF file
#  @param infokeys a vector of characters to define to-be-extracted entries in GTF file's column 9
#
#' @importFrom data.table  data.table fread setnames
#' @importFrom rtracklayer readGFF
#'
setMethod(
    'initFromGTFFile',
    c('GTF', 'character', 'vector'),
    function(obj, fgtf, infokeys, ...) {
        if ( ! file.exists(fgtf) ) stop(paste0('Cannot find :', fgtf, "\n"))

        obj = GTF()
        obj@fgtf     = fgtf
        obj@infokeys = infokeys

        ecl = list(...)
        obj@origin = ifelse(is.null(ecl$origin), 'UNKNOWN', ecl$origin)

        gtf_cols=c('seqid', 'type', 'start', 'end', 'strand')
        if ( length(infokeys) > 0 ) {
            dt = data.table(readGFF(fgtf, columns=gtf_cols, tags=infokeys))
        } else {
            dt = data.table(readGFF(fgtf, columns=gtf_cols))
        }
        setnames(dt, c('seqid', 'type'), c('chrom', 'feature'))

       #lines = readLines(fgtf)
       #nskip = 0
       #for ( line in lines ) {
       #    if ( ! grepl('^#', line, perl=T) ) break
       #    nskip = nskip + 1
       #}

       #dt = data.table()
       #if ( (length(lines) > 0) & (nskip < length(lines)) ) {
       #    dt = fread(fgtf, header=F, sep="\t", showProgress=F, skip=nskip,
       #               colClasses=c('character', 'NULL', 'character',
       #                            rep('integer', 2), 'NULL', 'character',
       #                            'NULL', 'character'))
       #    setnames(dt, 1:6,
       #             c('chrom', 'feature', 'start', 'end', 'strand', 'info'))

       #    if ( length(infokeys) > 0 ) {
       #        for ( infokey in infokeys ) {
       #            dt[, eval(infokey) :=
       #                gsub(paste0('.*', infokey, ' ([^;]+);.*'),
       #                     '\\1', info, perl=T) ]
       #            dt[, eval(infokey) := gsub('"', '', get(infokey), fixed=T)]
       #        }
       #    }
       #    dt[, info := NULL ]
       #}

        obj@grangedt = dt

        return(obj)
    }
)


#  construct a GTF object from a GenomicRanges object
#
#  @param obj a GTF object to be initialized
#  @param grs a GenomicRanges object to define ranges in GTF file
#
#' @import data.table
#' @importFrom GenomicRanges as.data.frame
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


#  construct a GTF object from a data.table object
#
#  @param obj  a GTF object to be initialized
#  @param dt   a data.table object to define genomic ranges
#  @param infokeys  a vector of characters to define to-be-extracted entries
#                   in GTF file's column 9
#
#
setMethod('initFromDataTable', c('GTF', 'data.table', 'vector'),
    function(obj, dt, infokeys, ...) {
        obj = GTF()
        obj@grangedt = dt
        infokeys(obj) = infokeys

        ecl = list(...)
        obj@origin  = ifelse(is.null(ecl$origin), 'UNKNOWN', ecl$origin)

        return(obj)
    }
)


#  write a GTF object to a GTF file
#
#  @param x a GTF object
#  @param fout a character object of GTF file name
#  @param append a boolean to inidicate if to append to a GTF file or not
#
#' @importFrom utils write.table
#
setMethod('writeGTF',
    c('GTF', 'character', 'logical'),
    function(x, fout, append) {
        grdt = grangedt(x)
        outdt = data.table()
        if ( nrow(grdt) > 0 ) {
            outdt = copy(grdt)
            if ( ! ('feature' %in% names(outdt)) ) {
                outdt[, feature := 'unknown']
            }
            ori_field = ifelse(length(origin(x)) == 0, 'UNKNOWN', origin(x))
            outdt[, `:=`( source = ori_field,
                          score  = '.',
                          frame  = '.'     ) ]

            info_keys = infokeys(x)
            for ( infokey in info_keys ) {
                outdt[, eval(infokey) := paste0(infokey, ' "', get(infokey),
                                                '"')]
            }
            outdt[, irow := .I]
            outdt[, attr := paste0(paste0(.SD, collapse='; '), ';'),
                  by=irow, .SDcols=info_keys]
            outdt = outdt[, list(chrom, source, feature, start, end, score,
                                 strand, frame, attr)]
        } else {
            warning(paste0('emptry GRange in ', fgtf(x), "\n"))
        }

        outdir = dirname(fout)
        if ( ! file.exists(outdir) ) dir.create(outdir, recursive=T)
        write.table(outdt, fout, quote=F, sep="\t", col.names=F, row.names=F,
                    append=append)
    }
)


#' parse a GTF file
#'
#' @param  fgtf       input GTF file
#' @param  info_keys  a vector of characters for attribute to be extracted
#'                    from GTF file's 9th column
#'                    Default: c('transcript_id', 'gene_id')
#'
#' @return a data.table object
#'
#' @export
#'
readGTF <- function(fgtf, info_keys=c('transcript_id', 'gene_id')) {
    gtf = new('GTF')
    gtf = initFromGTFFile(gtf, fgtf, info_keys)
    outdt = grangedt(gtf)

    return(outdt)
}
