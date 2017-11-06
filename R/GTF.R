#' @export
GTF = setClass( 'GTF',
                slots = list( fgtf     = 'character',
                              origin   = 'character',
                              infokeys = 'vector',
                              exon     = 'data.table' )
)


setGeneric('getFgtf',     function(x) standardGeneric('getFgtf'))
setGeneric('getOrigin',   function(x) standardGeneric('getOrigin'))
setGeneric('getInfokeys', function(x) standardGeneric('getInfokeys'))
setGeneric('getExon',     function(x) standardGeneric('getExon'))
setGeneric('writeGTF', function(x, fout, to_append) standardGeneric('writeGTF'))
## 2nd argument to be named as 'value'
setGeneric('fgtf<-', function(x, value) standardGeneric('fgtf<-'))

setMethod('getFgtf',     'GTF', function(x) x@fgtf)
setMethod('getOrigin',   'GTF', function(x) x@origin)
setMethod('getInfokeys', 'GTF', function(x) x@infokeys)
setMethod('getExon',     'GTF', function(x) x@exon)
setReplaceMethod('fgtf', 'GTF', function(x, value) {x@fgtf = value; x})


setMethod('show', 'GTF',
    function(object) {
        cat('fgtf:',     getFgtf(object),     "\n")
        cat('infokeys:', getInfokeys(object), "\n")
        cat('origin:',   getOrigin(object),   "\n")
        cat("exon:\n")
        print(getExon(object))
    }
)


setMethod('initialize', 'GTF',
    function(.Object, fgtf=fgtf, infokeys=infokeys, ...) {
        .Object@fgtf     = fgtf
        .Object@infokeys = infokeys

        ecl = list(...)
        .Object@origin = ifelse(is.null(ecl$origin), 'UNKNOWN', ecl$origin)

        lines = readLines(fgtf)
        nskip = 0
        for ( i in seq_along(lines) ) {
            if ( ! grepl('^#', lines[i], perl=T) ) break
            nskip = nskip + 1
        }

        indt = fread(fgtf, header=F, sep="\t", colClasses=c('character',
                     'NULL', 'character', rep('integer', 2), 'NULL',
                     'character', 'NULL', 'character'), skip=nskip)
        setnames(indt, 1:6, c('chrom', 'ftr', 'start', 'end', 'strand', 'info'))
        dt = subset(indt, ftr == 'exon')

        for ( infokey in infokeys ) {
            dt[, eval(infokey) := gsub(paste0('.*', infokey, ' ([^;]+);.*'),
                                       '\\1', info, perl=T) ]
            dt[, eval(infokey) := gsub('"', '', get(infokey), fixed=T)]
        }
        dt[, `:=`( ftr=NULL, info=NULL )]

        .Object@exon = dt
        return(.Object)
    }
)


#' @importFrom GenomicRanges as.data.frame
#setMethod(
#   'initFromGRanges',
#   c('GTF', 'GRanges'),
#   function(.Object, grs) {
#       .Object@
#   }
#)


#' @importFrom data.table copy
setMethod('writeGTF', c('GTF', 'character', 'logical'),
    function(x, fout, to_append) {
        outdt = copy(getExon(x))
        ori_field = ifelse(length(getOrigin(x)) == 0, 'UNKNOWN', getOrigin(x))
        outdt[, `:=`( source  = ori_field,
                      feature = 'exon',
                      score   = '.',
                      frame   = '.'     ) ]

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
