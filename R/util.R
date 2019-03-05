#  get operation system info
#
getOS <- function() {
    os = toupper(.Platform$OS.type)
    sysinf = Sys.info()
    if ( !is.null(sysinf) ) {
        os = sysinf['sysname']
        if ( os == 'Darwin' ) {
            os = "OSX"
        } else if ( os == 'Linux' ) {
            os = "LINUX"
        } else {
            os = "WINDOWS"
        }
    } else {
        if ( grepl("^darwin", R.version$os, perl=TRUE) ) {
            os = "OSX"
        } else if ( grepl("linux-gnu", R.version$os, perl=TRUE) ) {
            os = "LINUX"
        } else {
            os = "WINDOWS"
        }
    }

    return(os)
}


# convert ori 2 strand
#
convertOri2Strand <- function(ori) {
    if ( length(setdiff(unique(ori), c( '+', '-' ))) > 0  ) {
        msg = paste0("Ori must be one of '+' or '-'\n")
        stop(msg)
    }

    strand = ifelse(ori == '+', 'plus', ifelse( ori == '-', 'minus', NA))

    return(strand)
}


# convert strand 2 ori
#
convertStrand2Ori <- function(strand) {
    if ( length(setdiff(unique(strand), c( 'plus', 'minus' ))) > 0 ) {
        msg = paste0("Strand must be one of 'plus' or 'minus'\n")
        stop(msg)
    }

    ori = ifelse(strand == 'plus', '+', ifelse(strand == 'minus', '-', NA))

    return(ori)
}


#  Get a data.table of unique combinations of chrom and ori from input GRanges
#
#' @importFrom  GenomeInfoDb  seqnames
#' @importFrom  BiocGenerics  strand
#'
getUniChromOriDt <- function(grs) {
    dt = data.table( 
        chrom = as.character(seqnames(grs)),
        ori   = as.character(strand(grs)) )
    setkey(dt, NULL)
    unidt = unique(dt, by=c('chrom', 'ori'))

    return(unidt)
}


#' @importFrom  utils download.file untar
#'
#downloadAndUntar <- function(url, dldir) {
#   tgz_name = basename(url)
#   ftgz_dest = paste0(dldir, '/', tgz_name)
#   download.file(url, ftgz_dest, quiet=FALSE)

#   untar(ftgz_dest, exdir=dldir)
#   ex_dir = paste0(dldir, '/', gsub('.tar.gz', '', tgz_name, fixed=TRUE))

#   return(ex_dir)
#}


## file.exists cannot give T/F for character(0)
fileExists <- function(file) {
    is_existed = ifelse( identical(file, character(0)), FALSE,
        ifelse( ! file.exists(file), FALSE, TRUE))
    return(is_existed)
}


#' @importFrom rtracklayer readGFF
#' @importFrom data.table data.table
#'
getDTFromGTFFile <- function(fgtf, tags) {
    dt = data.table(readGFF(
        filepath = fgtf,
        columns  = c('seqid', 'source', 'type', 'start', 'end', 'strand'),
        tags     = tags ))
    setnames(dt, c('seqid', 'type'), c('chrom', 'feature'))

    return(dt)
}


#' @importFrom utils write.table
#'
writeDT2GTFFile <- function(grdt, fout, tags=NULL, append=FALSE) {
    feature = irow = chrom = score = frame = NULL
    outdt = data.table()
    if ( nrow(grdt) > 0 ) {
        outdt = copy(grdt)
        if ( ! ('feature' %in% names(outdt)) ) {
            outdt[, feature := 'unknown']
        }
        ori_field = ifelse('source' %in% names(grdt), grdt$source[1], 'UNKNOWN')
        outdt[, `:=`(   source = ori_field,
                        score  = '.',
                        frame  = '.'     ) ]

        if ( is.null(tags) ) {
            tags = setdiff(names(outdt), c('chrom', 'source', 'feature', 
                'start', 'end', 'score', 'strand', 'frame'))
        }
        for ( infokey in tags ) {
            outdt[, eval(infokey) := paste0(infokey, ' "', get(infokey), '"')]
        }
        outdt[, irow := .I]
        outdt[, attr := paste0(paste0(.SD, collapse='; '), ';'),
                by=irow, .SDcols=tags]
        outdt = outdt[, list(chrom, source, feature, start, end, score,
                            strand, frame, attr)]
    } else {
        warning("empty data.table to write to GTF file\n")
    }

    outdir = dirname(fout)
    if ( ! file.exists(outdir) ) dir.create(outdir, recursive=TRUE)
    write.table(outdt, fout, quote=FALSE, sep="\t", col.names=FALSE, 
                row.names=FALSE, append=append)
}
