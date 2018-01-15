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
        }
    } else {
        if ( grepl("^darwin", R.version$os, perl=T) ) {
            os = "OSX"
        } else if ( grepl("linux-gnu", R.version$os, perl=T) ) {
            os = "LINUX"
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
    dt = data.table( chrom = as.character(seqnames(grs)),
                     ori   = as.character(strand(grs)) )
    setkey(dt, NULL)
    unidt = unique(dt, by=c('chrom', 'ori'))

    return(unidt)
}
