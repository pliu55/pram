#' @title Define intergenic genomic regions
#'
#' @param  in_gtf  An input GTF file for defining genomic coordinates of
#'                 existing
#'                 genes.  Required to have `gene_id` in the attribute column
#'                 (column 9)
#'
#' @param chromgrs  A GRanges object defining chromosome sizes.
#'
#' @param  genome  Version of the genome. Will be used when `chromgrs` is
#'                 missing. Currently supported ones are:
#'                 \itemize{
#'                     \item hg19
#'                     \item hg38
#'                     \item mm9
#'                     \item mm10
#'                 }
#'                 All the above genomes have sizes for all chromosomes
#'                 including random and alt ones.
#'                 Default: NULL
#'
#' @param  fchromsize  Name of a file defining chromosome sizes. Will be used
#'                     when `chromgrs` and `genome` are missing.
#'                     It can be downloaded from
#'                     UCSC, e.g. for hg19, http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/chromInfo.txt.gz
#'                     Required to have at least two tab-delimited columns
#'                     without any header:
#'                     \enumerate{
#'                        \item chromosome name, e.g. chr1
#'                        \item chromosome length, e.g. 249250621
#'                     }
#'                     Both uncompressed and gzipped files are supported.
#'                     Default: NULL
#'
#' @param  radius  Region length (bp) of gene's upstream and downstream to be
#'                 excluded from intergenic region.
#'                 Default: 10,000
#'
#' @param  feat  Feature in the GTF file (column 3) to-be-used for defining
#'               genic region.
#'               Default: exon
#'
#' @param  chroms  A vector of chromosomes names to define intergenic regions.
#'                 e.g. c('chr10', 'chr11')
#'                 Default: NULL
#'
#' @return  a GRanges object of intergenic regions
#'
#' @importFrom  GenomicRanges  reduce
#'
#' @export
#'
#' @examples
#' fgtf = system.file('extdata/gtf/defIgRanges_in.gtf', package='pram')
#'
#' defIgRanges(fgtf, genome='hg38')
#'
defIgRanges <- function(in_gtf, chromgrs, genome=NULL, fchromsize=NULL,
                        radius=1e+4, feat='exon', chroms=NULL){
    feature = chrom = gene_id = gn_start = gn_end = NULL
    fgtf = in_gtf
    if ( missing(chromgrs) ) {
        chromgrs = getChromGRanges(genome, fchromsize, chroms)
    }

    gtf = new('GTF')
    gtf = initFromGTFFile(gtf, fgtf, infokeys=c('gene_id'))
    grdt = grangedt(gtf)
    seldt = grdt[ feature == feat ]
    if ( ! is.null(chroms) ) {
        seldt = seldt[ chrom %in% chroms ]
    }
    gndt = seldt[, list( chrom    = unique(chrom),
                         gn_start = min(start),
                         gn_end   = max(end) ), by=gene_id]

    gndt[, `:=`( start  = gn_start - radius,
                 end    = gn_end + radius,
                 strand = '*' )]

    gngrs = makeGRangesFromDataFrame(gndt, keep.extra.columns=F)
    gngrs = reduce(gngrs)

    iggrs = setdiff(chromgrs, gngrs)

    return(iggrs)
}


getChromGRanges <- function(genome, fchromsize, chroms) {
    chrom = NULL
    chromdt = NULL
    avail_genomes = c('hg19', 'hg38', 'mm10', 'mm9')

    if ( is.null(genome) & is.null(fchromsize) ) {
        stop("either [genome] or [fchromsize] needs to be defined\n")
    } else if ( ! is.null(genome) ) {
        genome = tolower(genome)
        if ( genome %in% avail_genomes ) {
            fin = system.file(paste0('extdata/chromsize/', genome, '.tsv.gz'),
                              package='pram')
            chromdt = readChromSize(fin)
        } else  {
            msg = paste0('defIg: genome ', genome, " is not implemented.\n",
                         "Supported genomes are ",
                         paste(avail_genomes, collapse=','), "\n")
            stop(msg)
        }
    } else if ( is.null(genome) & (! is.null(fchromsize)) ) {
        chromdt = readChromSize(fchromsize)
    }

    outdt = chromdt
    if ( ! is.null(chroms) ) {
        outdt = chromdt[ chrom %in% chroms ]
    }
    outdt[, chrom := as.character(chrom) ] ## avoid chrom w/ 0 length
    outgrs = makeGRangesFromDataFrame(outdt, keep.extra.columns=F)

    return(outgrs)
}


#' @importFrom  data.table  data.table fread
#' @importFrom  utils       read.table
readChromSize <- function(fin) {
    V1 = V2 = NULL
    `.` = function(...) NULL
    dt = data.table(read.table(fin, header=F, sep="\t"))

    if ( nrow(dt) == 0 ) {
        stop(paste0('Fail to read or file is empty: ', fin, "\n"))
    }
    outdt = dt[, .(V1, V2)]
    setnames(outdt, c('chrom', 'end'))
    outdt[, start := 1]

    return(outdt)
}
