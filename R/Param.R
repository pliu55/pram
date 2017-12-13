#' save all parameters and binary urls
#'
Param = setClass('Param',
    slots = list(
        STARURL           = 'character',
        CUFFLINKSLINUXURL = 'character',
        CUFFLINKSOSXURL   = 'character',
        STRINGTIELINUXURL = 'character',
        STRINGTIEOSXURL   = 'character',
        RSEMURL           = 'character',

        STARBIN      = 'character',
        CUFFLINKSBIN = 'character',
        STRINGTIEBIN = 'character',
        RSEMBINREF   = 'character',
        RSEMBINEXPR  = 'character',

        TEMPDIR = 'character',

        MAX_YIELD_SIZE = 'numeric',

        MAX_UNI_N_DUP_ALN = 'numeric',
        MAX_MUL_N_DUP_ALN = 'numeric',

        FR1STSTRAND2MATE2FLAG = 'list',

        OS = 'character'
    ),

    prototype = list(
        STARURL = 'https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz',

        CUFFLINKSLINUXURL = 'http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz',

        CUFFLINKSOSXURL = 'http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.OSX_x86_64.tar.gz',

        STRINGTIELINUXURL = 'http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz',

        STRINGTIEOSXURL = 'http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.OSX_x86_64.tar.gz',

        RSEMURL = 'https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz',


        TEMPDIR = paste0(tempdir(), '/'),

        MAX_YIELD_SIZE = 200000000, ## overwrite filterBam's default 1M

        MAX_UNI_N_DUP_ALN = 10, ## max # of duplicated alignments on uni-frags
        MAX_MUL_N_DUP_ALN = 10, ## max # of duplicated alignments on mul-frags

        ## for fr-1ststrand
        FR1STSTRAND2MATE2FLAG = list(
            'plus' = list(
                '1stmate' = Rsamtools::scanBamFlag( isProperPair     = T,
                                                    hasUnmappedMate  = F,
                                                    isMinusStrand    = T,
                                                    isFirstMateRead  = T,
                                                    isSecondMateRead = F  ),

                '2ndmate' = Rsamtools::scanBamFlag( isProperPair     = T,
                                                    hasUnmappedMate  = F,
                                                    isMinusStrand    = F,
                                                    isFirstMateRead  = F,
                                                    isSecondMateRead = T  ) ),

            'minus' = list(
                '1stmate' = Rsamtools::scanBamFlag( isProperPair     = T,
                                                    hasUnmappedMate  = F,
                                                    isMinusStrand    = F,
                                                    isFirstMateRead  = T,
                                                    isSecondMateRead = F  ),

                '2ndmate' = Rsamtools::scanBamFlag( isProperPair     = T,
                                                    hasUnmappedMate  = F,
                                                    isMinusStrand    = T,
                                                    isFirstMateRead  = F,
                                                    isSecondMateRead = T  ) )
        )
    )
)


setGeneric('StarBin<-', function(x, value) standardGeneric('StarBin<-'))
setGeneric('getStarBin',         function(x) standardGeneric('getStarBin'))
setGeneric('getMaxYieldSize',    function(x) standardGeneric('getMaxYieldSize'))
setGeneric('getTempDir',         function(x) standardGeneric('getTempDir'))
setGeneric('getMaxUniNDupAln',  function(x) standardGeneric('getMaxUniNDupAln'))
setGeneric('getMaxMulNDupAln',  function(x) standardGeneric('getMaxMulNDupAln'))
setGeneric('getFR1stStrand2Mate2Flag',
           function(x) standardGeneric('getFR1stStrand2Mate2Flag'))

setReplaceMethod('StarBin', 'Param', function(x, value) {x@STARBIN = value; x})
setMethod('getStarBin',         'Param', function(x) x@STARBIN)
setMethod('getMaxYieldSize',    'Param', function(x) x@MAX_YIELD_SIZE)
setMethod('getTempDir',         'Param', function(x) x@TEMPDIR)
setMethod('getMaxUniNDupAln',    'Param', function(x) x@MAX_UNI_N_DUP_ALN)
setMethod('getMaxMulNDupAln',    'Param', function(x) x@MAX_MUL_N_DUP_ALN)
setMethod('getFR1stStrand2Mate2Flag', 'Param',
          function(x) x@FR1STSTRAND2MATE2FLAG)


#' @importFrom Rsamtools scanBamFlag
#'
setMethod('initialize',
    'Param',
    function(.Object) {
        return(.Object)
    }
)
