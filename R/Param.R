#' save all parameters and binary urls
#'
Param = setClass('Param',
    slots = list(
        STAR_URL            = 'character',
        CUFFLINKS_LINUX_URL = 'character',
        CUFFLINKS_OSX_URL   = 'character',
        STRINGTIE_LINUX_URL = 'character',
        STRINGTIE_OSX_URL   = 'character',
        RSEM_URL            = 'character',

        STAR_BIN      = 'character',
        CUFFLINKS_BIN = 'character',
        STRINGTIE_BIN = 'character',
        RSEM_BIN_REF  = 'character',
        RSEM_BIN_EXPR = 'character',

        OUT_DIR  = 'character',
        TMP_DIR = 'character',

        NTHREADS = 'numeric',

        LIB_TYPE = 'character',
        MIN_ISOFORM_FRACTION    = 'numeric',
        MAX_MULTIREAD_FRACTION  = 'numeric',
        MIN_FRAGS_PER_TRANSFRAG = 'numeric',
        MIN_TR_FPKM_TO_INCLUDE  = 'numeric',
        MIN_TR_TPM_TO_INCLUDE   = 'numeric',


        MAX_YIELD_SIZE = 'numeric',

        MAX_UNI_N_DUP_ALN = 'numeric',
        MAX_MUL_N_DUP_ALN = 'numeric',

        FR1STSTRAND2MATE2FLAG = 'list',

        OS = 'character'
    ),

    prototype = list(
        STAR_URL = 'https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz',

        CUFFLINKS_LINUX_URL = 'http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz',

        CUFFLINKS_OSX_URL = 'http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.OSX_x86_64.tar.gz',

        STRINGTIE_LINUX_URL = 'http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz',

        STRINGTIE_OSX_URL = 'http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.OSX_x86_64.tar.gz',

        RSEM_URL = 'https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz',


        NTHREADS = 1,


        ### for model building
        LIB_TYPE = 'fr-firststrand', ## paired-end, ENCODE RNA-seq
                                    ## fr read1 -> <- read2
                                    ## ff read1 -> -> read2
        MIN_ISOFORM_FRACTION    = 0.1, # abundance cutoff to suppress Tr
        MAX_MULTIREAD_FRACTION  = 1.0, # max frac of allowed multiread per tr
        MIN_FRAGS_PER_TRANSFRAG = 1,   # min # of fragment to report a transfrag
        MIN_TR_FPKM_TO_INCLUDE  = 0.0,  ## StringTie-merge
        MIN_TR_TPM_TO_INCLUDE   = 0.0,  ## StringTie-merge, TACO


        MAX_YIELD_SIZE = 200000000, ## overwrite filterBam's default 1M

        ### for extracting alignments
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
                                                    isSecondMateRead = T  ) ) )
    )
)


## need to be 'value', other names won't work
setGeneric('starbin<-',  function(x, value) standardGeneric('starbin<-'))
setGeneric('nthreads<-', function(x, value) standardGeneric('nthreads<-'))
setGeneric('outdir<-',   function(x, value) standardGeneric('outdir<-'))
setGeneric('tmpdir<-',   function(x, value) standardGeneric('tmpdir<-'))
setGeneric('cufflinks<-',   function(x, value) standardGeneric('cufflinks<-'))
setGeneric('stringtie<-',   function(x, value) standardGeneric('stringtie<-'))
setGeneric('starbin',       function(x) standardGeneric('starbin'))
setGeneric('cufflinks',     function(x) standardGeneric('cufflinks'))
setGeneric('stringtie',     function(x) standardGeneric('stringtie'))
setGeneric('maxyieldsize',  function(x) standardGeneric('maxyieldsize'))
setGeneric('outdir',        function(x) standardGeneric('outdir'))
setGeneric('tmpdir',        function(x) standardGeneric('tmpdir'))
setGeneric('nthreads',      function(x) standardGeneric('nthreads'))
setGeneric('maxunindupaln', function(x) standardGeneric('maxunindupaln'))
setGeneric('maxmulndupaln', function(x) standardGeneric('maxmulndupaln'))
setGeneric('fr1ststrand2mate2flag',
           function(x) standardGeneric('fr1ststrand2mate2flag'))
setGeneric('libtype', function(x) standardGeneric('libtype'))
setGeneric('minisoformfraction',
           function(x) standardGeneric('minisoformfraction'))
setGeneric('maxmultireadfraction',
           function(x) standardGeneric('maxmultireadfraction'))
setGeneric('minfragspertransfrag',
           function(x) standardGeneric('minfragspertransfrag'))
setGeneric('mintrfpkmtoinclude',
           function(x) standardGeneric('mintrfpkmtoinclude'))
setGeneric('mintrtpmtoinclude',
           function(x) standardGeneric('mintrtpmtoinclude'))


setReplaceMethod('starbin',  'Param', function(x, value) {x@STAR_BIN=value; x})
setReplaceMethod('nthreads', 'Param', function(x, value) {x@NTHREADS=value; x})
setReplaceMethod('outdir',   'Param', function(x, value) {x@OUT_DIR=value; x})
setReplaceMethod('tmpdir',   'Param', function(x, value) {x@TMP_DIR=value; x})
setReplaceMethod('cufflinks', 'Param',
                 function(x, value) {x@CUFFLINKS_BIN=value; x})
setReplaceMethod('stringtie', 'Param',
                 function(x, value) {x@STRINGTIE_BIN=value; x})
setMethod('starbin',       'Param', function(x) x@STAR_BIN)
setMethod('cufflinks',     'Param', function(x) x@CUFFLINKS_BIN)
setMethod('maxyieldsize',  'Param', function(x) x@MAX_YIELD_SIZE)
setMethod('outdir',        'Param', function(x) x@OUT_DIR)
setMethod('tmpdir',        'Param', function(x) x@TMP_DIR)
setMethod('nthreads',      'Param', function(x) x@NTHREADS)
setMethod('maxunindupaln', 'Param', function(x) x@MAX_UNI_N_DUP_ALN)
setMethod('maxmulndupaln', 'Param', function(x) x@MAX_MUL_N_DUP_ALN)
setMethod('fr1ststrand2mate2flag', 'Param',
          function(x) x@FR1STSTRAND2MATE2FLAG)
setMethod('libtype', 'Param', function(x) x@LIB_TYPE)
setMethod('minisoformfraction', 'Param', function(x) x@MIN_ISOFORM_FRACTION)
setMethod('maxmultireadfraction', 'Param',
          function(x) x@MAX_MULTIREAD_FRACTION)
setMethod('minfragspertransfrag', 'Param',
          function(x) x@MIN_FRAGS_PER_TRANSFRAG)
setMethod('mintrfpkmtoinclude', 'Param',
          function(x) x@MIN_TR_FPKM_TO_INCLUDE)
setMethod('mintrtpmtoinclude', 'Param', function(x) x@MIN_TR_TPM_TO_INCLUDE)


#' @importFrom Rsamtools scanBamFlag
#'
setMethod('initialize',
    'Param',
    function(.Object) {
        return(.Object)
    }
)
