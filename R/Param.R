#' save all parameters and binary urls
#'
Param = setClass('Param',
    slots = list(
        OS2CUFFLINKS_URL = 'list',
        OS2STRINGTIE_URL = 'list',
        OS2TACO_URL      = 'list',

        STAR_BIN      = 'character',
        CUFFLINKS_BIN = 'character',
        CUFFMERGE_BIN = 'character',
        STRINGTIE_BIN = 'character',
        TACO_BIN      = 'character',

        OUT_DIR  = 'character',
        TMP_DIR = 'character',

        CHROM_ORI_DT = 'data.table',

        GTF_INFO_KEYS = 'vector',

        NTHREADS = 'numeric',

        CUFFLINKS_LIB_TYPE = 'character',
        STRINGTIE_LIB_TYPE = 'character',

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
       #STAR_URL = 'https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz',

       #RSEM_URL = 'https://github.com/deweylab/RSEM/archive/v1.3.0.tar.gz',

        OS2CUFFLINKS_URL = list(
            'LINUX' = 'http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz',

            'OSX' = 'http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.OSX_x86_64.tar.gz' ),

        OS2STRINGTIE_URL = list(
            'LINUX' = 'http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz',

            'OSX' = 'http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.OSX_x86_64.tar.gz' ),

        OS2TACO_URL = list(
            'LINUX' = 'https://github.com/tacorna/taco/releases/download/v0.7.0/taco-v0.7.0.Linux_x86_64.tar.gz',

            'OSX' = 'https://github.com/tacorna/taco/releases/download/v0.7.0/taco-v0.7.0.OSX_x86_64.tar.gz' ),


        GTF_INFO_KEYS = c('gene_id', 'transcript_id'),

        NTHREADS = 1,


        ### for model building
        CUFFLINKS_LIB_TYPE = 'fr-firststrand', ## paired-end, ENCODE RNA-seq
                                    ## fr read1 -> <- read2
                                    ## ff read1 -> -> read2
        STRINGTIE_LIB_TYPE = '--rf',


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
setGeneric('nthreads<-',    function(x, value) standardGeneric('nthreads<-'))
setGeneric('outdir<-',      function(x, value) standardGeneric('outdir<-'))
setGeneric('tmpdir<-',      function(x, value) standardGeneric('tmpdir<-'))
setGeneric('chromoridt<-',  function(x, value) standardGeneric('chromoridt<-'))
setGeneric('cufflinks<-',   function(x, value) standardGeneric('cufflinks<-'))
setGeneric('cuffmerge<-',   function(x, value) standardGeneric('cuffmerge<-'))
setGeneric('stringtie<-',   function(x, value) standardGeneric('stringtie<-'))
setGeneric('taco<-',        function(x, value) standardGeneric('taco<-'))
setGeneric('cufflinks',     function(x) standardGeneric('cufflinks'))
setGeneric('cuffmerge',     function(x) standardGeneric('cuffmerge'))
setGeneric('stringtie',     function(x) standardGeneric('stringtie'))
setGeneric('taco',          function(x) standardGeneric('taco'))
setGeneric('maxyieldsize',  function(x) standardGeneric('maxyieldsize'))
setGeneric('outdir',        function(x) standardGeneric('outdir'))
setGeneric('tmpdir',        function(x) standardGeneric('tmpdir'))
setGeneric('chromoridt',    function(x) standardGeneric('chromoridt'))
setGeneric('gtfinfokeys',   function(x) standardGeneric('gtfinfokeys'))
setGeneric('nthreads',      function(x) standardGeneric('nthreads'))
setGeneric('maxunindupaln', function(x) standardGeneric('maxunindupaln'))
setGeneric('maxmulndupaln', function(x) standardGeneric('maxmulndupaln'))
setGeneric('fr1ststrand2mate2flag',
           function(x) standardGeneric('fr1ststrand2mate2flag'))
setGeneric('cufflinkslibtype', function(x) standardGeneric('cufflinkslibtype'))
setGeneric('stringtielibtype', function(x) standardGeneric('stringtielibtype'))
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
setGeneric('os2cufflinks_url', function(x) standardGeneric('os2cufflinks_url'))
setGeneric('os2stringtie_url', function(x) standardGeneric('os2stringtie_url'))
setGeneric('os2taco_url', function(x) standardGeneric('os2taco_url'))


setReplaceMethod('nthreads', 'Param', function(x, value) {x@NTHREADS=value; x})
setReplaceMethod('outdir',   'Param', function(x, value) {x@OUT_DIR=value; x})
setReplaceMethod('tmpdir',   'Param', function(x, value) {x@TMP_DIR=value; x})
setReplaceMethod('chromoridt', 'Param',
                 function(x, value) {x@CHROM_ORI_DT=value; x})
setReplaceMethod('cufflinks', 'Param',
                 function(x, value) {x@CUFFLINKS_BIN=value; x})
setReplaceMethod('stringtie', 'Param',
                 function(x, value) {x@STRINGTIE_BIN=value; x})
setMethod('cufflinks',     'Param', function(x) x@CUFFLINKS_BIN)
setMethod('maxyieldsize',  'Param', function(x) x@MAX_YIELD_SIZE)
setMethod('outdir',        'Param', function(x) x@OUT_DIR)
setMethod('tmpdir',        'Param', function(x) x@TMP_DIR)
setMethod('chromoridt',    'Param', function(x) x@CHROM_ORI_DT)
setMethod('gtfinfokeys',   'Param', function(x) x@GTF_INFO_KEYS)
setMethod('nthreads',      'Param', function(x) x@NTHREADS)
setMethod('maxunindupaln', 'Param', function(x) x@MAX_UNI_N_DUP_ALN)
setMethod('maxmulndupaln', 'Param', function(x) x@MAX_MUL_N_DUP_ALN)
setMethod('fr1ststrand2mate2flag', 'Param',
          function(x) x@FR1STSTRAND2MATE2FLAG)
setMethod('cufflinkslibtype', 'Param', function(x) x@CUFFLINKS_LIB_TYPE)
setMethod('stringtielibtype', 'Param', function(x) x@STRINGTIE_LIB_TYPE)
setMethod('minisoformfraction',   'Param', function(x) x@MIN_ISOFORM_FRACTION)
setMethod('maxmultireadfraction', 'Param',
          function(x) x@MAX_MULTIREAD_FRACTION)
setMethod('minfragspertransfrag', 'Param',
          function(x) x@MIN_FRAGS_PER_TRANSFRAG)
setMethod('mintrfpkmtoinclude', 'Param',
          function(x) x@MIN_TR_FPKM_TO_INCLUDE)
setMethod('mintrtpmtoinclude', 'Param', function(x) x@MIN_TR_TPM_TO_INCLUDE)
setMethod('os2cufflinks_url',  'Param', function(x) x@OS2CUFFLINKS_URL)
setMethod('os2stringtie_url',  'Param', function(x) x@OS2STRINGTIE_URL)
setMethod('os2taco_url',       'Param', function(x) x@OS2TACO_URL)


#' @importFrom Rsamtools scanBamFlag
#'
setMethod('initialize',
    'Param',
    function(.Object) {
        return(.Object)
    }
)


setGeneric('checkCufflinksBin',
           function(cufflinks_bin, prm) standardGeneric('checkCufflinksBin'))

setMethod('checkCufflinksBin',
    c('character', 'Param'),
    function(cufflinks_bin, prm) {
        if ( ! file.exists(cufflinks_bin) ) {
            url = os2cufflinks_url(prm)[[getOS()]]
            msg = paste0('Cufflinks not found: ', cufflinks_bin, "\n",
                         'It can be downloaded at ', url, "\n")
            stop(msg)
        } else {
            cufflinks(prm) = cufflinks_bin
        }
    }
)


setGeneric('checkCuffmergeBin',
           function(cuffmerge_bin, prm) standardGeneric('checkCuffmergeBin'))

setMethod('checkCuffmergeBin',
    c('character', 'Param'),
    function(cuffmerge_bin, prm) {
        if ( ! file.exists(cuffmerge_bin) ) {
            url = os2cufflinks_url(prm)[[getOS()]]
            msg = paste0('Cuffmerge not found: ', cuffmerge_bin, "\n",
                         'It can be downloaded at ', url, "\n")
            stop(msg)
        } else {
            cuffmerge(prm) = cuffmerge_bin
        }
    }
)


setGeneric('checkStringTieBin',
           function(stringtie_bin, prm) standardGeneric('checkStringTieBin'))

setMethod('checkStringTieBin',
    c('character', 'Param'),
    function(stringtie_bin, prm) {
        if ( ! file.exists(stringtie_bin) ) {
            url = os2stringtie_url(prm)[[getOS()]]
            msg = paste0('StringTie not found: ', stringtie_bin, "\n",
                         'It can be downloaded at ', url, "\n")
            stop(msg)
        } else {
            stringtie(prm) = stringtie_bin
        }
    }
)


setGeneric('checkTacoBin',
           function(taco_bin, prm) standardGeneric('checkTacoBin'))

setMethod('checkTacoBin',
    c('character', 'Param'),
    function(taco_bin, prm) {
        if ( ! file.exists(taco_bin) ) {
            url = os2taco_url(prm)[[getOS()]]
            msg = paste0('TACO not found: ', taco_bin, "\n",
                         'It can be downloaded at ', url, "\n")
            stop(msg)
        } else {
            taco(prm) = taco_bin
        }
    }
)
