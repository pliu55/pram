#  save all parameters and binary urls
#
Param = setClass('Param',
    slots = list(
        OS2CUFFLINKS_URL = 'list',
        OS2STRINGTIE_URL = 'list',
        OS2TACO_URL      = 'list',

        STAR_BIN      = 'character',
        CUFFLINKS_BIN = 'character',
        STRINGTIE_BIN = 'character',
        TACO_BIN      = 'character',

        RANDOM_SEED = 'numeric',

        FUSERBAMS = 'vector',
        IGGRS     = 'GRanges',
        MODE      = 'character',
        FOUTGTF   = 'character',
        NTHREADS  = 'numeric',
        CUFFLINKS_REF_FA = 'character',

        TMP_DIR = 'character',

        MANAGER_DT   = 'data.table',
        CHROM_ORI_DT = 'data.table',

        GTF_INFO_KEYS = 'vector',


        CUFFLINKS_LIB_TYPE = 'character',
        STRINGTIE_LIB_TYPE = 'character',

        MIN_ISOFORM_FRACTION    = 'numeric',
        MAX_MULTIREAD_FRACTION  = 'numeric',
        MIN_FRAGS_PER_TRANSFRAG = 'numeric',
        MIN_TR_FPKM_TO_INCLUDE  = 'numeric',
        MIN_TR_TPM_TO_INCLUDE   = 'numeric',

        MAX_YIELD_SIZE = 'numeric',
        MAX_CHROM_LEN  = 'numeric',

        MAX_UNI_N_DUP_ALN = 'numeric',
        MAX_MUL_N_DUP_ALN = 'numeric',

        FR1STSTRAND2MATE2FLAG = 'list',


        EXPR_MIN_TPM          = 'numeric', ## 1
        CHIPSEQ_MAX_N_DUP_ALN = 'numeric', ## 5
        TSS_TES_EXT_WIDTH     = 'numeric', ## 1000
        CV_N_FOLDS            = 'numeric', ## 10


        OS = 'character'
    ),

    prototype = list(
        OS2CUFFLINKS_URL = list(
            'LINUX' = paste0('http://cole-trapnell-lab.github.io/cufflinks/',
                'assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz'),

            'OSX' = paste0('http://cole-trapnell-lab.github.io/cufflinks/',
                'assets/downloads/cufflinks-2.1.1.OSX_x86_64.tar.gz'),

            'WINDOWS' = 'http://cole-trapnell-lab.github.io/cufflinks/'
        ),

        OS2STRINGTIE_URL = list(
            'LINUX' = paste0('http://ccb.jhu.edu/software/stringtie/dl/',
                'stringtie-1.3.3b.Linux_x86_64.tar.gz'),

            'OSX' = paste0('http://ccb.jhu.edu/software/stringtie/dl/', 
                'stringtie-1.3.3b.OSX_x86_64.tar.gz'),

            'WINDOWS' = 'https://ccb.jhu.edu/software/stringtie/'
        ),

        OS2TACO_URL = list(
            'LINUX' = paste0('https://github.com/tacorna/taco/releases/',
                'download/v0.7.0/taco-v0.7.0.Linux_x86_64.tar.gz'),

            'OSX' = paste0('https://github.com/tacorna/taco/releases/', 
                'download/v0.7.0/taco-v0.7.0.OSX_x86_64.tar.gz'),

            'WINDOWS' = 'https://tacorna.github.io'
        ),


        RANDOM_SEED = 88888888,

        GTF_INFO_KEYS = c('gene_id', 'transcript_id'),

        NTHREADS = 1,

        TMP_DIR = tempdir(),

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


        MAX_YIELD_SIZE = 600000000, ## overwrite filterBam's default 1M
        MAX_CHROM_LEN = 536870912,  ## max allowed by scanBam, can be applied
                                    ## to hg38, hg19, mm10

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
                                                    isSecondMateRead = T  ) ) ),

        EXPR_MIN_TPM          = 1.0,
        CHIPSEQ_MAX_N_DUP_ALN = 5,
        TSS_TES_EXT_WIDTH     = 1000,
        CV_N_FOLDS            = 10
    )
)


## need to be 'value', other names won't work
setGeneric('cufflinks<-',   function(x, value) standardGeneric('cufflinks<-'))
setGeneric('stringtie<-',   function(x, value) standardGeneric('stringtie<-'))
setGeneric('taco<-',        function(x, value) standardGeneric('taco<-'))
setGeneric('nthreads<-',    function(x, value) standardGeneric('nthreads<-'))
setGeneric('fuserbams<-',   function(x, value) standardGeneric('fuserbams<-'))
setGeneric('iggrs<-',       function(x, value) standardGeneric('iggrs<-'))
setGeneric('mode<-',        function(x, value) standardGeneric('mode<-'))
setGeneric('foutgtf<-',     function(x, value) standardGeneric('foutgtf<-'))
setGeneric('tmpdir<-',      function(x, value) standardGeneric('tmpdir<-'))
setGeneric('cufflinksreffa<-', 
    function(x, value) standardGeneric('cufflinksreffa<-'))
setGeneric('managerdt<-',   function(x, value) standardGeneric('managerdt<-'))
setGeneric('chromoridt<-',  function(x, value) standardGeneric('chromoridt<-'))
setGeneric('exprmintpm<-',  function(x, value) standardGeneric('exprmintpm<-'))
setGeneric('cvnfolds<-',    function(x, value) standardGeneric('cvnfolds<-'))
setGeneric('maxunindupaln<-', 
    function(x, value) standardGeneric('maxunindupaln<-'))
setGeneric('maxmulndupaln<-', 
    function(x, value) standardGeneric('maxmulndupaln<-'))

setGeneric('cufflinks',      function(x) standardGeneric('cufflinks'))
setGeneric('stringtie',      function(x) standardGeneric('stringtie'))
setGeneric('taco',           function(x) standardGeneric('taco'))
setGeneric('randomseed',     function(x) standardGeneric('randomseed'))
setGeneric('maxyieldsize',   function(x) standardGeneric('maxyieldsize'))
setGeneric('maxchromlen',    function(x) standardGeneric('maxchromlen'))
setGeneric('fuserbams',      function(x) standardGeneric('fuserbams'))
setGeneric('iggrs',          function(x) standardGeneric('iggrs'))
setGeneric('mode',           function(x) standardGeneric('mode'))
setGeneric('foutgtf',        function(x) standardGeneric('foutgtf'))
setGeneric('tmpdir',         function(x) standardGeneric('tmpdir'))
setGeneric('cufflinksreffa', function(x) standardGeneric('cufflinksreffa'))
setGeneric('managerdt',      function(x) standardGeneric('managerdt'))
setGeneric('chromoridt',     function(x) standardGeneric('chromoridt'))
setGeneric('gtfinfokeys',    function(x) standardGeneric('gtfinfokeys'))
setGeneric('nthreads',       function(x) standardGeneric('nthreads'))
setGeneric('maxunindupaln',  function(x) standardGeneric('maxunindupaln'))
setGeneric('maxmulndupaln',  function(x) standardGeneric('maxmulndupaln'))
setGeneric(
    'fr1ststrand2mate2flag',
    function(x) standardGeneric('fr1ststrand2mate2flag'))
setGeneric('cufflinkslibtype', function(x) standardGeneric('cufflinkslibtype'))
setGeneric('stringtielibtype', function(x) standardGeneric('stringtielibtype'))
setGeneric(
    'minisoformfraction', function(x) standardGeneric('minisoformfraction'))
setGeneric(
    'maxmultireadfraction', function(x) standardGeneric('maxmultireadfraction'))
setGeneric(
    'minfragspertransfrag', function(x) standardGeneric('minfragspertransfrag'))
setGeneric(
    'mintrfpkmtoinclude', function(x) standardGeneric('mintrfpkmtoinclude'))
setGeneric(
    'mintrtpmtoinclude', function(x) standardGeneric('mintrtpmtoinclude'))
setGeneric('os2cufflinks_url', function(x) standardGeneric('os2cufflinks_url'))
setGeneric('os2stringtie_url', function(x) standardGeneric('os2stringtie_url'))
setGeneric('os2taco_url',      function(x) standardGeneric('os2taco_url'))
setGeneric('exprmintpm',       function(x) standardGeneric('exprmintpm'))
setGeneric(
    'chipseqmaxndupaln', function(x) standardGeneric('chipseqmaxndupaln'))
setGeneric('tsstesextwidth',   function(x) standardGeneric('tsstesextwidth'))
setGeneric('cvnfolds',         function(x) standardGeneric('cvnfolds'))

setReplaceMethod('nthreads',  'Param', function(x, value) {x@NTHREADS=value; x})
setReplaceMethod('fuserbams', 'Param', function(x, value) {x@FUSERBAMS=value;x})
setReplaceMethod('iggrs',     'Param', function(x, value) {x@IGGRS=value; x})
setReplaceMethod('mode',      'Param', function(x, value) {x@MODE=value; x})
setReplaceMethod('foutgtf',   'Param', function(x, value) {x@FOUTGTF=value; x})
setReplaceMethod('tmpdir',    'Param', function(x, value) {x@TMP_DIR=value; x})
setReplaceMethod('cufflinksreffa', 'Param', 
    function(x, value) {x@CUFFLINKS_REF_FA=value; x})
setReplaceMethod('managerdt', 'Param',
    function(x, value) {x@MANAGER_DT=value; x})
setReplaceMethod('chromoridt', 'Param',
    function(x, value) {x@CHROM_ORI_DT=value; x})
setReplaceMethod('cufflinks', 'Param',
    function(x, value) {x@CUFFLINKS_BIN=value; x})
setReplaceMethod('stringtie', 'Param',
    function(x, value) {x@STRINGTIE_BIN=value; x})
setReplaceMethod('taco',       'Param',
    function(x, value) {x@TACO_BIN=value; x})
setReplaceMethod('exprmintpm', 'Param',
    function(x, value) {x@EXPR_MIN_TPM=value; x})
setReplaceMethod('cvnfolds', 'Param',
    function(x, value) {x@CV_N_FOLDS=value; x})
setReplaceMethod('maxunindupaln', 'Param',
    function(x, value) {x@MAX_UNI_N_DUP_ALN=value; x})
setReplaceMethod('maxmulndupaln', 'Param',
    function(x, value) {x@MAX_MUL_N_DUP_ALN=value; x})

setMethod('cufflinks',      'Param', function(x) x@CUFFLINKS_BIN)
setMethod('stringtie',      'Param', function(x) x@STRINGTIE_BIN)
setMethod('taco',           'Param', function(x) x@TACO_BIN)
setMethod('randomseed',     'Param', function(x) x@RANDOM_SEED)
setMethod('maxyieldsize',   'Param', function(x) x@MAX_YIELD_SIZE)
setMethod('maxchromlen',    'Param', function(x) x@MAX_CHROM_LEN)
setMethod('fuserbams',      'Param', function(x) x@FUSERBAMS)
setMethod('iggrs',          'Param', function(x) x@IGGRS)
setMethod('mode',           'Param', function(x) x@MODE)
setMethod('foutgtf',        'Param', function(x) x@FOUTGTF)
setMethod('tmpdir',         'Param', function(x) x@TMP_DIR)
setMethod('cufflinksreffa', 'Param', function(x) x@CUFFLINKS_REF_FA)
setMethod('managerdt',      'Param', function(x) x@MANAGER_DT)
setMethod('chromoridt',     'Param', function(x) x@CHROM_ORI_DT)
setMethod('gtfinfokeys',    'Param', function(x) x@GTF_INFO_KEYS)
setMethod('nthreads',       'Param', function(x) x@NTHREADS)
setMethod('maxunindupaln',  'Param', function(x) x@MAX_UNI_N_DUP_ALN)
setMethod('maxmulndupaln',  'Param', function(x) x@MAX_MUL_N_DUP_ALN)
setMethod('fr1ststrand2mate2flag', 'Param',
    function(x) x@FR1STSTRAND2MATE2FLAG)
setMethod('cufflinkslibtype',     'Param', function(x) x@CUFFLINKS_LIB_TYPE)
setMethod('stringtielibtype',     'Param', function(x) x@STRINGTIE_LIB_TYPE)
setMethod('minisoformfraction',   'Param', function(x) x@MIN_ISOFORM_FRACTION)
setMethod('maxmultireadfraction', 'Param',
    function(x) x@MAX_MULTIREAD_FRACTION)
setMethod('minfragspertransfrag', 'Param',
    function(x) x@MIN_FRAGS_PER_TRANSFRAG)
setMethod('mintrfpkmtoinclude',   'Param',
    function(x) x@MIN_TR_FPKM_TO_INCLUDE)
setMethod('mintrtpmtoinclude', 'Param', function(x) x@MIN_TR_TPM_TO_INCLUDE)
setMethod('os2cufflinks_url',  'Param', function(x) x@OS2CUFFLINKS_URL)
setMethod('os2stringtie_url',  'Param', function(x) x@OS2STRINGTIE_URL)
setMethod('os2taco_url',       'Param', function(x) x@OS2TACO_URL)
setMethod('exprmintpm',        'Param', function(x) x@EXPR_MIN_TPM)
setMethod('chipseqmaxndupaln', 'Param', function(x) x@CHIPSEQ_MAX_N_DUP_ALN)
setMethod('tsstesextwidth',    'Param', function(x) x@TSS_TES_EXT_WIDTH)
setMethod('cvnfolds',          'Param', function(x) x@CV_N_FOLDS)


#' @importFrom Rsamtools scanBamFlag
#'
setMethod('initialize',
    'Param',
    function(.Object) {
        return(.Object)
    }
)


setGeneric('checkCufflinksBin',
    function(prm) standardGeneric('checkCufflinksBin'))

setMethod('checkCufflinksBin', 'Param',
    function(prm) {
        cufflinks_bin = cufflinks(prm)
        if ( ! fileExists(cufflinks_bin) ) {
            url = os2cufflinks_url(prm)[[getOS()]]

            ## check if bin exists in tempdir()
            temp_cufflinks_bin = paste0(tempdir(), '/',
                                        gsub('.tar.gz', '', basename(url)), '/',
                                        'cufflinks')
            if ( ! file.exists(temp_cufflinks_bin) ) {
                stop(paste0("\nCufflinks not found. You may download it from ",
                    url, "\n"))
            } else {
                cufflinks(prm) = temp_cufflinks_bin
            }
        }
        return(prm)
    }
)


setGeneric('checkCuffmergeRequiredBins',
    function(prm) standardGeneric('checkCuffmergeRequiredBins'))

setMethod('checkCuffmergeRequiredBins', 'Param',
    function(prm) {
        prm = checkCufflinksBin(prm)
        to_dl = F
        if ( fileExists( dirname(cufflinks(prm)) ) ) {
            cufflinksdir = paste0(dirname(cufflinks(prm)), '/')
            cuffmerge_bin   = paste0(cufflinksdir, 'cuffmerge')
            cuffcompare_bin = paste0(cufflinksdir, 'cuffcompare')
            gtftosam_bin    = paste0(cufflinksdir, 'gtf_to_sam')

            if (( ! fileExists(cuffmerge_bin)   ) |
                ( ! fileExists(cuffcompare_bin) ) |
                ( ! fileExists(gtftosam_bin)    ) ) {
                    to_dl = T
            }
        } else {
            to_dl = T
        }

        if ( to_dl ) {
            url = os2cufflinks_url(prm)[[getOS()]]
            stop(paste0("\nCuffmerge not found. You may download it from ",
                url, "\n"))
        }

        return(prm)
    }
)


setGeneric('checkStringTieBin',
    function(prm) standardGeneric('checkStringTieBin'))

setMethod('checkStringTieBin', 'Param',
    function(prm) {
        stringtie_bin = stringtie(prm)
        if ( ! fileExists(stringtie_bin) ) {
            url = os2stringtie_url(prm)[[getOS()]]
            ## check if bin exists in tempdir()
            temp_stringtie_bin = paste0(tempdir(), '/',
                                        gsub('.tar.gz', '', basename(url)), '/',
                                        'stringtie')
            if ( ! file.exists(temp_stringtie_bin) ) {
                stop(paste0("\nStringTie not found. You may download it from ",
                    url, "\n"))
            } else {
                stringtie(prm) = temp_stringtie_bin
            }
        }
        return(prm)
    }
)


setGeneric('checkTacoBin',
    function(prm) standardGeneric('checkTacoBin'))

setMethod('checkTacoBin', 'Param',
    function(prm) {
        taco_bin = taco(prm)
        if ( ! fileExists(taco_bin) ) {
            url = os2taco_url(prm)[[getOS()]]
            ## check if bin exists in tempdir()
            temp_taco_bin = paste0(
                tempdir(), '/', gsub('.tar.gz', '', basename(url)), '/',
                'taco_run')
            if ( ! file.exists(temp_taco_bin) ) {
                stop(paste0("\nTACO not found. You may download it from ", url,
                    "\n"))
            } else {
                taco(prm) = temp_taco_bin
            }
        }
        return(prm)
    }
)
