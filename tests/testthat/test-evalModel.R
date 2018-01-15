library(data.table)
suppressMessages(library(GenomicRanges))

main <- function() {
    context('evalModel')

    testBenchmark()
}


testBenchmark <- function() {
    ftgt = system.file('extdata/benchmark/tgt.tsv', package='pram')
    tgtdt = fread(ftgt, header=T, sep="\t")

    method2results = list(
                  ##--- indi jnc ---##--- tr jnc ---##-------- nuc --------##
                  ## TP    FN   FP  ## TP   FN   FP ##   TP       FN     FP
        'cfmg' = c( 2637, 414, 549,   969, 287, 530,  1549929, 282989, 6929 ),
        'cfpl' = c( 2889, 162, 192,  1138, 118, 252,  1723581, 109337, 8424 ),
        'stmg' = c( 2739, 312,   0,  1005, 251,  89,  1735326,  97592,  505 ),
        'stpl' = c( 2864, 187,   0,  1098, 158,  70,  1748793,  84125,   35 ),
        'taco' = c( 2723, 328, 251,  1038, 218, 417,  1533130, 299788, 9128 )
    )

    nthr = 4

    mdldtlist = mclapply(names(method2results), readModel, mc.cores=nthr)
    names(mdldtlist) = names(method2results)

    mclapply(names(method2results), testEvalModelByDT, method2results,
             mdldtlist, tgtdt, mc.cores=nthr)

    mclapply(names(method2results), testEvalModelByGR, method2results,
             mdldtlist, tgtdt, mc.cores=nthr)

    mclapply(names(method2results), testEvalModelByGTF, method2results,
             mdldtlist, tgtdt, mc.cores=nthr)
}


testEvalModelByGTF <- function(method, method2results, mdldtlist, tgtdt) {
    mdldt = mdldtlist[[method]]
    mdlgtf = new('GTF')
    tgtgtf = new('GTF')

    if ( ! 'transcript_id' %in% names(mdldt) ) {
        setnames(mdldt, 'trid', 'transcript_id')
    }

    if ( ! 'transcript_id' %in% names(tgtdt) ) {
        setnames(tgtdt, 'trid', 'transcript_id')
    }

    mdldt[, feature := 'exon']
    tgtdt[, feature := 'exon']
    mdlgtf = initFromDataTable(mdlgtf, mdldt, c('transcript_id'))
    tgtgtf = initFromDataTable(tgtgtf, tgtdt, c('transcript_id'))
    fmdlgtf = paste0(tempdir(), '/', method, '_mdl.gtf')
    ftgtgtf = paste0(tempdir(), '/', method, '_tgt.gtf')
    writeGTF(mdlgtf, fmdlgtf, append=F)
    writeGTF(tgtgtf, ftgtgtf, append=F)

    evaldt = evalModel(fmdlgtf, ftgtgtf)

    results = method2results[[method]]
    testResults(results, evaldt, 'testEvalModelByGTF', method)
}


testEvalModelByGR <- function(method, method2results, mdldtlist, tgtdt) {
    mdldt = mdldtlist[[method]]
    mdlgrs = makeGRangesFromDataFrame(mdldt, keep.extra.columns=T)
    tgtgrs = makeGRangesFromDataFrame(tgtdt, keep.extra.columns=T)

    evaldt = evalModel(mdlgrs, tgtgrs)

    results = method2results[[method]]
    testResults(results, evaldt, 'testEvalModelByGR', method)
}


testEvalModelByDT <- function(method, method2results, mdldtlist, tgtdt) {
    mdldt = mdldtlist[[method]]
    evaldt = evalModel(mdldt, tgtdt)

    results = method2results[[method]]
    testResults(results, evaldt, 'testEvalModelByDT', method)
}


readModel <- function(method) {
    fmdl = system.file(paste0('extdata/benchmark/', method, '.tsv'),
                       package='pram')
    mdldt = fread(fmdl, header=T, sep="\t")

    return(mdldt)
}


testResults <- function(results, evaldt, func_name, method) {
    indi_jnc_tp = results[1]
    indi_jnc_fn = results[2]
    indi_jnc_fp = results[3]

    tr_jnc_tp   = results[4]
    tr_jnc_fn   = results[5]
    tr_jnc_fp   = results[6]

    nuc_tp      = results[7]
    nuc_fn      = results[8]
    nuc_fp      = results[9]

    indi_jnc_precision = indi_jnc_tp/(indi_jnc_tp + indi_jnc_fp)
    indi_jnc_recall    = indi_jnc_tp/(indi_jnc_tp + indi_jnc_fn)

    tr_jnc_precision   = tr_jnc_tp/(tr_jnc_tp + tr_jnc_fp)
    tr_jnc_recall      = tr_jnc_tp/(tr_jnc_tp + tr_jnc_fn)

    nuc_precision      = nuc_tp/(nuc_tp + nuc_fp)
    nuc_recall         = nuc_tp/(nuc_tp + nuc_fn)

    test_that(paste0('evalModel::', func_name, method), {
        expect_equal(indi_jnc_tp, evaldt[feat=='indi_jnc', ntp])
        expect_equal(indi_jnc_fn, evaldt[feat=='indi_jnc', nfn])
        expect_equal(indi_jnc_fp, evaldt[feat=='indi_jnc', nfp])

        expect_equal(tr_jnc_tp,   evaldt[feat=='tr_jnc', ntp])
        expect_equal(tr_jnc_fn,   evaldt[feat=='tr_jnc', nfn])
        expect_equal(tr_jnc_fp,   evaldt[feat=='tr_jnc', nfp])

        expect_equal(nuc_tp,      evaldt[feat=='exon_nuc', ntp])
        expect_equal(nuc_fn,      evaldt[feat=='exon_nuc', nfn])
        expect_equal(nuc_fp,      evaldt[feat=='exon_nuc', nfp])

        expect_equal(indi_jnc_precision, evaldt[feat=='indi_jnc', precision ])
        expect_equal(indi_jnc_recall,    evaldt[feat=='indi_jnc', recall    ])

        expect_equal(tr_jnc_precision,   evaldt[feat=='tr_jnc', precision   ])
        expect_equal(tr_jnc_recall,      evaldt[feat=='tr_jnc', recall      ])

        expect_equal(nuc_precision,      evaldt[feat=='exon_nuc', precision ])
        expect_equal(nuc_recall,         evaldt[feat=='exon_nuc', recall    ])
    })
}


main()
