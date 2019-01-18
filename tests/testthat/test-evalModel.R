library(data.table)
suppressMessages(library(GenomicRanges))

main <- function() {
    context('evalModel')

    testBenchmark()
}


testBenchmark <- function() {
    ftgt = system.file('extdata/benchmark/tgt.tsv', package='pram')
    tgtdt = fread(ftgt, header=TRUE, sep="\t")

    mode2results = list(
                  ##--- indi jnc ---##--- tr jnc ---##-------- nuc --------##
                  ## TP    FN   FP  ## TP   FN   FP ##   TP       FN     FP
        'plcf' = c( 2889, 162, 192,  1138, 118, 252,  1723581, 109337, 8424 )
      # 'plst' = c( 2864, 187,   0,  1098, 158,  70,  1748793,  84125,   35 ),
      # 'cfmg' = c( 2637, 414, 549,   969, 287, 530,  1549929, 282989, 6929 ),
      # 'stmg' = c( 2739, 312,   0,  1005, 251,  89,  1735326,  97592,  505 ),
      # 'cftc' = c( 2723, 328, 251,  1038, 218, 417,  1533130, 299788, 9128 )
    )

    mdldtlist = lapply(names(mode2results), readModel)
    names(mdldtlist) = names(mode2results)

    lapply(names(mode2results), testEvalModelByDT, mode2results,
           mdldtlist, tgtdt)

    testEvalModelByGR('plcf', mode2results, mdldtlist, tgtdt)

    testEvalModelByGTF('plcf', mode2results, mdldtlist, tgtdt)
}


testEvalModelByGTF <- function(mode, mode2results, mdldtlist, tgtdt) {
    mdldt = mdldtlist[[mode]]
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
    fmdlgtf = paste0(tempdir(), '/', mode, '_mdl.gtf')
    ftgtgtf = paste0(tempdir(), '/', mode, '_tgt.gtf')
    writeGTF(mdlgtf, fmdlgtf, append=FALSE)
    writeGTF(tgtgtf, ftgtgtf, append=FALSE)

    evaldt = evalModel(fmdlgtf, ftgtgtf)

    results = mode2results[[mode]]
    testResults(results, evaldt, 'testEvalModelByGTF', mode)
}


testEvalModelByGR <- function(mode, mode2results, mdldtlist, tgtdt) {
    mdldt = mdldtlist[[mode]]
    mdlgrs = makeGRangesFromDataFrame(mdldt, keep.extra.columns=TRUE)
    tgtgrs = makeGRangesFromDataFrame(tgtdt, keep.extra.columns=TRUE)

    evaldt = evalModel(mdlgrs, tgtgrs)

    results = mode2results[[mode]]
    testResults(results, evaldt, 'testEvalModelByGR', mode)
}


testEvalModelByDT <- function(mode, mode2results, mdldtlist, tgtdt) {
    mdldt = mdldtlist[[mode]]
    evaldt = evalModel(mdldt, tgtdt)

    results = mode2results[[mode]]
    testResults(results, evaldt, 'testEvalModelByDT', mode)
}


readModel <- function(mode) {
    fmdl = system.file(paste0('extdata/benchmark/', mode, '.tsv'),
                       package='pram')
    mdldt = fread(fmdl, header=TRUE, sep="\t")

    return(mdldt)
}


testResults <- function(results, evaldt, func_name, mode) {
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

    test_that(paste0('evalModel::', func_name, mode), {
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
