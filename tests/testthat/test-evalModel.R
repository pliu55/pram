library(data.table)

main <- function() {
    context('evalModel')

    testBenchmark()
}


testBenchmark <- function() {
    ftgt = system.file('extdata/benchmark/tgt.tsv', package='pram')
    tgtdt = fread(ftgt, header=T, sep="\t")
    tgttr = Transcript(tgtdt)

    method2results = list(
                  ##--- indi jnc ---##--- tr jnc ---##-------- nuc --------##
                  ## TP    FN   FP  ## TP   FN   FP ##   TP       FN     FP
        'cfmg' = c( 2637, 414, 549,   969, 287, 530,  1549929, 282989, 6929 ),
        'cfpl' = c( 2889, 162, 192,  1138, 118, 252,  1723581, 109337, 8424 ),
        'stmg' = c( 2739, 312,   0,  1005, 251,  89,  1735326,  97592,  505 ),
        'stpl' = c( 2864, 187,   0,  1098, 158,  70,  1748793,  84125,   35 ),
        'taco' = c( 2723, 328, 251,  1038, 218, 417,  1533130, 299788, 9128 )
    )

    lapply(names(method2results), testBenchmarkMethod, method2results, tgttr)
}


testBenchmarkMethod <- function(method, method2results, tgttr) {
    fmdl = system.file(paste0('extdata/benchmark/', method, '.tsv'),
                       package='pram')
    mdldt = fread(fmdl, header=T, sep="\t")
    mdltr = Transcript(mdldt)
    evaldt = evalModel(mdltr, tgttr)

    indi_jnc_tp = method2results[[method]][1]
    indi_jnc_fn = method2results[[method]][2]
    indi_jnc_fp = method2results[[method]][3]

    tr_jnc_tp   = method2results[[method]][4]
    tr_jnc_fn   = method2results[[method]][5]
    tr_jnc_fp   = method2results[[method]][6]

    nuc_tp      = method2results[[method]][7]
    nuc_fn      = method2results[[method]][8]
    nuc_fp      = method2results[[method]][9]

    indi_jnc_precision = indi_jnc_tp/(indi_jnc_tp + indi_jnc_fp)
    indi_jnc_recall    = indi_jnc_tp/(indi_jnc_tp + indi_jnc_fn)

    tr_jnc_precision   = tr_jnc_tp/(tr_jnc_tp + tr_jnc_fp)
    tr_jnc_recall      = tr_jnc_tp/(tr_jnc_tp + tr_jnc_fn)

    nuc_precision      = nuc_tp/(nuc_tp + nuc_fp)
    nuc_recall         = nuc_tp/(nuc_tp + nuc_fn)

    test_that(paste0('evalModel::testBenchmarkMethod::', method), {
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
