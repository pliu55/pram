library(data.table)

main <- function() {
    context('evalMdl')

   #testRda()

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
    evaldt = evalMdl(mdltr, tgttr)

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

    test_that(paste0('evalMdl::testBenchmarkMethod::', method), {
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


testRda <- function() {
    data(tgtexondt)
    data(mdlexondt)

    tgttr = Transcript(tgtexondt)
    mdltr = Transcript(mdlexondt)
    dt = evalMdl(mdltr, tgttr)

    exon_nuc = subset(dt, feat == 'exon_nuc')
    indi_jnc = subset(dt, feat == 'indi_jnc')
    tr_jnc   = subset(dt, feat == 'tr_jnc')

    exon_nuc_ntp = 1723581
    exon_nuc_nfn = 109337
    exon_nuc_nfp = 8424
    exon_nuc_precision = exon_nuc_ntp/(exon_nuc_ntp + exon_nuc_nfp)
    exon_nuc_recall    = exon_nuc_ntp/(exon_nuc_ntp + exon_nuc_nfn)

    indi_jnc_ntp = 2889
    indi_jnc_nfn =  162
    indi_jnc_nfp =  192
    indi_jnc_precision = indi_jnc_ntp/(indi_jnc_ntp + indi_jnc_nfp)
    indi_jnc_recall    = indi_jnc_ntp/(indi_jnc_ntp + indi_jnc_nfn)

    tr_jnc_ntp = 1138
    tr_jnc_nfn =  118
    tr_jnc_nfp =  252
    tr_jnc_precision = tr_jnc_ntp/(tr_jnc_ntp + tr_jnc_nfp)
    tr_jnc_recall    = tr_jnc_ntp/(tr_jnc_ntp + tr_jnc_nfn)

    test_that('evalMdl::testRda', {
              expect_equal(exon_nuc[, ntp],       exon_nuc_ntp)
              expect_equal(exon_nuc[, nfn],       exon_nuc_nfn)
              expect_equal(exon_nuc[, nfp],       exon_nuc_nfp)
              expect_equal(exon_nuc[, precision], exon_nuc_precision)
              expect_equal(exon_nuc[, recall],    exon_nuc_recall)

              expect_equal(indi_jnc[, ntp],       indi_jnc_ntp)
              expect_equal(indi_jnc[, nfn],       indi_jnc_nfn)
              expect_equal(indi_jnc[, nfp],       indi_jnc_nfp)
              expect_equal(indi_jnc[, precision], indi_jnc_precision)
              expect_equal(indi_jnc[, recall],    indi_jnc_recall)

              expect_equal(tr_jnc[, ntp], 1138)
              expect_equal(tr_jnc[, nfn],  118)
              expect_equal(tr_jnc[, nfp],  252)
              } )
}


main()
