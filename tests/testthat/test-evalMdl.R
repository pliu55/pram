context('evalMdl')

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

test_that('evalMdl', {
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
