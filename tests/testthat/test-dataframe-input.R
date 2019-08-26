context("data.frame ranks/DGE input to multiGSEA")

xdf <- exampleDgeResult()
scores <- setNames(xdf$logFC, xdf$feature_id)

gdb <- getMSigGeneSetDb("h", "human", "ensembl")
gdb <- conform(gdb, xdf$feature_id)
gs.idx <- as.list(gdb, active.only = TRUE, value = "x.idx")

test_that("ranks-based gsea works", {
  mg <- multiGSEA(gdb, scores, methods = c("cameraPR"), xmeta. = xdf)
  mgres <- result(mg)
  mgres$key <- encode_gskey(mgres)

  cpr <- limma::cameraPR(scores, gs.idx, sort = FALSE)

  expect_equal(mgres$key, rownames(cpr))
  expect_equal(mgres$n, cpr$NGenes)
  expect_equal(mgres$Direction, cpr$Direction)
  expect_equal(mgres$pval, cpr$PValue)
  expect_equal(mgres$padj, cpr$FDR)
})

test_that("enrichment-based methods work", {
  fbias <- setNames(xdf$effective_length, xdf$feature_id)
  mg <- multiGSEA(gdb, scores, methods = "goseq", xmeta. = xdf,
                  feature.bias = fbias)
  mgres <- result(mg, "goseq")
  mgres$key <- encode_gskey(mgres)

  gseq <- expect_warning({
    multiGSEA::goseq(
      gdb,
      selected = subset(xdf, significant)$feature_id,
      universe = xdf$feature_id,
      feature.bias = fbias)
  }, "initial point")

  expect_equal(mgres$key, gseq$category)
  expect_equal(mgres$pval, gseq$over_represented_pvalue)
  expect_equal(mgres$pval.under, gseq$under_represented_pvalue)
  expect_equal(mgres$n.drawn, gseq$numDEInCat)
  expect_equal(mgres$n, gseq$numInCat)
})
