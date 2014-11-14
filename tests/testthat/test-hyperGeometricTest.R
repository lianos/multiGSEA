context("Hyper Geometric Test")

test_that('do.hyperGeometricTest performs like standard hyperG test', {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gst <- GeneSetTable(gsets.lol, vm)

  ## Calculate expected resulted from hyperG test
  tt <- calculateIndividualLogFC(vm, vm$design, ncol(vm$design),
                                 provide='table')
  univ <- rownames(vm)
  selected <- rownames(subset(tt, abs(logFC) >= 1 & padj <= 0.1))

#   expected <- lapply(gst@table$membership, function(idxs) {
#     membership <- featureNames(gst)[idxs]
#     p <- new("HyperGParams", geneIds=rownames(vm)[gst$])
#   })
#
#   expect_equal(hyperGresults, my.results)
})
