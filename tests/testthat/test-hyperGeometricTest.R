context("Hyper Geometric Test")

##test_that("do.hyperGeometricTest performs like standard hyperG test", {
##  vm <- exampleExpressionSet(do.voom=TRUE)
##  gsl <- exampleGeneSets()
##  gsd <- conform(GeneSetDb(gsl), vm)
##
##  min.logFC <- 1
##  max.padj <- 0.1
##
##  my.hg <- multiGSEA:::do.hyperGeometricTest(gsd,vm,vm$design,ncol(vm$design))
##
##  ## Calculate expected resulted from hyperG test
##  tt <- calculateIndividualLogFC(vm, vm$design, ncol(vm$design),
##                                 min.logFC=min.logFC, max.padj=max.padj)
##  tt[, hyperG.selected := abs(logFC) >= min.logFC & padj <= max.padj]
##})

test_that("hyperGeometricTest performs like do.hyperGeometricTest performs", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), vm)

  min.logFC <- 1
  max.padj <- 0.1

  mg <- multiGSEA(gsd, vm, vm$design, methods='hyperGeometricTest',
                  feature.min.logFC=min.logFC, feature.max.padj=max.padj,
                  really.use.hyperGeometricTest=TRUE)
  res <- results(mg)

  selected <- subset(logFC(mg), significant)
  selected <- unique(selected$featureId)

  hg <- hyperGeometricTest(gsd, selected, rownames(vm))

  expect_equal(hg$pval, res$pval)
})
