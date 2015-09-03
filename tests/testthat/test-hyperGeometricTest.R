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
                  feature.min.logFC=min.logFC, feature.max.padj=max.padj)
  res <- results(mg)

  selected <- subset(logFC(mg), abs(logFC) >= min.logFC & padj <= max.padj)
  selected <- unique(selected$featureId)

  hg <- hyperGeometricTest(gsd, selected, rownames(vm))

  expect_equal(hg$pval, res$pval)
})

test_that("do.hyperGeometricTest sets incoming logFC correctly", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), vm)

  min.logFC <- 1
  max.padj <- 0.1

  tt <- calculateIndividualLogFC(vm, vm$design, ncol(vm$design))
  expect_true(is.null(tt$hypuerG.selected))

  my.hg <-
    multiGSEA:::do.hyperGeometricTest(gsd, vm, vm$design,
                                      ncol(vm$design), min.logFC=min.logFC,
                                      max.padj=max.padj, logFC=tt)
  expect_true(is.logical(tt$hyperG.selected))

  tt[, expect.selected := abs(logFC) >= min.logFC & padj <= max.padj]
  expect_equal(tt$hyperG.selected, tt$expect.selected)
})
