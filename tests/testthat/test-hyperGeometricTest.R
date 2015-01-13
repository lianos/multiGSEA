context("Hyper Geometric Test")

test_that("do.hyperGeometricTest performs like standard hyperG test", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), vm)

  min.logFC <- 1
  max.padj <- 0.1

  my.hg <- multiGSEA:::do.hyperGeometricTest(gsd, vm, vm$design,
                                             ncol(vm$design))

  ## Calculate expected resulted from hyperG test
  tt <- calculateIndividualLogFC(vm, vm$design, ncol(vm$design),
                                 min.logFC=min.logFC, max.padj=max.padj)
  tt[, hyperG.selected := abs(logFC) >= min.logFC & padj <= max.padj]
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
  expect_true(all(tt$hyperG.selected == tt$expect.selected))
})
