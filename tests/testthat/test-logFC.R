context("Pass-through logFC GSEA method")

test_that("logFC pass through generates expected gene-set stats", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)

  mgc <- multiGSEA(gsd, vm, vm$design, methods='camera')
  mgf <- multiGSEA(gsd, vm, vm$design)
  expect_equal(geneSets(mgc), geneSets(mgf))
})

test_that("t-stats and logFCs match full design when only stats passed", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)

  mgf <- suppressWarnings(multiGSEA(gsd, vm, vm$design))
  x <- logFC(mgf)

  tstats <- setNames(x$t, x$featureId)
  lfc <- setNames(x$logFC, x$featureId)

  mgt <- multiGSEA(gsd, tstats)
  mgl <- multiGSEA(gsd, lfc)

  ro <- suppressWarnings(results(mgf))
  rt <- suppressWarnings(results(mgt))
  rl <- suppressWarnings(results(mgl))

  expect_equal(ro$mean.t, rt$mean.t)
  expect_equal(ro$mean.t.trim, rt$mean.t.trim)

  expect_equal(ro$mean.logFC, rl$mean.logFC)
  expect_equal(ro$mean.logFC.trim, rl$mean.logFC.trim)
})
