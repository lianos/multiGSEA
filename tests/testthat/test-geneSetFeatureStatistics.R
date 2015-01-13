context("geneSetFeatureStatistics")

test_that("geneSetFeatureStatistics", {

  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), vm)

  mg <- multiGSEA(gsd, vm, vm$design, methods='geneSetTest', ranks.only=TRUE,
                  use.cache=FALSE)

  trim <- 0.10
  min.logFC <- 1
  max.padj <- 0.10
  gs.stats <- geneSetFeatureStatistics(mg, min.logFC, max.padj, trim)

  ## calculate expected
  istats <- calculateIndividualLogFC(vm, vm$design, ncol(vm$design))
  setkeyv(istats, 'featureId')

  ## The code below is the same used in do.geneSetScores -- let's think of a
  ## more manual way to do this.
  expected <- geneSets(gsd)[, {
    ids <- featureIds(gsd, .BY[[1]], .BY[[2]])
    lfc <- istats[J(ids)]
    list(mean.logFC=mean(lfc$logFC, na.rm=TRUE),
         mean.logFC.trim=mean(lfc$logFC, na.rm=TRUE, trim=trim),
         mean.t=mean(lfc$t, na.rm=TRUE),
         mean.t.trim=mean(lfc$t, na.rm=TRUE, trim=trim))
  }, by=c('collection', 'name')]

  expect_equal(expected, gs.stats[, names(expected), with=FALSE])
})
