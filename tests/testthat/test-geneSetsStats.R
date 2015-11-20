context("geneSetsStats")

test_that("geneSetsStats", {

  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), vm)

  mg <- multiGSEA(gsd, vm, vm$design, methods='geneSetTest', ranks.only=TRUE)

  trim <- 0.10
  min.logFC <- 1
  max.padj <- 0.10
  gs.stats <- geneSetsStats(mg, min.logFC, max.padj, trim, .external=FALSE)

  ## calculate expected
  istats <- calculateIndividualLogFC(vm, vm$design, ncol(vm$design),
                                     .external=FALSE)
  data.table::setkeyv(istats, 'featureId')

  ## The code below is the same used in do.geneSetScores -- let's think of a
  ## more manual way to do this.
  expected <- geneSets(gsd, .external=FALSE)[, {
    ids <- featureIds(gsd, .BY[[1]], .BY[[2]])
    lfc <- istats[J(ids)]
    t.nona <- lfc$t[!is.na(lfc$t)]
    list(JG=sum(t.nona) / sqrt(length(t.nona)),
         mean.logFC=mean(lfc$logFC, na.rm=TRUE),
         mean.logFC.trim=mean(lfc$logFC, na.rm=TRUE, trim=trim),
         mean.t=mean(lfc$t, na.rm=TRUE),
         mean.t.trim=mean(lfc$t, na.rm=TRUE, trim=trim))
  }, by=c('collection', 'name')]

  expect_equal(expected, gs.stats[, names(expected), with=FALSE])
})
