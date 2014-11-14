context("Gene Set Scores")

test_that("geneSetScores match re-implememtation", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gst <- GeneSetTable(exampleGeneSets('lol'), vm)
  trim <- 0.10

  stats <- calculateIndividualLogFC(vm, vm$design, ncol(vm$design))
  stats <- as.data.table(stats)[, featureId := rownames(stats)]
  setkeyv(stats, 'featureId')

  gs.scores <- do.geneSetScores(vm, gst, vm$design, ncol(vm$design), trim=trim)

  expected <- gst@table[, {
    ids <- featureIds(gst, .BY[[1]], .BY[[2]])
    lfc <- stats[J(ids)]
    list(mean.logFC=mean(lfc$logFC, na.rm=TRUE),
         tmean.logFC=mean(lfc$logFC, na.rm=TRUE, trim=trim),
         mean.t=mean(lfc$t, na.rm=TRUE),
         tmean.t=mean(lfc$t, na.rm=TRUE, trim=trim))
  }, by=c('group', 'id')]

  expect_equal(expected, gs.scores)
})
