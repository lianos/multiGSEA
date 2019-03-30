context("Using DGEList's as input")

test_that("edgeR::quasiLikelihood pipeline run for logFC's of DGEList input", {
  library(edgeR)
  gdb <- GeneSetDb(exampleGeneSets())
  y <- exampleExpressionSet(do.voom = FALSE)
  d <-  model.matrix(~ y$samples$Cancer_Status)
  colnames(d)[2] <- 'tumor'

  y <- estimateDisp(y, d, robust = TRUE)

  fit <- glmQLFit(y, d, robust = TRUE)
  res <- glmQLFTest(fit, 2)
  tte <- as.data.frame(topTags(res, Inf, sort.by='none'))

  mge <- multiGSEA(gdb, y, d, 2, use.treat=FALSE)
  lfc <- logFC(mge)

  ## multiGSEA calls does some reordering of output
  tte <- tte[lfc$featureId,]
  expect_equal(lfc$pval, tte$PValue)
  expect_equal(lfc$AveExpr, tte$logCPM)
})
