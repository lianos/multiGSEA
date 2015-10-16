context("Using DGEList's as input")

test_that("edgeR::quasiLikelihood pipeline run for logFC's of DGEList input", {
  library(edgeR)
  gdb <- GeneSetDb(exampleGeneSets())
  es <- exampleExpressionSet(do.voom=FALSE)
  d <-  model.matrix(~ es$Cancer_Status)
  colnames(d)[2] <- 'tumor'

  y <- DGEList(exprs(es), group=es$Cancer_Status, genes=fData(es))
  y <- calcNormFactors(y)
  y <- estimateDisp(y, d, robust=TRUE)

  fit <- glmQLFit(y, d, robust=TRUE)
  res <- glmQLFTest(fit, 2)
  tte <- as.data.frame(topTags(res, Inf, sort.by='none'))

  mge <- multiGSEA(gdb, y, d, 2, use.treat=FALSE)
  lfc <- logFC(mge)

  ## multiGSEA calls does some reordering of output
  tte <- tte[lfc$featureId,]
  expect_equal(lfc$pval, tte$PValue)
  expect_equal(lfc$AveExpr, tte$logCPM)

  if (FALSE) {
    ## compare to voom?
    vm <- voom(y, d)
    mgv <- multiGSEA(gdb, vm, d, 2)
    vfc <- logFC(mgv)
    plot(-log10(lfc$pval), -log10(vfc$pval), pch=16, col="#3f3f3f33")
    abline(0, 1, col='red')

    not.voom <- subset(vfc, lfc$significant & !significant)
  }
})
