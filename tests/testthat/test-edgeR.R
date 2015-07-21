context("Using DGEList's as input")

test_that("logFC calculations are similar between voom and DGEList input", {
  gdb <- GeneSetDb(exampleGeneSets())
  es <- exampleExpressionSet(do.voom=FALSE)
  d <-  model.matrix(~ es$Cancer_Status)
  colnames(d)[2] <- 'tumor'

  y <- DGEList(exprs(es), group=es$Cancer_Status, genes=fData(es))
  y <- calcNormFactors(y)
  y <- estimateDisp(y, d)

  vm <- voom(y, d)

  mgv <- multiGSEA(gdb, vm, d, 2)
  mge <- multiGSEA(gdb, y, d, 2)

  expect_equal(logFC(mge), logFC(mgv))
})
