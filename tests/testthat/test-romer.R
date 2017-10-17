context("romer")

## TODO: Test why the mixied pvalues are so significant always!!?!

test_that('romer runs equivalently from do.romer vs direct call', {
  es <- exampleExpressionSet(do.voom=FALSE)
  gdb <- GeneSetDb(exampleGeneSets())

  d <-  model.matrix(~ es$Cancer_Status)
  colnames(d)[2] <- 'tumor'

  y <- edgeR::DGEList(exprs(es), group=es$Cancer_Status, genes=fData(es))
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, d)

  gdb <- conform(gdb, y)

  ## We have to ensure that the genesets are tested in the same order as they
  ## are tested from the GeneSetDb for the pvalues to be equivalent given
  ## the same random seed.
  gsd.idxs <- as.list(gdb, value='x.idx')

  ## nrot <- 10000
  nrot <- 250
  seed <- 123
  set.seed(seed)
  expected <- limma::romer(y, gsd.idxs, d, ncol(d), nrot=nrot)

  set.seed(seed)
  my <- multiGSEA:::do.romer(gdb, y, d, ncol(d), nrot=nrot)

  ## Test that inernal call matches direct limma call
  expect_true(setequal(rownames(my), rownames(expected)))
  expected <- expected[rownames(my),,drop=FALSE]
  expect_equal(my, expected, check.attributes=FALSE)

  ## order of geneset should be the same as the GeneSetDb
  expect_equal(rownames(my), encode_gskey(geneSets(gdb)))
  expect_equal(my[, 'NGenes'], geneSets(gdb)$n, check.attributes=FALSE)

  ## multiGSEA pass through & result call matches raw result
  set.seed(seed)
  mg <- multiGSEA(gdb, y, d, ncol(d), method='romer', nrot=nrot)
  res <- result(mg, 'romer')
  res$key <- encode_gskey(res)

  expect_equal(res$key, rownames(my))
  expect_equal(res$n, my[, 'NGenes'], check.attributes=FALSE)
  expect_equal(res$pval.up, my[, 'Up'], check.attributes=FALSE)
  expect_equal(res$pval.down, my[, 'Down'], check.attributes=FALSE)
  expect_equal(res$pval, my[, 'Mixed'], check.attributes=FALSE)
})
