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
  my <- multiGSEA:::do.romer(gdb, y, d, ncol(d), nrot=nrot, use.cache=FALSE)

  ## order of geneset should be the same as gsd
  expect_equal(geneSets(gdb, as.dt=TRUE)[, list(collection, name)],
               my[, list(collection, name)])
  my[, n := geneSets(gdb)$n]

  ## Columns of camera output are NGenes, Correlation, Direction, PValue, FDR
  ## make `my` look like that, and test for equality
  comp <- local({
    out <- my[, list(n, pval.up, pval.down, pval)]
    data.table::setnames(out, colnames(expected))
    out <- as.matrix(out)
    rownames(out) <- paste(my$collection, my$name, sep=';;')
    out[rownames(expected),]
  })

  expect_equal(comp, expected)
})

test_that("multiGSEA romer fails when not given a DGEList", {

})
