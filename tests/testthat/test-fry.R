context("fry")

test_that('fry runs equivalently from do.roast vs direct call', {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsi <- exampleGeneSets(vm)
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), vm)

  ## We have to ensure that the genesets are tested in the same order as they
  ## are tested from the GeneSetDb for the pvalues to be equivalent given
  ## the same random seed.
  gsd.idxs <- as.list(gsd, value='x.idx')
  gsi <- gsi[names(gsd.idxs)]

  fried <- limma::fry(vm, gsi, vm$design, ncol(vm$design), sort=FALSE)
  my <- multiGSEA:::do.fry(gsd, vm, vm$design, ncol(vm$design))

  ## order of geneset should be the same as gsd
  expect_equal(geneSets(gsd, as.dt=TRUE)[, list(collection, name)],
               my[, list(collection, name)])
  my[, NGenes := geneSets(gsd, as.dt=TRUE)$n]

  ## Columns of camera output are NGenes, Correlation, Direction, PValue, FDR
  ## make `my` look like that, and test for equality
  comp <- local({
    out <- my[, !names(my) %in% c('collection', 'name'), with=FALSE]
    data.table::setnames(out,
                         c('pval', 'padj', 'pval.mixed', 'padj.mixed'),
                         c('PValue', 'FDR', 'PValue.Mixed', 'FDR.Mixed'))
    out <- as.data.frame(out)[, names(fried)]
    rownames(out) <- paste(my$collection, my$name, sep=';;')
    out[rownames(fried),]
  })

  expect_equal(fried, comp)
})

test_that("fry runs through multiGSEA wrapper", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsi <- exampleGeneSets(vm)
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), vm)

  ## We have to ensure that the genesets are tested in the same order as they
  ## are tested from the GeneSetDb for the pvalues to be equivalent given
  ## the same random seed.
  gsd.idxs <- as.list(gsd, value='x.idx')
  gsi <- gsi[names(gsd.idxs)]

  fried <- limma::fry(vm, gsi, vm$design, ncol(vm$design), sort=FALSE)
  my <- multiGSEA(gsd, vm, vm$design, ncol(vm$design), methods='fry')

})
