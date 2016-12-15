context("roast")

## TODO: Test why the mixied pvalues are so significant always!!?!

test_that('roast runs equivalently from do.roast vs direct call', {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsi <- exampleGeneSets(vm)
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), vm)

  ## We have to ensure that the genesets are tested in the same order as they
  ## are tested from the GeneSetDb for the pvalues to be equivalent given
  ## the same random seed.
  gsd.idxs <- as.list(gsd, value='x.idx')
  gsi <- gsi[names(gsd.idxs)]
  ## nrot <- 10000
  nrot <- 250
  seed <- 123
  set.seed(seed)
  roasted <- limma::mroast(vm, gsi, vm$design, ncol(vm$design), nrot=nrot,
                           sort='none')

  set.seed(seed)
  my <- multiGSEA:::do.roast(gsd, vm, vm$design, ncol(vm$design), nrot=nrot,
                             use.cache=FALSE)

  ## order of geneset should be the same as gsd
  expect_equal(geneSets(gsd, .external=FALSE)[, list(collection, name)],
               my[, list(collection, name)])
  my[, n := geneSets(gsd, .external=FALSE)$n]

  ## Columns of camera output are NGenes, Correlation, Direction, PValue, FDR
  ## make `my` look like that, and test for equality
  comp <- local({
    out <- my[, list(n, PropDown, PropUp, Direction, pval, padj,
                     pval.mixed, padj.mixed)]
    data.table::setnames(out, names(roasted))
    out <- as.data.frame(out)
    rownames(out) <- paste(my$collection, my$name, sep=';;')
    out[rownames(roasted),]
  })

  expect_equal(roasted, comp)
})

