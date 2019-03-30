context("fry")

test_that('fry runs equivalently from do.roast vs direct call', {
  vm <- exampleExpressionSet()
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

  expect_equal(my, fried, check.attributes=FALSE)

  mgr <- multiGSEA:::mgres.fry(my, gsd)
  gs.tuple <- split_gskey(rownames(my))

  expect_equal(mgr$collection, gs.tuple$collection)
  expect_equal(mgr$name, gs.tuple$name)
})

test_that("fry runs through multiGSEA wrapper", {
  vm <- exampleExpressionSet()
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
  mgres <- transform(result(my, 'fry'),
                     key=encode_gskey(collection, name),
                     stringsAsFactors=FALSE)

  expect_true(setequal(rownames(fried), mgres$key))
  comp <- fried[mgres$key,]

  expect_equal(mgres$n, comp$NGenes)
  expect_equal(mgres$pval, comp$PValue)
})
