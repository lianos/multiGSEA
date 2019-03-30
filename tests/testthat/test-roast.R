context("roast")

## TODO: Test why the mixied pvalues are so significant always!!?!

test_that('roast runs equivalently from do.roast vs direct call', {
  vm <- exampleExpressionSet()
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

  ## Internal result should match external call
  expect_true(setequal(rownames(my), rownames(roasted)))
  roasted <- roasted[rownames(my),,drop=FALSE]
  expect_equal(my, roasted, check.attributes=FALSE)

  ## order of geneset should be the same as the GeneSetDb
  expect_equal(rownames(my), encode_gskey(geneSets(gsd)))
  expect_equal(my$NGenes, geneSets(gsd)$n)
})

