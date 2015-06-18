context("geneSetTest")

test_that("geneSetTest matches re-implememtation", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsi <- exampleGeneSets(vm)
  gsd <- conform(GeneSetDb(gsl), vm)

  ## We have to ensure that the genesets are tested in the same order as they
  ## are tested from the GeneSetDb for the pvalues to be equivalent given
  ## the same random seed.
  gsd.idxs <- multiGSEA:::as.expression.indexes(gsd, value='x.idx')
  gsi <- gsi[names(gsd.idxs)]

  stats <- calculateIndividualLogFC(vm, vm$design, ncol(vm$design))

  nsim <- 250
  seed <- 123
  set.seed(seed)
  expected <- sapply(gsi, geneSetTest, stats$t, nsim=nsim)

  set.seed(seed)
  my.gsd <- multiGSEA:::do.geneSetTest(gsd, vm, logFC=stats, nsim=nsim,
                                        score.by='t')
  comp <- with(my.gsd, {
    out <- setNames(pval, paste(collection, name, sep=';;'))
    out[names(expected)]
  })
  expect_equal(expected, comp)
})
