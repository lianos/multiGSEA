context("camera")

test_that('camera runs equivalently from do.camera vs direct call', {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsi <- exampleGeneSets(vm)
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), vm)

  photo <- camera(vm, gsi, vm$design, ncol(vm$design))
  my <- multiGSEA:::do.camera(gsd, vm, vm$design, ncol(vm$design))

  ## order of geneset should be the same as gsd
  expect_equal(geneSets(gsd)[, list(collection, name)],
               my[, list(collection, name)])
  my[, n := geneSets(gsd)$n]

  ## Columns of camera output are NGenes, Correlation, Direction, PValue, FDR
  ## make `my` look like that, and test for equality
  comp <- local({
    out <- my[, list(n, Correlation, Direction, pval, padj)]
    setnames(out, names(photo))
    out <- as.data.frame(out)
    rownames(out) <- paste(my$collection, my$name, sep='.')
    out[rownames(photo),]
  })

  expect_equal(photo, comp)
})

