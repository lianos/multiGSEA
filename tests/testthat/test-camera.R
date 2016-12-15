context("camera")

test_that('camera runs equivalently from do.camera vs direct call', {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsd <- conform(exampleGeneSetDb(), vm)
  gsi <- as.list(gsd, value='x.idx')

  photo <- limma::camera(vm, gsi, vm$design, ncol(vm$design))
  my <- multiGSEA:::do.camera(gsd, vm, vm$design, ncol(vm$design))

  ## order of geneset should be the same as gsd
  expect_equal(geneSets(gsd, .external=FALSE)[, list(collection, name)],
               my[, list(collection, name)])
  my[, n := geneSets(gsd, .external=FALSE)$n]

  ## Columns of camera output are NGenes, Correlation, Direction, PValue, FDR
  ## make `my` look like that, and test for equality
  comp <- local({
    ## NOTE: As of Bioc3.3 camera has a default inter.gene.cor of 0.01. If this
    ##       parameter isn't explicitly set to NULL, then it is used and the
    ##       Correlation column of camera's output is dropped
    if ('Correlation' %in% names(photo)) {
      out <- my[, list(n, Correlation, Direction, pval, padj)]
    } else {
      out <- my[, list(n, Direction, pval, padj)]
    }
    data.table::setnames(out, names(photo))
    out <- as.data.frame(out)
    rownames(out) <- paste(my$collection, my$name, sep=';;')
    out[rownames(photo),]
  })

  expect_equal(photo, comp)
})

