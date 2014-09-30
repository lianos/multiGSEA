context("camera")

test_that('camera runs equivalently from do.camera vs direct call', {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gsets <- exampleGeneSets('limma')

  gst <- GeneSetTable(gsets.lol, vm)

  photo <- camera(vm, gsets, vm$design, ncol(vm$design))
  my <- do.camera(vm, gst, vm$design, ncol(vm$design))

  ## Columns of camera output are NGenes, Correlation, Direction, PValue, FDR
  ## make `my` look like that, and test for equality
  comp <- local({
    out <- my[, list(n, Correlation, Direction, pval, padj)]
    setnames(out, names(photo))
    out <- as.data.frame(out)
    rownames(out) <- paste(my$group, my$id, sep='.')
    out
  })

  expect_equal(photo, comp)
})

