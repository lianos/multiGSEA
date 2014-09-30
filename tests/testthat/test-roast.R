context("roast")

test_that('roast runs equivalently from do.roast vs direct call', {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gsets <- exampleGeneSets('limma')

  gst <- GeneSetTable(gsets.lol, vm)

  set.seed(123)
  roasted <- mroast(vm, gsets, vm$design, ncol(vm$design), nrot=500,
                    sort='none')

  my <- do.roast(vm, gst, vm$design, ncol(vm$design), nrot=500, .seed=123)

  ## Columns of camera output are NGenes, Correlation, Direction, PValue, FDR
  ## make `my` look like that, and test for equality
  comp <- local({
    out <- my[, list(n, PropDown, PropUp, Direction, pval, padj,
                     pval.mixed, padj.mixed)]
    setnames(out, names(roasted))
    out <- as.data.frame(out)
    rownames(out) <- paste(my$group, my$id, sep='.')
    out[rownames(roasted),]
  })

  expect_equal(roasted, comp)
})

