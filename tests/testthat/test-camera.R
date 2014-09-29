context("camera")

test_that('camera runs equivalently from do.camera vs direct call', {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gsets <- exampleGeneSets('limma')

  gst <- GeneSetTable(vm, gsets.lol)
  multi.cres <- do.camera(vm, gst, vm$design, ncol(vm$design))

  cres <- camera(vm, gsets, vm$design, ncol(vm$design))

  expect_equal(cres$FDR, multi.cres$padj)
})

