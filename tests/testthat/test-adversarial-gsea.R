context("Testing GSEA results in adversarial situations")

## How do our methods behave when we have a geneset that has
test_that("Gene Sets with logFCs in oppossite extremes", {
  ## Do the diferential expression analysis on and create some custom genesets
  ## that have members that are:
  ## 1. All on one side of the expression or the other
  ## 2. On mixed sides of the expresion spectrum
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets <- exampleGeneSets('lol')
  lfc <- calculateIndividualLogFC(vm, vm$design, 'tumor')

  ## TODO: get custom gene sets and add them to gsets$
  ## Add custom
  gsets <- exampleGeneSets('limma')

  gst <- GeneSetTable(gsets.lol, vm)

  photo <- camera(vm, gsets, vm$design, ncol(vm$design))
  my <- do.camera(vm, gst, vm$design, ncol(vm$design))



})

