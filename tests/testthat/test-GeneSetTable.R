context("GeneSetTable")

test_that("GeneSetTable constructor minimally works", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gsets <- exampleGeneSets('limma')

  gst <- GeneSetTable(vm, gsets.lol)
  expect_is(gst, 'GeneSetTable')
  expect_equal(nrow(vm), length(gst@table$membership[[1]]))

  ## Built a GeneSetTabe that subsets genesets base on size
  gst50 <- GeneSetTable(vm, gsets.lol, min.gs.size=50)
  expect_true(all(sapply(gst50@table$membership, sum) >= 50))
})

test_that("GeneSetTable mapping works", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gsets <- unlist(gsets.lol, recursive=FALSE)

  gst <- GeneSetTable(vm, gsets.lol, min.gs.size=1)
  gt <- gst@table
  gt[, gid := paste(group, id, sep='.')]

  for (i in seq(gsets)) {
    gs.id <- names(gsets)[i]
    gs.entrez <- gsets[[i]]
    gs.intersect <- intersect(gs.entrez, rownames(vm))

    gst.entrez <- rownames(vm)[gt$membership[[i]]]

    ## The mapping is totally hosed -- how!?
    expect_true(setequal(gs.intersect, gst.entrez), info=gs.id)
  }

})
