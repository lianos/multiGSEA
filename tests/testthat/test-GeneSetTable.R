context("GeneSetTable")

test_that("GeneSetTable constructor minimally works", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gsets <- exampleGeneSets('limma')

  gst <- GeneSetTable(gsets.lol, vm)
  expect_is(gst, 'GeneSetTable')
  expect_equal(nrow(vm), length(gst@table$membership[[1]]))

  ## Built a GeneSetTabe that subsets genesets base on size
  gst50 <- GeneSetTable(gsets.lol, vm, min.gs.size=50)
  expect_true(all(sapply(gst50@table$membership, sum) >= 50))
})

test_that("GeneSetTable mapping works", {
  es <- exampleExpressionSet(do.voom=FALSE)
  gsets.lol <- exampleGeneSets('lol')

  gst <- GeneSetTable(gsets.lol, es, min.gs.size=1)
  gt <- gst@table

  ## This is the name of the geneset in `gsets`
  gt[, gid := paste(group, id, sep='.')]
  gt[, feature.id := lapply(membership, function(mm) rownames(es)[mm])]

  ## Test that feature.id's are correct
  for (i in 1:nrow(gt)) {
    id <- gt$id[[i]]
    group <- gt$group[[i]]
    ifeatures <- gt$feature.id[[i]]

    ## expected features: we reverse engineer these from the original feature
    ## ids in gsets.lol and keep the ones that are found in es
    ## in gsets and the rownames of es
    efeatures <- intersect(gsets.lol[[group]][[id]], rownames(es))
    expect_true(setequal(efeatures, ifeatures),
                info=sprintf('%s (intersected entrez)', id))
  }
})

test_that("conform,GeneSetTable works", {
  es.1 <- exampleExpressionSet(do.voom=FALSE)
  es.2 <- es.1[sample(nrow(es.1))]

  gsets.lol <- exampleGeneSets('lol')

  gst.1 <- GeneSetTable(gsets.lol, es.1, min.gs.size=1)
  gst.2 <- GeneSetTable(gsets.lol, es.2, min.gs.size=1)
  gst.c <- conform(gst.1, es.2)

  expect_equal(gst.2@table, gst.c@table)
  expect_equal(gst.2@feature.lookup, gst.c@feature.lookup)
})
