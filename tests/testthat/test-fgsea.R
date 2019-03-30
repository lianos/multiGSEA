context("fgsea")

test_that("multiGSEA calculate t and preranked t match fgsea results", {
  nperm <- 1000
  gseaParam <- 1

  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()

  # Since Bioc 3.5, running fgsea warns about ties in preranked stats
  expect_warning({
    mgt <- multiGSEA(gdb, vm, vm$design, 'tumor', 'fgsea', nperm=nperm,
                     gseaParam=gseaParam, score.by='t')
  }, "ties")


  gs.idxs <- as.list(geneSetDb(mgt), active.only=TRUE, value='x.id')
  min.max <- range(sapply(gs.idxs, length))

  lfc <- logFC(mgt)
  ranks.lfc <- setNames(lfc[['logFC']], lfc[['featureId']])
  ranks.t <- setNames(lfc[['t']], lfc[['featureId']])

  expect_warning({
    rest <- fgsea::fgsea(gs.idxs, ranks.t, nperm, min.max[1], min.max[2],
                         gseaParam=gseaParam)
  }, "ties")

  mgres <- mgt %>%
    result("fgsea") %>%
    mutate(pathway = encode_gskey(collection, name))
  expect_equal(nrow(mgres), nrow(rest))
  expect_equal(mgres$pathway, rest$pathway)
  expect_equal(rest$size, mgres$n)
  same.sign <- sign(mgres$ES) == sign(rest$ES)
  expect_true(all(same.sign))

  ## passing in a preranked vector gives same results
  expect_warning({
    mgpre <- multiGSEA(gdb, ranks.t, methods='fgsea', nperm=nperm,
                       gseaParam=gseaParam, score.by='t')
  }, "ties")

  rpre <- result(mgpre, 'fgsea')
  comp.cols <- c('collection', 'name', 'N', 'n', 'size')
  expect_equal(rpre[, comp.cols], mgres[, comp.cols])
  expect_equal(sign(rpre$ES), sign(mgres$ES))
})
