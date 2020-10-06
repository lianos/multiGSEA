context("fgsea")

test_that("multiGSEA calculate t and preranked t match fgsea results", {
  gseaParam <- 1

  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()

  # Since Bioc 3.5, running fgsea warns about ties in preranked stats
  set.seed(123)
  expect_warning({
    mgt <- multiGSEA(gdb, vm, vm$design, 'tumor', 'fgsea', score.by = 't',
                     nPermSimple = nperm, gseaParam = gseaParam)
  }, "ties")
  mgres <- mgt %>%
    result("fgsea") %>%
    mutate(pathway = encode_gskey(collection, name))


  gs.idxs <- as.list(geneSetDb(mgt), active.only=TRUE, value='x.id')
  min.max <- range(sapply(gs.idxs, length))

  lfc <- logFC(mgt)
  ranks.lfc <- setNames(lfc[['logFC']], lfc[['feature_id']])
  ranks.t <- setNames(lfc[['t']], lfc[['feature_id']])

  set.seed(123)
  expect_warning({
    rest <- fgsea::fgsea(gs.idxs, ranks.t,
                         minSize = min.max[1], maxSize = min.max[2],
                         gseaParam = gseaParam)
  }, "ties")

  expect_equal(nrow(mgres), nrow(rest))
  expect_equal(mgres$pathway, rest$pathway)
  expect_equal(rest$size, mgres$n)
  expect_equal(mgres$pval, rest$pval)
  same.sign <- sign(mgres$ES) == sign(rest$ES)
  expect_true(all(same.sign))

  # passing in a preranked vector gives same results ---------------------------
  set.seed(123)
  expect_warning({
    mgpre <- multiGSEA(gdb, ranks.t, methods = "fgsea", nperm = nperm,
                       gseaParam = gseaParam, score.by = "t")
  }, "ties")

  rpre <- result(mgpre, 'fgsea')
  comp.cols <- c('collection', 'name', 'N', 'n', 'size', 'pval', 'ES')
  expect_equal(rpre[, comp.cols], mgres[, comp.cols])

  # Passing in data.frame works, too -------------------------------------------
  set.seed(123)
  expect_warning({
    mgdf <- multiGSEA(gdb, lfc, methods = "fgsea", nperm = nperm,
                      rank_by = "t", rank_order = "descending",
                      gseaParam = gseaParam)
  }, "ties")
  res.df <- result(mgdf)
  expect_equal(res.df[, comp.cols], mgres[, comp.cols])
})

