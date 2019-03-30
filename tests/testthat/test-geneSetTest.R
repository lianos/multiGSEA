context("geneSetTest")

test_that("geneSetTest matches re-implememtation", {
  nsim <- 250
  seed <- 123
  vm <- exampleExpressionSet()
  gdb <- conform(exampleGeneSetDb(), vm)

  stats <- calculateIndividualLogFC(vm, vm$design, ncol(vm$design),
                                    as.dt = TRUE)
  tstats <- setNames(stats$t, stats$featureId)
  gdb <- conform(gdb, names(tstats))
  gsi <- as.list(gdb, value = "x.idx")

  set.seed(seed)
  expected <- sapply(gsi, limma::geneSetTest, tstats, nsim = nsim)

  set.seed(seed)
  mg <- multiGSEA(gdb, vm, vm$design, ncol(vm$design), methods = "geneSetTest",
                  score.by = "t")
  res <- result(mg, "geneSetTest")
  my <- setNames(res$pval, encode_gskey(res))
  my <- my[names(expected)]
  expect_equal(my, expected)
})
