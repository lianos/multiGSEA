context("limma::kegga wrapper ('enrichtest')")

dfinput <- exampleDgeResult("human", "ensembl")
gdb. <- getMSigGeneSetDb("h", "human", "ensembl")

test_that("basic 'enrichtest' works", {
  # control for lenght bias
  lbias <- enrichtest(gdb., dfinput, selected = "selected",
                      feature.bias = "effective_length",
                      plot.bias = FALSE)
  if (FALSE) {
    plot_enrichtest_bias(dfinput, "selected", "effective_length")
  }

  # control for expression bias
  ebias <- enrichtest(gdb., dfinput, selected = "selected",
                      feature.bias = "AveExpr",
                      plot.bias = FALSE)
  # expect_true(all(nbias$P.all <= ebias$P.all))
  if (FALSE) {
    plot(-log10(lbias$P.all), -log10(ebias$P.all))
  }


  # no bias
  nbias <- enrichtest(gdb., dfinput, selected = "selected",
                      feature.bias = NULL)
  expect_equal(lbias$Pathway, ebias$Pathway)
  expect_equal(nbias$Pathway, ebias$Pathway)
  expect_equal(lbias$N, ebias$N)
  expect_equal(lbias$all, ebias$all)
  expect_equal(lbias$up, ebias$up)
  expect_equal(lbias$down, ebias$down)
})

test_that("induced length associattion to significance is accounted for", {
  biased <- dfinput %>%
    arrange(pval) %>%
    mutate(effective_length = sort(effective_length, decreasing = TRUE))

  nbias <- enrichtest(gdb., biased, selected = "selected",
                      feature.bias = NULL)
  lbias <- enrichtest(gdb., biased, selected = "selected",
                      feature.bias = "effective_length")
  if (FALSE) {
    plot_enrichtest_bias(biased, "selected", "effective_length")
  }
  expect_equal(nbias$Pathway, lbias$Pathway)
  if (FALSE) {
    plot(-log10(nbias$P.all), -log10(lbias$P.all))
    abline(0,1,col = "red")
  }
  # majority of pvalues when corrected for effective_length should be penalized,
  # which is to say: higher.
  frac.less <- mean(lbias$P.all > nbias$P.all)
  expect_true(frac.less > 0.70)

  # ranomizing length should negate penalty
  set.seed(0xBEEF)
  rando <- mutate(biased, effective_length = sample(effective_length))
  rbias <- enrichtest(gdb., rando, selected = "selected",
                      feature.bias = "effective_length")
  if (FALSE) {
    plot_enrichtest_bias(rando, "selected", "effective_length")
  }
})


test_that("enrichment,groups variable accepts column or list", {
  group.list <- split(dfinput$feature_id, dfinput$direction)
  group.cols <- paste0("P.", c("all", "down", "up"))

  g1 <- enrichtest(gdb, dfinput, selected = "selected", groups = "direction")
  g2 <- enrichtest(gdb, dfinput, selected = "selected", groups = group.list)
  expect_equal(g1$Pathway, g2$Pathway)
  for (pname in group.cols) {
    expect_numeric(g1[[pname]], info = pname)
    expect_equal(g1[[pname]], g2[[pname]], info = pname)
  }
})
