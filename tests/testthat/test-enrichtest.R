context("limma::kegga wrapper ('enrichtest')")

gdb. <- getMSigGeneSetDb("h", "human", "ensembl")

test_that("induced length associattion to significance is accounted for", {
  biased <- exampleDgeResult("human", "ensembl",
                             induce.bias = "effective_length")

  nbias <- enrichtest(gdb., biased, selected = "selected",
                      feature.bias = NULL)
  lbias <- enrichtest(gdb., biased, selected = "selected",
                      feature.bias = "effective_length")
  if (FALSE) {
    plot_enrichtest_bias(biased, "selected", "effective_length")
  }
  expect_equal(nbias$Pathway, lbias$Pathway)
  expect_equal(nbias$N, lbias$N)
  expect_equal(nbias$all, lbias$all)
  expect_equal(nbias$up, lbias$up)
  expect_equal(nbias$down, lbias$down)

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
  expect_equal(rbias$P.all, nbias$P.all, tolerance = 0.005)

  if (FALSE) {
    plot_enrichtest_bias(rando, "selected", "effective_length")
  }
})

test_that("enrichtest,groups variable accepts column or list", {
  dfinput <- exampleDgeResult("human", "ensembl")
  group.list <- split(dfinput$feature_id, dfinput$direction)
  p.cols <- paste0("P.", c("all", "down", "up"))

  g1 <- enrichtest(gdb., dfinput, selected = "selected", groups = "direction")
  g2 <- enrichtest(gdb., dfinput, selected = "selected", groups = group.list)
  expect_equal(g1$Pathway, g2$Pathway)
  for (pname in p.cols) {
    expect_numeric(g1[[pname]], info = pname)
    expect_equal(g1[[pname]], g2[[pname]], info = pname)
  }

  # This is some metadata that is required for correctly processing these
  # results inside the 'multiGSEA' pipeline
  expect_true(attr(g1, "mgunlist"))
  expect_setequal(attr(g1, "groups"), c("all", names(group.list)))
  expect_true(attr(g1, "rawresult"))
})

test_that("enrichtest and goseq give probably approximately correct answers", {
  dfinput <- exampleDgeResult("human", "ensembl",
                              induce.bias = "effective_length")

  # no bias correction
  e1 <- enrichtest(gdb., dfinput, selected = "selected", groups = "direction")
  g1 <- expect_warning({
    multiGSEA::goseq(
      gdb.,
      subset(dfinput, selected)$feature_id,
      dfinput$feature_id,
      setNames(dfinput$effective_length, dfinput$feature_id),
      method = "Hypergeometric")
  }, "initial point")

  expect_equal(e1$Pathway, g1$category)
  expect_equal(e1$P.all, g1$over_represented_pvalue)

  # length correction
  e2 <- enrichtest(gdb., dfinput, selected = "selected", groups = "direction",
                   feature.bias = "effective_length")
  g2 <- expect_warning({
    multiGSEA::goseq(
      gdb.,
      subset(dfinput, selected)$feature_id,
      dfinput$feature_id,
      setNames(dfinput$effective_length, dfinput$feature_id),
      method = "Wallenius")
  }, "initial point")

  if (FALSE) {
    par(mfrow = c(1, 2))
    plot(-log10(e1$P.all), -log10(e2$P.all),
         main = "Uncorrected vs corrected enrichtest",
         xlab = "Uncorrected", ylab = "Corrected")
    abline(0, 1, col = "red")
    plot(-log10(g1$over_represented_pvalue), -log10(g2$over_represented_pvalue),
         main = "Uncorrected vs corrected goseq",
         xlab = "Uncorrected", ylab = "Corrected")
    abline(0, 1, col = "red")

    par(mfrow = c(1,1))

    # enrichtest is a bit more conservative
    plot(-log10(e2$P.all), -log10(g2$over_represented_pvalue),
         main = "Corrected enrichtest vs goseq",
         xlab = "enrichtest", ylab = "goseq")
    abline(0, 1, col = "red")
  }
  # test that average difference is less than a threshold
  pval.diff <- abs(e2$P.all - g2$over_represented_pvalue)
  expect_lt(mean(pval.diff), 0.025)
})

test_that("'naked' enrichtest call vs multiGSEA pipeline are equivalent", {
  dfinput <- exampleDgeResult("human", "ensembl",
                              induce.bias = "effective_length")
  nres <- enrichtest(gdb., dfinput, selected = "selected", groups = "direction",
                     feature.bias = "effective_length")
  mres <- multiGSEA(gdb., setNames(dfinput$t, dfinput$feature_id),
                    methods = "enrichtest", feature.bias = "effective_length",
                    xmeta. = dfinput)

  groups <- c(all = "enrichtest", up = "enrichtest.up",
              down = "enrichtest.down")
  expect_setequal(resultNames(mres), groups)

  for (i in seq(groups)) {
    # Call enrichtest direct
    ename <- names(groups)[i]
    pcol <- paste0("P.", ename)
    cmp <- nres[, c("Pathway", "N", ename, pcol)]

    # Pull out of MultiGSEAResult object
    mname <- groups[i]
    mg <- result(mres, mname)
    mg <- mg[, c("collection", "name", "N", "n", "n.drawn", "pval")]

    expect_equal(mg$name, sub(".*;;", "", nres$Pathway), info = ename)
    expect_equal(mg$n, cmp$N, info = ename)
    expect_equal(mg$n.drawn, cmp[[ename]], info = ename)
    expect_equal(mg$pval, cmp[[pcol]], info = ename)
  }
})
