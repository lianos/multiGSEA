context("goseq")

test_that("internal goseq mimics goseq package", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  gsd <- conform(gsd, vm)

  mg <- multiGSEA(gsd, vm, vm$design)
  lfc <- logFC(mg)
  selected <- subset(lfc, significant)$featureId
  universe <- rownames(vm)

  ## lens <- goseq::getlength(universe, 'hg19', 'knownGene')
  mylens <- multiGSEA:::create.glength.vector(vm)

  degenes <- setNames(integer(length(universe)), universe)
  degenes[selected] <- 1L

  my.res <- multiGSEA::goseq(gsd, selected, universe, method='Wallenius',
                             use_genes_without_cat=TRUE)
  ## pwf <- attr(my.res, 'pwf')

  ## expected
  g2c <- transform(as.data.frame(gsd), category=name)
  g2c <- g2c[, c('category', 'featureId')]
  pwf <- goseq::nullp(degenes, bias.data=mylens)
  goseq.res <- goseq::goseq(pwf, gene2cat=g2c, method='Wallenius',
                            use_genes_without_cat=TRUE)

  ## Match up and compare
  goseq.res <- goseq.res[match(my.res$name, goseq.res$category),]
  my <- setnames(copy(my.res),
                 c("name", "n", "n.drawn"),
                 c("category", "numInCat", "numDEInCat"))
  my <- my[, names(goseq.res), with=FALSE]
  setDF(my)
  expect_equal(my, goseq.res, check.attributes=FALSE)

  mgs <- multiGSEA(gsd, vm, vm$design, methods='goseq')

  r <- result(mgs, 'goseq')
  rc <- r[, list(collection, name, n, n.drawn=n.sig,
                 over_represented_pvalue=pval,
                 under_represented_pvalue=pval.under)]
  all.equal(rc, my.res, info="multiGSEA(method='goseq')", check.attributes=FALSE)
})

test_that("goseq hypergeometric test is like do.hyperGeometricTest", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  gsd <- conform(gsd, vm)

  mg <- multiGSEA(gsd, vm, vm$design)
  lfc <- logFC(mg)
  selected <- subset(lfc, significant)$featureId
  universe <- rownames(vm)

  mygoh <- multiGSEA::goseq(gsd, selected, universe, method='Hypergeometric')

  myhyp <- multiGSEA::hyperGeometricTest(gsd, selected, universe)
  if (FALSE) {
    plot(-log10(mygoh$over_represented_pvalue), -log10(myhyp$pval))
    abline(0,1,col='red')
  }

  expect_equal(mygoh$over_represented_pvalue, myhyp$pval, tolerance=0.015)
})
